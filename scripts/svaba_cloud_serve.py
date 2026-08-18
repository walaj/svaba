#!/usr/bin/env python3
"""
svaba_cloud_serve.py -- tiny local backend for docs/svaba_cloud_launcher.html.

Serves the launcher and a few JSON endpoints that shell out to your *existing*
gcloud/gsutil auth, turning the static page into a live dashboard:
  GET  /                      -> the launcher (token injected so it knows backend is live)
  GET  /api/images?filter=    -> live worker-image list           (gcloud images list)
  GET  /api/progress?ws=gs..  -> {done:[..]} from .done_part markers (gsutil ls)
  GET  /api/vms?prefix=       -> running worker VMs + status       (gcloud instances list)
  GET  /api/joblog?id=        -> live log + running/rc for a launched job
  POST /api/run {cmd}         -> run a generated command in the repo, stream to a log

Security: binds 127.0.0.1 ONLY (no remote access). Mutating calls (/api/run)
require a per-session token that is injected into the page it serves, plus a
same-origin check -- so another website in your browser can't trigger commands.
It runs as you, with your shell + gcloud creds: only ever point it at commands
you'd run yourself.

Run:
  python3 scripts/svaba_cloud_serve.py [--port 8787] [--project osteosarc]
then open the printed http://127.0.0.1:PORT/ URL.
"""
import argparse, http.server, json, os, re, secrets, socketserver, subprocess, urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
HTML = os.path.join(REPO, "docs", "svaba_cloud_launcher.html")

ap = argparse.ArgumentParser()
ap.add_argument("--port", type=int, default=int(os.environ.get("PORT", "8787")))
ap.add_argument("--project", default=os.environ.get("PROJECT", ""))
ARGS = ap.parse_args()
TOKEN = secrets.token_urlsafe(18)
JOBS = {}   # id -> {proc, logpath, cmd}


def run(args, timeout=40):
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
        return p.stdout, p.returncode, p.stderr
    except Exception as e:
        return "", 1, str(e)


def proj(args):
    return args + ([f"--project={ARGS.project}"] if ARGS.project else [])


class H(http.server.BaseHTTPRequestHandler):
    def log_message(self, *a):  # quiet
        pass

    def _send(self, code, body, ctype="application/json"):
        b = body if isinstance(body, bytes) else body.encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(b)))
        self.end_headers()
        self.wfile.write(b)

    def _json(self, code, obj):
        self._send(code, json.dumps(obj))

    def do_GET(self):
        u = urllib.parse.urlparse(self.path)
        q = urllib.parse.parse_qs(u.query)
        if u.path in ("/", "/index.html"):
            try:
                html = open(HTML).read()
                inj = ("<script>window.SVABA_BACKEND=true;window.SVABA_TOKEN="
                       + json.dumps(TOKEN) + ";</script></head>")
                return self._send(200, html.replace("</head>", inj, 1), "text/html")
            except Exception as e:
                return self._send(500, str(e), "text/plain")

        if u.path == "/api/images":
            filt = (q.get("filter") or ["svaba"])[0]
            out, rc, err = run(proj(["gcloud", "compute", "images", "list",
                                     "--no-standard-images", "--filter=name~" + filt,
                                     "--format=value(name)"]))
            imgs = [l.strip() for l in out.splitlines() if l.strip()]
            return self._json(200, {"ok": rc == 0, "images": imgs, "err": err.strip()})

        if u.path == "/api/progress":
            ws = (q.get("ws") or [""])[0].rstrip("/")
            if not ws:
                return self._json(400, {"ok": False, "err": "ws required"})
            out, rc, err = run(["gsutil", "ls", ws + "/.done_part*"])
            done = sorted({int(m.group(1)) for m in re.finditer(r"\.done_part(\d+)", out)})
            return self._json(200, {"ok": True, "done": done})

        if u.path == "/api/vms":
            pref = (q.get("prefix") or ["svaba"])[0]
            out, rc, err = run(proj(["gcloud", "compute", "instances", "list",
                                     "--filter=name~" + pref,
                                     "--format=value(name,status,zone.basename())"]))
            vms = []
            for l in out.splitlines():
                if not l.strip():
                    continue
                p = l.split("\t")
                vms.append({"name": p[0], "status": p[1] if len(p) > 1 else "",
                            "zone": p[2] if len(p) > 2 else ""})
            return self._json(200, {"ok": rc == 0, "vms": vms, "err": err.strip()})

        if u.path == "/api/serial":
            vm = (q.get("vm") or [""])[0]
            zone = (q.get("zone") or [""])[0]
            if not vm:
                return self._json(400, {"ok": False, "err": "vm required"})
            # 1) preferred: tail the worker's log file over SSH (carries live svaba stderr)
            ssh = ["gcloud", "compute", "ssh", vm,
                   "--command", "sudo tail -c 120000 /var/log/svaba_startup.log 2>/dev/null",
                   "--ssh-flag=-oConnectTimeout=10", "--ssh-flag=-oStrictHostKeyChecking=no"]
            if zone:
                ssh.append("--zone=" + zone)
            out, rc, err = run(proj(ssh), timeout=55)
            source = "ssh /var/log/svaba_startup.log"
            if not out.strip():
                # 2) fallback: serial console, filtered to startup-script lines
                ser = ["gcloud", "compute", "instances", "get-serial-port-output", vm]
                if zone:
                    ser.append("--zone=" + zone)
                sout, src, serr = run(proj(ser), timeout=60)
                keep = [m.group(1) for ln in sout.splitlines()
                        for m in [re.search(r'command\("/bin/bash"\):\s?(.*)', ln)] if m]
                out = "\n".join(keep) if keep else sout
                err = serr
                source = "serial console"
            out = re.sub(r"\x1b\[[0-9;?]*[ -/]*[@-~]", "", out)[-100000:]
            return self._json(200, {"ok": True, "log": out, "source": source, "err": err.strip()[-200:]})

        if u.path == "/api/joblog":
            jid = (q.get("id") or [""])[0]
            j = JOBS.get(jid)
            if not j:
                return self._json(404, {"ok": False, "err": "no such job"})
            try:
                log = open(j["logpath"]).read()[-40000:]
                log = re.sub(r"\x1b\[[0-9;?]*[ -/]*[@-~]", "", log)  # strip ANSI/control seqs
            except Exception:
                log = ""
            running = j["proc"].poll() is None
            return self._json(200, {"ok": True, "running": running,
                                    "rc": j["proc"].returncode, "cmd": j["cmd"], "log": log})

        return self._send(404, "not found", "text/plain")

    def do_POST(self):
        # token + same-origin guard so a stray website can't trigger commands
        if self.headers.get("X-Svaba-Token") != TOKEN:
            return self._json(403, {"ok": False, "err": "bad/missing token"})
        origin = self.headers.get("Origin", "")
        if origin and not re.match(r"^http://(127\.0\.0\.1|localhost):%d$" % ARGS.port, origin):
            return self._json(403, {"ok": False, "err": "bad origin"})
        ln = int(self.headers.get("Content-Length", "0"))
        try:
            data = json.loads(self.rfile.read(ln) or b"{}")
        except Exception:
            data = {}
        u = urllib.parse.urlparse(self.path)

        if u.path == "/api/run":
            cmd = (data.get("cmd") or "").strip()
            if not cmd:
                return self._json(400, {"ok": False, "err": "empty cmd"})
            jid = secrets.token_hex(4)
            logp = os.path.join("/tmp", "svaba_job_" + jid + ".log")
            lf = open(logp, "w")
            proc = subprocess.Popen(["/bin/bash", "-lc", cmd], cwd=REPO,
                                    stdout=lf, stderr=subprocess.STDOUT)
            JOBS[jid] = {"proc": proc, "logpath": logp, "cmd": cmd}
            return self._json(200, {"ok": True, "id": jid})

        return self._json(404, {"ok": False, "err": "not found"})


class Server(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True


if __name__ == "__main__":
    httpd = Server(("127.0.0.1", ARGS.port), H)
    print("=" * 60)
    print("  svaba cloud launcher backend  ->  http://127.0.0.1:%d/" % ARGS.port)
    print("  project: %s" % (ARGS.project or "(gcloud default)"))
    print("  localhost-only; auth = your gcloud; Ctrl-C to stop")
    print("=" * 60)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nbye")
