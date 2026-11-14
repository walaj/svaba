/* ===========================================================================
 * truth_store.js -- local, file://-friendly "ground truth" labeling for the
 * svaba variant viewers (bps_explorer.html, comparison.html, ...).
 *
 * Click a label (Valid / Artifact / Germline / Ambiguous) while reviewing a
 * variant; it is saved silently + immediately to the browser (localStorage) and
 * -- if you link a file once -- AUTO-SAVED to a real ground_truth.json on disk
 * (e.g. checked into this repo) on every change, via the File System Access API.
 *
 * Assumptions (per request): a single fixed test BAM, so a variant is keyed
 * purely by breakend coordinates + coarse type (NOT svaba version / bp_id /
 * cname). Keys are orientation- and strand-agnostic (A/B order canonicalized)
 * with a +/-POS_TOL fuzzy fallback so minor cross-version coord shifts still
 * resolve to the same label.
 *
 * Persistence layers:
 *   1. localStorage  -- always on; instant; survives sessions on this machine.
 *   2. linked file   -- optional; Chrome/Edge; auto-saves ground_truth.json on
 *                       disk so it can be committed. Handle persists in
 *                       IndexedDB; one permission re-grant click per session.
 *   3. Export/Import -- manual JSON download/upload; works in every browser.
 *
 * No server, no build deps. Include before the viewer's main script; use
 * window.SvabaTruth.
 * =========================================================================== */
window.SvabaTruth = (function () {
  "use strict";
  var LSKEY = "svaba_truth_v1";
  var POS_TOL = 10;
  var LABELS = [
    { id: "valid",     name: "Valid",     color: "#16a34a" },
    { id: "artifact",  name: "Artifact",  color: "#dc2626" },
    { id: "germline",  name: "Germline",  color: "#2563eb" },
    { id: "ambiguous", name: "Ambiguous", color: "#d97706" }
  ];
  var COLOR = {}; LABELS.forEach(function (L) { COLOR[L.id] = L.color; });

  var store = load();
  var idx = buildIndex();
  var notify = null;             // toolbar refresh hook (set by makeToolbar)

  function load() { try { return JSON.parse(localStorage.getItem(LSKEY)) || {}; } catch (e) { return {}; } }
  function persist() { try { localStorage.setItem(LSKEY, JSON.stringify(store)); } catch (e) { console.warn("SvabaTruth: localStorage write failed", e); } }
  function typeClass(t) {
    t = (t || "").toUpperCase();
    if (t.indexOf("DEL") === 0) return "DEL";
    if (t.indexOf("DUP") === 0) return "DUP";
    if (t.indexOf("INV") === 0) return "INV";
    if (t.indexOf("INS") === 0) return "INS";
    if (t.indexOf("BND") === 0 || t.indexOf("TRA") === 0) return "BND";
    if (t.indexOf("IND") === 0) return "INDEL";
    return t || "NA";
  }
  function key(c1, p1, c2, p2, type) {
    var a = c1 + ":" + (+p1), b = c2 + ":" + (+p2);
    var x = (a <= b) ? a : b, y = (a <= b) ? b : a;
    return x + "__" + y + "__" + typeClass(type);
  }
  function buildIndex() {
    var ix = {};
    for (var k in store) {
      if (!store.hasOwnProperty(k)) continue;
      var first = k.split("__")[0], sp = first.split(":");
      (ix[sp[0]] = ix[sp[0]] || []).push([+sp[1], k]);
    }
    for (var cc in ix) ix[cc].sort(function (a, b) { return a[0] - b[0]; });
    return ix;
  }
  function get(c1, p1, c2, p2, type) {
    var k = key(c1, p1, c2, p2, type);
    if (store[k]) return { key: k, label: store[k].label, exact: true };
    var tc = typeClass(type), ends = [[c1, +p1], [c2, +p2]];
    for (var i = 0; i < ends.length; i++) {
      var arr = idx[ends[i][0]]; if (!arr) continue;
      for (var j = 0; j < arr.length; j++) {
        var pos = arr[j][0], kk = arr[j][1];
        if ((Math.abs(pos - (+p1)) <= POS_TOL || Math.abs(pos - (+p2)) <= POS_TOL)
            && kk.lastIndexOf("__" + tc) === kk.length - tc.length - 2) {
          return { key: k, label: store[kk].label, exact: false, matchedKey: kk };
        }
      }
    }
    return { key: k, label: null };
  }
  function set(c1, p1, c2, p2, type, label, meta) {
    var k = key(c1, p1, c2, p2, type);
    if (!label) { if (store[k]) delete store[k]; }
    else { var rec = { label: label, ts: Date.now() }; if (meta) for (var m in meta) rec[m] = meta[m]; store[k] = rec; }
    persist(); idx = buildIndex(); scheduleSave();
    return k;
  }
  function summary() {
    var s = { total: 0 }; LABELS.forEach(function (L) { s[L.id] = 0; });
    for (var k in store) { if (!store.hasOwnProperty(k)) continue; s.total++; var l = store[k].label; if (l in s) s[l]++; }
    return s;
  }
  function jsonText() {
    return JSON.stringify({
      version: "svaba_truth_v1",
      _note: "svaba ground-truth labels for the shared test BAM. Auto-saved by docs/truth_store.js. Keyed by breakend coords + coarse type; orientation/strand-agnostic. Commit to share.",
      exported: new Date().toISOString(),
      labels: store
    }, null, 1);
  }
  function exportJSON(filename) {
    var blob = new Blob([jsonText()], { type: "application/json" });
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = filename || "ground_truth.json";
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(function () { URL.revokeObjectURL(a.href); }, 2000);
  }
  // newest-ts-wins merge (used by import + file link)
  function mergeNewer(incoming) {
    var n = 0;
    for (var k in incoming) {
      if (!incoming.hasOwnProperty(k)) continue;
      var v = incoming[k], rec = (typeof v === "string") ? { label: v, ts: 0 } : v;
      var cur = store[k];
      if (!cur || (rec.ts || 0) >= (cur.ts || 0)) { store[k] = rec; n++; }
    }
    persist(); idx = buildIndex(); return n;
  }
  function importJSON(obj, mode) {
    var incoming = (obj && obj.labels) ? obj.labels : (obj || {});
    if (mode === "replace") store = {};
    var n = mergeNewer(incoming); scheduleSave(); return n;
  }
  function clearAll() { store = {}; persist(); idx = {}; scheduleSave(); }
  function colorOf(label) { return COLOR[label] || "#888"; }
  function labelOf(c1, p1, c2, p2, type) { return get(c1, p1, c2, p2, type).label; }

  /* ---- File System Access auto-save ------------------------------------- */
  var fsHandle = null, fsAuto = false, fsTimer = null, fsHasStored = false;
  function fsSupported() { return typeof window !== "undefined" && "showOpenFilePicker" in window; }

  function idbOpen() {
    return new Promise(function (res, rej) {
      var r = indexedDB.open("svaba_truth_db", 1);
      r.onupgradeneeded = function () { r.result.createObjectStore("kv"); };
      r.onsuccess = function () { res(r.result); };
      r.onerror = function () { rej(r.error); };
    });
  }
  function idbPut(k, v) { return idbOpen().then(function (db) { return new Promise(function (res, rej) { var tx = db.transaction("kv", "readwrite"); tx.objectStore("kv").put(v, k); tx.oncomplete = function () { db.close(); res(); }; tx.onerror = function () { rej(tx.error); }; }); }); }
  function idbGet(k) { return idbOpen().then(function (db) { return new Promise(function (res, rej) { var tx = db.transaction("kv", "readonly"); var rq = tx.objectStore("kv").get(k); rq.onsuccess = function () { db.close(); res(rq.result); }; rq.onerror = function () { rej(rq.error); }; }); }); }
  function idbDel(k) { return idbOpen().then(function (db) { return new Promise(function (res) { var tx = db.transaction("kv", "readwrite"); tx.objectStore("kv").delete(k); tx.oncomplete = function () { db.close(); res(); }; }); }); }

  function ensurePerm(handle, mode) {
    var opts = { mode: mode || "readwrite" };
    return Promise.resolve(handle.queryPermission(opts)).then(function (p) {
      if (p === "granted") return true;
      return Promise.resolve(handle.requestPermission(opts)).then(function (q) { return q === "granted"; });
    });
  }
  function writeFile() {
    if (!fsHandle) return Promise.resolve(false);
    return ensurePerm(fsHandle, "readwrite").then(function (ok) {
      if (!ok) return false;
      return fsHandle.createWritable().then(function (w) {
        return w.write(jsonText()).then(function () { return w.close(); }).then(function () { return true; });
      });
    }).catch(function (e) { console.warn("SvabaTruth: auto-save failed", e); return false; });
  }
  function scheduleSave() {
    if (!fsAuto || !fsHandle) return;
    if (fsTimer) clearTimeout(fsTimer);
    fsTimer = setTimeout(function () { fsTimer = null; writeFile().then(function () { if (notify) notify(); }); }, 400);
  }
  function readHandleFile(handle) {
    return handle.getFile().then(function (f) { return f.text(); }).then(function (t) {
      if (!t || !t.trim()) return 0;
      var obj = JSON.parse(t); return mergeNewer(obj.labels || obj);
    });
  }
  // Link an EXISTING file (pick ground_truth.json): load+merge, then auto-save.
  function linkOpen() {
    if (!fsSupported()) return Promise.reject(new Error("File System Access API not available in this browser; use Export/Import."));
    return window.showOpenFilePicker({ multiple: false, types: [{ description: "JSON", accept: { "application/json": [".json"] } }] })
      .then(function (picks) { fsHandle = picks[0]; return ensurePerm(fsHandle, "readwrite"); })
      .then(function (ok) { if (!ok) throw new Error("permission denied"); return readHandleFile(fsHandle); })
      .then(function (n) { return idbPut("gt_handle", fsHandle).then(function () { fsAuto = true; fsHasStored = true; return writeFile().then(function () { return n; }); }); });
  }
  // Reconnect a previously-linked file after a reload (one click / session).
  function reconnect() {
    return idbGet("gt_handle").then(function (h) {
      if (!h) throw new Error("no stored file");
      fsHandle = h; return ensurePerm(fsHandle, "readwrite");
    }).then(function (ok) { if (!ok) throw new Error("permission denied"); return readHandleFile(fsHandle); })
      .then(function (n) { fsAuto = true; return writeFile().then(function () { return n; }); });
  }
  function unlink() { fsAuto = false; fsHandle = null; return idbDel("gt_handle"); }
  function fsState() { return fsAuto ? "linked" : (fsHasStored ? "reconnect" : "off"); }
  // detect a previously linked file on load (so the toolbar offers Reconnect)
  if (fsSupported()) { idbGet("gt_handle").then(function (h) { if (h) { fsHasStored = true; if (notify) notify(); } }).catch(function () {}); }

  /* ---- UI --------------------------------------------------------------- */
  function injectCSS() {
    if (document.getElementById("svaba-truth-css")) return;
    var css =
      ".truth-widget{display:flex;align-items:center;gap:6px;flex-wrap:wrap;margin:2px 0 9px}" +
      ".truth-lbl{font-size:.7rem;opacity:.7;margin-right:2px;font-weight:600}" +
      ".truth-btn{font:inherit;font-size:.7rem;line-height:1;padding:4px 10px;border:1px solid #999;border-radius:6px;background:transparent;color:inherit;cursor:pointer;transition:all .08s}" +
      ".truth-btn:hover{filter:brightness(1.2)}" +
      ".truth-bar{display:flex;align-items:center;gap:8px;flex-wrap:wrap;font-size:.72rem}" +
      ".truth-dot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:3px;vertical-align:middle}" +
      ".truth-fuzzy{font-size:.62rem;opacity:.6;margin-left:4px}" +
      ".truth-sync{font:inherit;font-size:.7rem;padding:4px 10px;border:1px solid #999;border-radius:6px;background:transparent;color:inherit;cursor:pointer}" +
      ".truth-sync.on{border-color:#16a34a;color:#16a34a}.truth-sync.recon{border-color:#d97706;color:#d97706}";
    var s = document.createElement("style"); s.id = "svaba-truth-css"; s.textContent = css;
    (document.head || document.documentElement).appendChild(s);
  }
  injectCSS();

  function makeWidget(c1, p1, c2, p2, type, onSet) {
    var info = get(c1, p1, c2, p2, type), cur = info.label;
    var wrap = document.createElement("div"); wrap.className = "truth-widget";
    var lab = document.createElement("span"); lab.className = "truth-lbl"; lab.textContent = "Ground truth:"; wrap.appendChild(lab);
    var btns = [];
    function paint() {
      btns.forEach(function (o) {
        var on = (cur === o.L.id);
        o.b.style.background = on ? o.L.color : "transparent";
        o.b.style.color = on ? "#fff" : "";
        o.b.style.borderColor = on ? o.L.color : "#999";
        o.b.style.fontWeight = on ? "700" : "400";
      });
    }
    LABELS.forEach(function (L) {
      var b = document.createElement("button");
      b.type = "button"; b.className = "truth-btn"; b.textContent = L.name;
      b.title = "Mark as " + L.name + " (click again to clear). Saved instantly.";
      b.onclick = function () { cur = (cur === L.id) ? null : L.id; set(c1, p1, c2, p2, type, cur); paint(); if (onSet) onSet(cur); };
      btns.push({ b: b, L: L }); wrap.appendChild(b);
    });
    paint();
    if (info.label && !info.exact) {
      var f = document.createElement("span"); f.className = "truth-fuzzy";
      f.textContent = "(~match ±" + POS_TOL + "bp)";
      f.title = "resolved from a nearby labeled variant; clicking re-saves at exact coords";
      wrap.appendChild(f);
    }
    return wrap;
  }

  function makeToolbar(onChange) {
    var bar = document.createElement("div"); bar.className = "truth-bar";
    var txt = document.createElement("span");
    var sync = document.createElement("button"); sync.className = "truth-sync";
    var exp = document.createElement("button"); exp.className = "truth-btn"; exp.textContent = "Export JSON"; exp.onclick = function () { exportJSON(); };
    var imp = document.createElement("button"); imp.className = "truth-btn"; imp.textContent = "Import JSON";
    var fi = document.createElement("input"); fi.type = "file"; fi.accept = ".json,application/json"; fi.style.display = "none";
    imp.onclick = function () { fi.click(); };
    fi.onchange = function (e) {
      var f = e.target.files[0]; if (!f) return; var rd = new FileReader();
      rd.onload = function () { try { var n = importJSON(JSON.parse(rd.result), "merge"); refresh(); if (onChange) onChange(); alert("Imported " + n + " labels (merged into " + summary().total + " total)."); } catch (err) { alert("Import failed: " + err); } };
      rd.readAsText(f); fi.value = "";
    };
    sync.onclick = function () {
      var st = fsState();
      var p = (st === "linked") ? writeFile() : (st === "reconnect" ? reconnect() : linkOpen());
      Promise.resolve(p).then(function () { refresh(); if (onChange) onChange(); })
        .catch(function (e) { if (e && e.name === "AbortError") return; alert("Auto-save link failed: " + (e && e.message || e)); });
    };
    bar.appendChild(txt); bar.appendChild(sync); bar.appendChild(exp); bar.appendChild(imp); bar.appendChild(fi);
    function refresh() {
      var s = summary();
      txt.innerHTML = '<b>Ground truth:</b> ' + s.total + ' labeled &nbsp; ' +
        LABELS.map(function (L) { return '<span class="truth-dot" style="background:' + L.color + '"></span>' + s[L.id]; }).join(" &nbsp; ");
      if (!fsSupported()) { sync.style.display = "none"; }
      else {
        sync.className = "truth-sync";
        var st = fsState();
        if (st === "linked") { sync.classList.add("on"); sync.textContent = "● auto-save → " + (fsHandle ? fsHandle.name : "file"); sync.title = "Auto-saving on every change. Click to flush a save now."; }
        else if (st === "reconnect") { sync.classList.add("recon"); sync.textContent = "⟳ reconnect ground_truth.json"; sync.title = "Re-grant write access to resume auto-save into your linked file"; }
        else { sync.textContent = "○ auto-save to file…"; sync.title = "Pick ground_truth.json (e.g. in the svaba repo) to auto-save into it on every label"; }
      }
    }
    notify = refresh; refresh(); bar._refresh = refresh;
    return bar;
  }

  return {
    LABELS: LABELS, key: key, get: get, set: set, summary: summary,
    exportJSON: exportJSON, importJSON: importJSON, clearAll: clearAll,
    colorOf: colorOf, labelOf: labelOf, makeWidget: makeWidget, makeToolbar: makeToolbar,
    linkOpen: linkOpen, reconnect: reconnect, unlink: unlink, fsState: fsState,
    fsSupported: fsSupported, all: function () { return store; }, POS_TOL: POS_TOL
  };
})();
