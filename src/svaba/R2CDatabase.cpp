#include "R2CDatabase.h"
#include <stdexcept>
#include <string>

// SQLite3 is an OPTIONAL dependency (see R2CDatabase.h and CMakeLists.txt).
// When -DSVABA_HAS_SQLITE3 is defined we pull in sqlite3.h and provide
// the real implementation. Otherwise this translation unit compiles
// without sqlite3 — the class still exists so other files can reference
// it, but every method is a stub. The ctor throws so any caller that
// bypasses the available() gate fails loudly instead of silently
// dropping data; SvabaOptions emits a single startup warning at
// --dump-reads parse time so users know up front.

#ifdef SVABA_HAS_SQLITE3

#include "sqlite3.h"

bool R2CDatabase::available() { return true; }

R2CDatabase::R2CDatabase(const std::string& path) {
  int rc = sqlite3_open(path.c_str(), &db_);
  if (rc != SQLITE_OK) {
    std::string msg = "R2CDatabase: failed to open " + path + ": " +
                      sqlite3_errmsg(db_);
    sqlite3_close(db_);
    db_ = nullptr;
    throw std::runtime_error(msg);
  }
  open_ = true;

  // Write-speed pragmas (write-once database, crash safety not needed).
  exec("PRAGMA page_size=8192;");
  exec("PRAGMA journal_mode=WAL;");
  exec("PRAGMA synchronous=OFF;");
  exec("PRAGMA cache_size=-64000;");

  // Create tables.
  exec(
    "CREATE TABLE IF NOT EXISTS contigs ("
    "  cname      TEXT PRIMARY KEY,"
    "  contig_len INT,"
    "  seq        TEXT,"
    "  frags      TEXT,"
    "  bps        TEXT,"
    "  n_reads    INT"
    ");"
  );

  exec(
    "CREATE TABLE IF NOT EXISTS reads ("
    "  cname        TEXT,"
    "  read_id      TEXT,"
    "  chrom        TEXT,"
    "  pos          INT,"
    "  flag         INT,"
    "  r2c_cigar    TEXT,"
    "  r2c_start    INT,"
    "  r2c_rc       INT,"
    "  r2c_nm       INT,"
    "  support      TEXT,"
    "  split_bps    TEXT,"
    "  disc_bps     TEXT,"
    "  r2c_score    REAL,"  // alignment score under r2c CIGAR (svaba::readAlignmentScore)
    "  native_score REAL,"  // alignment score under BAM/corrected CIGAR
    "  seq          TEXT"
    ");"
  );

  exec("CREATE INDEX IF NOT EXISTS idx_reads_cname ON reads(cname);");

  // Prepare insert statements.
  const char* sql_contig =
    "INSERT INTO contigs (cname, contig_len, seq, frags, bps, n_reads) "
    "VALUES (?1, ?2, ?3, ?4, ?5, ?6);";

  rc = sqlite3_prepare_v2(db_, sql_contig, -1, &stmt_contig_, nullptr);
  if (rc != SQLITE_OK) {
    throw std::runtime_error(
      std::string("R2CDatabase: prepare contig insert: ") + sqlite3_errmsg(db_));
  }

  const char* sql_read =
    "INSERT INTO reads (cname, read_id, chrom, pos, flag, r2c_cigar, "
    "r2c_start, r2c_rc, r2c_nm, support, split_bps, disc_bps, "
    "r2c_score, native_score, seq) "
    "VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15);";

  rc = sqlite3_prepare_v2(db_, sql_read, -1, &stmt_read_, nullptr);
  if (rc != SQLITE_OK) {
    throw std::runtime_error(
      std::string("R2CDatabase: prepare read insert: ") + sqlite3_errmsg(db_));
  }

  // Begin the write transaction.
  exec("BEGIN TRANSACTION;");
}

R2CDatabase::~R2CDatabase() {
  close();
}

void R2CDatabase::insert_contig(const std::string& cname,
                                int contig_len,
                                const std::string& seq,
                                const std::string& frags,
                                const std::string& bps,
                                int n_reads) {
  sqlite3_reset(stmt_contig_);
  sqlite3_bind_text(stmt_contig_, 1, cname.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_int(stmt_contig_, 2, contig_len);
  sqlite3_bind_text(stmt_contig_, 3, seq.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt_contig_, 4, frags.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt_contig_, 5, bps.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_int(stmt_contig_, 6, n_reads);

  int rc = sqlite3_step(stmt_contig_);
  if (rc != SQLITE_DONE) {
    throw std::runtime_error(
      std::string("R2CDatabase: insert_contig failed: ") + sqlite3_errmsg(db_));
  }
}

void R2CDatabase::insert_read(const std::string& cname,
                              const std::string& read_id,
                              const std::string& chrom,
                              int pos,
                              int flag,
                              const std::string& r2c_cigar,
                              int r2c_start,
                              int r2c_rc,
                              int r2c_nm,
                              const std::string& support,
                              const std::string& split_bps,
                              const std::string& disc_bps,
                              double r2c_score,
                              double native_score,
                              const std::string& seq) {
  sqlite3_reset(stmt_read_);
  sqlite3_bind_text(stmt_read_, 1, cname.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt_read_, 2, read_id.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt_read_, 3, chrom.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_int(stmt_read_, 4, pos);
  sqlite3_bind_int(stmt_read_, 5, flag);
  sqlite3_bind_text(stmt_read_, 6, r2c_cigar.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_int(stmt_read_, 7, r2c_start);
  sqlite3_bind_int(stmt_read_, 8, r2c_rc);
  sqlite3_bind_int(stmt_read_, 9, r2c_nm);
  sqlite3_bind_text(stmt_read_, 10, support.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt_read_, 11, split_bps.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt_read_, 12, disc_bps.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_double(stmt_read_, 13, r2c_score);
  sqlite3_bind_double(stmt_read_, 14, native_score);
  sqlite3_bind_text(stmt_read_, 15, seq.c_str(), -1, SQLITE_TRANSIENT);

  int rc = sqlite3_step(stmt_read_);
  if (rc != SQLITE_DONE) {
    throw std::runtime_error(
      std::string("R2CDatabase: insert_read failed: ") + sqlite3_errmsg(db_));
  }
}

void R2CDatabase::commit() {
  if (open_ && db_) {
    exec("COMMIT;");
  }
}

void R2CDatabase::merge_from(const std::string& other_path) {
  if (!open_ || !db_)
    throw std::runtime_error("R2CDatabase::merge_from: database not open");

  // ATTACH the other database under a fixed alias `src`. We always
  // detach in a finally-style cleanup so a partial INSERT failure
  // doesn't leave stray attachments on the connection.
  std::string attach = "ATTACH DATABASE '";
  // Single-quote escaping for SQL string literal: double any embedded '.
  for (char c : other_path) {
    attach.push_back(c);
    if (c == '\'') attach.push_back('\'');
  }
  attach += "' AS src;";
  exec(attach.c_str());

  bool inserted_ok = false;
  try {
    // INSERT OR IGNORE on contigs because cname is the PRIMARY KEY:
    // duplicates only arise from accidental re-runs and dropping them
    // is safe (each contig is owned by exactly one worker).
    exec("INSERT OR IGNORE INTO main.contigs SELECT * FROM src.contigs;");
    // reads has no uniqueness constraint, so a plain INSERT is fine.
    exec("INSERT INTO main.reads SELECT * FROM src.reads;");
    inserted_ok = true;
  } catch (...) {
    // Fall through to DETACH; rethrow after cleanup so the caller still
    // gets the original error.
  }

  // Always detach, even on failure.
  try { exec("DETACH DATABASE src;"); } catch (...) { /* ignore */ }

  if (!inserted_ok)
    throw std::runtime_error("R2CDatabase::merge_from: insert failed for " +
                             other_path);
}

void R2CDatabase::close() {
  if (!open_) return;

  if (stmt_contig_) {
    sqlite3_finalize(stmt_contig_);
    stmt_contig_ = nullptr;
  }
  if (stmt_read_) {
    sqlite3_finalize(stmt_read_);
    stmt_read_ = nullptr;
  }
  if (db_) {
    sqlite3_close(db_);
    db_ = nullptr;
  }
  open_ = false;
}

void R2CDatabase::exec(const char* sql) {
  char* err = nullptr;
  int rc = sqlite3_exec(db_, sql, nullptr, nullptr, &err);
  if (rc != SQLITE_OK) {
    std::string msg = std::string("R2CDatabase: exec failed: ") +
                      (err ? err : "unknown error") + " [SQL: " + sql + "]";
    sqlite3_free(err);
    throw std::runtime_error(msg);
  }
}

void R2CDatabase::stamp_pass_cnames(
    const std::string& db_path,
    const std::unordered_set<std::string>& pass_cnames,
    const std::unordered_set<std::string>& som_cnames) {

  // We bypass the R2CDatabase ctor here because (a) it's tailored for
  // the contigs/reads insert paths and (b) we want a fresh transaction
  // around the pass_cnames table only — opening the existing file via
  // the ctor would re-CREATE the contigs/reads tables (harmless) and
  // re-prepare insert statements (wasted setup). Direct sqlite3 keeps
  // intent local.
  sqlite3* db = nullptr;
  if (sqlite3_open(db_path.c_str(), &db) != SQLITE_OK) {
    const std::string msg = std::string("stamp_pass_cnames: open failed: ") +
                            (db ? sqlite3_errmsg(db) : "unknown");
    sqlite3_close(db);
    throw std::runtime_error(msg);
  }

  auto run = [&](const char* sql) {
    char* err = nullptr;
    if (sqlite3_exec(db, sql, nullptr, nullptr, &err) != SQLITE_OK) {
      const std::string msg = err ? err : "unknown error";
      sqlite3_free(err);
      sqlite3_close(db);
      throw std::runtime_error(
        std::string("stamp_pass_cnames: ") + msg + " [" + sql + "]");
    }
  };

  // Match the write-speed pragmas R2CDatabase uses on the per-thread
  // db files; this is a small bulk insert (typically <100k rows even
  // on noisy WGS) but no reason not to be fast.
  run("PRAGMA journal_mode=WAL;");
  run("PRAGMA synchronous=OFF;");
  run("DROP TABLE IF EXISTS pass_cnames;");
  run(
    "CREATE TABLE pass_cnames ("
    "  cname    TEXT PRIMARY KEY,"
    "  somatic  INT NOT NULL DEFAULT 0"
    ");"
  );
  run("BEGIN TRANSACTION;");

  sqlite3_stmt* ins = nullptr;
  if (sqlite3_prepare_v2(db,
        "INSERT INTO pass_cnames (cname, somatic) VALUES (?1, ?2);",
        -1, &ins, nullptr) != SQLITE_OK) {
    sqlite3_close(db);
    throw std::runtime_error("stamp_pass_cnames: prepare insert failed");
  }

  for (const auto& cname : pass_cnames) {
    sqlite3_reset(ins);
    sqlite3_bind_text(ins, 1, cname.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int (ins, 2, som_cnames.count(cname) ? 1 : 0);
    if (sqlite3_step(ins) != SQLITE_DONE) {
      const std::string msg = sqlite3_errmsg(db);
      sqlite3_finalize(ins);
      sqlite3_close(db);
      throw std::runtime_error("stamp_pass_cnames: insert failed: " + msg);
    }
  }
  sqlite3_finalize(ins);

  run("COMMIT;");
  sqlite3_close(db);
}

#else  // !SVABA_HAS_SQLITE3 — stub implementations -----------------------------

bool R2CDatabase::available() { return false; }

R2CDatabase::R2CDatabase(const std::string& path) {
  // svaba was built without sqlite3 support. The expected flow is:
  //   1. SvabaOptions parses --dump-reads
  //   2. emits a one-line warning when !available()
  //   3. SvabaThreadUnit::ctor gates allocation on available()
  //   4. so this branch never runs in practice
  // We throw rather than silently no-op to make any future caller that
  // bypasses the gate fail visibly with a real error.
  throw std::runtime_error(
    "R2CDatabase: svaba was built without sqlite3 support; rebuild "
    "with libsqlite3-dev / sqlite-devel / `brew install sqlite` to "
    "use --dump-reads (path: " + path + ")");
}

R2CDatabase::~R2CDatabase() = default;

void R2CDatabase::insert_contig(const std::string&, int, const std::string&,
                                const std::string&, const std::string&, int) {}
void R2CDatabase::insert_read(const std::string&, const std::string&,
                              const std::string&, int, int, const std::string&,
                              int, int, int, const std::string&,
                              const std::string&, const std::string&,
                              double, double, const std::string&) {}
void R2CDatabase::commit() {}
void R2CDatabase::merge_from(const std::string&) {}
void R2CDatabase::close() {}
void R2CDatabase::exec(const char*) {}

void R2CDatabase::stamp_pass_cnames(
    const std::string&,
    const std::unordered_set<std::string>&,
    const std::unordered_set<std::string>&) {
  // No-op in stub mode. Postprocess gates this on
  // R2CDatabase::available() so it never reaches us; defending just
  // in case a future caller forgets.
}

#endif // SVABA_HAS_SQLITE3
