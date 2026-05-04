#include "R2CDatabase.h"
#include <stdexcept>
#include <string>

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
    "  r2c_score    INT,"
    "  native_score INT,"
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
                              int r2c_score,
                              int native_score,
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
  sqlite3_bind_int(stmt_read_, 13, r2c_score);
  sqlite3_bind_int(stmt_read_, 14, native_score);
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
