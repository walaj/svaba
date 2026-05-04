#ifndef SVABA_R2C_DATABASE_H
#define SVABA_R2C_DATABASE_H

#include <string>
#include "sqlite3.h"

// Write-only SQLite emitter for the r2c alignment data.
// Replaces the gzipped TSV with a queryable database that supports
// O(1) contig lookup and efficient batch filtering by cname set.
//
// Usage:
//   R2CDatabase db("output.r2c.db");
//   db.insert_contig(...);
//   db.insert_read(...);
//   db.commit();
//   // destructor calls close()

class R2CDatabase {
public:
  // Opens (or creates) the database at `path`, creates tables if they
  // don't exist, sets write-speed pragmas, and begins a transaction.
  explicit R2CDatabase(const std::string& path);

  ~R2CDatabase();

  // Non-copyable, non-movable (owns sqlite handles).
  R2CDatabase(const R2CDatabase&) = delete;
  R2CDatabase& operator=(const R2CDatabase&) = delete;

  void insert_contig(const std::string& cname,
                     int contig_len,
                     const std::string& seq,
                     const std::string& frags,
                     const std::string& bps,
                     int n_reads);

  void insert_read(const std::string& cname,
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
                   const std::string& seq);

  // Commits the open transaction.
  void commit();

  // Finalizes prepared statements and closes the database.
  // Safe to call multiple times.
  void close();

private:
  void exec(const char* sql);

  sqlite3* db_ = nullptr;
  sqlite3_stmt* stmt_contig_ = nullptr;
  sqlite3_stmt* stmt_read_ = nullptr;
  bool open_ = false;
};

#endif // SVABA_R2C_DATABASE_H
