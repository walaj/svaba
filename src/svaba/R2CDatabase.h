#ifndef SVABA_R2C_DATABASE_H
#define SVABA_R2C_DATABASE_H

#include <string>
#include <unordered_set>

// Write-only SQLite emitter for the r2c alignment data.
// Replaces the gzipped TSV with a queryable database that supports
// O(1) contig lookup and efficient batch filtering by cname set.
//
// SQLite is an OPTIONAL build dependency. When CMake's
// find_package(SQLite3) succeeds, svaba is compiled with the
// SVABA_HAS_SQLITE3 macro defined and this class is fully functional.
// When sqlite3 isn't available, this header still compiles (so
// downstream files can include it unconditionally and reference the
// type), but the implementations are stubs:
//   - available()  returns false
//   - the ctor throws std::runtime_error
//   - all other methods are no-ops
// Callers are expected to gate construction on `available()` so the
// throw never fires; the gate also lets the option parser emit a
// single startup warning when --dump-reads is used in a no-sqlite
// build, instead of N per-thread warnings inside ctors.
//
// Usage (in a sqlite3-enabled build):
//   if (R2CDatabase::available()) {
//     R2CDatabase db("output.r2c.db");
//     db.insert_contig(...);
//     db.insert_read(...);
//     db.commit();
//     // destructor calls close()
//   }

// Forward-declare the sqlite types so the header doesn't need to
// include sqlite3.h. In a no-sqlite build these never get instantiated;
// the fields below still take up space (a few pointers per instance),
// but the class is only ever constructed once per worker thread, so
// the overhead is negligible.
struct sqlite3;
struct sqlite3_stmt;

class R2CDatabase {
public:
  // Returns true iff svaba was compiled with sqlite3 support
  // (i.e. -DSVABA_HAS_SQLITE3). Call this before constructing an
  // R2CDatabase — if false, construction will throw.
  static bool available();

  // Stamp a `pass_cnames(cname TEXT PK, somatic INT)` lookup table
  // into an existing .r2c.db file. Used by `svaba postprocess` after
  // Phase 3 identifies PASS / PASS+somatic cname sets from
  // bps.txt.gz. Idempotent (DROP TABLE IF EXISTS first).
  //
  // No-op if available() is false. Throws std::runtime_error on any
  // SQLite error in the available() == true path.
  static void stamp_pass_cnames(
    const std::string& db_path,
    const std::unordered_set<std::string>& pass_cnames,
    const std::unordered_set<std::string>& som_cnames);

  // Opens (or creates) the database at `path`, creates tables if they
  // don't exist, sets write-speed pragmas, and begins a transaction.
  // Throws std::runtime_error in a no-sqlite build (callers must gate
  // on available()).
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
                   double r2c_score,
                   double native_score,
                   const std::string& seq);

  // Commits the open transaction.
  void commit();

  // Merge another R2CDatabase file (path) into this one. Used by
  // postprocess to coalesce the per-thread .r2c.db files into the final
  // ${ID}.r2c.db. Implementation: ATTACH DATABASE the other file, INSERT
  // INTO ... SELECT * FROM the attached tables, then DETACH. Caller is
  // responsible for ensuring this database is in a transaction (the
  // ctor begins one); merge_from does NOT commit, so multiple merges
  // accumulate in a single fast transaction.
  //
  // Duplicate cnames in the contigs table are silently ignored
  // (INSERT OR IGNORE) — each contig is owned by one worker thread, so
  // duplicates only arise from accidental double-runs and dropping
  // them is safe.
  void merge_from(const std::string& other_path);

  // Finalizes prepared statements and closes the database.
  // Safe to call multiple times.
  void close();

private:
  void exec(const char* sql);

  sqlite3*      db_          = nullptr;
  sqlite3_stmt* stmt_contig_ = nullptr;
  sqlite3_stmt* stmt_read_   = nullptr;
  bool          open_        = false;
};

#endif // SVABA_R2C_DATABASE_H
