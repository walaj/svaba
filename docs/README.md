# SvABA BPS Viewer

This folder contains a standalone static viewer for SvABA `bps.txt` and `bps.txt.gz` files.

## Files

- `index.html`: app shell
- `styles.css`: layout and visual styling
- `app.js`: parser, filters, table, and detail pane
- `svaba/am.bps.txt.gz`: bundled example file already present in this repo

## Usage

You can open `index.html` directly and use the file picker for local files.

If you want the bundled example button to work, serve this folder over HTTP:

```bash
cd /Users/jeremiahwala/git/svaba/viewer
python3 -m http.server 8000
```

Then open [http://localhost:8000](http://localhost:8000).

## Current filter set

The UI supports:

- free-text search across locus, contig, alleles, confidence, DBSNP, and sample fields
- numeric filtering for `somlod`, `maxlod`, `qual`, split support, ALT support, and discordant support
- button filters for derived event type, raw evidence (`type`), confidence (`conf`), and somatic state
- export of the currently filtered rows as TSV

## Notes

- The raw header contains duplicate `alt` names. In the viewer these are treated as `alt_allele` and `alt_count`.
- The raw `type` column is evidence (`INDEL`, `ASDIS`, `DSCRD`, `ASSMB`). The viewer also derives a separate event class to make filtering more usable.
