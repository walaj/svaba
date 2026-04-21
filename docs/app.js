(function () {
  "use strict";

  const BASE_COLUMNS = [
    "chr1",
    "pos1",
    "strand1",
    "chr2",
    "pos2",
    "strand2",
    "ref",
    "alt_allele",
    "span",
    "split",
    "alt_count",
    "cov",
    "cigar",
    "cigar_near",
    "dmq1",
    "dmq2",
    "dcn",
    "dct",
    "mapq1",
    "mapq2",
    "nm1",
    "nm2",
    "as1",
    "as2",
    "sub1",
    "sub2",
    "homol",
    "insert",
    "repeat",
    "contig_and_region",
    "naln",
    "conf",
    "evidence",
    "qual",
    "secondary",
    "somatic",
    "somlod",
    "maxlod",
    "dbsnp"
  ];

  const NUMERIC_FIELDS = new Set([
    "pos1",
    "pos2",
    "span",
    "split",
    "alt_count",
    "cov",
    "cigar",
    "cigar_near",
    "dmq1",
    "dmq2",
    "dcn",
    "dct",
    "mapq1",
    "mapq2",
    "nm1",
    "nm2",
    "as1",
    "as2",
    "sub1",
    "sub2",
    "naln",
    "qual",
    "secondary",
    "somlod",
    "maxlod"
  ]);

  const SORT_DEFAULTS = {
    key: "maxlod",
    direction: "desc"
  };

  const state = {
    fileName: "",
    originalHeader: [],
    sampleHeaders: [],
    rows: [],
    filteredRows: [],
    selectedId: null,
    page: 1,
    pageSize: 100,
    sort: { ...SORT_DEFAULTS },
    filters: {
      search: "",
      somlodMin: null,
      somlodMax: null,
      maxlodMin: null,
      maxlodMax: null,
      qualMin: null,
      splitMin: null,
      altMin: null,
      discordantMin: null,
      dbsnp: "all",
      eventTypes: new Set(),
      evidence: new Set(),
      confidence: new Set(),
      somatic: new Set()
    }
  };

  const els = {
    fileInput: document.getElementById("file-input"),
    loadExample: document.getElementById("load-example"),
    exportButton: document.getElementById("export-button"),
    resetFilters: document.getElementById("reset-filters"),
    searchInput: document.getElementById("search-input"),
    somlodMin: document.getElementById("somlod-min"),
    somlodMax: document.getElementById("somlod-max"),
    maxlodMin: document.getElementById("maxlod-min"),
    maxlodMax: document.getElementById("maxlod-max"),
    qualMin: document.getElementById("qual-min"),
    splitMin: document.getElementById("split-min"),
    altMin: document.getElementById("alt-min"),
    discordantMin: document.getElementById("discordant-min"),
    dbsnpFilter: document.getElementById("dbsnp-filter"),
    pageSize: document.getElementById("page-size"),
    tableBody: document.getElementById("table-body"),
    detailContent: document.getElementById("detail-content"),
    statusText: document.getElementById("status-text"),
    resultsTitle: document.getElementById("results-title"),
    resultsSummary: document.getElementById("results-summary"),
    prevPage: document.getElementById("prev-page"),
    nextPage: document.getElementById("next-page"),
    pageSummary: document.getElementById("page-summary"),
    statTotal: document.getElementById("stat-total"),
    statShown: document.getElementById("stat-shown"),
    statPass: document.getElementById("stat-pass"),
    statTopEvidence: document.getElementById("stat-top-evidence"),
    eventTypeChips: document.getElementById("event-type-chips"),
    evidenceChips: document.getElementById("evidence-chips"),
    confidenceChips: document.getElementById("confidence-chips"),
    somaticChips: document.getElementById("somatic-chips")
  };

  document.querySelectorAll(".sort-button").forEach((button) => {
    button.addEventListener("click", () => {
      const key = button.dataset.sort;
      if (state.sort.key === key) {
        state.sort.direction = state.sort.direction === "asc" ? "desc" : "asc";
      } else {
        state.sort.key = key;
        state.sort.direction = key === "locus" || key === "eventType" || key === "conf" || key === "evidence" || key === "contig_and_region" ? "asc" : "desc";
      }
      render();
    });
  });

  document.querySelectorAll(".preset-button").forEach((button) => {
    button.addEventListener("click", () => {
      const target = document.getElementById(button.dataset.target);
      target.value = button.dataset.value;
      syncFilterState();
      render();
    });
  });

  els.fileInput.addEventListener("change", async (event) => {
    const file = event.target.files && event.target.files[0];
    if (!file) {
      return;
    }
    await loadFromBlob(file, file.name);
  });

  els.loadExample.addEventListener("click", async () => {
    await loadBundledExample();
  });

  els.exportButton.addEventListener("click", () => {
    exportFilteredRows();
  });

  els.resetFilters.addEventListener("click", () => {
    resetFilters();
    render();
  });

  els.prevPage.addEventListener("click", () => {
    if (state.page > 1) {
      state.page -= 1;
      renderTable();
    }
  });

  els.nextPage.addEventListener("click", () => {
    const totalPages = Math.max(1, Math.ceil(state.filteredRows.length / state.pageSize));
    if (state.page < totalPages) {
      state.page += 1;
      renderTable();
    }
  });

  [
    "searchInput",
    "somlodMin",
    "somlodMax",
    "maxlodMin",
    "maxlodMax",
    "qualMin",
    "splitMin",
    "altMin",
    "discordantMin",
    "dbsnpFilter",
    "pageSize"
  ].forEach((key) => {
    els[key].addEventListener("input", () => {
      syncFilterState();
      render();
    });
    els[key].addEventListener("change", () => {
      syncFilterState();
      render();
    });
  });

  resetFilters();
  render();

  function resetFilters() {
    state.pageSize = 100;
    state.filters.search = "";
    state.filters.somlodMin = null;
    state.filters.somlodMax = null;
    state.filters.maxlodMin = null;
    state.filters.maxlodMax = null;
    state.filters.qualMin = null;
    state.filters.splitMin = null;
    state.filters.altMin = null;
    state.filters.discordantMin = null;
    state.filters.dbsnp = "all";
    state.filters.eventTypes = new Set();
    state.filters.evidence = new Set();
    state.filters.confidence = new Set();
    state.filters.somatic = new Set();
    state.page = 1;

    els.searchInput.value = "";
    els.somlodMin.value = "";
    els.somlodMax.value = "";
    els.maxlodMin.value = "";
    els.maxlodMax.value = "";
    els.qualMin.value = "";
    els.splitMin.value = "";
    els.altMin.value = "";
    els.discordantMin.value = "";
    els.dbsnpFilter.value = "all";
    els.pageSize.value = String(state.pageSize);
  }

  function syncFilterState() {
    state.filters.search = els.searchInput.value.trim().toLowerCase();
    state.filters.somlodMin = parseNullableNumber(els.somlodMin.value);
    state.filters.somlodMax = parseNullableNumber(els.somlodMax.value);
    state.filters.maxlodMin = parseNullableNumber(els.maxlodMin.value);
    state.filters.maxlodMax = parseNullableNumber(els.maxlodMax.value);
    state.filters.qualMin = parseNullableNumber(els.qualMin.value);
    state.filters.splitMin = parseNullableNumber(els.splitMin.value);
    state.filters.altMin = parseNullableNumber(els.altMin.value);
    state.filters.discordantMin = parseNullableNumber(els.discordantMin.value);
    state.filters.dbsnp = els.dbsnpFilter.value;
    state.pageSize = Number(els.pageSize.value) || 100;
    state.page = 1;
  }

  async function loadBundledExample() {
    setStatus("Loading bundled example from ./svaba/am.bps.txt.gz ...");
    try {
      const response = await fetch("./svaba/am.bps.txt.gz");
      if (!response.ok) {
        throw new Error("Could not fetch ./svaba/am.bps.txt.gz");
      }
      const text = await readResponseText(response, "am.bps.txt.gz");
      ingestText(text, "svaba/am.bps.txt.gz");
    } catch (error) {
      setStatus(
        "Bundled example load failed. If you opened index.html directly as a file, use the file picker or serve viewer/ over HTTP."
      );
      console.error(error);
    }
  }

  async function loadFromBlob(blob, fileName) {
    setStatus("Reading " + fileName + " ...");
    try {
      const text = await readBlobText(blob, fileName);
      ingestText(text, fileName);
    } catch (error) {
      setStatus("Could not read " + fileName + ". " + error.message);
      console.error(error);
    } finally {
      els.fileInput.value = "";
    }
  }

  async function readBlobText(blob, fileName) {
    if (fileName.endsWith(".gz")) {
      if (typeof DecompressionStream === "undefined") {
        throw new Error("This browser does not support gzip decompression. Use a plain .txt file here.");
      }
      const stream = blob.stream().pipeThrough(new DecompressionStream("gzip"));
      return await new Response(stream).text();
    }
    return await blob.text();
  }

  async function readResponseText(response, fileName) {
    if (fileName.endsWith(".gz")) {
      if (typeof DecompressionStream === "undefined") {
        throw new Error("This browser does not support gzip decompression.");
      }
      const stream = response.body.pipeThrough(new DecompressionStream("gzip"));
      return await new Response(stream).text();
    }
    return await response.text();
  }

  function ingestText(text, fileName) {
    const lines = text.replace(/\r\n?/g, "\n").split("\n").filter((line) => line.length > 0);
    if (lines.length < 2) {
      throw new Error("Expected a header and at least one row.");
    }

    const originalHeader = lines[0].split("\t");
    if (originalHeader.length < BASE_COLUMNS.length) {
      throw new Error("Unexpected header width: " + originalHeader.length);
    }

    const sampleHeaders = originalHeader.slice(BASE_COLUMNS.length);
    const rows = [];

    for (let index = 1; index < lines.length; index += 1) {
      const line = lines[index];
      const cells = line.split("\t");
      if (cells.length !== originalHeader.length) {
        continue;
      }
      rows.push(makeRow(index, cells, sampleHeaders));
    }

    state.fileName = fileName;
    state.originalHeader = originalHeader;
    state.sampleHeaders = sampleHeaders;
    state.rows = rows;
    resetFilters();
    state.selectedId = rows.length ? rows[0].id : null;
    state.page = 1;
    state.sort = { ...SORT_DEFAULTS };

    render();
    setStatus("Loaded " + rows.length.toLocaleString() + " calls from " + fileName + ".");
  }

  function makeRow(index, cells, sampleHeaders) {
    const row = { id: index, originalCells: cells };

    BASE_COLUMNS.forEach((field, fieldIndex) => {
      const value = cells[fieldIndex];
      row[field] = NUMERIC_FIELDS.has(field) ? parseMaybeNumber(value) : value;
    });

    row.chr1 = stripLeadingHash(row.chr1);
    row.locus = row.chr1 + ":" + row.pos1 + " -> " + row.chr2 + ":" + row.pos2;
    row.sampleInfo = sampleHeaders.map((header, offset) => parseSampleField(header, cells[BASE_COLUMNS.length + offset], row.evidence));
    row.eventType = deriveEventType(row);
    row.spanAbs = row.span === null || row.span === -1 ? Infinity : Math.abs(row.span);
    row.supportPeak = Math.max(
      zeroIfNull(row.split),
      zeroIfNull(row.alt_count),
      zeroIfNull(row.cigar),
      zeroIfNull(row.dcn),
      zeroIfNull(row.dct)
    );
    row.hasDbsnp = Boolean(row.dbsnp && row.dbsnp !== "x");
    row.searchText = [
      row.locus,
      row.eventType,
      row.evidence,
      row.conf,
      row.ref,
      row.alt_allele,
      row.contig_and_region,
      row.insert,
      row.repeat,
      row.homol,
      row.dbsnp,
      row.somatic
    ]
      .concat(row.sampleInfo.map((sample) => sample.label + " " + sample.raw))
      .join(" ")
      .toLowerCase();
    return row;
  }

  function deriveEventType(row) {
    if (row.evidence === "INDEL") {
      const refLength = (row.ref || "").length;
      const altLength = (row.alt_allele || "").length;
      if (altLength > refLength) {
        return "Insertion";
      }
      if (altLength < refLength) {
        return "Deletion";
      }
      return "Substitution";
    }

    if (row.chr1 !== row.chr2) {
      return "Translocation";
    }

    if (row.strand1 === row.strand2) {
      return "Inversion";
    }

    if (row.strand1 === "-" && row.strand2 === "+") {
      return "Duplication-like";
    }

    if (row.strand1 === "+" && row.strand2 === "-") {
      return "Deletion-like";
    }

    return "Rearrangement";
  }

  function parseSampleField(header, raw, evidence) {
    const label = header.split("_")[0];
    const parts = raw.split(":");
    const isIndel = evidence === "INDEL";
    return {
      label,
      header,
      raw,
      genotype: parts[0] || "",
      altSupport: parseMaybeNumber(parts[1]),
      depth: parseMaybeNumber(parts[2]),
      splitSupport: parseMaybeNumber(parts[3]),
      pairedSupport: parseMaybeNumber(parts[4]),
      supportLabel: isIndel ? "cigar" : "discordant",
      gq: parseMaybeNumber(parts[5]),
      pl: parts[6] || "",
      lod: parseMaybeNumber(parts[7]),
      lodNormal: parseMaybeNumber(parts[8])
    };
  }

  function rebuildChipGroups() {
    renderChipList(els.eventTypeChips, countBy(state.rows, "eventType"), "eventTypes");
    renderChipList(els.evidenceChips, countBy(state.rows, "evidence"), "evidence");
    renderChipList(els.confidenceChips, countBy(state.rows, "conf"), "confidence");
    renderChipList(els.somaticChips, countBy(state.rows, "somatic"), "somatic");
  }

  function renderChipList(container, counts, filterKey) {
    const entries = Array.from(counts.entries()).sort((left, right) => {
      if (right[1] !== left[1]) {
        return right[1] - left[1];
      }
      return left[0].localeCompare(right[0]);
    });

    if (!entries.length) {
      container.innerHTML = '<span class="mini-note">Load a file to populate this filter.</span>';
      return;
    }

    container.innerHTML = "";
    entries.forEach(([value, count]) => {
      const button = document.createElement("button");
      button.type = "button";
      button.className = "chip";
      button.textContent = value + " (" + count.toLocaleString() + ")";
      if (state.filters[filterKey].has(value)) {
        button.classList.add("active");
      }
      button.addEventListener("click", () => {
        toggleSetValue(state.filters[filterKey], value);
        state.page = 1;
        render();
      });
      container.appendChild(button);
    });
  }

  function render() {
    state.filteredRows = applyFilters(state.rows);
    sortRows(state.filteredRows);
    reconcileSelection();
    updateSummary();
    renderTable();
    renderDetail();
    rebuildChipGroups();
  }

  function applyFilters(rows) {
    return rows.filter((row) => {
      if (state.filters.search && !row.searchText.includes(state.filters.search)) {
        return false;
      }
      if (!passesNumericRange(row.somlod, state.filters.somlodMin, state.filters.somlodMax)) {
        return false;
      }
      if (!passesNumericRange(row.maxlod, state.filters.maxlodMin, state.filters.maxlodMax)) {
        return false;
      }
      if (!passesNumericMinimum(row.qual, state.filters.qualMin)) {
        return false;
      }
      if (!passesNumericMinimum(row.split, state.filters.splitMin)) {
        return false;
      }
      if (!passesNumericMinimum(row.alt_count, state.filters.altMin)) {
        return false;
      }
      if (!passesNumericMinimum(Math.max(zeroIfNull(row.dcn), zeroIfNull(row.dct)), state.filters.discordantMin)) {
        return false;
      }
      if (state.filters.dbsnp === "present" && !row.hasDbsnp) {
        return false;
      }
      if (state.filters.dbsnp === "absent" && row.hasDbsnp) {
        return false;
      }
      if (state.filters.eventTypes.size && !state.filters.eventTypes.has(row.eventType)) {
        return false;
      }
      if (state.filters.evidence.size && !state.filters.evidence.has(row.evidence)) {
        return false;
      }
      if (state.filters.confidence.size && !state.filters.confidence.has(row.conf)) {
        return false;
      }
      if (state.filters.somatic.size && !state.filters.somatic.has(row.somatic)) {
        return false;
      }
      return true;
    });
  }

  function sortRows(rows) {
    const direction = state.sort.direction === "asc" ? 1 : -1;
    const key = state.sort.key;
    rows.sort((left, right) => {
      const leftValue = sortableValue(left, key);
      const rightValue = sortableValue(right, key);

      if (typeof leftValue === "number" || typeof rightValue === "number") {
        const leftNumber = normalizeSortNumber(leftValue);
        const rightNumber = normalizeSortNumber(rightValue);
        if (leftNumber !== rightNumber) {
          return (leftNumber - rightNumber) * direction;
        }
      } else {
        const compare = String(leftValue).localeCompare(String(rightValue));
        if (compare !== 0) {
          return compare * direction;
        }
      }
      return left.id - right.id;
    });
  }

  function sortableValue(row, key) {
    if (key === "locus") {
      return row.locus;
    }
    if (key === "spanAbs") {
      return row.spanAbs;
    }
    return row[key];
  }

  function updateSummary() {
    const totalRows = state.rows.length;
    const shownRows = state.filteredRows.length;
    const passCount = state.filteredRows.filter((row) => row.conf === "PASS").length;
    const topEvidence = topCountValue(state.filteredRows, "evidence");
    const selected = findSelectedRow();

    els.statTotal.textContent = totalRows.toLocaleString();
    els.statShown.textContent = shownRows.toLocaleString();
    els.statPass.textContent = passCount.toLocaleString();
    els.statTopEvidence.textContent = topEvidence || "-";
    els.resultsTitle.textContent = totalRows
      ? state.fileName + " (" + totalRows.toLocaleString() + " calls)"
      : "No variants loaded";
    els.resultsSummary.textContent = totalRows
      ? shownRows.toLocaleString() + " rows after filters. Sorted by " + state.sort.key + " (" + state.sort.direction + ")."
      : "Load a file to start browsing calls.";
    els.exportButton.disabled = shownRows === 0;
    if (selected && !state.filteredRows.length) {
      state.selectedId = null;
    }
  }

  function renderTable() {
    const totalPages = Math.max(1, Math.ceil(state.filteredRows.length / state.pageSize));
    if (state.page > totalPages) {
      state.page = totalPages;
    }

    const start = (state.page - 1) * state.pageSize;
    const end = start + state.pageSize;
    const rows = state.filteredRows.slice(start, end);

    if (!rows.length) {
      els.tableBody.innerHTML =
        '<tr><td colspan="10" class="empty-table">No rows match the current filter set.</td></tr>';
    } else {
      const html = rows
        .map((row) => {
          const isSelected = row.id === state.selectedId ? " selected" : "";
          const spanLabel = row.span === -1 ? "interchrom" : formatNumber(row.span);
          const support = [
            "split " + formatNumber(row.split),
            "alt " + formatNumber(row.alt_count),
            "disc " + formatNumber(Math.max(zeroIfNull(row.dcn), zeroIfNull(row.dct)))
          ].join(" / ");
          return (
            '<tr class="' +
            isSelected.trim() +
            '" data-row-id="' +
            row.id +
            '">' +
            '<td><div class="locus-cell"><strong>' +
            escapeHtml(row.locus) +
            '</strong><span>' +
            escapeHtml(row.ref + " -> " + row.alt_allele) +
            "</span></div></td>" +
            "<td>" +
            escapeHtml(row.eventType) +
            "</td>" +
            "<td><span class=\"badge\">" +
            escapeHtml(row.evidence) +
            "</span></td>" +
            "<td>" +
            escapeHtml(row.conf) +
            "</td>" +
            "<td>" +
            formatNumber(row.qual) +
            "</td>" +
            "<td>" +
            formatNumber(row.somlod) +
            "</td>" +
            "<td>" +
            formatNumber(row.maxlod) +
            "</td>" +
            "<td>" +
            escapeHtml(support) +
            "</td>" +
            "<td>" +
            escapeHtml(spanLabel) +
            "</td>" +
            "<td>" +
            escapeHtml(row.contig_and_region) +
            "</td>" +
            "</tr>"
          );
        })
        .join("");
      els.tableBody.innerHTML = html;

      els.tableBody.querySelectorAll("tr[data-row-id]").forEach((rowElement) => {
        rowElement.addEventListener("click", () => {
          state.selectedId = Number(rowElement.dataset.rowId);
          renderTable();
          renderDetail();
        });
      });
    }

    els.prevPage.disabled = state.page <= 1 || !state.filteredRows.length;
    els.nextPage.disabled = state.page >= totalPages || !state.filteredRows.length;
    els.pageSummary.textContent = "Page " + state.page + " of " + totalPages;
  }

  function renderDetail() {
    const row = findSelectedRow();
    if (!row) {
      els.detailContent.className = "detail-content detail-empty";
      els.detailContent.textContent =
        "Select a row to inspect the raw breakpoint call, support metrics, and sample-level values.";
      return;
    }

    const confidenceClass = row.conf === "PASS" ? "pass" : row.conf === "LOWLOD" || row.conf === "LOWMAPQ" ? "warn" : "low";
    const sampleCards = row.sampleInfo
      .map((sample) => {
        return (
          '<article class="sample-card">' +
          "<strong>" +
          escapeHtml(sample.label) +
          "</strong>" +
          '<div class="mini-note">' +
          escapeHtml(sample.header) +
          "</div>" +
          "<dl>" +
          detailPair("GT", sample.genotype) +
          detailPair("AD", formatNumber(sample.altSupport)) +
          detailPair("DP", formatNumber(sample.depth)) +
          detailPair("SR", formatNumber(sample.splitSupport)) +
          detailPair(sample.supportLabel, formatNumber(sample.pairedSupport)) +
          detailPair("GQ", formatNumber(sample.gq)) +
          detailPair("LO", formatNumber(sample.lod)) +
          detailPair("LO_n", formatNumber(sample.lodNormal)) +
          "</dl>" +
          "</article>"
        );
      })
      .join("");

    els.detailContent.className = "detail-content";
    els.detailContent.innerHTML =
      '<div class="detail-header">' +
      "<h3>" +
      escapeHtml(row.locus) +
      "</h3>" +
      '<p class="detail-subtitle">' +
      escapeHtml(row.eventType + " / " + row.evidence + " / " + row.conf) +
      "</p>" +
      '<div class="detail-tags">' +
      '<span class="tag ' +
      confidenceClass +
      '">' +
      escapeHtml("conf " + row.conf) +
      "</span>" +
      '<span class="tag">' +
      escapeHtml("somatic " + row.somatic) +
      "</span>" +
      '<span class="tag">' +
      escapeHtml("somlod " + formatNumber(row.somlod)) +
      "</span>" +
      '<span class="tag">' +
      escapeHtml("maxlod " + formatNumber(row.maxlod)) +
      "</span>" +
      "</div>" +
      "</div>" +
      '<div class="detail-grid">' +
      detailCard("Reference", row.ref) +
      detailCard("Alternate", row.alt_allele) +
      detailCard("Breakends", row.chr1 + ":" + row.pos1 + " " + row.strand1 + " -> " + row.chr2 + ":" + row.pos2 + " " + row.strand2) +
      detailCard("Support", "split " + formatNumber(row.split) + ", alt " + formatNumber(row.alt_count) + ", cigar " + formatNumber(row.cigar) + ", dcn " + formatNumber(row.dcn) + ", dct " + formatNumber(row.dct)) +
      detailCard("Mapping and alignment", "mapq " + formatNumber(row.mapq1) + "/" + formatNumber(row.mapq2) + ", AS " + formatNumber(row.as1) + "/" + formatNumber(row.as2) + ", NM " + formatNumber(row.nm1) + "/" + formatNumber(row.nm2)) +
      detailCard("Context", "homol " + presentOrX(row.homol) + ", insert " + presentOrX(row.insert) + ", repeat " + presentOrX(row.repeat)) +
      detailCard("Contig", row.contig_and_region) +
      detailCard("Flags", "qual " + formatNumber(row.qual) + ", secondary " + formatNumber(row.secondary) + ", span " + (row.span === -1 ? "interchrom" : formatNumber(row.span)) + ", dbsnp " + presentOrX(row.dbsnp)) +
      "</div>" +
      '<h4 class="section-title">Per-sample fields</h4>' +
      '<div class="sample-grid">' +
      sampleCards +
      "</div>" +
      '<h4 class="section-title">Raw row</h4>' +
      '<pre class="raw-block">' +
      escapeHtml(state.originalHeader.join("\t")) +
      "\n" +
      escapeHtml(row.originalCells.join("\t")) +
      "</pre>";
  }

  function detailCard(label, value) {
    return "<div><dt>" + escapeHtml(label) + "</dt><dd>" + escapeHtml(value) + "</dd></div>";
  }

  function detailPair(label, value) {
    return "<dt>" + escapeHtml(label) + "</dt><dd>" + escapeHtml(value) + "</dd>";
  }

  function reconcileSelection() {
    const selectedStillVisible = state.filteredRows.some((row) => row.id === state.selectedId);
    if (!selectedStillVisible) {
      state.selectedId = state.filteredRows.length ? state.filteredRows[0].id : null;
    }
  }

  function findSelectedRow() {
    return state.filteredRows.find((row) => row.id === state.selectedId) || state.rows.find((row) => row.id === state.selectedId) || null;
  }

  function exportFilteredRows() {
    if (!state.filteredRows.length) {
      return;
    }
    const lines = [state.originalHeader.join("\t")].concat(state.filteredRows.map((row) => row.originalCells.join("\t")));
    const blob = new Blob([lines.join("\n")], { type: "text/tab-separated-values;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement("a");
    const baseName = state.fileName.replace(/\.gz$/i, "").replace(/\.txt$/i, "");
    anchor.href = url;
    anchor.download = baseName + ".filtered.tsv";
    anchor.click();
    URL.revokeObjectURL(url);
  }

  function countBy(rows, key) {
    const map = new Map();
    rows.forEach((row) => {
      const value = row[key];
      map.set(value, (map.get(value) || 0) + 1);
    });
    return map;
  }

  function topCountValue(rows, key) {
    const counts = countBy(rows, key);
    let bestValue = "";
    let bestCount = -1;
    counts.forEach((count, value) => {
      if (count > bestCount) {
        bestCount = count;
        bestValue = value;
      }
    });
    return bestValue;
  }

  function toggleSetValue(set, value) {
    if (set.has(value)) {
      set.delete(value);
    } else {
      set.add(value);
    }
  }

  function setStatus(message) {
    els.statusText.textContent = message;
  }

  function passesNumericRange(value, min, max) {
    if (min === null && max === null) {
      return true;
    }
    if (value === null || !Number.isFinite(value)) {
      return false;
    }
    if (min !== null && value < min) {
      return false;
    }
    if (max !== null && value > max) {
      return false;
    }
    return true;
  }

  function passesNumericMinimum(value, min) {
    if (min === null) {
      return true;
    }
    if (value === null || !Number.isFinite(value)) {
      return false;
    }
    return value >= min;
  }

  function parseNullableNumber(raw) {
    if (raw === "" || raw === null || raw === undefined) {
      return null;
    }
    const value = Number(raw);
    return Number.isFinite(value) ? value : null;
  }

  function parseMaybeNumber(raw) {
    if (raw === "" || raw === "NA" || raw === "x" || raw === "NOTSET") {
      return null;
    }
    const value = Number(raw);
    return Number.isFinite(value) ? value : null;
  }

  function formatNumber(value) {
    if (value === null || value === undefined || value === Infinity) {
      return "-";
    }
    if (!Number.isFinite(value)) {
      return String(value);
    }
    return Math.abs(value) >= 100 || Number.isInteger(value)
      ? value.toLocaleString()
      : value.toFixed(3).replace(/\.?0+$/, "");
  }

  function presentOrX(value) {
    return value && value !== "x" ? value : "x";
  }

  function normalizeSortNumber(value) {
    if (value === Infinity) {
      return Number.MAX_SAFE_INTEGER;
    }
    if (value === null || value === undefined || Number.isNaN(value)) {
      return -Number.MAX_SAFE_INTEGER;
    }
    return Number(value);
  }

  function zeroIfNull(value) {
    return value === null || value === undefined ? 0 : value;
  }

  function stripLeadingHash(value) {
    return typeof value === "string" ? value.replace(/^#/, "") : value;
  }

  function escapeHtml(value) {
    return String(value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }
})();
