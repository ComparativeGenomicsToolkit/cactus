# Universal Column Coordinate Index (`.tui`) — spec

Status: implemented in `taffy` (`taffy index -u` build, `taffy view -U`
query), validated on real rodent (TAF) + apes (MAF) universal alignments.
This doc reflects the implementation at the C1-redesign commit; it supersedes
all earlier drafts (the older "reconstruct Index A from Index B + drive the
`.tai` per row-0 piece" query model was found to silently drop blocks and is
abandoned — see Index X §).

`cactus-hal2maf --universal` produces a MAF/TAF whose alignment columns, in
file order, form one shared coordinate system: the entire top-node ancestral
genome plus each descendant ancestor's novel columns, in canonical order
(ancestor pre-order → sequence → ascending start; maximal/merged). The
universal **column id** is the k-th alignment column in file order.

The `.tui` is a single provenanced ONEcode container, written next to the
MAF/TAF (`<file>.tui`), holding three things — and it is the **complete**
runtime input for `taffy view -U`: **no `.tai` is needed for the universal
query path.**

- **Index B** — any genome's coords → universal column(s).
- **Index X** — universal column → file offset (the random-access engine
  that replaces the `.tai` for `-U`).
- **Index A** — universal column → row-0 (ancestral) coordinate. A recorded
  reference track (the materialized canonical `.bed`); used for *reporting
  which ancestor a column belongs to*, **not** for extraction.

---

## The universal column space

- Column id = k-th alignment column in MAF/TAF file order. `maxRefGap = 0`
  ⇒ every block has exactly `blocklen` columns, so the running count is a
  stable global integer. Row-0 (the per-block ancestral reference) is
  gap-free.
- `T` = total columns = Σ blocklen = Σ row-0 lengths; space `[0, T)`; fits
  in 64-bit (T ≈ 2.7e9 rodent subtree, 3.3e9 apes).
- **Exactly-once** (`--noRefDupes` + reference-rooted `--universal`): every
  covered base of every genome (ancestral or leaf) maps to 0 or 1 columns.
  0 = leaf/lineage-specific insertion, intentionally uncovered.
- **File order == column order** — the property the whole query engine rests
  on. Crucially the row-0 *ancestor* is NOT a single monotone coordinate:
  it changes per block and a given ancestor sequence recurs non-monotonically
  (root full-export interleaved with descendant `--novel` blocks). The
  *universal column*, by contrast, is a single globally-monotone key.

---

# Index B — any genome → universal column

Entry point from any genome's own coordinates into the shared column space
(`genome.seq:pos → column`), so e.g. `hg38` and `catshark` cross-lift via the
common space though they never co-occur in a reference-projected MAF.

- Per `genome.sequence`: a `t_start`-sorted run list `(t_start, g_start,
  length, strand)`. Query `seq:p`: bsearch the run with `t_start ≤ p <
  t_start+length` → `g_start + (p − t_start)` (`+`) or
  `g_start + (t_start+length−1 − p)` (`−`). Gaps between runs = uncovered →
  unmapped. Partial function, 0-or-1 column/base.
- Range query `seq:[s,e)` → a sorted, **merged** set of universal-column
  intervals (`tui_query`).

### Build (two-phase, memory-bounded for any #genomes / size)

**Phase 1 — one streaming C scan, O(1) RAM.** State = running column counter
`k` + current block. For every row (incl. row-0) split at internal gaps into
within-block runs and append `(genome, seq, t_start, g_start, length,
strand)` to that genome's spill file (in column order). Never Python over the
MAF (vertebrate scale ⇒ days).

**Phase 2 — per-genome finalize.** Sort each spill by `(seq, t_start)`,
colinear-merge runs across block boundaries, write the genome's `S`/`R`
objects. One genome at a time; genomes independent ⇒ parallelizable.

### Run serialization (`R`)

ONElib's built-in `INT_LIST` Huffman is poor on absolute `(t,g,len)` triples
(measured ~13.7 B/run, 1.16 GiB on the rodent subtree after the merge). So
`R` is an explicit per-sequence codec: merge colinear runs, then a
structure-of-arrays blob — three concatenated LEB128-varint streams
`gap | gsk | lenc`, `gap = t-(prevT+prevLen)` (≈0 for ~99 % of splits → ≈
all-zero stream), `gsk = g-(prevG+prevLen)` (the irreducible signal),
`lenc = len<<1|strand` — zlib-deflated. `R = (inflatedLen INT, nMeta…,
deflated STRING)`. Measured 253 MiB on the rodent subtree (4.4× vs
merged-absolute, 6.6× vs unmerged). zlib (already linked) chosen over
zstd/xz; PForDelta worse; lzma ≈15 % smaller but a new dependency.

---

# Index X — universal column → file offset (the query engine)

This replaces the `.tai` for the `-U` path and is the fix for the review's
Critical **C1**.

### Why the `.tai` cannot drive a universal query

`.tai`/`tai_iterator` is keyed on the row-0 **sequence coordinate** and
assumes each contig's blocks are file-position-monotone in that contig's
coords (true for canonical single-reference TAF). In a universal MAF the
row-0 ancestor changes per block and recurs **non-monotonically**, so
`tai_iterator`'s monotone-per-contig seek/stop **silently drops blocks**
(reproduced: `AncA→AncB→AncA` emitted 2 of 3 blocks, exit 0, no error). The
fix is to key random access on the **universal column** — a single globally
monotone axis (file order == column order) — so that failure is
*structurally impossible*.

### What is stored (`X`)

A sparse `(universal_column → file_offset)` index, sampled every
`TUI_IDX_BLOCK = 10000` universal columns at a **resync-able block start**.
Stored as `X = (inflatedLen INT, nRec INT, deflated STRING)`: two
concatenated monotone-varint streams (column-delta, file_pos-delta), zlib-
deflated (rodent +1.4 MB, apes +1.8 MB on top of the rest of the `.tui`).
Anchors are strictly increasing in both column and file offset (the build
asserts this; bgzf virtual offsets are monotone in read order).

`file_pos` = `LI_get_position(li)` = start of the *peeked* block-first line
(== the `tai_create_taf` `LI_tell`-after-consuming convention; a block-
oriented builder needs `li->pos`, not `li->prev_pos`, hence the new
`LI_get_position`). For TAF, a line is anchor-eligible only if it carries a
coordinate for **every** row — exactly `tai_create_taf`'s rule, reused via
the newly exposed `tai_taf_line_is_anchor` (= `parse_coordinates_line !=
NULL`). MAF blocks are self-contained ⇒ every block start is eligible, so
`idxCol[0] == 0` unconditionally and anchors are dense.

### Query: `taffy view -U -r ANYGENOME.seq:s-e`

1. Index B (`tui_query`) → sorted/merged universal-column intervals.
2. Binary-search `X` for the anchor ≤ the first interval start; `LI_seek`
   there; `LI_get_next_line`; for TAF `tai_resync_taf_line` (the `s→i`
   resync) — both **reused from `.tai`** (newly exposed
   `tai_resync_taf_line` / `tai_maf_read_block`, one implementation shared
   with `.tai`).
3. Single forward column-ordered scan: read blocks (maintaining the running
   universal column = anchor col + Σ column_number, and the TAF-delta
   context), **emit each block whose column span overlaps an interval,
   clipped to the covered universal-column run(s)** — one physical block can
   yield several sub-blocks where the queried genome has internal gaps;
   non-overlapping intermediate blocks are read (for context) and discarded.
   Stop past the last interval. Output is whole alignment blocks (all
   species), in column order, each exactly once — **no dedup, no `.tai`, no
   row-0 / ancestor-name heuristics, C1 impossible**.

Clipping reuses `clip_alignment`'s start/length convention (`start +=`
removed-left non-gap; `length =` kept non-gap; identical for `+`/`-`), so a
partial `-U` query is byte-identical to the equivalent plain `taffy view
-r row0seq:…` extraction.

For TAF output the emitted blocks are generally not file-adjacent (skipped
intermediates), so they are written self-contained (`taf_write_block2(NULL,
…)`) — never delta-encoded against a non-adjacent predecessor.

---

# Index A — universal column → row-0 ancestor (recorded, not the query path)

A reference track, in universal-column order, one record per maximal row-0
segment: `(colStart, sOrd, row0Start, len)`, strand always `+` (asserted).
It is the materialized canonical `.bed`: `colStart` strictly increasing,
segments tile `[0,T)` exactly (every block has exactly one row-0). Stored
deflated SoA (`A = (inflatedLen INT, nSeg INT, deflated STRING)`); `colStart`
is **not** stored (= prefix sum of `len`, since the segments tile). Built in
the Phase-1 scan (row-0 is `aln->row`), colinear-merged on the fly.

Recorded so a column can be reported as a row-0 `(ancestor.seq, pos)` (e.g.
"which ancestor does universal column c belong to", debug / future
`colToAnc`). It is **not** used by `-U` extraction (Index X drives that).
The old "use Index A to produce row-0 pieces and run `tai_iterator` per
piece" model is the C1 bug and is gone. `tui_col_range_to_ref` remains
available for the column→ancestor mapping but is off the extraction path.

---

# Serialization: ONEcode container

One provenanced ONEcode binary file (`ONElib`, BSD-3-style; vendored from its
own upstream). Binary form has a persisted per-object-type footer index
(`oneGoto(of,type,i)` = O(1) `fseek`); ASCII (`ONEview`) is debug-only.
Provenance/reference records tie it to the source MAF/TAF, `T`, cactus + hal
commits.

### Schema (as built — `taffy/impl/tui.c` `TUI_SCHEMA`)

```
P 3 tui                              universal column index

D t 1 3 INT                          total columns T (global)
D A 3 3 INT 3 INT 6 STRING           Index A: inflatedLen, nSeg, deflate(SoA)
D X 3 3 INT 3 INT 6 STRING           Index X: inflatedLen, nRec, deflate(SoA)
D d 4 6 STRING 3 INT 3 INT 3 INT     dir: seqName, S-ordinal, seqLen, isRef
O S 2 6 STRING 3 INT                 sequence object: seqName, seqLen
D R 2 3 INT 6 STRING                 runs: inflatedLen, deflate(SoA delta blob)
```

Write order: `t`, `A`, `X` (front-of-file, cheap whole-load), then the `d`
directory, then per-sequence `S`/`R`. Genome is derived from the seq name
(`genome_of`), so no genome group/name list. `S` is the only indexed object
type (one footer index; `oneGoto(of,'S',k)` — the reader currently re-scans,
a known baseline, not a correctness issue).

### Load / access

- Load `t`, `A`, `X`, `d` (front of file) → `T`, the Index-A segment array,
  the Index-X anchor arrays, `seqName ↔ ordinal`. Bounded RAM (X ≈ T/10000
  anchors ≈ tens of MB at vertebrate scale; A ≈ #merged row-0 segments).
- `-U` query: Index B (`R` of that seq, located via the directory) →
  universal intervals → Index X anchor + forward column scan (above). No
  `.tai`.
- Bulk lift: iterate a genome's `S` objects sequentially (ONEcode-native).

---

# Validation (performed, real data)

- **Rodent (TAF, `.taf.gz`)**: `-U` root-genome region == plain `taffy view
  -r` **byte-identical** at scale; mouse `GCF…:1000000-1000100` exact;
  **C1 multi-ancestor** `MuridaeAnc4refChr183:0-40000` → 428 blocks, row-0
  = 321 MuridaeAnc3 + 107 MuridaeAnc4 (the non-monotone pattern), coverage
  exactly `[0,40000)`, **0 blocks dropped** (old design dropped exactly the
  MuridaeAnc3-row-0 blocks here).
- **Apes (MAF, `.maf.gz`)**: root-identity byte-identical;
  `hg38.chr1:1000000-1000100` exact.
- **Synthetics (MAF + TAF, incl. RLE)**: `AncA→AncB→AncA` → all 3 blocks;
  reverse-strand leaf; internal-gap blocks → multi-run sub-blocks; intervals
  spanning / abutting block boundaries; TUI_IDX_BLOCK boundary; all-gap rows;
  empty regions; `-m/-p/-c/-a/-n`. Every result byte-identical to plain `-r`
  ground truth. Suite 21/21. valgrind: zero per-block leak / invalid access
  over a 70k-block scan.

Self-checks (taffy builds with asserts on): every `R`, `A`, `X` blob is
decode-round-tripped against its in-RAM source at build time; the X-anchor
column monotonicity is asserted; spill write/close failures (disk full)
abort loudly rather than emit a silently-truncated index.

---

# Invariants / open items

- `T` is the canonical column count; Index-A segments tile `[0,T)`; Index-X
  anchors strictly increasing in (column, file_pos); `idxCol[0] == 0` for
  MAF and any valid TAF.
- Format note: the `.tui` schema (`A`, `X`, deflate-`R`) changed during this
  work — any pre-existing `.tui` must be rebuilt (`taffy index -u`);
  independent of MAF generation.
- Known low-priority follow-ups (do not affect lift correctness): clipped
  sub-blocks drop per-column `@` tags (N/A for cactus universal raw
  MAF — no tags); a disk-full idx-spill abort leaves spills (consistent with
  the deliberate fail-loud spill design); a pathologically sparse-anchor
  *TAF* (huge/zero `repeat_coordinates_every_n_columns`) makes a `-U` scan
  start far from the target (perf only, correctness holds; N/A for MAF
  input); `-n` + `-U` keys `tui_query` on the post-name-map name (pre-existing
  `-r` semantic).
- This spec file's home is temporary (`cactus.chain/src/cactus/maf/`) — must
  be relocated to a proper docs location before any merge.
