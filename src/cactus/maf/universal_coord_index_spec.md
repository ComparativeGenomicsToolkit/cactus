# Universal Column Coordinate Index (`.tui`) — spec

Status: implemented in `taffy` (`taffy index -u` build, `taffy view -U`
query, `taffy lift` source-side range + reverse leaf lift), validated on
rodent (TAF), apes (MAF), and fish (MAF) universal alignments at vertebrate
scale.  This doc has been re-synced to the implementation after the
chunked-`R` + per-genome lift redesigns; it supersedes all earlier drafts.
Previous drafts that described "Index A" as a recorded reference track
or "reconstruct Index A from Index B + drive the `.tai` per row-0 piece"
are obsolete — Index A is not recorded and `.tai` is not consulted on the
`-U` path.

`cactus-hal2maf --universal` produces a MAF/TAF whose alignment columns, in
file order, form one shared coordinate system: the entire top-node ancestral
genome plus each descendant ancestor's novel columns, in canonical order
(ancestor pre-order → sequence → ascending start; maximal/merged). The
universal **column id** is the k-th alignment column in file order.

The `.tui` is a single provenanced ONEcode container, written next to the
MAF/TAF (`<file>.tui`), holding three things — and it is the **complete**
runtime input for `taffy view -U` and `taffy lift`: **no `.tai` is needed
for the universal query path.**

- **Index B** — any genome's coords → universal column(s).  Source-side
  range queries (`tui_query`, `tui_load_seq_runs`).
- **Index L** — universal column → a target genome's coords + paralog set.
  Reverse leaf lift (`tui_genome_lift_load` + `tui_genome_lift_column`).
- **Index X** — universal column → file offset (the random-access engine
  that replaces the `.tai` for `-U`).

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

- Per `genome.sequence`: a chunked run list `(t_start, g_start, length,
  strand)`.  Query `seq:p`: bsearch the chunk + run with `t_start ≤ p <
  t_start+length` → `g_start + (p − t_start)` (`+`) or
  `g_start + (t_start+length−1 − p)` (`−`). Gaps between runs = uncovered
  → unmapped. Partial function, 0-or-1 column/base.
- Range query `seq:[s,e)` → a sorted, **merged** set of universal-column
  intervals (`tui_query`).  Per-chunk t-range metadata lets the reader skip
  chunks that don't overlap `[s,e)` without paying the zlib decompress.
- Batch query "all runs of one seq" returned t-sorted for in-RAM bsearch
  by callers that do many lookups against the same source seq
  (`tui_load_seq_runs`).  Skips deliberately not applied here -- the
  contract is "return ALL runs."

### Build (two-phase)

**Phase 1 — one streaming MAF/TAF scan.** State = running column counter
`k` + current block.  For every row (incl. row-0) split at internal gaps
into within-block runs and append `(genome, seq, t_start, g_start, length,
strand)` to that genome's spill file (in column order).  Never Python
over the MAF (vertebrate scale ⇒ days).

Per-genome spill file descriptors are managed via a bounded LRU pool so
the writer never blows past `RLIMIT_NOFILE` at 1k+-genome scale.  See
`taffy/impl/tui.c` `Phase1::spill_ents` for the eviction logic.

**Phase 2 — per-genome finalize.** Sort each spill by `(seq, t_start)`,
colinear-merge runs across block boundaries, sort each chunk's runs by
`g_start`, write the genome's `S` + (`C`, `R`)+ chunk objects.  One genome
at a time; genomes independent ⇒ parallelizable.

Memory: bounded PER GENOME by that genome's spill size, not absolutely.
At vertebrate scale the largest leaf's `Run[]` can be ~1.2 GB; the .tui
build host budgets accordingly.

### Run serialization (`R` blob, per chunk)

Each sequence's runs are split into chunks of `TUI_CHUNK_RUNS` (= 65536)
runs at write time.  Within a chunk, runs are sorted by `g_start` (tight
per-chunk `[g_min, g_max)` for the column-keyed lookups in Index L);
between chunks, runs are still seq-major / t-major in the spill.  Each
chunk pairs one `C` header object with one `R` payload data record.

`R` is an explicit codec, not ONElib's built-in `INT_LIST` Huffman (which
was poor on absolute `(t,g,len)` triples — measured ~13.7 B/run, 1.16 GiB
pre-chunked on the rodent subtree).  A structure-of-arrays blob: three
concatenated LEB128-varint streams `gap | gsk | lenc`, `gap = t -
(prevT+prevLen)` (≈0 for ~99 % of splits → ≈ all-zero stream),
`gsk = g - (prevG+prevLen)` (the irreducible signal),
`lenc = len<<1|strand` — zlib-deflated.  `R = (inflatedLen INT, deflated
STRING)`.  Per-chunk metadata (g and t ranges) lives in the `C` header,
not the `R` payload, so the reader can skip whole chunks without
decompressing.

### Chunk metadata (`C` header)

```
O C 6 INT INT INT INT INT INT     # g_min, g_max, parent S-ord, self c_ord,
                                  # t_min, t_max
```

`g_min, g_max` are the chunk's universal-column range (used by Index L's
column-keyed lookup to early-out via `chunk_max_end[]`); `t_min, t_max`
are the source-coord range (used by `tui_query` to skip non-overlapping
chunks before decompress); `parent S-ord` ties the chunk back to its
sequence; `self c_ord` is the file-order `C` ordinal for `oneGoto(C,
c_ord)` during lazy `R` decode.

---

# Index L — universal column → target genome's coords (reverse leaf lift)

Per-target-genome reverse-direction lookup: given a universal column
`c`, return every base of the target genome aligned at that column
(typically 0-1, but multiple under paralog).

Drives `taffy lift -g <target>` (BED-in / wig-in → BED-out / wig-out in
target coords).  Replaces the prior "scan the whole .tui for every
column" approach with O(log n_chunks) lookups via the chunked C+R
layout.

### Load (`tui_genome_lift_load`)

For target `G`: walk the `d` directory entries with prefix `"G."`,
collect each matching seq's S-ordinal, jump to each S, read every C
header (NOT the R payloads) into an in-memory `TGLChunk[]` array
sorted by chunk `g_min`.  Build a running `max(g_max)` prefix array so
the column-keyed lookup can early-exit when no earlier chunk could
contain the queried column.  R payloads stay on disk; chunks are
decoded LAZILY on first column hit (each `TGLChunk` keeps its file
ordinal for `oneGoto(C, c_ord)` + `oneReadLine R`).

Memory at load: only chunk metadata (small).  Worst case after all
chunks have been lazily decoded: the full target's lift table in RAM
(~640 MB for a 3-Gb mammal).  The lift access pattern is a sequential
sweep, so LRU eviction would thrash; we accept the steady-state ceiling.

### Query (`tui_genome_lift_column`)

For column `c`: binary-search the chunks for the largest chunk with
`g_min ≤ c`, then walk backward through chunks whose
`chunk_max_end[j-1] > c` (overlapping candidates).  For each candidate
chunk, decode its R payload on first hit, then for each run in the
chunk whose `[g_start, g_start+length)` contains `c`, emit one
`(seq, pos, strand)` match.  Multiple matches per call = paralogs (the
same target genome aligns multiple times at this column).

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

## (Removed: Index A)

Earlier drafts of this spec described an "Index A" — a recorded column-
ordered row-0 ancestor track + a `tui_col_range_to_ref` query helper.
That track is NOT recorded in the on-disk format and the helper does
not exist.  The column-scan extractor (Index X + `tui_extract_*`) drives
all extraction; if a column→ancestor reporting feature comes back, it
would be added as a new line type.

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
D X 3 3 INT 3 INT 6 STRING           Index X: inflatedLen, nRec, deflate(SoA)
O d 3 6 STRING 3 INT 3 INT           dir: seqName, S-ordinal, seqLen
O S 2 6 STRING 3 INT                 sequence object: seqName, seqLen
O C 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT  chunk header:
                                     g_min, g_max,
                                     parent S-ord, self c_ord,
                                     t_min, t_max
D R 2 3 INT 6 STRING                 runs (one per chunk):
                                     inflatedLen, deflate(SoA delta blob)
```

Write order: `t`, `X` (front-of-file, cheap whole-load), then the `d`
directory in NAME-SORTED order (so the reader can binary-search it by name
via `oneGoto(of,'d',mid)`), then per-sequence `S` followed by (`C`, `R`)+
chunk pairs in genome-major order.  Genome is derived from the seq name
at build (`genome_of`), used only to group per-genome spills; the reader
does no genome resolution.

`d`, `S`, and `C` are all indexed object types (each gets a footer
index).  `oneGoto(of,'d',i)` jumps to the i-th name-sorted directory
entry; `oneGoto(of,'S',k)` jumps to the k-th genome-major sequence;
`oneGoto(of,'C',c_ord)` jumps to a specific chunk header by file-order
ordinal (used by Index L for lazy `R` decode).

### Load / access

- `tui_load`: read `t` and `X` only (front of file).  Bounded RAM (X ≈
  T/10000 anchors ≈ tens of MB at vertebrate scale).  Directory size
  `n_d` is read from the ONElib footer.  The directory and per-sequence
  chunks stay on disk.
- `-U` query path: binary-search the `d` directory by name (~O(log
  n_seqs) `oneGoto`s) → S-ordinal → `oneGoto(of,'S',ord+1)` → walk
  this seq's `C` headers, decompressing only the chunks whose
  `[t_min, t_max)` overlaps the query → universal-column intervals →
  Index X anchor + forward column scan.  No `.tai`.
- Reverse leaf lift (`taffy lift -g`): `tui_genome_lift_load(G)` walks
  the `d` directory for the `"G."` prefix, builds the per-target chunk
  index (metadata only); subsequent `tui_genome_lift_column(c)` calls
  binary-search by `g_min` + walk back through overlapping chunks,
  lazy-decoding `R` payloads on first hit per chunk.

### What this format is NOT

- No explicit `Index A` reference track (column-ordered row-0 ancestor
  segments).  The column-scan extractor (Index X + `tui_extract_*`)
  handles all extraction; if a column→ancestor reporting feature comes
  back, it would be added as a new line type (probably reusing the same
  SoA delta-codec).
- No `tui_col_range_to_ref` helper.  Earlier drafts mentioned this; it
  was never implemented.
- No `.tai` consultation on the universal query path.  The `-U` query is
  end-to-end driven by Index B + Index X.

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

- `T` is the canonical column count; Index-X anchors strictly increasing
  in `(column, file_pos)`; `idxCol[0] == 0` for MAF and any valid TAF.
- Per-chunk `(g_min, g_max)` is the tight range of `g` over the chunk's
  runs after the writer's g-sort; per-chunk `(t_min, t_max)` is the tight
  range of `t`.  These are written verbatim from the in-memory chunk; the
  build does not assert tightness post hoc.
- Format note: the `.tui` schema added per-chunk t-range fields (`O C`
  went from 4 INTs to 6) and the chunked `(C, R)+` layout in the same
  re-sync.  Older `.tui` files have a 4-field `C` and no t-range skip;
  the reader checks `info['C']->nField` at open and falls back to the
  decode-every-chunk path on the source side.  Newly-built `.tui` files
  transparently get the skip + reverse-lift speedups.
- Validation numbers (rodent 1.16 GiB → 253 MiB after delta+deflate)
  predate the chunked layout; the chunked format adds a tiny per-chunk
  overhead (16 B for the two new t-range INTs × ~6k-100M chunks at
  vertebrate scale = 0.02-1% size growth measured), with sublinear
  decompress cost per query because of the skip.
- Known low-priority follow-ups (do not affect lift correctness): clipped
  sub-blocks drop per-column `@` tags (N/A for cactus universal raw
  MAF — no tags); a disk-full idx-spill abort leaves spills (consistent
  with the deliberate fail-loud spill design); a pathologically sparse-
  anchor *TAF* (huge/zero `repeat_coordinates_every_n_columns`) makes a
  `-U` scan start far from the target (perf only, correctness holds;
  N/A for MAF input); `-n` + `-U` keys `tui_query` on the post-name-map
  name (pre-existing `-r` semantic).
- This spec file's home is temporary (`cactus.chain/src/cactus/maf/`) —
  must be relocated to a proper docs location before any merge.
