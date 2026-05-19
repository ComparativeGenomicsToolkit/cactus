# Universal Column Coordinate Index — spec

Status: design.

`cactus-hal2maf --universal` produces a MAF whose row-0 (block-reference)
intervals form a single shared coordinate system: the entire top-node
ancestral genome plus each descendant ancestor's novel intervals, in canonical
order (ancestor pre-order → sequence → ascending start; maximal/merged;
identical to MAF block order and to the `.tai`).

Two logical indexes over that space:

- **Index B** — any genome's coordinates → universal column. Built by one
  streaming scan of the MAF over **all rows** (including row-0). Per-`(genome,
  seq)` run lists.
- **Index A** — universal column ↔ ancestral coordinate. An **explicit
  reference track** recorded in the container at build time (the canonical
  `.bed`, materialized): per maximal row-0 segment, in column order,
  `(colStart, sOrd, row0Start, len)`.  NOT reconstructed from Index B's
  per-genome runs — that only works under `--noAncestors` (see Index A §).

Consequence: the runtime query system needs only **`.tai` + the ONEcode Index
B container**. The `.bed` is *not* a runtime input — it reduces to an optional
human-readable / debug / independent cross-check artifact (still emitted
cheaply by cactus-hal2maf's awk piggyback, so keeping it costs nothing).

Tool-agnostic; home TBD — current lean: Index B (and thus the whole index) in
`taffy`. The `.bed`→Index-A standalone path can live anywhere.

---

## The universal column space (definitions)

- The universal **column id** = the k-th alignment column in MAF file order.
  Because blocks are in canonical order and `maxRefGap = 0` makes every block
  have exactly `blocklen` columns, this running count is a stable global
  integer. Row-0 (the ancestral reference) is gap-free, so a block's row-0
  spans ancestral `[start, start+len)` ↔ columns `[k, k+len)` contiguously.
- `T` = total columns = Σ over blocks of `blocklen` = Σ row-0 lengths. Column
  space is `[0, T)`. `T` fits in 64-bit.
- Exactly-once guarantee (from the `--noRefDupes` + reference-rooted
  `--universal` design, validated on fly): every covered base of **every**
  genome (ancestral or leaf) maps to exactly **0 or 1** columns. No multimap,
  no policy. 0 = leaf/lineage-specific insertion, intentionally uncovered
  (out of current scope to cover; future relaxation noted elsewhere).

---

# Index B — any genome → universal column

## Purpose

Entry point from any genome's own coordinates into the shared column space
(`genome.seq:pos → column`), so e.g. `hg38` and `catshark` can be cross-lifted
via the common space even though they never co-occur in a reference-projected
MAF. Also the *source* of Index A (its reference-genome subset).

## Definitions

- Build needs **only a sequential column scan** of the MAF — not the `.bed`,
  not the `.tai`, not Index A. The column id is the running file-order column
  count; that is identical, by construction, to the `.bed`'s `G` numbering.
- Scan records **every row, including row-0**, into Index B. Row-0 is *also*
  captured separately, in column order, as the explicit **reference track**
  (Index A) — see Index A §. Row-0 is gap-free ⇒ one segment per block,
  colinear-merged on the fly (cheap, BED-sized).
- Per `genome.sequence`: a `t_start`-sorted run list
  `runs[] = (t_start, g_start, length, strand)`. Query `seq:p` →
  binary-search the run with `t_start ≤ p < t_start+length` →
  `g_start + (p − t_start)` (`+`) or `g_start + (t_start+length−1 − p)` (`−`).
  Gaps between runs = uncovered bases → unmapped.
- Partial function, 0-or-1 column per base (exactly-once).

## Cost reality

Inherently **O(all alignment rows)** — something must touch every row once.
Design goal: a single **streaming C** scan, never Python over the MAF
(vertebrate scale ⇒ a Python pass = days).

## Serialized form

The run list *is* the colocality compression (a colinear stretch of any length
= one record; no per-base storage). Serialized as the ONEcode container below;
logical shape per `(genome,seq)`:

```
record:  t_start  g_start  length  strand        (run; ≈16 B fixed-width equiv)
per (genome,seq): records sorted by t_start
query:   locate (genome,seq) object -> its runs -> bsearch t_start   # O(log #runs)
```

## Build pipeline (two-phase; memory-bounded regardless of #genomes / size)

**Phase 1 — single streaming C MAF scan, O(1) memory.** State = running
column counter `k` + current block only. For **each row (incl. row-0)**, split
at internal gaps into within-block runs and *append*
`(genome, seqId, t_start, g_start, length, strand)` to that genome's **spill
file**. No per-genome accumulation in RAM — spilled as produced, in column
(`g`) order. One scan; bounded memory for any #genomes / alignment size.

**Phase 2 — per-genome finalize, bounded memory, parallelizable.** Per genome:
external-sort its spill by `(seqId, t_start)`; merge colinear runs contiguous
across block boundaries (heals maxBlockLen / chunk splits); write its objects
into the ONEcode container. Memory = one genome at a time; genomes independent
⇒ parallel.

The `t_start` sort is intrinsic (single scan is column-ordered; queries are
genome-coordinate-ordered). Spills are per-genome and ≪ MAF.

⇒ Memory is never the constraint: Phase 1 spills to disk (O(1) RAM); Phase 2
is one genome at a time.

## Query API / tool surface

```
taffy lift-index -i universal.maf [-o out.ucx]   # Phase 1 + Phase 2 -> ONEcode container
taffy lift  <genome.seq:pos | BED>               # -> universal column(s)
   composes with the explicit Index-A track: column -> (anc.seq,pos) -> .tai -> full column
```

Depends only on **taffy** (already a streaming, column-oriented MAF/TAF
processor tracking per-row coords/strands/names; owns `.tai`) + a vendored
ONElib. No Python over the MAF.

## Phasing (simplest-first)

1. ONEcode container per `(genome,seq)`.  ONElib's built-in `INT_LIST`
   Huffman was *measured* and is poor on absolute (t,g,len) triples
   (~13.7 B/run on the 577-way: 1.16 GiB after the colinear merge).  So the
   run-line instead carries an explicit per-sequence codec: merge colinear
   runs, then a structure-of-arrays blob — three concatenated LEB128-varint
   streams `gap | gsk | lenc` where `gap = t-(prevT+prevLen)` (≈0 for ~99% of
   splits: sequence forward-contiguous, only OTHER lineages' columns
   intervened → stream ≈ all zeros), `gsk = g-(prevG+prevLen)` (intervening
   universal columns; the irreducible signal), `lenc = len<<1|strand` —
   zlib-`deflate`d; stored as `R = (inflatedLen INT, deflated bytes STRING)`,
   small header carries `m` and the two stream lengths.  Measured 253 MiB on
   the 577-way (4.4× vs merged-absolute, 6.6× vs original unmerged 1.74 GiB).
   zlib (already linked) over zstd/xz; PForDelta on `gsk` measured worse
   (heavy-tailed); lzma ≈15% smaller but new dependency + slow builds.
2. Further shrink only if a real consumer needs it (lzma; byte-plane
   transpose of `gsk` ≈8%).

## Open knobs

- Spill granularity: one file per genome (simple; ~hundreds of FDs for
  vertebrates — batch-flush if thousands) vs single partitioned spill.
- Phase 2 home: inside `taffy lift-index` vs separate finalizer (taffy doing
  both is cleanest).
- Strand: pack into the high bit of `length`.
- Which genomes: all by default (free RAM-wise — disk spill).

## Test plan

- Sampled `(genome,seq,p)`: `liftB(p)=g`; cross-check via the Index-A track
  → `.tai` → the column actually contains `genome.seq:p`.
- Round-trip small case: every base's run-list column == its block's running
  column id.
- Uncovered (leaf-insertion) bases fall in run gaps → unmapped.
- Chunk-invariance: result independent of chunkSize / batch count.
- Memory bound: Phase 1 RSS flat vs alignment size / #genomes.

---

# Index A — column ↔ ancestral (a view over Index B)

## What it is

The bijection between a global integer **column** and an **ancestral
coordinate** `(genome.sequence, position)`. Semantics:

- `colToAnc(col)`, `0 ≤ col < T` → `(anc.seq, pos)`.
- `ancToCol(anc.seq, pos)` → `col`, or `None` if not covered.
- Range forms `colRangeToAnc(a,b)` / `ancRangeToCol(anc.seq,s,e)` —
  first-class (the "later work" is interval-based).

Exactly-once ⇒ a clean bijection on the covered set; no multimap.

## Source: an EXPLICIT reference track in the container (not reconstructed)

Earlier drafts said Index A is a reconstructed view over Index B's
reference-genome runs. **That is fragile and is abandoned.** It only works
under `cactus-hal2maf --noAncestors`:

- `.tai` (the extraction engine we reuse) is keyed strictly on the literal
  **row-0 sequence+pos** per block.
- Index B records *every* genome's column coverage, with only a per-`(genome,
  seq)` `isRef` boolean — it does **not** mark *which* coverage was row-0.
- Without `--noAncestors`, an internal ancestor also appears as a **non-row-0
  row** in other blocks. Its Index-B runs then cover columns where it is *not*
  the `.tai` key, and ≥2 ancestors' runs cover the same column. "Which
  ancestor is the row-0 `.tai` indexed here" becomes unrecoverable from Index
  B → silently wrong / empty extraction (no error). Verified at scale
  (apes 26.4M blocks: row-0 spans Anc0..Anc6; the clean property held only
  because those files used `--noAncestors`).

The universal↔row-0 map is **known a priori** — it is exactly the canonical
`.bed` `cactus-hal2maf --universal` emits, and exactly what the index builder
sees as it scans blocks in column order (`aln->row` *is* row-0). So **store
it explicitly**; never infer it. This dissolves the `--noAncestors`
dependency: the track records the *actual* row-0 per block regardless of
whether other ancestors are also present as non-row-0 rows.

### Reference track (= Index A, = the `.bed`, materialized)

In universal-column order, one record per maximal row-0 segment:

    (colStart, sOrd, row0Start, len)        # all INT; strand always '+'

- `colStart` — first universal column of the segment (strictly increasing;
  segments tile `[0,T)` exactly: every block has exactly one row-0).
- `sOrd` — directory S-ordinal of the row-0 sequence (resolve to the name via
  the `d` directory; avoids repeating long ancestor names).
- `row0Start` — start of the segment in that row-0 sequence's own coords.
- `len` — column count = row-0 base count (row-0 is gap-free, maxRefGap==0).
- Strand: row-0 is always `+` (and `tai_create` itself rejects a `-` row-0),
  so it is asserted at build, not stored. `col∈[colStart,colStart+len)` ⇒
  `row0pos = row0Start + (col − colStart)` (affine, no RC wrinkle).

Schema addition (one list line, written right after `t`, before `d`):

    D A 1 8 INT_LIST          # [colStart, sOrd, row0Start, len] * nSeg, column-ordered

Start simple: a single `A` INT_LIST. The track is small (merged row-0
segmentation ≪ the leaf runs; e.g. rodent root `refChr0` = one segment for
3642 columns). If a real file shows it large, upgrade `A` to the same
deflated delta+varint blob form as `R` (follow-up, only if measured).

### Builder change (Phase 1, no Phase-2 sort needed)

The Phase-1 scan is already in column order and already visits `aln->row`
(row-0) per block. Maintain one *open* segment; per block:

- `seg = (colStart = running T, sOrd = dir-ordinal of row0 seq,
  row0Start = aln->row->start, len = column_number)`;
  assert `colStart == running T` (columns globally sequential) and
  `aln->row->strand == '+'`.
- If it extends the open segment colinearly — same sOrd, same strand,
  `open.row0Start + open.len == seg.row0Start` (column contiguity is
  automatic) — grow `open.len += column_number`; else flush `open`, start a
  new one. Flush the last open segment at end of scan.

On-the-fly merge bounds RAM/size to the *merged* segment count. Emit the
(already column-ordered, already merged) track as the `A` line in Phase 2.

### Reader / lookup

- `tui_load` additionally reads `A` into a RAM array `RefSeg[]` (sorted by
  `colStart` by construction) and an `sOrd → name` table from `d`.
- `colRangeToRef(a, b)`: binary-search `RefSeg` by `colStart`; split `[a,b)`
  at segment boundaries; per piece emit
  `(name(sOrd), row0Start + (pieceStart − colStart), pieceLen)`.
  One universal interval can yield several pieces on **different** ancestors
  (e.g. `hg38.chr1:1-100 → AncR…  + AncC…`) — handled naturally.

### Step-2 in `taffy view` (replaces the `-U` print-only mode)

`-r ANYGENOME.seq:s-e` →
1. Index B (step 1, done): → universal intervals.
2. `colRangeToRef` → `(row0Seq, refStart, refLen)` pieces; merge/dedup
   per `row0Seq` (a leaf query can produce many tiny universal intervals that
   fall in the same row-0 block — collapse so `.tai` extracts each block once).
3. For each piece: existing `tai_iterator` on the existing `.tai`; emit blocks
   in column order. No ancestor-name heuristics anywhere.

### Validation

- The `A` track must equal the canonical `.bed` `cactus-hal2maf --universal`
  emitted (same column ids by construction — cross-check on the rodent file).
- Root-genome region query == plain `taffy view -r root…` (affine identity
  where row-0 == the queried genome).
- Rodent + apes: leaf region → step1→step2→extract; the queried genome's row
  in the output must cover exactly the queried interval.

### Format note

Adding the `A` line is a `.tui` format change ⇒ existing `.tui` files must be
rebuilt (`taffy index -u`) to gain the track. The builder edit is small; the
rebuild is the existing ~tens-of-minutes index pass. Independent of MAF
generation, so it does not block the cluster fish-MAF run.

## Optional independent path (debug / cross-check, NOT runtime)

The same Index A can be built standalone from the canonical `.bed` (one pass,
prefix-sum `G`, in-memory arrays). It is **identical by construction** to the
reconstructed-from-B view (both derive from the same gap-free row-0 column
stream). Use it only as an independent verification artifact or for tooling
that has the `.bed` but not the container. Not required at runtime.

## Composition with `.tai`

`colToAnc(col)` → `(anc.seq,pos)` → `taffy view -r anc.seq:pos-… -i maf` (via
`.tai`) → the full column (all species). Reconstructed Index A + `.tai` give
random access to any column's complete contents; Index B is the entry point
from a non-ancestral genome's coordinates.

## Sizes

#reference-genome runs ≈ #`.bed` intervals ≈ 1e6–1e7 → ≤ ~200 MB reconstructed
in RAM. One-time `g`-sort at load is fast.

## Invariants / guards

- Numbering is defined by *the* canonical column stream (≡ `.bed` order ≡
  `.tai` order). Container stores `T`; assert it matches the `.tai` column
  extent and (if present) the `.bed` length sum.
- `colToAnc` total on `[0,T)`; `ancToCol` `None` on uncovered input (no
  exceptions).
- Reconstructed-from-B Index A and `.bed`-built Index A must be identical — a
  cheap CI cross-check (this is exactly why the `.bed` is worth keeping as a
  debug artifact).

---

# Serialization: ONEcode container

One provenanced ONEcode file (`ONElib`, as used by FASTGA) holds **every
genome's run lists** (ancestral genomes included — Index A is the reference
subset). `.tai` + this container = the complete runtime input.

## Confirmations (verified against `ONElib.c` / `FASTGA/LICENSE`)

- **License:** `FASTGA/LICENSE` is permissive **BSD-3-Clause-style** (source +
  binary redistribution with retained notice + no-endorsement). Vendor-OK.
  ONElib is a *separate upstream project* (Durbin/Myers, 2019–) with its own
  copyright — for a clean vendor, take `ONElib.{h,c}` from its own upstream
  repo with that repo's LICENSE + attribution. No legal blocker.
- **Persisted object index:** binary ONE files write a **footer** with a
  per-object-type byte-offset index (schema `D - … offset of footer`,
  `D & … li->index`); loaded on open (`memcpy li->index`); `oneGoto(of,type,i)`
  is a single `fseek(index[i])` — O(1) seek to the i-th object, **no scan**.
- **Caveats baked in:**
  - Footer/index exist only in the **binary** form (ASCII has none →
    `oneGoto` fails on ASCII). Write/query binary; ASCII (`ONEview`) is
    debug-only. Duality holds; the seekable index is a binary property.
  - `oneGoto` is by **object ordinal**, not key. Objects are `(genome,seq)`;
    a point query is directory → ordinal → `oneGoto` (O(1) fseek) → read that
    seq object's runs → in-RAM bsearch by `t_start`. Bulk lift = sequential
    object read (ONEcode-native, ideal; matches bulk-liftover).
- **Integrity tie:** source MAF + `.tai` paths/sizes, `T`, cactus + hal
  commits → ONEcode **provenance/reference** header records
  (`oneAddProvenance`/`oneInheritReference`). Self-describing; no rolled
  checksum.

## Schema (`P ucx`) — unified, one model for all genomes

ONEcode type-name encoding is length-prefixed (`3 INT`, `6 STRING`,
`8 INT_LIST`, `4 CHAR`). Provenance/reference are standard header records via
the API.

As built (`taffy/impl/tui.c` `TUI_SCHEMA`), with the new Index-A line `A`:

```
P 3 tui                              universal column index

D t 1 3 INT                          total columns T (global)
D A 1 8 INT_LIST                     Index-A reference track:
                                       [colStart, sOrd, row0Start, len]*nSeg
                                       column-ordered, segments tile [0,T)
D d 4 6 STRING 3 INT 3 INT 3 INT     dir: seqName, S-ordinal, seqLen, isRef
O S 2 6 STRING 3 INT                 sequence object: seqName, seqLen
                                       (the oneGoto unit)
D R 2 3 INT 6 STRING                 runs: inflatedLen, deflate(SoA delta blob)
```

(Earlier drafts had `P 3 ucx`, a `D x` reference-name list, a `D g` group,
and `D R … INT_LIST`. Superseded: genome is derived from the seq name
(`genome_of`); ONElib's `G` group carries object-count semantics we don't
need; `R` is a per-sequence zigzag-delta+varint+zlib blob, not `INT_LIST`
(measured ~5× smaller). `x` is gone — `isRef` lives per `d` entry, and Index A
is now the **explicit `A` track**, not a reconstruction.)

Notes:
- **`A` is the materialized `.bed`.** `colStart` strictly increasing;
  segments partition `[0,T)` (every block has exactly one row-0). `sOrd`
  resolves to a name via `d`. Strand always `+` (asserted at build).
  Written once, right after `t`, before `d` (front-of-file, cheap load).
- One deflated SoA blob `R` per `S` object (see the encoding section).
- `S` is the only indexed object type → one footer index;
  `oneGoto(of,'S',k)` seeks the k-th sequence object; `d` maps
  `seqName → k`. (Reader currently re-scans instead of `oneGoto`; baseline.)

## Access model

- **Load:** read `t`, `A`, `d` (front of file) → `T`, the Index-A segment
  array (already column-ordered), `seqName ↔ ordinal`. No reconstruction.
- **Point/range lift** `genome.seq:s-e`: Index B (`R` of that seq) →
  universal intervals → binary-search `A` (`colRangeToRef`) →
  `(row0Seq, refStart, refLen)` pieces (possibly several ancestors) →
  existing `.tai` per piece → full columns.
- **Bulk lift:** iterate a genome's `S` objects sequentially (ONEcode-native).

## Open knobs

- Index home: `taffy` (lean) vs standalone — gates whether ONElib is vendored
  in-tree.
- `.bed` retention: keep as optional debug / independent Index-A cross-check
  artifact (free from cactus-hal2maf), but it is **not** a runtime input.
- Single container (all genomes incl. ancestral) is the design; no split
  needed since Index A is just the reference subset.
