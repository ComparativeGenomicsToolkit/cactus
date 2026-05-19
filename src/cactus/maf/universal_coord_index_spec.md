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
- **Index A** — universal column ↔ ancestral coordinate. **Not a separately
  built structure**: it is a reconstructed in-memory *view* over the
  **reference-genome** (row-0/ancestral) run lists already inside Index B.

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
- Scan records **every row, including row-0**. The ancestral reference is just
  another genome; its runs are exactly the Index-A intervals (see Index A).
  Row-0 is gap-free ⇒ exactly one run per block for it (cheap).
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
   composes with reconstructed Index A: column -> (anc.seq,pos) -> .tai -> full column
```

Depends only on **taffy** (already a streaming, column-oriented MAF/TAF
processor tracking per-row coords/strands/names; owns `.tai`) + a vendored
ONElib. No Python over the MAF.

## Phasing (simplest-first)

1. ONEcode container per `(genome,seq)`.  ONElib's built-in `INT_LIST`
   Huffman was *measured* and is poor on absolute (t,g,len) triples
   (~13.7 B/run on the 577-way: 1.16 GiB).  So the run-line instead carries
   an explicit per-sequence codec: consecutive colinear runs merged, then
   zigzag-delta (t,g carried as running absolutes) + LEB128 varint, then
   zlib `deflate`; stored as `R = (inflatedLen INT, deflated bytes STRING)`.
   Measured ~4.2 B/run (≈3.3× smaller, ~355 MB) — zlib (already linked)
   chosen over zstd/xz (≈10–20 % more) to avoid a new dependency.
2. Further shrink (PForDelta/bit-packing, ~150–200 MB) only if measured-need.

## Open knobs

- Spill granularity: one file per genome (simple; ~hundreds of FDs for
  vertebrates — batch-flush if thousands) vs single partitioned spill.
- Phase 2 home: inside `taffy lift-index` vs separate finalizer (taffy doing
  both is cleanest).
- Strand: pack into the high bit of `length`.
- Which genomes: all by default (free RAM-wise — disk spill).

## Test plan

- Sampled `(genome,seq,p)`: `liftB(p)=g`; cross-check via reconstructed Index A
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

## Source: reconstructed from Index B's reference-genome runs

Index A is **not built separately**. The reference (row-0/ancestral) genomes
are recorded in Index B like any other genome; their runs *are* the Index-A
intervals (row-0 gap-free ⇒ one run per block; Phase-2 merge yields exactly the
maximal `.bed` intervals; column ids identical by construction).

Reconstruction at load (Index-A-sized = #intervals; cheap):

- `ancToCol`: directly the reference genomes' per-`(seq)` `t_start`-sorted run
  lists already in Index B (binary search by `pos`).
- `colToAnc`: concatenate the reference genomes' runs, sort by `g_start` once
  in memory (one-time, BED-sized), binary search by `col`.

Ancestral row-0 is always `+` strand, so no reverse-complement wrinkle.

Identifying the reference genomes: an explicit **reference-genome marker** in
the container (see schema). Robust fallback: the reference genomes are exactly
those whose runs tile `[0,T)` with no gaps/overlaps — but the marker is
cleaner and authoritative.

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

```
P 3 ucx                              universal column coordinate index

. ---- globals / directory (cheap front-of-file load) ----
D t 1 3 INT                          total columns T
D x 1 11 STRING_LIST                 reference (ancestral row-0) genome names
D d 4 6 STRING 6 STRING 3 INT 3 INT  dir entry: genome, seq, S-ordinal, seqLen

. ---- one uniform model for ALL genomes (ancestral + leaf) ----
G g 1 6 STRING                       genome name : groups the S objects below
O S 2 6 STRING 3 INT                 (genome) sequence: name, seqLen ; indexed object (oneGoto unit)
D R 1 8 INT_LIST                     run list: [t,g,len,strandBit]*nRuns, sorted by t (compressed)
```

Notes:
- **No separate Index-A object type.** Index A = the `S`/`R` objects whose
  genome is in `x` (the reference marker), reinterpreted (`ancToCol` direct;
  `colToAnc` via one-time `g`-sort of those runs).
- One compressed `INT_LIST` `R` line per `S` object — minimal records; ONEcode
  list compression is the deferred delta/varint phase, for free. Expanded
  per-run ASCII via `ONEview` for debugging.
- `strandBit` packed into `len`'s high bit.
- `S` is the only indexed object type → one footer index; `oneGoto(of,'S',k)`
  seeks the k-th sequence object; `g` group + `d` directory map
  `(genome,seq) → k`; `x` marks which genomes reconstruct Index A.

## Access model

- **Load:** read `t`, `x`, `d` (front of file) → `T`, reference set,
  name→ordinal. Reconstruct Index A in RAM from the reference genomes' `S`/`R`
  runs (incl. the one-time `g`-sort for `colToAnc`).
- **Point lift** `genome.seq:p`: dir → ordinal → `oneGoto('S',k)` → its `R`
  list → bsearch `t` → column → `colToAnc` (reconstructed Index A) → `.tai` →
  full column.
- **Bulk lift:** iterate a genome's `S` objects sequentially (ONEcode-native).

## Open knobs

- Index home: `taffy` (lean) vs standalone — gates whether ONElib is vendored
  in-tree.
- `.bed` retention: keep as optional debug / independent Index-A cross-check
  artifact (free from cactus-hal2maf), but it is **not** a runtime input.
- Single container (all genomes incl. ancestral) is the design; no split
  needed since Index A is just the reference subset.
