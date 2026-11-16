---
name: determinism
description: Cactus must stay deterministic - the same input and options give the same alignment or graph. Use when changing anything that could affect output content or order - C code in caf, bar, reference or the pinch/cactus graphs (hashing, sorting, iteration order, threads, random numbers), how work is split into jobs or chunks, or settings derived from cores, memory or the batch system.
---

# Determinism

Cactus is deterministic: the same input and the same options give the same output. Keep it that way.

- **Equivalent, not necessarily bit-identical.** Files may differ where they hold timestamps or other metadata (HDF5 does), but the alignment (HAL, MAF) or the graph (GFA, GBZ, VCF) must be the same.
- **Not affected by the run's circumstances:** the number of cores or threads, the batch system (single machine, Slurm, ...), system resources (memory, disk, the machine itself), job order, retries, or the working directory.
- **So output must never depend on** memory addresses or pointer order, thread scheduling, ties in an unstable sort, unseeded or time-seeded randomness, or how much work a resource limit happens to put in one job.
- **Two easy ones to miss in Python:** iterating a `set` of strings (string hashing is randomized in every process, so each Toil job can see a different order), and directory listings (`glob`, `os.listdir`: their order is the filesystem's and differs between machines). Iterate a list, or sort. Order matters even when the contents are the same: cactus_consolidated names sequences in the order they are given, and CAF and BAR break ties by name.
- **To check a change,** run the same alignment twice, and again with a different `--maxCores`, and compare the output by content, not file checksums: `hal2fasta` of each genome, `halStats`, `hal2maf`.
