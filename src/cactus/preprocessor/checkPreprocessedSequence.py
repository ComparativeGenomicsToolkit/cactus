#!/usr/bin/env python3
"""Verify that a preprocessing step did not alter the underlying sequence.

Preprocessors are allowed to rename sequences, change the case of bases (this is
how soft-masking is represented), re-wrap the fasta at a different line length,
and reorder sequences (cactus_redPrefilter moves the filtered contigs to the end
of the file).  Anything else -- a base that changed, a sequence that was
truncated, dropped or added -- is data corruption, and is fatal.

The one sanctioned exception is hard-masking, which turns bases into N.  Pass
allow_hardmask=True for the preprocessors configured to do that.  That check
needs the sequences themselves rather than a checksum, so it is much slower; it
only runs when the fast path has already found a difference.

Clipping preprocessors (dna-brnn / maskFile with action="clip") remove sequence
on purpose, so the caller skips them entirely.

Each sequence is summarised by its length and a checksum of its upper-cased
bases.  Truncation -- the failure this exists to catch -- is caught exactly by
the length; the checksum only has to catch substitutions that preserve length,
so crc32 is used in preference to a cryptographic hash, which is twice as slow
for no practical gain here.  Swap _checksum for hashlib.md5 if that tradeoff
ever needs revisiting.
"""

import re
import zlib
from collections import Counter

# maps a-z to A-Z and leaves everything else alone.  Combined with a deletion
# set, bytes.translate does the upper-casing and the newline stripping in a
# single pass over the data.
_UPPER = bytes(range(256)).upper()
_STRIP = b'\r\n'

# small enough that each block stays in cache across the translate and the
# checksum; measured ~650 MB/s here, against ~170 MB/s for megabyte blocks.
_BLOCK = 1 << 15

_checksum = zlib.crc32


def preprocessed_fasta_id(seq_id):
    """Pull the fasta out of whatever a preprocessing step returned.

    dna-brnn and maskFile return (fasta, bed, merged bed); every other
    preprocessor returns the fasta on its own.  Those two run last in the chain,
    so before this check existed only the caller at the very end of the pipeline
    had to know about it.
    """
    return seq_id[0] if isinstance(seq_id, (tuple, list)) else seq_id


def fasta_digest(path):
    """Return [(name, length, checksum)], one entry per record, in file order.

    Length and checksum are over the upper-cased sequence with all whitespace
    removed, so both are invariant to soft-masking and to line wrapping.
    """
    records = []
    name = None
    crc = 0
    length = 0
    leftover = b''
    with open(path, 'rb') as f:
        while True:
            block = f.read(_BLOCK)
            if not block:
                break
            data = leftover + block if leftover else block
            cut = data.rfind(b'\n')
            if cut < 0:
                leftover = data
                continue
            leftover = data[cut + 1:]
            if cut + 1 != len(data):
                data = data[:cut + 1]
            # headers are rare, so the whole block is usually pure sequence
            if b'>' not in data:
                if name is None:
                    raise RuntimeError('{}: sequence data before the first header'.format(path))
                bases = data.translate(_UPPER, _STRIP)
                crc = _checksum(bases, crc)
                length += len(bases)
                continue
            for line in data.splitlines():
                if line.startswith(b'>'):
                    if name is not None:
                        records.append((name, length, crc))
                    fields = line[1:].split()
                    name = fields[0].decode('utf-8', 'replace') if fields else ''
                    crc = 0
                    length = 0
                elif line:
                    if name is None:
                        raise RuntimeError('{}: sequence data before the first header'.format(path))
                    bases = line.translate(_UPPER, _STRIP)
                    crc = _checksum(bases, crc)
                    length += len(bases)
    if leftover:
        # final line of a file that does not end in a newline
        if leftover.startswith(b'>'):
            if name is not None:
                records.append((name, length, crc))
            fields = leftover[1:].split()
            name = fields[0].decode('utf-8', 'replace') if fields else ''
            crc = 0
            length = 0
        elif name is not None:
            bases = leftover.translate(_UPPER, _STRIP)
            crc = _checksum(bases, crc)
            length += len(bases)
    if name is not None:
        records.append((name, length, crc))
    return records


def iter_fasta(path):
    """Yield (name, upper-cased sequence) one record at a time.

    Only used by the hard-masking path, which needs the bases themselves.
    """
    name = None
    chunks = []
    with open(path, 'rb') as f:
        for line in f:
            if line.startswith(b'>'):
                if name is not None:
                    yield name, b''.join(chunks)
                fields = line[1:].split()
                name = fields[0].decode('utf-8', 'replace') if fields else ''
                chunks = []
            elif name is not None:
                chunks.append(line.translate(_UPPER, _STRIP))
    if name is not None:
        yield name, b''.join(chunks)


def hardmask_ok(in_seq, out_seq):
    """True if out_seq is in_seq with some bases replaced by N, and nothing else.

    in_seq must already be upper-cased.  Only the stretches between N runs are
    compared, so no second copy of the sequence is built.
    """
    if len(in_seq) != len(out_seq):
        return False
    pos = 0
    for match in re.finditer(b'N+', out_seq):
        start, end = match.span()
        if in_seq[pos:start] != out_seq[pos:start]:
            return False
        pos = end
    return in_seq[pos:] == out_seq[pos:]


def _describe(records, limit=5):
    total = sum(length for _, length, _ in records)
    head = ', '.join('{} ({} bp)'.format(n, l) for n, l, _ in records[:limit])
    if len(records) > limit:
        head += ', ...'
    return '{} sequences totalling {} bp [{}]'.format(len(records), total, head)


def check_sequence_preserved(in_path, out_path, event_name=None, step_name=None,
                             allow_hardmask=False):
    """Raise RuntimeError if out_path is not in_path modulo names, case and wrapping.

    With allow_hardmask, bases are additionally permitted to become N.
    """
    context = 'preprocessor "{}"'.format(step_name) if step_name else 'preprocessor'
    if event_name:
        context += ' on genome "{}"'.format(event_name)

    in_records = fasta_digest(in_path)
    out_records = fasta_digest(out_path)

    # fast path: the same multiset of sequences, whatever the names or the order
    if (Counter((l, c) for _, l, c in in_records) ==
            Counter((l, c) for _, l, c in out_records)):
        return

    in_total = sum(length for _, length, _ in in_records)
    out_total = sum(length for _, length, _ in out_records)

    if len(in_records) != len(out_records):
        raise RuntimeError(
            '{} changed the number of sequences: input had {}, output has {}. '
            'Input: {}. Output: {}.'.format(
                context, len(in_records), len(out_records),
                _describe(in_records), _describe(out_records)))

    # pair up: by name when the names survived, otherwise by position
    in_names = [n for n, _, _ in in_records]
    out_names = [n for n, _, _ in out_records]
    if sorted(in_names) == sorted(out_names) and len(set(in_names)) == len(in_names):
        by_name = {n: (l, c) for n, l, c in out_records}
        pairs = [(n, l, c) + by_name[n] for n, l, c in in_records]
        paired_by = 'name'
    else:
        pairs = [(i_n, i_l, i_c, o_l, o_c)
                 for (i_n, i_l, i_c), (_, o_l, o_c) in zip(in_records, out_records)]
        paired_by = 'position'

    truncated = [(n, i_l, o_l) for n, i_l, _, o_l, _ in pairs if i_l != o_l]
    if truncated:
        detail = '; '.join('{}: {} bp -> {} bp ({:+d})'.format(n, i, o, o - i)
                           for n, i, o in truncated[:10])
        if len(truncated) > 10:
            detail += '; ... {} more'.format(len(truncated) - 10)
        raise RuntimeError(
            '{} changed the length of {} sequence(s), total {} bp -> {} bp ({:+d}). '
            'Sequences paired by {}. {}'.format(
                context, len(truncated), in_total, out_total, out_total - in_total,
                paired_by, detail))

    changed = [n for n, _, i_c, _, o_c in pairs if i_c != o_c]
    if not changed:
        return

    if allow_hardmask:
        # lengths already match, so the only legal difference left is base -> N.
        # this needs the bases themselves, and relies on the order being
        # unchanged, which holds for every preprocessor that can hard-mask.
        bad = []
        for (in_name, in_seq), (_, out_seq) in zip(iter_fasta(in_path),
                                                   iter_fasta(out_path)):
            if not hardmask_ok(in_seq, out_seq):
                bad.append(in_name)
                if len(bad) >= 10:
                    break
        if not bad:
            return
        raise RuntimeError(
            '{} altered {} sequence(s) in a way that is not hard-masking: {}'.format(
                context, len(bad), ', '.join(bad)))

    detail = ', '.join(changed[:10])
    if len(changed) > 10:
        detail += ', ... {} more'.format(len(changed) - 10)
    raise RuntimeError(
        '{} changed the bases of {} sequence(s); the lengths are unchanged, so this '
        'is not truncation. Sequences paired by {}. Affected: {}'.format(
            context, len(changed), paired_by, detail))


def main():
    import argparse
    import sys
    parser = argparse.ArgumentParser(
        description='Check that a preprocessed fasta has the same sequence as its input. '
                    'Renaming, case changes, re-wrapping and reordering are allowed.')
    parser.add_argument('input_fasta')
    parser.add_argument('output_fasta')
    parser.add_argument('--allow-hardmask', action='store_true',
                        help='also allow bases to become N')
    parser.add_argument('--event', default=None, help='genome name, for the error message')
    options = parser.parse_args()
    try:
        check_sequence_preserved(options.input_fasta, options.output_fasta,
                                 event_name=options.event,
                                 allow_hardmask=options.allow_hardmask)
    except RuntimeError as e:
        sys.stderr.write('SEQUENCE CHECK FAILED: {}\n'.format(e))
        sys.exit(1)
    sys.stderr.write('sequence preserved\n')


if __name__ == '__main__':
    main()
