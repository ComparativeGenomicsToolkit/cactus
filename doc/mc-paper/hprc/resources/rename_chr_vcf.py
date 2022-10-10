import argparse
import gzip

# remove a prefix from the chromosome/seqnames.
# for example the prefix added in the pangenome "GRCh38.chr1" -> "chr1"
parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input gzipped VCF', required=True)
parser.add_argument('-p', help='Prefix to remove', default='GRCh38.')
args = parser.parse_args()

for line in gzip.open(args.i, 'rb'):
    line = line.rstrip().decode('ascii')
    if line[0] == '#':
        if 'contig' in line:
            line = line.replace(args.p, '')
        print(line)
        continue
    line = line.split('\t')
    line[0] = line[0].replace(args.p, '')
    print('\t'.join(line))
