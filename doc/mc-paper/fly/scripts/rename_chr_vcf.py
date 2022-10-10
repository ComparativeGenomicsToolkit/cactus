import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input gzipped VCF', required=True)
parser.add_argument('-r', help='Prefix to remove', default='')
parser.add_argument('-a', help='Prefix to add', default='')
args = parser.parse_args()

if args.r == '' and args.a == '':
    print('either -r or -a arguments must be specified')
    exit()

pref_rm = False
def rename_chr(contig):
    return(args.a + contig)
if args.r != '':
    pref_rm = True
    def rename_chr(contig):
        return(contig.replace(args.r, ''))
    
for line in gzip.open(args.i, 'rb'):
    line = line.rstrip().decode('ascii')
    if line[0] == '#':
        if 'contig' in line:
            line = line.replace('##contig=<ID=', '')
            line = '##contig=<ID=' + rename_chr(line)
        print(line)
        continue
    line = line.split('\t')
    line[0] = rename_chr(line[0])
    print('\t'.join(line))
