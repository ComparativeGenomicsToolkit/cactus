import fileinput


# parse a TSV-like stream created by jq from the json-converted GAM records
# and count how many reads
# in each profile (mapping quality x perfect alignement x aligment score)
#    mapping quality: number. '-1' for unmapped reads
#    perfect: boolean to specify if the reads aligned perfectly
#    alignment score: value of the AS tag
# the input tsv columns are: name, identity, score, query_position, sequence

# to tally the number of reads/records in each profile
records = {}

# to deal with secondary (skip them) and split reads (flag them)
cur_reads = ''
cur_pos = ''
cur_len = ''
cur_score = ''
cur_split = False
cur_identity = ''

for line in fileinput.input():
    line = line.rstrip().split('\t')
    # only use information from the first record for each read
    # (assuming the secondary/splits are right next to the primary)
    if line[0] == cur_reads:
        # check if overlap in read space
        # if no, flag as a "split read"
        if (not cur_split):
            if (int(line[3]) > cur_pos + cur_len) or \
               (int(line[3]) + len(line[4]) < cur_pos):
                cur_split = True
        continue
    else:
        # it's a new read, so write the previous one
        if cur_reads != '':
            # increment records
            rid = '{}\t{}\t{}'.format(cur_identity, cur_score, cur_split)
            if rid not in records:
                records[rid] = 1
            else:
                records[rid] += 1
            cur_split = False
    cur_reads = line[0]
    cur_identity = int(100*float(line[1]))
    cur_score = line[2]
    cur_pos = int(line[3])
    cur_len = len(line[4])

# print the counts for each rid/profiles
for rec in records.keys():
    print('{}\t{}'.format(records[rec], rec))
