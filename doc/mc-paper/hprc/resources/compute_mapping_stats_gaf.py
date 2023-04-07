import fileinput

# parse the GAF record and count how many reads
# in each profile (mapping quality x perfect alignement x aligment score)
#    mapping quality: number. '-1' for unmapped reads
#    perfect: boolean to specify if the reads aligned perfectly
#    alignment score: value of the AS tag

# to tally the number of reads/records in each profile
records = {}

for line in fileinput.input():
    line = line.rstrip().split('\t')
    ## handle unmapped reads
    if line[9] == "*":
        line[11] = '-1'
    ## extract AS tag
    as_tag = False
    for ii in range(11, len(line)):
        tag = line[ii].split(':')
        if tag[0] == 'AS':
            as_tag = tag[2]
            break
    ## check if perfectly aligned
    perfect = False
    if line[9] != "*" and line[9] == line[10]:
        perfect = True
    ## increment records
    rid = '{}\t{}\t{}'.format(line[11], perfect, as_tag)
    if rid not in records:
        records[rid] = 1
    else:
        records[rid] += 1

# print the counts for each rid/profiles
for rec in records.keys():
    print('{}\t{}'.format(records[rec], rec))
