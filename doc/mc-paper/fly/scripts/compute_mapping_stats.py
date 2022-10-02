import fileinput

records = {}

for line in fileinput.input():
    if line[0] == "@":
        continue
    line = line.rstrip().split('\t')
    ## skip if secondary or supplementary alignment
    if int(int(line[1])/256)%2 == 1 or int(int(line[1])/2048)%2 == 1:
        continue
    ## handle unmapped reads
    if line[5] == "*":
        line[4] = '-1'
    ## check if perfectly aligned
    perfect = False
    if line[5] == str(len(line[9])) + "M":
        ## check MD tag
        for ii in range(11, len(line)):
            tag = line[ii].split(':')
            if tag[0] == 'MD':
                if line[ii] == "MD:Z:" + str(len(line[9])):
                    perfect = True
                break
    ## increment records
    rid = '{}\t{}'.format(line[4], perfect)
    if rid not in records:
        records[rid] = 1
    else:
        records[rid] += 1

for rec in records.keys():
    print('{}\t{}'.format(records[rec], rec))
