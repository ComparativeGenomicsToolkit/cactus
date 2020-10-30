#todo: consider making extract_single_mappings clever-er so that it doesn't have to scan every possible "single_mapping_region" when its determining what regions to extract.

import copy
import operator
from toil.job import Job
from cigar import Cigar
import collections as col

def parse_large_mappings(job, paf, min_size_mapping, min_mapq):
    mappings = dict()
    with open(job.fileStore.readGlobalFile(paf)) as inf:
        for line in inf:
            parsed = line.split()

            for i in range(1,4):
                parsed[i] = int(parsed[i])
            for i in range(6,12):
                parsed[i] = int(parsed[i])
                
            if (parsed[11] != 255 and parsed[11] >= min_mapq) and parsed[10] >= min_size_mapping:
                if parsed[5] not in mappings:
                    # add the reference chrom to the mappings dict:
                    mappings[parsed[5]] = list()
                mappings[parsed[5]].append(parsed)
    return mappings

def get_single_mapping_regions(mappings):
    single_mapping_regions = dict()
    for chrom, mapping_list in mappings.items():
        # get points, in format (ref_position, start_bool). 
        # Start_bool indicates whether point is a start of mapping, or end of mapping.
        mapping_points = list()
        for mapping in mapping_list:
            mapping_points.append((int(mapping[7]), True))
            mapping_points.append((int(mapping[8]), False))

        # sort the points by position on the reference.
        mapping_points.sort(key=operator.itemgetter(0))
        
        # get regions without overlap:
        single_mapping_regions[chrom] = list()
        cur_start = 0
        cur_depth = 0
        debug_i = 0
        for point in mapping_points:
            if cur_depth == 1:
                # then we're about to conclude a region on the reference which had exactly one mapping. Record the region!
                if (point[0] - cur_start) != 0:
                    single_mapping_regions[chrom].append((cur_start, point[0]))
                cur_start = None

            if point[1]:
                cur_depth += 1
            else:
                cur_depth -= 1

            if cur_depth == 1:
                # then we've begun a region on the ref with exactly one mapping. Save start position.
                cur_start = point[0]

    # print("single_mapping_regions", single_mapping_regions)

    return single_mapping_regions

def drop_unadjusted_fields(mapping):
    """Drops all fields beyond mapping[0:12], except for the cigar and alignment type field.

    Args:
        mapping ([type]): [description]
    """
    # print("mapping before drop:", mapping)
    # Find the cigar and mapping-type fields.
    keep_fields = list()
    for field in range(len(mapping[12:])):
        field += 12
        if mapping[field][:5] == "cg:Z:" or mapping[field][:5] == "tp:A:":
            keep_fields.append(field)
    
    fixed_mapping = mapping[:12]
    for kept_field in keep_fields:
        fixed_mapping.append(mapping[kept_field])
    # print("mapping after drop:", fixed_mapping)
    
    # drop all fields after field 11, except for the cigar.
    return fixed_mapping
    
def adjust_mapping(mapping, overlap_region):
    # don't adjust the original mapping, adjust a copy.
    mapping = copy.deepcopy(mapping)
    
    debug_original_mapping = mapping[:]
    # find the cigar & its field from the SAM-liked typed key-value pairs at the end.
    cig = None
    cig_field = None
    for field in range(len(mapping[12:])):
        field += 12
        if mapping[field][:5] == "cg:Z:":
            cig = col.deque(Cigar(mapping[field][5:]).items())
            cig_field = field
            break
    # print ("adjusting mapping", cig, "overlap region", overlap_region)
        
    ## Note: mapping[X], X=2,3 is query start & end; X=7,8 is target start and end; 
    #        9 is # of M bases; 10 is # of M+I+D bases
    
    if mapping[7] < overlap_region[0]:
        #adjust start of the mapping
        adjust_amount = overlap_region[0] - mapping[7]
        while adjust_amount > 0:
            if cig[0][1] == "M":
                if cig[0][0] <= adjust_amount:
                    adjust_amount -= cig[0][0]
                    mapping[2] += cig[0][0]
                    mapping[7] += cig[0][0]
                    # mapping[9] -= cig[0][0]
                    mapping[10] -= cig[0][0]
                    cig.popleft()
                else:
                    # we only need to slightly adjust this cig field, because the adjust amount is less than the cig field.
                    cig[0] = (cig[0][0] - adjust_amount, cig[0][1])
                    mapping[2] += adjust_amount
                    mapping[7] += adjust_amount
                    # mapping[9] -= adjust_amount
                    mapping[10] -= adjust_amount
                    adjust_amount = 0
            elif cig[0][1] == "D":
                if cig[0][0] <= adjust_amount:
                    adjust_amount -= cig[0][0]
                    mapping[7] += cig[0][0]
                    mapping[10] -= cig[0][0]
                    cig.popleft()
                else:
                    # we only need to slightly adjust this cig field, because the adjust amount is less than the cig field.
                    cig[0] = (cig[0][0] - adjust_amount, cig[0][1])
                    mapping[7] += adjust_amount
                    mapping[10] -= adjust_amount
                    adjust_amount = 0
            elif cig[0][1] == "I":
                mapping[2] += cig[0][0]
                mapping[10] -= cig[0][0]
                cig.popleft()
            elif cig[0][1] in "HS": 
                cig.popleft()
            else:
                raise ValueError('Inappropiate Cigar value ' + cig[0][1] + ' in cigar ' + mapping[cig_field])

    # print ("mapping with adj start", cig, "overlap region", overlap_region)

    if mapping[8] > overlap_region[1]:
        #adjust end of the mapping
        adjust_amount = mapping[8] - overlap_region[1]
        while adjust_amount > 0:
            if cig[-1][1] == "M":
                if cig[-1][0] <= adjust_amount:
                    adjust_amount -= cig[-1][0]
                    mapping[3] -= cig[-1][0]
                    mapping[8] -= cig[-1][0]
                    # mapping[9] -= cig[-1][0]
                    mapping[10] -= cig[-1][0]
                    cig.pop()
                else:
                    # we only need to slightly adjust this cig field, because the adjust amount is less than the cig field.
                    cig[-1] = (cig[-1][0] - adjust_amount, cig[-1][1])
                    mapping[3] -= adjust_amount
                    mapping[8] -= adjust_amount
                    mapping[9] -= adjust_amount
                    mapping[10] -= adjust_amount
                    adjust_amount = 0
            elif cig[-1][1] == "D":
                if cig[-1][0] <= adjust_amount:
                    adjust_amount -= cig[-1][0]
                    mapping[8] -= cig[-1][0]
                    mapping[10] -= cig[-1][0]
                    cig.pop()
                else:
                    # we only need to slightly adjust this cig field, because the adjust amount is less than the cig field.
                    cig[-1] = (cig[-1][0] - adjust_amount, cig[-1][1])
                    mapping[8] -= adjust_amount
                    mapping[10] -= adjust_amount
                    adjust_amount = 0
            elif cig[-1][1] == "I":
                mapping[3] -= cig[-1][0]
                mapping[10] -= cig[-1][0]
                cig.pop()
            elif cig[-1][1] in "HS": 
                cig.pop()
            else:
                raise ValueError('Inappropiate Cigar value ' + cig[-1][1] + ' in cigar ' + mapping[cig_field])

    # print ("mapping with adj end", cig, "overlap region", overlap_region)

    # if mapping[7] == 28596469 or mapping[7] == 17555049:
    #     print("mapping of interest before final adjustment:", mapping)
    #     print("cigar", cig)

    #todo: consider if there may be no mapping left - just a deletion field. Is it possible to accidentally find an "overlap region" that just shows unmapped reference, i.e. a deletion? This would mean below while-loop would crash eventually. Need to catch that eventuality and have it return a None mapping, instead.
    while cig[-1][1] in "ID":
        #while the cigar ends on an insertion/deletion, (which doesn't make biological sense), remove it.
        if cig[-1][1] == "D":
            mapping[8] -= cig[-1][0]
            mapping[10] -= cig[-1][0]
            cig.pop()
        elif cig[-1][1] == "I":
            mapping[3] -= cig[-1][0]
            mapping[10] -= cig[-1][0]
            cig.pop()

    # print ("dropped all in 'ID' in end.", cig, "overlap region", overlap_region)

            
    while cig[0][1] in "ID":
        #while the cigar starts on an insertion/deletion, (which doesn't make biological sense), remove it.
        if cig[0][1] == "D":
            mapping[7] += cig[0][0]
            mapping[10] -= cig[0][0]
            cig.popleft()
        elif cig[0][1] == "I":
            mapping[2] += cig[0][0]
            mapping[10] -= cig[0][0]
            cig.popleft()

    # if mapping[7] == 28596469 or mapping[7] == 17555049:
    #     print("mapping of interest after final adjustment:", mapping)
    #     print("cigar", cig)

    # adjust the "match count" field.
    mapping[9] = 0
    for i in cig:
        if i[1] == "M":
            mapping[9] += i[0]

    adjusted_cigar = ""
    for i in cig:
        adjusted_cigar += str(i[0]) + i[1]

    mapping[cig_field] = "cg:Z:" + adjusted_cigar

    mapping = drop_unadjusted_fields(mapping)

    return mapping

def extract_single_mappings(parsed_mappings, single_mapping_regions, min_var_len):
    """Extracts any segments of mappings >=min_var_len that don't overlap with any other
     mappings in parsed_mappings.

    Args:
        parsed_mappings ([type]): [description]
        single_mapping_regions ([type]): [description]

    Raises:
        ZeroDivisionError: [description]
    """
    print(single_mapping_regions)
    extracted_mappings = dict()
    for chrom, mappings in parsed_mappings.items():
        extracted_mappings[chrom] = list()
        for mapping in mappings:
            # if mapping[8] == 28596469:
            #     print("we're at the mapping of interest.", mapping[0:11])
                
            # i = 0
            # print(single_mapping_regions, mapping, single_mapping_regions)
            # print(single_mapping_regions, mapping[8], single_mapping_regions[chrom][0][0])
            # print(int(mapping[7]), (i < len(single_mapping_regions)), (single_mapping_regions[chrom][i][0] < int(mapping[8])))
            # print("looking at mapping:", mapping[7], mapping[8])

            if mapping[10] >= min_var_len:
                # if mapping[8] == 28596469:
                #     print("len(single_mapping_regions[chrom])", len(single_mapping_regions[chrom]))
                #     print("list(range(len(single_mapping_regions[chrom])))", list(range(len(single_mapping_regions[chrom]))))
                for i in range(len(single_mapping_regions[chrom])):
                    # if mapping[8] == 28596469:
                    #     print("in the for loop. with mapping of interest at int i", i)
                        # print("single_mapping_regions", single_mapping_regions)

                # while (i < len(single_mapping_regions[chrom])) and (single_mapping_regions[chrom][i][0] < mapping[8]) and (mapping[10] >= min_var_len):
                    # while there is still mapping regions to look at, and we're still in mapping regions that potentially overlap the mapping:

                    # overlap_region = (max(single_mapping_regions[chrom][i][0], mapping[7]), min(single_mapping_regions[chrom][i][1], mapping[8])+1)
                    overlap_region = (max(single_mapping_regions[chrom][i][0], mapping[7]), min(single_mapping_regions[chrom][i][1], mapping[8]))

                    
                    #debug
                    # if mapping[8] == 28596469 and single_mapping_regions[chrom][i][1] in [28625484, 28626723]:
                    #     print("single_mapping_region at", i, "is", single_mapping_regions[chrom][i])                
                    #     print("overlap_region at", i, "is", overlap_region)

                    # print(overlap_region)
                    if overlap_region[1] > overlap_region[0]:
                        #then we have found an overlap.
                        adj_mapping = adjust_mapping(mapping, overlap_region)
                        if adj_mapping is not None:
                            extracted_mappings[chrom].append(adj_mapping)
                        
                    # if mapping[8] == 28596469:
                    #     print("end of for loop. with mapping of interest at int i", i)
                    #     print("len(single_mapping_regions[chrom])", len(single_mapping_regions[chrom]))

                    # i += 1

    return extracted_mappings


def apply_dipcall_bed_filter(job, paf, min_var_len=50000, min_size_mapping=10000, min_mapq=5):
# def apply_dipcall_bed_filter(job, paf, min_var_len=5000, min_size_mapping=1000, min_mapq=5):
    """Built to imitate the filter used to produce the bedfile in dipcall.
    Guarantees that there will be no overlapping mappings.
    * includes mappings >=min_var_len in size.
    * BUT: will exclude regions of these mappings >=min_size_mapping in size.
    * all mappings considered for inclusion or overlap must have >= 5 mapQ.

    Args:
        job ([type]): [description]
        paf ([type]): [description]

    Returns:
        [type]: [description]
    """
    # First, compile a sorted list of all the mappings >=min_size_mapping in size, with mapq>=min_mapq
    # key of mappings is reference chrom, value is above sorted list for mappings to that ref.
    parsed_mappings = parse_large_mappings(job, paf, min_size_mapping, min_mapq)

    # Find the start, stop pairs that are regions on the ref covered by exactly one mapping.
    single_mapping_regions = get_single_mapping_regions(parsed_mappings)

    # Extract the parts of mappings that are involved in only a single mapping. Only include mappings with original length  >= min_var_len.
    parsed_single_mappings = extract_single_mappings(parsed_mappings, single_mapping_regions, min_var_len)

    #todo: planned testing stages:
    # 3rd I will see if I can cut those mappings to produce new valid mappings where there is overlap with min_size_mapping. 
    #TODO: very important! Make some unit tests to make sure I get out reasonable mappings with permissible cigars (that match up with the rest of the mapping data), that truly are only 1-depth.
    
    paf_filtered = job.fileStore.getLocalTempFile()
    with open(paf_filtered, "w") as outf:
        for mapping_list in parsed_single_mappings.values():
            for mapping in mapping_list:

                # convert numerical fields back to str.
                for i in range(1,4):
                    mapping[i] = str(mapping[i])
                for i in range(6,12):
                    mapping[i] = str(mapping[i])
                
                outf.write('\t'.join(mapping) + '\n')
                
    return job.fileStore.writeGlobalFile(paf_filtered)