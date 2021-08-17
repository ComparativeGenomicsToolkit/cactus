from matplotlib import pyplot
from matplotlib import cm
import argparse
import os


def parse_cigar_as_tuples(cigar_string):
    operations = list()

    length_string = ""
    prev_is_numeric = True

    for i,c in enumerate(cigar_string):
        if c.isalpha():
            if i > 0 and not prev_is_numeric:
                exit("ERROR: cigar string contains impossible sequence of numeric and alphabetic characters")

            operations.append((c, int(length_string)))
            length_string = ""
            prev_is_numeric = False

        if c.isnumeric():
            length_string += c
            prev_is_numeric = True

    return operations


def is_reference_move(cigar_type):
    if cigar_type == 'M':
        return True
    elif cigar_type == 'I':
        return False
    elif cigar_type == 'D':
        return True
    elif cigar_type == '=':
        return True
    elif cigar_type == 'X':
        return True
    elif cigar_type == 'S':
        return False
    elif cigar_type == 'H':
        return False
    else:
        exit("ERROR: unrecognized cigar type: " + cigar_type)


def is_query_move(cigar_type):
    if cigar_type == 'M':
        return True
    elif cigar_type == 'I':
        return True
    elif cigar_type == 'D':
        return False
    elif cigar_type == '=':
        return True
    elif cigar_type == 'X':
        return True
    elif cigar_type == 'S':
        return True
    elif cigar_type == 'H':
        return False
    else:
        exit("ERROR: unrecognized cigar type: " + cigar_type)


def get_ref_alignment_length(cigar_operations):
    l = 0
    for item in cigar_operations:
        l += item[1]*is_reference_move(item[0])

    return l


class PafElement:
    def __init__(self, paf_line, store_cigar=False):
        tokens = paf_line.strip().split('\t')

        self.query_name = tokens[0]
        self.ref_name = tokens[5]
        self.query_length = int(tokens[1])
        self.ref_length = int(tokens[6])
        self.query_start = int(tokens[2])
        self.query_stop = int(tokens[3])
        self.ref_start = int(tokens[7])
        self.ref_stop = int(tokens[8])
        self.map_quality = int(tokens[11])
        self.reversal = (tokens[4] == '-')
        self.cigar_operations = None

        if store_cigar:
            for token in tokens[12:]:
                if token.startswith("cg:Z:"):
                    self.cigar_operations = parse_cigar_as_tuples(token[5:])


def plot_abridged_alignment(paf_element, axes):
    viridis = cm.get_cmap('viridis', 256)

    x1 = paf_element.ref_start
    x2 = paf_element.ref_stop
    y1 = paf_element.query_start
    y2 = paf_element.query_stop

    color = viridis((60 - paf_element.map_quality)/60)
    axes.plot([x1,x2], [y1,y2], color=color, linewidth=0.5)


def plot_full_alignment(paf_element, axes):
    colors = {'M':"blue",
              'I':"orange",
              'D':"orange",
              '=':"blue",
              'X':"purple",
              'S':"black",
              'H':"black"}

    print(paf_element.query_name, paf_element.ref_name, paf_element.reversal, paf_element.ref_start)
    print(paf_element.cigar_operations[:10])

    if paf_element.reversal:
        paf_element.cigar_operations = reversed(paf_element.cigar_operations)

    ref_index = paf_element.ref_start
    query_index = paf_element.query_start

    for o,operation in enumerate(paf_element.cigar_operations):
        x1 = ref_index
        x2 = x1 + int(is_reference_move(operation[0]))*(1-2*int(paf_element.reversal))*operation[1]
        y1 = query_index
        y2 = y1 + int(is_query_move(operation[0]))*(1-2*int(paf_element.reversal))*operation[1]

        # print(operation, is_reference_move(operation[0]), is_query_move(operation[0]), "(%d,%d) -> (%d,%d)" % (x1,y1,x2,y2))

        ref_index = x2
        query_index = y2

        if operation[0] != "S" and operation[0] != "H":
            axes.plot([x1,x2], [y1,y2], color=colors[operation[0]], linewidth=0.5)


def dotplot_from_paf(paf_path, use_full_alignment):
    figure = pyplot.figure()
    axes = pyplot.axes()

    ref_names = set()
    ref_length = 0

    with open(paf_path, 'r') as file:
        for l,line in enumerate(file):
            paf_element = PafElement(paf_line=line, store_cigar=use_full_alignment)

            ref_names.add(paf_element.ref_name)
            ref_length = paf_element.ref_length

            if len(ref_names) > 1:
                exit("ERROR: more than one reference sequence in PAF file")

            if paf_element.map_quality < 1:
                continue

            if use_full_alignment:
                plot_full_alignment(paf_element=paf_element, axes=axes)
            else:
                plot_abridged_alignment(paf_element=paf_element, axes=axes)

    axes.set_aspect('equal')
    axes.set_ylim([0,ref_length])
    axes.set_xlim([0,ref_length])

    axes.set_xlabel("Reference")
    axes.set_ylabel("Query")

    # pyplot.savefig(os.path.join(output_directory, title+"_self_alignment.png"), dpi=200)

    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        type=str,
        required=True,
        help="Path of PAF file"
    )
    parser.add_argument(
        "--use_cigar","-c",
        dest="use_cigar",
        required=False,
        action="store_true",
        help="Add this boolean flag to optionally parse and plot the full base-level alignment (can be prohibitively slow)"
    )

    args = parser.parse_args()

    dotplot_from_paf(args.i, use_full_alignment=args.use_cigar)
