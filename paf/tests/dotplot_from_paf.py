from mpl_toolkits.axes_grid1 import make_axes_locatable
from collections import defaultdict
from matplotlib import collections
from matplotlib import pyplot
from matplotlib import colors
from matplotlib import cm
import matplotlib
import argparse
import random
import numpy
import sys

from zlib import crc32


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap


HSV_COLORMAP = cm.get_cmap('hsv', 256)
GNUPLOT_COLORMAP = truncate_colormap(cm.get_cmap('gnuplot2', 256), minval=0.12, maxval=0.80)
VIRIDIS_COLORMAP = cm.get_cmap('viridis', 256)


"""
https://stackoverflow.com/a/42909410
"""
def bytes_to_float(b):
    return float(crc32(b) & 0xffffffff) / 2**32


def str_to_float(s, encoding="utf-8"):
    return bytes_to_float(s.encode(encoding))


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
    def __init__(self, paf_line, store_tags=False):
        if store_tags:
            self.tokens = paf_line.strip().split()
        else:
            self.tokens = paf_line.strip().split()[:12]

        self.tokens[1] = int(self.tokens[1])
        self.tokens[6] = int(self.tokens[6])
        self.tokens[2] = int(self.tokens[2])
        self.tokens[3] = int(self.tokens[3])
        self.tokens[7] = int(self.tokens[7])
        self.tokens[8] = int(self.tokens[8])
        self.tokens[11] = int(self.tokens[11])
        self.tokens[4] = (self.tokens[4] == '-')

        self.cigar_index = None
        self.cigar = None

    def __str__(self):
        s = "query_name: " + str(self.get_query_name()) + '\n' + \
            "ref_name: " + str(self.get_ref_name()) + '\n' + \
            "query_length: " + str(self.get_query_length()) + '\n' + \
            "ref_length: " + str(self.get_ref_length()) + '\n' + \
            "query_start: " + str(self.get_query_start()) + '\n' + \
            "query_stop: " + str(self.get_query_stop()) + '\n' + \
            "ref_start: " + str(self.get_ref_start()) + '\n' + \
            "ref_stop: " + str(self.get_ref_stop()) + '\n' + \
            "map_quality: " + str(self.get_map_quality()) + '\n' + \
            "reversal:" + str(self.get_reversal())

        return s

    def get_data_by_column(self, c):
        data = self.tokens[c]

        # Parse differently if index indicates a tag field
        if c > 11:
            data = data.split(":")

            if len(data) < 2:
                exit("ERROR: specified a tag that has no colon delimiter")

            # Only keep what follows the last ":" delimiter
            data = data[-1]

        return data

    def get_query_name(self):
        return self.tokens[0]

    def get_ref_name(self):
        return self.tokens[5]

    def get_query_length(self):
        return self.tokens[1]

    def get_ref_length(self):
        return self.tokens[6]

    def get_query_start(self):
        return self.tokens[2]

    def get_query_stop(self):
        return self.tokens[3]

    def get_ref_start(self):
        return self.tokens[7]

    def get_ref_stop(self):
        return self.tokens[8]

    def get_map_quality(self):
        return self.tokens[11]

    def get_reversal(self):
        return self.tokens[4]

    def get_cigar(self):
        if self.cigar_index is None:
            self.find_cigar_index()

        if self.cigar is None:
            self.cigar = parse_cigar_as_tuples(self.tokens[self.cigar_index].split("cg:Z:")[-1])
            self.tokens[self.cigar_index] = ""

        return self.cigar

    def find_cigar_index(self):
        self.cigar_index = self.find_tag_index("cg:Z:")

    def find_tag_index(self, tag_prefix):
        for t,token in enumerate(self.tokens[12:]):
            if token.startswith(tag_prefix):
                return 12 + t

    def find_tag(self, tag_substring):
        for token in self.tokens[12:]:
            if token.startswith(tag_substring):
                return token


def plot_abridged_alignment(paf_element, plot_data, use_random_color, use_endpoints, color_index=11, color_scale_max=60):
    if paf_element.get_reversal():
        x1 = paf_element.get_ref_stop()
        x2 = paf_element.get_ref_start()

    else:
        x1 = paf_element.get_ref_start()
        x2 = paf_element.get_ref_stop()

    y1 = paf_element.get_query_start()
    y2 = paf_element.get_query_stop()

    plot_data.lines.append([(x1,y1), (x2,y2)])

    if use_random_color:
        color = HSV_COLORMAP(random.uniform(0,1))
        plot_data.colors.append(color)
    else:
        data = paf_element.get_data_by_column(color_index)

        if color_scale_max is not None:
            color = GNUPLOT_COLORMAP((float(data) + 1e-9)/color_scale_max)
        else:
            color = HSV_COLORMAP(str_to_float(str(data)))

        plot_data.colors.append(color)

    if use_endpoints:
        plot_data.dots_x.append([x1,x2])
        plot_data.dots_y.append([y1,y2])


def plot_full_alignment(paf_element, plot_data, color_index, color_scale_max, use_random_color, use_endpoints):
    color_key = {'M':"blue",
                 'I':"orange",
                 'D':"orange",
                 '=':"blue",
                 'X':"purple",
                 'S':"black",
                 'H':"black"}

    print(paf_element.get_query_name(), paf_element.get_ref_name(), paf_element.get_reversal(), paf_element.get_ref_start())
    print(paf_element.get_cigar()[:10])

    ref_index = paf_element.get_ref_start()
    query_index = paf_element.get_query_start()

    color = None
    gray = (0.5,0.5,0.5,0.1)

    if use_random_color:
        hsv = cm.get_cmap('hsv', 256)
        color = hsv(random.uniform(0,1))
    elif color_index is not None:
        data = paf_element.get_data_by_column(color_index)

        color = None
        if color_scale_max is not None:
            # Interpret as a scalar if specified by user
            color = GNUPLOT_COLORMAP((float(data) + 1e-9)/color_scale_max)
        else:
            # Interpret as a deterministic function of string (by default)
            color = HSV_COLORMAP(str_to_float(str(data)))

    cigar_operations = paf_element.get_cigar()
    if paf_element.get_reversal():
        cigar_operations = reversed(cigar_operations)
        ref_index = paf_element.get_ref_stop()
        query_index = paf_element.get_query_start()

    for o,operation in enumerate(cigar_operations):
        x1 = ref_index
        x2 = x1 + int(is_reference_move(operation[0]))*(1-2*int(paf_element.get_reversal()))*operation[1]
        y1 = query_index
        y2 = y1 + int(is_query_move(operation[0]))*operation[1]

        # print(operation, is_reference_move(operation[0]), is_query_move(operation[0]), "(%d,%d) -> (%d,%d)" % (x1,y1,x2,y2))

        ref_index = x2
        query_index = y2

        if operation[0] != "S" and operation[0] != "H":
            plot_data.lines.append([(x1,y1), (x2,y2)])

            if use_random_color:
                if operation[0] == 'M' or operation[0] == '=' or operation[0] == 'X':
                    plot_data.colors.append(color)
                else:
                    plot_data.colors.append(gray)
            elif color_index is not None:
                plot_data.colors.append(color)

            else:
                plot_data.colors.append(color_key[operation[0]])

    if use_endpoints:
        if paf_element.get_reversal():
            x1 = paf_element.get_ref_stop()
            x2 = paf_element.get_ref_start()
        else:
            x1 = paf_element.get_ref_start()
            x2 = paf_element.get_ref_stop()

        y1 = paf_element.get_query_start()
        y2 = paf_element.get_query_stop()

        plot_data.dots_x.append([x1,x2])
        plot_data.dots_y.append([y1,y2])


class PlotData:
    def __init__(self):
        self.lines = list()
        self.dots_x = list()
        self.dots_y = list()
        self.colors = list()


def dotplot_from_paf(paf_path,
                     min_mapq,
                     min_length,
                     color_by,
                     color_scale_max,
                     use_full_alignment,
                     use_random_color,
                     use_endpoints):

    plot_data_per_alignment_pair = defaultdict(lambda: PlotData())
    figures_per_alignment_pair = defaultdict(lambda: pyplot.figure())
    lengths_per_alignment_pair = dict()

    # ref_names = set()
    # lines = list()
    # dots_x = list()
    # dots_y = list()
    # colors = list()

    color_index = None
    store_tags = use_full_alignment or (color_by is not None)

    with open(paf_path, 'r') as file:
        for l,line in enumerate(file):
            paf_element = PafElement(paf_line=line, store_tags=store_tags)
            pair_identifier = paf_element.get_query_name() + "_VS_" + paf_element.get_ref_name()

            if color_by is not None:
                if color_index is None:
                    if color_by.isdigit():
                        color_index = int(color_by)
                    elif color_by.isalnum():
                        color_index = paf_element.find_tag_index(color_by)

                    data = paf_element.get_data_by_column(color_index)

                    # Check if this field is compatible with user's choice for arbitrary vs numeric
                    if color_scale_max is not None:
                        if type(data) == str and not data.isnumeric():
                            exit("ERROR: cannot interpret data as numeric: \n\t" + data)

            lengths_per_alignment_pair[pair_identifier] = max(paf_element.get_ref_length(), paf_element.get_query_length())

            if paf_element.get_map_quality() < min_mapq:
                continue

            if abs(paf_element.get_ref_stop() - paf_element.get_ref_start()) < min_length:
                continue

            if use_full_alignment:
                plot_full_alignment(
                    paf_element=paf_element,
                    plot_data=plot_data_per_alignment_pair[pair_identifier],
                    color_index=color_index,
                    color_scale_max=color_scale_max,
                    use_random_color=use_random_color,
                    use_endpoints=use_endpoints)
            else:
                plot_abridged_alignment(
                    paf_element=paf_element,
                    plot_data=plot_data_per_alignment_pair[pair_identifier],
                    color_index=color_index,
                    color_scale_max=color_scale_max,
                    use_random_color=use_random_color,
                    use_endpoints=use_endpoints)

    for pair_identifier,data in plot_data_per_alignment_pair.items():
        figure = figures_per_alignment_pair[pair_identifier]
        axes = figure.add_subplot()

        line_collection = collections.LineCollection(data.lines, colors=data.colors, linewidths=1.2)
        axes.add_collection(line_collection)

        axes.scatter(data.dots_x, data.dots_y, color="black", s=0.3, zorder=sys.maxsize)

        axes.set_aspect('equal')
        axes.set_ylim([0,lengths_per_alignment_pair[pair_identifier]])
        axes.set_xlim([0,lengths_per_alignment_pair[pair_identifier]])

        axes.set_xlabel("Reference")
        axes.set_ylabel("Query")

        axes.set_title(pair_identifier)

        if not use_random_color and not color_scale_max is None:
            # create an axes on the right side of ax. The width of cax will be 5%
            # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            divider = make_axes_locatable(axes)
            cax = divider.append_axes("right", size="5%", pad=0.05)

            colorbar_tick_labels = [round(x,3) for x in numpy.linspace(0, color_scale_max, 10)]
            colorbar_ticks = numpy.linspace(0, 1, 10)
            print(colorbar_ticks)

            colorbar = matplotlib.colorbar.ColorbarBase(cax, cmap=GNUPLOT_COLORMAP)
            colorbar.set_label(color_by)
            colorbar.set_ticks(colorbar_ticks)
            colorbar.set_ticklabels(colorbar_tick_labels)

            figure.tight_layout()

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
        "--color_by",
        type=str,
        required=False,
        default="0",
        help="Which data in the PAF to color the alignments with. If a string, it must indicate a tag (e.g. 'cm') "
             "and that tag must be found past column 11 (0-based) and be followed by a colon. If a number (integer), "
             "it must indicate the (0-based) column in which the data is expected to be found."
    )
    parser.add_argument(
        "--scalar_color",
        type=float,
        required=False,
        default=None,
        help="Intended for use in combination with 'color_by' parameter. For integers/floats, this parameter indicates"
             "that the alignment should be interpreted as a quantitative value, and the provided value indicates the "
             "maximum of the color range. If this parameter is unused, the data specified by the 'color_by' param "
             "will be interpreted as a string identifier and assigned a random (deterministic) color."
    )
    parser.add_argument(
        "--min_mapq","-q",
        type=int,
        required=False,
        default=0,
        help="Minimum map quality to plot"
    )
    parser.add_argument(
        "--min_length","-l",
        type=int,
        required=False,
        default=0,
        help="Minimum length alignment (in ref coord) to plot"
    )
    parser.add_argument(
        "--use_cigar","-c",
        dest="use_cigar",
        required=False,
        action="store_true",
        help="Add this boolean flag to optionally parse and plot the full base-level alignment (can be prohibitively slow)"
    )
    parser.add_argument(
        "--random_color","-r",
        dest="random_color",
        required=False,
        action="store_true",
        help="Add this boolean flag to optionally plot all alignments with a random individual color (as default matplotlib color cycle)"
    )
    parser.add_argument(
        "--endpoints","-e",
        dest="endpoints",
        required=False,
        action="store_true",
        help="Add this boolean flag to optionally plot black circle delimiters on the ends of alignment lines"
    )

    args = parser.parse_args()

    dotplot_from_paf(args.i,
                     min_mapq=args.min_mapq,
                     min_length=args.min_length,
                     color_by=args.color_by,
                     color_scale_max=args.scalar_color,
                     use_full_alignment=args.use_cigar,
                     use_random_color=args.random_color,
                     use_endpoints=args.endpoints)