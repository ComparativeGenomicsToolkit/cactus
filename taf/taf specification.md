# The TAF File Format

This is a specification for a "transposed alignment format", or maybe "terrific alignment format" (.taf). The idea 
is to describe the alignment as a series of columns, one per line, with bases in the
columns optionally run-length encoded, and row coordinates given as needed after the bases in each column.
Where coordinates are not given it is assumed the coordinates continue from the previous column.
User defined tags, as key:value pairs, are also included to allow the alignemnt
columns to be annotated.

The format supports line based indexing for rapid retrieval of any column or contiguous sequence of columns from 
the file. The format is intentionally simple but should prove quite space efficient for large alignments. 

Its key potential benefits over the MAF format are that it is easier to parse (no blocks), does not suffer the same
issue with fragmentation as the number of sequences grows (there are no blocks!), is often less verbose (particularly 
for large alignments), is very easy to index as each column of the alignment is a single line, and supports 
extensible column annotations.

The format is composed of a sequence of a header and some columns. 
Tokens are separated by white-space. The syntax is defined as follows:

    .taf -> header '\n' columns
    
    header -> tag_string
    
    tag_string -> '#' tags
    
    tags -> tag tags
         -> tag
         -> empty_string   
    
    tag -> key':'value
    
    key -> alphanumerical string
    value -> alphanumerical string

    empty_string -> zero length string, e.g. ""

(The header is therefore just a possibly empty sequence of key:value pairs)
    
    columns -> column '\n' columns
               column 'n'
    
    column -> bases coordinates_column tag_string 
           -> bases coordinates_column
           -> bases tag_string
           -> bases

(Columns are encoded, one per line. Each column encodes an alignment of the
bases in the column (see below) and then some optional coordinates and tags.)

    bases -> run_length_encoded_bases
          -> bases

    run_length_encoded_bases -> base count run_length_encoded_bases
                             -> base count
    
    count -> integer > 0
    
    base -> alphabet character or '-' or '*' ([A-Z,a-z,-,*])

    bases -> string or alphabet characters or '-' or '*' ([A-Z,a-z,-,*]+)

(Aligned bases are either a run-length encoded
representation or simply as a sequence of aligned characters. To specify
which format to use we use the tag "run_length_encoded_bases=1" in the header.
If "run_length_encoded_bases=0" or the tag is not specified the format is
to NOT use run length encoding.)

    
    coordinates_column -> ';' coordinates
    
    coordinates -> coordinate_operation coordinates 
                -> coordinate_operation
                -> empty_string
    
    coordinate_operation -> 'i' row coordinate
                         -> 'd' row
                         -> 's' row coordinate
                         -> 'g' row gap_length
                         -> 'G' row gap_string
                     

(The 'i' stands for insertion, the 'd' for deletion, the 's' for substitution, 'g' and 'G' for gap. These operations
allow us to update the coordinates of the sequences as we go, and work as their name suggests. Rows are indexed from zero. 
The operations are affected in order, so, for example, inserting a row at position i will shift all remaining column indices by one place
for any operations that are specified after an insertion.
The substitute operation can be used to change the coordinates for a sequence in a row, or it can be used to periodically repeat
coordinates to prevent needing to scan back more than N rows to find the coordinates of a row.
Using the 'g' gap operation can be used to increment the coordinate of a sequence by a specified amount and is useful to inserting
long substrings. Using the 'G' instead of 'g' allows one to instead specify the gap substring.)

    row -> integer >= 0
    
    gap_length -> integer >= 0
    
    gap_string -> alphabet string ([A-Z,a-z]*)
    
    coordinate -> sequence_name offset strand
    
    sequence_name -> string without white space
    
    offset -> integer >= 0
    
    strand -> '+'
           -> '-'


The following shows by example the translation between the MAF format and TAF format.

The MAF file (602 bytes):

    ##maf version=1 scoring=N/A
    
    a
    s       simDog.chr6     437451  11      +       593897  CCCGTCAGTGT
    s       simHuman.chr6   446327  11      +       601863  TCCGCCAAGGT
    s       simMouse.chr6   460751  11      +       636262  TTCATCAGAGT
    s       simRat.chr6     470339  11      +       647215  TTCATTAGGGT
    
    a
    s       simCow.chr6     445326  8       +       602619  TTTTCCCA
    s       simDog.chr6     437462  8       +       593897  TT-TTCCG
    s       simHuman.chr6   446338  8       +       601863  TTCTTCCG
    s       simMouse.chr6   460762  8       +       636262  TTTTACCG
    s       simRat.chr6     470355  8       +       647215  TTTTACCG


The corresponding TAF file (262 bytes):

    # version:1 scoring:N/A
    CTTT ; i 0 simDog.chr6 437451 + i 1 simHuman.chr6 446327 + i 2 simMouse.chr6 460751 + i 3 simRat.chr6 470339 11
    CCTT
    CCCC
    GGAA
    CCCT
    AAAA
    GAGG
    TGAG
    GGGG
    TTTT
    TTTTT ; i 0 simCow.chr6 445326 + g 4 5
    TTTTT
    T-CTT
    TTTTT
    CTTAA
    CCCCC
    CCCCC
    AGGGG


