# The TAF File Format

This is a specification for a "transposed alignment format" or maybe "terrific alignment format" (.taf). The idea 
is to describe the alignment as a series of columns, one per line, with bases in the
columns run-length encoded, and row coordinates given as needed after the bases in each column.
Where coordinates are not given it is assumed the coordinates continue from the previous column.

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

tag -> key':'value

key -> alphanumerical string
value -> alphanumerical string

columns -> column '\n' columns
           column 'n'

column -> bases coordinates_column tag_string 
       -> bases coordinates_column
       -> bases tag_string
       -> bases

coordinates_column -> ';' coordinates

coordinates -> coordinate_operation coordinates 
            -> coordinate_operation

coordinate_operation -> 'i' row coordinate
                     -> 's' row coordinate
                     -> 'g' row gap_length
                     -> 'd' row 
                     -> 'a' coordinate

(The 'i' stands for insertion, the 's' for substitution, 'g' for gap, the 'd' for deletion and the 'a' for append. These operations
allow us to update the coordinates of the sequences as we go and work as their name suggests. Rows are indexed from zero. 
The operations are affected in order, so, for example, inserting a row at position i will shift all remaining column indices by one place
for any operations that are specified after an insertion.
The substitute operation can be used to change the coordinates for a sequence in a row, or it can be used to periodically repeat
coordinates to prevent needing to scan back more than N rows to find the coordinates of a row.
Appends always add rows at the end.
Using the 'g' gap operation can be used to increment the coordinate of a sequence by a specified amount and is useful to inserting
long substrings.)

row -> integer >= 0

gap_length -> integer >= 0

coordinate -> sequence_name offset strand

sequence_name -> string

offset -> integer >= 0

strand -> '+'
       -> '-'

bases -> base count bases
      -> base count

count -> integer > 0

base -> alphabet character or '-' or '*' ([A-Z,a-z,-,*])

The following shows by example the translation between the MAF format and TAF format.

The MAF file (602 bytes):
"
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
s       simRat.chr6     470350  8       +       647215  TTTTACCG
"

The corresponding TAF file (308 bytes):
"
# version=1 scoring=N/A
C 1 T 3 ; a simDog.chr6 437451 + a simHuman.chr6 446327 + a simMouse.chr6 460751 + a simRat.chr6 470339 11
C 2 T 2
C 4
G 2 A 2
C 3 T 1
A 4
G 1 A 1 G 2
T 1 G 1 A 1 G 1
G 4
T 4
T 5 ; i 0 simCow.chr6 445326 + s 1 simDog.chr6 437462 + 
T 5
T 1 - 1 C 1 T 2
T 5
C 1 T 2 A 2
C 5
C 5
A 1 G 4
"

