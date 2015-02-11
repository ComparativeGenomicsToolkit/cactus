# Rules for converting the txt documents to html.

DOCS = FAQ.html bisulfite.html last-map-probs.html		\
last-pair-probs.html last-postmask.html last-split.html		\
last-tutorial.html last.html maf-convert.html maf-cull.html

all: ${DOCS}

${DOCS}: last-doc.css Makefile

.SUFFIXES: .html .txt

# Below, we face 2 problems:
# 1. Figure out where html4css1.css is.
# 2. Work with rst2html version 0.4.1, which only allows one stylesheet-path.

# Ugh!  Is there a better way?
RST_CSS = `locate html4css1.css | head -n1`

.txt.html:
	cat ${RST_CSS} last-doc.css | rst2html --stylesheet-path=/dev/stdin --initial-header-level=2 --no-doc-info $< $@
