#-----------
# make-include.mak--
#	Defines variables used by all LASTZ Makefiles
#-----------

INSTALL =  install
ARCH    ?= $(shell uname -m)

installDir = ../../../bin
