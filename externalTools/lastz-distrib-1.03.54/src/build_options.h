//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: build_options.h
//
//----------
//
// This file contains a list of compile-time optional builds.  Normally these
// build options should be left commented out here, and enabled instead by
// use of the "make" command line, as described in the Makefile.
//
// We show them here for two reasons.  First, it provides a single place to
// describe the options, rather than scatter them over the .c files.  Second,
// it provides a common point to make such definitions and to "ensure" (to the
// extent that we can do so) that all modules will be built with the same
// settings.
//
//----------

#ifndef build_options_H			// (prevent multiple inclusion)
#define build_options_H

//#define allowBackToBackGaps	// if this is defined, gapped_extend.c is
								// .. modified to allow the opening of a delete
								// .. right after an insert, or vice versa


#undef global
#endif // build_options_H
