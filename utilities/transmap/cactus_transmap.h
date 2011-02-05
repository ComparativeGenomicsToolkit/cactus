/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef TRANSMAP_H_ 
#define TRANSMAP_H_ 

#include "cactus.h"

/*
 * The basic Transmap API functions.
 *
 * ------- Binary cap "connection", "order" and "orientation" -----
 *
 * Let A and B be two caps. Each cap can be thought of as having
 * complementarily paired copies, one one each strand. Let A and -A indicate the two oppositely stranded copies.
 *
 * A thread is a connected component of adjacency and segment edges that represents
 * a sequence of double stranded DNA. If we observe two caps in the same thread we say they are "connected".
 *
 * Let (A, B) represent a strand of a (sub)thread which starts at A and ends at B.
 * Note, if observe (A, B) then we also observe the 'mirror' (-A, -B).
 *
 * (A, B) and (B, A) (and their mirrors) have a different "order" of the caps,
 *  but keep the same strand, which we call maintaining "orientation".
 * Inversely, (A, B) and (A, -B) (and their mirrors) have a different orientation of the caps,
 * but maintain the same order of the caps.
 *
 * If two caps are contained in a thread, either the thread is a ring or a contig.
 *
 * If the thread is a ring and we observe (A, B) then we also observe (B, A) (and their mirrors).
 * Thus in a ring we can observe either (A, B) and (B, A) or (A, -B) and (-B, A) (and their mirrors).
 * Therefore we can not break the order of the caps in a ring, but we can break their orientation with
 * respect to one another.
 *
 * If the thread is a contig we can observe either (A, B), (A, -B), (B, A) or (-B, A) (and their mirrors), i.e.
 *  we can break order and orientation separately between the caps, giving us
 *  four possible independently observable combinations.
 *
 * -- Virtual threads --
 *
 * A face is "disrupted" if any adjacency edges are added to its bottom or intermediate nodes,
 * or new intermediate nodes are added to the face which are attached. (IS THIS SUFFICIENT?)
 *
 * The set of "virtual threads" for an AVG A is the set of threads that could be created
 * by extension of A without disrupting or creating any new non-trivial faces.
 *
 * -- Cactus Adjacency Recursion and Terminal Threads --
 *
 * An adjacency ADJ(A, B) between two caps A and B contains an interval of sequence between A and B.
 * If this is interval is empty we say ADJ(A, B) is terminal, else it is non-terminal.
 *
 * For a complete cactus-tree and flower N with non-terminal adjacency ADJ(A, B) there exists c(A) and c(B), copies of the
 * caps in a flower which is a child of N (see "Cactus Graphs for Genome Comparison", Recomb 2010).
 *
 * We say a thread is a "terminal thread" if it contains only terminal adjacencies. A non-terminal thread
 * has a "terminal replacement thread", if and only if each non-terminal adjacency ADJ(A, B) can be replaced with a terminal thread
 * starting with A and ending with B. The same terminology applies to virtual threads.
 *
 * -- Basic transmap functions --
 *
 * The following functions take two strand specific caps, A and B and an event E.
 *
 * Let e(A) be the event during which A existed. E must be or be an ancestor of the nearest common ancestor
 * of e(A) and e(B) else an assert error will be created.
 */

/*
 * Returns non-zero if and only if there exists a terminal thread or a virtual terminal thread containing
 * (A, B) and/or (B, A) (and their mirrors) at event E.
 */
bool transmap_connectivityWasPresentAtEvent(Event *E, Cap *A, Cap *B, int sample_size, int result_cutoff, int distance_multiplier);

/*
 * Returns non-zero if and only if there exists a terminal thread or a virtual terminal terminal thread containing
 * (A, B), (B, A), (A, -B), and/or (-B, A) (and their mirrors) at event E.
 */
bool transmap_connectivityAndOrderWasPresentAtEvent(Event *E, Cap *A, Cap *B, int sample_size, int result_cutoff, int distance_multiplier);

/*
 * Returns non-zero if and only if there exists a terminal thread or a virtual terminal thread containing
 * (A, B) and/or (A, -B) (and their mirrors) at event E.
 */
bool transmap_connectivityAndOrientationWasPresentAtEvent(Event *E, Cap *A, Cap *B, int sample_size, int result_cutoff, int distance_multiplier);

/*
 * Returns non-zero if and only if there exists a terminal thread or a virtual terminal thread containing
 * (A, B) at event E.
 *
 * All the above functions can be created by logic on this function.
 */
bool transmap_connectivityOrderAndOrientationWasPresentAtEvent(Event *E, Cap *A, Cap *B, int sample_size, int result_cutoff, int distance_multiplier);

#endif
