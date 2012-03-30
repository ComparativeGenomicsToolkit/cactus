/*
    PMimplementation.h

    Copyright 2008 Vladimir Kolmogorov (vnk@adastral.ucl.ac.uk)

    This software can be used for research purposes only. Commercial use is prohibited.
    Public redistribution of the code or its derivatives is prohibited.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef ASKHAKJSNTJAKSNBBAVASRA
#define ASKHAKJSNTJAKSNBBAVASRA

#ifdef _MSC_VER
#pragma warning(disable: 4311)
#pragma warning(disable: 4312)
#endif

#include "PerfectMatching.h"
#include "timer.h"
#include "PQ.h"
#include "LCA.h"


#define LCA_REPAIRS




#define IS_INT ( ((REAL)1 / 2) == 0 )
#define COST_FACTOR 2

#define PM_THRESHOLD ((REAL)1e-12)

struct PerfectMatching::Node
{
	unsigned int	is_outer : 1; // 0 - the node is contained in another blossom, 1 - it's the outermost node
	unsigned int	flag : 2; // 0 corresponds to +, 1 corresponds to -, 2 corresponds to a free node
	unsigned int	is_tree_root : 1;
	unsigned int	is_processed : 1;
	unsigned int	is_blossom : 1;
	unsigned int	is_marked : 1;
	unsigned int	is_removed : 1;

	Edge*		first[2];
	union
	{
		Arc*	match; // used if not a tree root (is_tree_root = 0) or it's an inner node (is_outer == 0)
		Node*	blossom_grandparent;
	};
	REAL		y;

	union
	{
		struct // used when is_outer = 0
		{
			Arc*	blossom_sibling;
			Node*	blossom_parent;
			union
			{
				Edge*	blossom_selfloops;
				Node*	blossom_ptr; // used in repairs
#ifdef LCA_REPAIRS
				int		lca_preorder; // used in repairs
#endif
			};
			REAL blossom_eps; // stores 'eps' of the tree at the moment when the node was shrunk into a blossom.
			                  // (it is used for determining slacks of self-loops)
		};
		struct // used when is_outer = 1
		{
			union
			{
				struct // used for "+" nodes (flag = 0)
				{
					Node*	first_tree_child;
					Node*	tree_sibling_prev; // circular list (with one exception: parent->first_tree_child->tree_sibling_prev->tree_sibling_next is NULL)
					Node*	tree_sibling_next;
				};
				Arc*	tree_parent; // used for "-" nodes (flag = 1)
			};
			union
			{
				Tree*	tree;
				Edge*	best_edge;  // used during InitGlobal() for non-tree nodes
#ifdef LCA_REPAIRS
				int		lca_size; // used in repairs
				LCATreeX*	lca;  // used in repairs
#endif
			};
		};
	};
};

struct PerfectMatching::Edge : PriorityQueue<REAL>::Item
{
	Node*	head[2];
	Node*	head0[2];
	Edge*	next[2];
	Edge*	prev[2];
};

typedef unsigned long POINTER_TYPE;
// if the declaration below fails, set POINTER_TYPE to be the appropriate integer type of the same length as (void*)
extern char dummy_array[2*(sizeof(void*)==sizeof(POINTER_TYPE))-1];

#define ARC_TO_EDGE_PTR(a)       ( (Edge*) ( ((POINTER_TYPE)(a)) & (~1)      ) )
#define ARC_TO_EDGE_DIR(a)       ( (int)   ( ((POINTER_TYPE)(a)) & 1         ) )
#define EDGE_DIR_TO_ARC(a, dir)  ( (Arc*)  ( (char*)(a) + (dir)) )

#define ARC_REV(a) ( (Arc*) ( ((POINTER_TYPE)(a)) ^ 1 ) )

#define ARC_TAIL(a)  (ARC_TO_EDGE_PTR(a)->head [1-ARC_TO_EDGE_DIR(a)])
#define ARC_TAIL0(a) (ARC_TO_EDGE_PTR(a)->head0[1-ARC_TO_EDGE_DIR(a)])
#define ARC_HEAD(a)  (ARC_TO_EDGE_PTR(a)->head [ARC_TO_EDGE_DIR(a)])
#define ARC_HEAD0(a) (ARC_TO_EDGE_PTR(a)->head0[ARC_TO_EDGE_DIR(a)])

struct PerfectMatching::PQPointers
{
	PriorityQueue<REAL> pq00; // plus-plus edges
	union
	{
		PriorityQueue<REAL> pq01[2]; // plus-minus, minus-plus edges. Used for tree edges.
		struct // used for trees.
		{
			PriorityQueue<REAL> pq0; // plus-free edges
			PriorityQueue<REAL> pq_blossoms;
		};
	};
};

struct PerfectMatching::Tree : PQPointers
{
	REAL		eps;
	TreeEdge*	first[2];
	Node*		root;

	PQPointers*	pq_current;
	int			dir_current;

	/////////////////////////////////////////
	// used while computing dual updates
	REAL		eps_delta;
	Tree*		next;
	union
	{
		int			id;
		TreeEdge*	dfs_parent;
	};
};

struct PerfectMatching::TreeEdge : PQPointers
{
	Tree*		head[2];
	TreeEdge*	next[2];
};



///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#define GET_PENULTIMATE_BLOSSOM(j)\
	{\
		Node* jtmp1 = j;\
		while ( 1 )\
		{\
			if (!j->blossom_grandparent->is_outer) j = j->blossom_grandparent;\
			else if (j->blossom_grandparent != j->blossom_parent) j->blossom_grandparent = j->blossom_parent;\
			else break;\
		}\
		Node* jtmp2;\
		for ( ; jtmp1!=j; jtmp1=jtmp2)\
		{\
			jtmp2 = jtmp1->blossom_grandparent;\
			jtmp1->blossom_grandparent = j;\
		}\
	}
#define GET_PENULTIMATE_BLOSSOM2(j)\
	{\
		Node* jtmp1 = j;\
		Node* jtmp_prev = NULL;\
		while ( 1 )\
		{\
			if (!j->blossom_grandparent->is_outer) { jtmp_prev = j; j = j->blossom_grandparent; }\
			else if (j->blossom_grandparent != j->blossom_parent) j->blossom_grandparent = j->blossom_parent;\
			else break;\
		}\
		if (jtmp_prev)\
		{\
			Node* jtmp2;\
			for ( ; jtmp1!=jtmp_prev; jtmp1=jtmp2)\
			{\
				jtmp2 = jtmp1->blossom_grandparent;\
				jtmp1->blossom_grandparent = jtmp_prev;\
			}\
		}\
	}


///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


struct PerfectMatching::EdgeIterator
{
	Edge*	a_last;
	int		start_flag;
};

#define FOR_ALL_EDGES(i, a, dir, I)\
	for ( dir = (i->first[0]) ? 0 : 1, I.a_last = a = i->first[dir], I.start_flag = (a) ? 0 : 1;\
		a != I.a_last || (I.start_flag ++ == 0) || (dir ++ == 0 && (I.a_last = a = i->first[1]));\
		a = a->next[dir] )

#define CONTINUE_FOR_ALL_EDGES(i, a, dir, I)\
	for ( a = a->next[dir];\
		a != I.a_last || (I.start_flag ++ == 0) || (dir ++ == 0 && (I.a_last = a = i->first[1]));\
		a = a->next[dir] )

#define REMOVE_EDGE(i, a, dir)\
	{\
		if ((a)->prev[dir]==(a)) (i)->first[dir] = NULL;\
		else\
		{\
			(a)->prev[dir]->next[dir] = (a)->next[dir];\
			(a)->next[dir]->prev[dir] = (a)->prev[dir];\
			(i)->first[dir] = (a)->next[dir];\
		}\
	}

#define ADD_EDGE(i, a, dir)\
	{\
		if ((i)->first[dir])\
		{\
			(a)->prev[dir] = (i)->first[dir]->prev[dir];\
			(a)->next[dir] = (i)->first[dir];\
			(i)->first[dir]->prev[dir]->next[dir] = (a);\
			(i)->first[dir]->prev[dir] = (a);\
		}\
		else (i)->first[dir] = (a)->prev[dir] = (a)->next[dir] = (a);\
		(a)->head[1-(dir)] = (i);\
	}

#define MOVE_EDGE(i_old, i_new, a, dir)\
	{\
		REMOVE_EDGE(i_old, a, dir);\
		ADD_EDGE(i_new, a, dir);\
	}

#define GET_OUTER_HEAD(a, dir, j)\
	{\
		j = (a)->head[dir];\
		if (!j->is_outer)\
		{\
			Node* j_orig = j;\
			GET_PENULTIMATE_BLOSSOM(j);\
			j = j->blossom_parent;\
			int dir_rev = 1 - (dir);\
			MOVE_EDGE(j_orig, j, a, dir_rev);\
		}\
	}

#define GET_TREE_PARENT(child, parent)\
	{\
		Arc* a = (child)->tree_parent;\
		Edge* e = ARC_TO_EDGE_PTR(a);\
		int dir = ARC_TO_EDGE_DIR(a);\
		GET_OUTER_HEAD(e, dir, parent);\
	}

///////////////////////////////////////////////////////////////////////////////////////////////

struct PerfectMatching::TreeEdgeIterator
{
	TreeEdge**	e_ptr;
};

#define FOR_ALL_TREE_EDGES(t, e, dir)\
	for ( dir = (t->first[0]) ? 0 : 1, e = t->first[dir];\
	      e || (dir ++ == 0 && (e = t->first[1]));\
	      e = e->next[dir] )

#define FOR_ALL_TREE_EDGES_X(t, e, dir, T)\
	for ( dir = (t->first[0]) ? 0 : 1, T.e_ptr = &t->first[dir], e = *T.e_ptr;\
	      e || (dir ++ == 0 && (e = *(T.e_ptr = &t->first[1])));\
	      e = *T.e_ptr )\
	if (e->head[dir] == NULL) { *T.e_ptr = e->next[dir]; tree_edges->Delete(e); }\
	else if ((T.e_ptr = &e->next[dir]))

///////////////////////////////////////////////////////////////////////////////////////////////

#define MOVE_NODE_IN_TREE(i)\
	{\
		if ((i)->first_tree_child) (i) = (i)->first_tree_child;\
		else\
		{\
			while (!(i)->is_tree_root && !(i)->tree_sibling_next) { (i) = ARC_HEAD((i)->match); GET_TREE_PARENT(i, i); }\
			if ((i)->is_tree_root) break;\
			(i) = (i)->tree_sibling_next;\
		}\
	}

// i=parent, j=child
#define ADD_TREE_CHILD(i, j)\
	{\
		(j)->flag = 0;\
		(j)->tree = (i)->tree;\
		(j)->first_tree_child = NULL;\
		(j)->tree_sibling_next = (i)->first_tree_child;\
		if ((i)->first_tree_child)\
		{\
			(j)->tree_sibling_prev = (i)->first_tree_child->tree_sibling_prev;\
			(i)->first_tree_child->tree_sibling_prev = j;\
		}\
		else\
		{\
			(j)->tree_sibling_prev = j;\
		}\
		(i)->first_tree_child = j;\
	}


#define REMOVE_FROM_TREE(i)\
	{\
		if ((i)->tree_sibling_next) (i)->tree_sibling_next->tree_sibling_prev = (i)->tree_sibling_prev;\
		else\
		{\
			Node* i_NEXT = ARC_HEAD((i)->match); i_NEXT = ARC_HEAD(i_NEXT->tree_parent); i_NEXT = i_NEXT->first_tree_child;\
			i_NEXT->tree_sibling_prev = (i)->tree_sibling_prev;\
		}\
		if ((i)->tree_sibling_prev->tree_sibling_next) (i)->tree_sibling_prev->tree_sibling_next = (i)->tree_sibling_next;\
		else\
		{\
			Node* i_PARENT = ARC_HEAD((i)->match); i_PARENT = ARC_HEAD(i_PARENT->tree_parent);\
			i_PARENT->first_tree_child = (i)->tree_sibling_next;\
		}\
	}


#endif

