/*

LCA.h - Imlementation of LCA data structure described in

O. Berkman, U. Vishkin: Recursive Star-Tree Parallel Data Structure. SIAM J. Comput. 22(2): 221-242 (1993)

It takes O(n \log n) space (and O(n \log n) time for construction). The memory requirement can be
reduced to O(n) by defining LCA_BLOCKS. The implementation then follows the approach described in

J. Fischer, V. Heun: Theoretical and Practical Improvements on the RMQ-Problem, with Applications to LCA and LCE.
	Proceedings of the 17th Annual Symposium on Combinatorial Pattern Matching
	(CPM'06), Lecture Notes in Computer Science 4009, 36-48, Springer-Verlag, 2006.

(although in-block queries are computed naively in O(K) worst-case time, where K is the size of the block,
rather than as described by Fischer and Heun).

Written by Vladimir Kolmogorov, vnk@adastral.ucl.ac.uk.

//  example for 
//      4
//     / \
//    2   3
//   / \
//  0   1
// (Note: the ordering of nodes is in the tree preorder!!!)

	LCATree* lca = new LCATree(5);
	char A0, A1, A2, A3, A4;

	lca->Add(&A0, &A2);
	lca->Add(&A1, &A2);
	lca->Add(&A2, &A4);
	lca->Add(&A3, &A4);
	lca->AddRoot(&A4);

	int result = lca->GetLCA(1, 3); // should be 4
	delete lca;


*/

#ifndef GNAKDLATHJSTHAJSRNAKSJDA
#define GNAKDLATHJSTHAJSRNAKSJDA

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define LCA_BLOCKS

class LCATree
{
public:
	typedef void* NodeId; // can be any type, e.g. int. (The code checks NodeId's only for equalities.)
	typedef int PreorderId;

	LCATree(int node_num_max);
	~LCATree();

	// construct tree. Nodes must be added in the tree preorder!!!
	// First call returns 0, second returns 1, and so on.
	PreorderId Add(NodeId i, NodeId i_parent);
	PreorderId AddRoot(NodeId i); // completes tree construction

	PreorderId GetLCA(PreorderId i, PreorderId j);
	// Let i0=i, j0=j be the input nodes, and let r = LCA(i0,j0).
	// This function sets i and j to be the immediate children of r 
	// such that i is a descendant of i0 and j is a descendant of j0.
	// There must hold i0!=r, j0!=r.
	void GetPenultimateNodes(PreorderId& i, PreorderId& j);


//////////////////////////////////////////////////////////////////////////
private:
	int n, n_max, K, k_max;
	int** array;

	NodeId* buf0;
	int* buf1;
	NodeId* parent_current; 
	int* child_current;

	int* parents;

	int _GetLCA(int i, int j); // same as GetLCA, but assumes that i<j
	int GetLCADirect(int i, int j);
};

inline LCATree::LCATree(int _n_max) : n_max(_n_max), array(NULL)
{
#ifdef LCA_BLOCKS
	K = -2;
	n = n_max;
	while (n > 0) { K ++; n /= 2; }
	if (K < 1) K = 1;
#else
	K = 1;
	n = 0;
#endif
	parents = new int[n_max];
	buf0 = new NodeId[n_max];
	buf1 = new int[n_max];
	parent_current = buf0;
	child_current = buf1;
}

inline LCATree::~LCATree()
{
	int k;
	delete [] parents;
	if (buf0) delete [] buf0;
	if (buf1) delete [] buf1;
	if (array)
	{
		for (k=1; k<=k_max; k++) delete [] array[k];
		delete [] array;
	}
}

inline LCATree::PreorderId LCATree::Add(NodeId i, NodeId i_parent)
{
	assert(n < n_max);

	if (n == 0)
	{
		*parent_current = i;
		*(++ parent_current) = i_parent;
		parents[0] = -1;
	}
	else
	{
		if (i == *parent_current)
		{
			int c = *child_current --;
			while ( 1 )
			{
				int c_next = parents[c];
				parents[c] = n;
				if (c_next < 0) break;
				c = c_next;
			}
			parent_current --;
		}
		if (i_parent == *parent_current) parents[n] = *child_current;
		else
		{
			*(++ parent_current) = i_parent;
			parents[n] = -1;
			child_current ++;
		}
	}
	*child_current = n;
	return n ++;
}


inline LCATree::PreorderId LCATree::AddRoot(NodeId i)
{
	assert(n < n_max);

	if (n > 0)
	{
		if (i != *parent_current || parent_current != buf0+1)
		{
			printf("Error in LCATree construction: wrong sequence of calls!\n");
			exit(1);
		}
		int c = *child_current --;
		while ( 1 )
		{
			int c_next = parents[c];
			parents[c] = n;
			if (c_next < 0) break;
			c = c_next;
		}
		child_current ++;
	}
	parents[n++] = -1;

	delete [] buf0;
	buf0 = NULL;
	delete [] buf1;
	buf1 = NULL;

	// initialize array
	int b, k = 1, block_num = (n-1)/K+1;
	if (block_num < 3) return n-1;
	int d = (block_num-1)/4;
	while (d) { k ++; d >>= 1; }
	k_max = k;

	array = new int*[k_max+1];
	array[0] = parents;
	for (k=1, d=2; k<=k_max; k++, d*=2)
	{
		array[k] = new int[block_num-d];
		if (k == 1)
		{
			for (b=0; b<block_num-d; b++) array[1][b] = GetLCADirect((b+1)*K-1, (b+d)*K);
		}
		else
		{
			for (b=0; b<block_num-d; b++)
			{
				int i = array[k-1][b];
				int j = array[k-1][b+d/2];
				if (i < j) i = j;
				j = array[1][b+d/2-1];
				array[k][b] = (i > j) ? i : j;
			}
		}
	}

	return n-1;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

inline int LCATree::GetLCADirect(int i, int j)
{
	while (i < j) i = parents[i];
	return i;
}


inline int LCATree::_GetLCA(int i, int j)
{
#ifdef LCA_BLOCKS

	int bi = i/K, bj = j/K;
	if (bi == bj) return GetLCADirect(i, j);
	int i_last = (bi+1)*K-1, j_first = bj*K;
	i = GetLCADirect(i, i_last);
	j = GetLCADirect(j_first, j);
	if (i < j) i = j;
	// set j = LCA(i_last, j_first)
	if (j_first - i_last == 1) j = parents[i_last];
	else
	{
		int k = 1, d = (bj-bi)/4;
		while (d) { k ++; d >>= 1; }
		int diff = 1<<k;
		//assert(bi+diff <= bj && bi+diff>bj-diff);
		j = (array[k][bi] > array[k][bj-diff]) ? array[k][bi] : array[k][bj-diff];
	}
	return (i > j) ? i : j;

#else

	if (j == i) return i;

	int k = 0, d = (j-i)/2;
	while (d) { k ++; d >>= 1; }
	int diff = 1<<k;
	//assert(i+diff <= j && i+diff>j-diff);
	return (array[k][i] > array[k][j-diff]) ? array[k][i] : array[k][j-diff];

#endif
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

inline LCATree::PreorderId LCATree::GetLCA(PreorderId i, PreorderId j)
{
	if (i > j) { PreorderId k = i; i = j; j = k; }
	return _GetLCA(i, j);
}

inline void LCATree::GetPenultimateNodes(PreorderId& _i, PreorderId& _j)
{
	int i, j, d, swap;
	if (_i < _j) { i = _i; j = _j; swap = 0; }
	else         { i = _j; j = _i; swap = 1; }
	int r = _GetLCA(i, j);
	assert(i!=r && j!=r);
	while (parents[i] != r)
	{
		int i0 = parents[i];
		d = (j - i0)/2;
		while ( (i=_GetLCA(i0, i0+d)) == r ) d /= 2;
	}
	while (parents[j] != r)
	{
		int j0 = parents[j];
		d = (r - j0)/2;
		while ( (j=_GetLCA(j0, j0+d)) == r ) d /= 2;
	}
	if (swap == 0) { _i = i; _j = j; }
	else           { _j = i; _i = j; }
}

#endif
