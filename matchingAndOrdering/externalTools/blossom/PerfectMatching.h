/*
    PerfectMatching.h - interface to min cost perfect matching code

    Copyright 2008 Vladimir Kolmogorov (vnk@adastral.ucl.ac.uk)

    This software can be used for research and evaluation purposes only. Commercial use is prohibited.
    Public redistribution of the code or its derivatives is prohibited.
    If you use this software for research purposes, you should cite the following paper in any resulting publication:
        Vladimir Kolmogorov. "Blossom V: A new implementation of a minimum cost perfect matching algorithm."
        In Mathematical Programming Computation (MPC), July 2009, 1(1):43-67.

    For commercial use of the software not covered by this agreement, please contact the author.

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

#ifndef HALSKDJDFHALSJASFDFASJGLA
#define HALSKDJDFHALSJASFDFASJGLA

#include <assert.h>
#include <limits.h>
#include "block.h"


// if defined, edge costs are of type 'double', otherwise 'int'
//#define PERFECT_MATCHING_DOUBLE

// Note: with floating point numbers polynomial complexity is not guaranteed;
// the code may even get stuck due to rounding errors. If the code terminates,
// the solution may not be optimal. It may be worth calling CheckPerfectMatchingOptimality()
// to see whether complementary slackness conditions are satisfied.
//
// Using single precision floating point numbers (float) is really not recommended.


class PerfectMatching
{
public:

#ifdef PERFECT_MATCHING_DOUBLE
	typedef double REAL; 
	#define PM_INFTY ((REAL)1e100)
#else
	typedef int REAL;
	#define PM_INFTY (INT_MAX/2)
#endif

	typedef int NodeId;
	typedef int EdgeId;

	PerfectMatching(int nodeNum, int edgeNumMax);
	~PerfectMatching();

	// first call returns 0, second 1, and so on. 
	EdgeId AddEdge(NodeId i, NodeId j, REAL cost);

	// Computes a perfect matching of minimum cost. 
	// NOTE: a perfect matching of finite cost must exist (otherwise the behaviour is not specified).
	// If finish is false, then the final matching is not computed (so GetSolution() cannot be called afterwards).
	void Solve(bool finish=true);

	///////////////////////////////////////////////////////////////
	// Read primal solution (can be called after Solve()).
	int GetSolution(EdgeId e); // returns 1 if e is in the matching, 0 otherwise
	NodeId GetMatch(NodeId i); // alternative way to get the result

	///////////////////////////////////////////////////////////////
	// Read dual solution (can be called after Solve()).
	// 'blossom_parents' and 'twice_y' must be arrays of size node_num+GetBlossomNum().
	// The function sets blossom_parent[i] to the parent of i (or to -1 for exterior nodes).
	void GetDualSolution(int* blossom_parents, REAL* twice_y);
	int GetBlossomNum();

	///////////////////////////////////////////////////////////////
	// Dynamic graph updates. After calling Solve() you may call //
	// StartUpdate(), ..., FinishUpdate() and then Solve() again //
	///////////////////////////////////////////////////////////////
	void StartUpdate();
	void FinishUpdate();

	// 3 functions below can be called only between StartUpdate() and FinishUpdate().
	REAL GetTwiceSum(NodeId i); // if 2*cost(i,j)>=GetTwiceSum(i)+GetTwiceSum(j) then adding new edge (i,j) is not necessary - optimal solution will not change
	EdgeId AddNewEdge(NodeId i, NodeId j, REAL cost, bool do_not_add_if_positive_slack=true); // if do_not_add_if_positive_slack is true and the slack of the edge turns out to be non-negative, then the edge will not be added and -1 will be returned
	void UpdateCost(EdgeId e, REAL delta_cost);


	// NOTE: with default options the dual vector is guaranteed to be half-integral
	// (if all input weights are integral). However, with some options there is
	// no such guarantee, in particular, if dual_greedy_update_option=2 or dual_LP_threshold>0.
	// These options can be used only if the code is compiled with REAL=double.
	struct Options
	{
		Options() : fractional_jumpstart(true),
		            dual_greedy_update_option(0),
		            dual_LP_threshold(0.00),
		            update_duals_before(false),
		            update_duals_after(false),
		            single_tree_threshold(1.00),
		            verbose(true)
		{}

		bool	fractional_jumpstart; // false: greedy, true: compute fractional matching

		int 	dual_greedy_update_option; // 0: compute connected components (as in Blossom IV)
				                           // 1: compute strongly connected components (discussed by Cook-Rohe, but not implemented)
				                           // 2: single eps for all trees (fixed eps approach)

		double	dual_LP_threshold; // if tree_num => dual_updates_threshold*node_num: greedy updates
		                           // if tree_num <  dual_updates_threshold*node_num: global updates (solve LP)

		bool	update_duals_before; // before tree growth
		bool	update_duals_after;  // after tree growth

		double	single_tree_threshold; // if (tree_num => single_tree_threshold*node_num && update_duals_after): try to grow a single tree as long as possible

		bool	verbose;
	} options;


	// save problem to a file. format=0 corresponds to DIMACS format,
	// format=1 corresponds to the format used by blossom4.
	// CANNOT BE CALLED AFTER Solve()!!!
	void Save(char* filename, int format=0); 

	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
private:
	struct Node;
	struct Arc; // no such struct, only Arc* is used (the pointer can be odd or even)
	struct Edge; // pointer Edge* is always even
	struct Tree;
	struct TreeEdge;
	struct PQPointers;
	struct EdgeIterator;
	struct TreeEdgeIterator;
	struct LCATreeX;

	Node*	nodes;
	Edge*	edges;
	char*	edges_orig;
	DBlock<Node>* blossoms;
	Tree*	trees;
	DBlock<TreeEdge>* tree_edges;
	struct ExpandTmpItem
	{
		Node*	i;
		Node*	blossom_parent;
		Node*	blossom_grandparent;
	};
	Block<ExpandTmpItem>* expand_tmp_list; // used inside Expand()

	int		node_num;
	int		edge_num, edge_num_max;
	int		tree_num, tree_num_max;

	Node*	removed_first;
	int		blossom_num;
	int		removed_num;

	void*	pq_buf;

	bool	first_solve;

	// stat
	struct Stat
	{
		int		shrink_count;
		int		expand_count;
		int		grow_count;
		double	shrink_time;
		double	expand_time;
		double	dual_time;
	} stat;

	////////////////////////////////////////////////////////////////////

	void InitGreedy(bool allocate_trees=true);

	void InitGlobal(); // compute fractional matching
	Node* FindBlossomRootInit(Edge* a0);
	void ShrinkInit(Edge* a0, Node* tree_root);
	void ExpandInit(Node* b);
	void AugmentBranchInit(Node* i0, Node* tree_root);

	void Finish(); // sets matching for inner nodes

	void ProcessNegativeEdge(Edge* a);

	void GetRealEndpoints(Edge* a, Node*& tail, Node*& head);
	Node* FindBlossomRoot(Edge* a0);
	void Shrink(Edge* a0);
	void Expand(Node* b);
	void Augment(Edge* a0);
	void AugmentBranch(Node* i0);
	void GrowNode(Node* i);
	void GrowTree(Node* root, bool new_subtree);
	bool ProcessEdge00(Edge* a, bool update_boundary_edge=true); // returns true if boundary edge, false otherwise
	void ProcessSelfloop(Node* b, Edge* a);

	void AddTreeEdge(Tree* t0, Tree* t1);

	void ComputeEpsSingle(); // called from UpdateDuals()
	void ComputeEpsCC(); // called from UpdateDuals()
	void ComputeEpsSCC(); // called from UpdateDuals()
	void ComputeEpsGlobal(); // called from UpdateDuals()
	bool UpdateDuals();

	void FreeRemoved();
	void CommitEps();

	void ReallocateEdges();

	void PrintAll();
};


// in the functions below, 'edges' is an array of size 2*edge_num (edge e = (edges[2*e],edges[2*e+1])), 
//                         'weights' is an array of size edge_num

// checks complementary slackness conditions. 
// returns 0 if success.
// returns 1 if complementary slackness conditions are violated (then the amount of violation is printed - could potentially happen for double's)
// returns 2 if the blossom tree structure is incorrect (or inconsistent with primal solution)
int CheckPerfectMatchingOptimality(int node_num, int edge_num, int* edges, int* weights, PerfectMatching* pm, PerfectMatching::REAL threshold=(PerfectMatching::REAL)(1e-10));


double ComputePerfectMatchingCost(int node_num, int edge_num, int* edges, int* weights, PerfectMatching* pm);

#endif
