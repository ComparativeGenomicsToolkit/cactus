/*
    MinCost.h - successive shortest path algorithm of Ford and Fulkerson for solving a minimum cost flow problem

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


#ifndef __MINCOST_H__
#define __MINCOST_H__

#include <string.h>
#include <assert.h>

// if GRAPH_ASSERT is defined then all calls to graph construction functions are assert'ed for correctness
// (e.g. that node_id's are valid id's and edge capacities are non-negative).
//#define GRAPH_ASSERT 

//#define MINCOST_DEBUG



template <typename FlowType, typename CostType> class MinCost
{
public:
	typedef int NodeId;
	typedef int EdgeId;

	MinCost(int NodeNum, int edgeNumMax, void (*err_function)(const char *) = NULL);

	// Destructor
	~MinCost();

	void AddNodeExcess(NodeId i, FlowType excess);

	// first call returns 0, second 1, and so on.
	// cap, rev_cap must be non-negative. 
	// cost can be negative.
	EdgeId AddEdge(NodeId i, NodeId j, FlowType cap, FlowType rev_cap, CostType cost);

	CostType Solve();

	///////////////////////////////////////////////////

	FlowType GetRCap(EdgeId e);
	void SetRCap(EdgeId e, FlowType new_rcap);
	FlowType GetReverseRCap(EdgeId e);
	void SetReverseRCap(EdgeId e, FlowType new_rcap);
	void PushFlow(EdgeId e, FlowType delta);
	void UpdateCost(EdgeId e, FlowType cap_orig, CostType delta);

	CostType GetDual(NodeId i) { return nodes[i].pi; }

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
	
protected:
	// internal variables and functions

	struct Node;
	struct Arc;

	struct Node
	{
		Arc			*firstNonsaturated;
		Arc			*firstSaturated;

		Arc			*parent;
		Node		*next; // list of nodes with positive excesses

		FlowType	excess;
		CostType	pi;
		int			flag;
		union
		{
			int		heap_ptr;
			Node*	next_permanent;
		};
#ifdef MINCOST_DEBUG
		int			id;
#endif
	};

	struct Arc
	{
		Node		*head;
		Arc			*prev;
		Arc			*next;
		Arc			*sister;	// reverse arc

		FlowType	r_cap;		// residual capacity
#ifdef MINCOST_DEBUG
		FlowType	cap_orig;
#endif
		CostType	cost;
		CostType GetRCost() { return cost + head->pi - sister->head->pi; }
	};

	int		nodeNum, edgeNum, edgeNumMax;
	Node	*nodes;
	Arc		*arcs;
	Node*	firstActive;
	int		counter;
	CostType cost;


	void	(*error_function)(const char *);	// this function is called if a error occurs,
										// with a corresponding error message
										// (or exit(1) is called if it's NULL)

	/////////////////////////////////////////////////////////////////////////

	struct PriorityQueue
	{
		PriorityQueue();
		~PriorityQueue();
		void Reset();
		CostType GetKey(Node* i);
		void Add(Node* i, CostType key);
		void DecreaseKey(Node* i, CostType key);
		Node* RemoveMin(CostType& key);

	private:
		struct Item
		{
			Node*		i;
			CostType	key;
		}* array;
		int N, arraySize;
		void Swap(int k1, int k2);
	};

	PriorityQueue queue;

	/////////////////////////////////////////////////////////////////////////

	void SetRCap(Arc* a, FlowType new_rcap);
	void PushFlow(Arc* a, FlowType delta);

	void Init();
	void DecreaseRCap(Arc* a, FlowType delta);
	void IncreaseRCap(Arc* a, FlowType delta);
	FlowType Augment(Node* start, Node* end);
	void Dijkstra(Node* start);

	void TestOptimality();
#ifdef MINCOST_DEBUG
	void TestCosts();
#endif
};


template <typename CostType> class DualMinCost : private MinCost<int, CostType>
{
public:
	typedef int NodeId;
	DualMinCost(int node_num, int constraint_num_max);
	~DualMinCost();

	void AddUnaryTerm(NodeId i, int objective_coef);
	void SetLowerBound(NodeId, CostType cmin);
	void SetUpperBound(NodeId, CostType cmax);
	void AddConstraint(NodeId i, NodeId j, CostType cmax); // xj - xi <= cmax

	void Solve();
	CostType GetSolution(NodeId i);

private:
	NodeId	source;
};









///////////////////////////////////////
// Implementation - inline functions //
///////////////////////////////////////



template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::AddNodeExcess(NodeId _i, FlowType excess)
{
	assert(_i>=0 && _i<nodeNum);
	nodes[_i].excess += excess;
	if (nodes[_i].excess > 0 && !nodes[_i].next)
	{
		nodes[_i].next = firstActive;
		firstActive = &nodes[_i];
	}
}

template <typename FlowType, typename CostType> 
	inline typename MinCost<FlowType, CostType>::EdgeId MinCost<FlowType, CostType>::AddEdge(NodeId _i, NodeId _j, FlowType cap, FlowType rev_cap, CostType cost)
{
	assert(_i>=0 && _i<nodeNum);
	assert(_j>=0 && _j<nodeNum);
	assert(_i!=_j && edgeNum<edgeNumMax);
	assert(cap >= 0);
	assert(rev_cap >= 0);

	Arc *a = &arcs[2*edgeNum];
	Arc *a_rev = a+1;
	edgeNum ++;

	Node* i = nodes + _i;
	Node* j = nodes + _j;

	a -> sister = a_rev;
	a_rev -> sister = a;
	if (cap > 0)
	{
		if (i->firstNonsaturated) i->firstNonsaturated->prev = a;
		a -> next = i -> firstNonsaturated;
		i -> firstNonsaturated = a;
	}
	else
	{
		if (i->firstSaturated) i->firstSaturated->prev = a;
		a -> next = i -> firstSaturated;
		i -> firstSaturated = a;
	}
	a->prev = NULL;
	if (rev_cap > 0)
	{
		if (j->firstNonsaturated) j->firstNonsaturated->prev = a_rev;
		a_rev -> next = j -> firstNonsaturated;
		j -> firstNonsaturated = a_rev;
	}
	else
	{
		if (j->firstSaturated) j->firstSaturated->prev = a_rev;
		a_rev -> next = j -> firstSaturated;
		j -> firstSaturated = a_rev;
	}
	a_rev->prev = NULL;

	a -> head = j;
	a_rev -> head = i;
	a -> r_cap = cap;
	a_rev -> r_cap = rev_cap;
	a -> cost = cost;
	a_rev -> cost = -cost;
#ifdef MINCOST_DEBUG
	a->cap_orig = cap;
	a_rev->cap_orig = rev_cap;
#endif

	if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
	if (a_rev->r_cap > 0 && a_rev->GetRCost() < 0) PushFlow(a_rev, a_rev->r_cap);

	return edgeNum-1;
}

///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::DecreaseRCap(Arc* a, FlowType delta)
{
	a->r_cap -= delta;
	if (a->r_cap == 0)
	{
		Node* i = a->sister->head;
		if (a->next) a->next->prev = a->prev;
		if (a->prev) a->prev->next = a->next;
		else         i->firstNonsaturated = a->next;
		a->next = i->firstSaturated;
		if (a->next) a->next->prev = a;
		a->prev = NULL;
		i->firstSaturated = a;
	}
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::IncreaseRCap(Arc* a, FlowType delta)
{
	if (a->r_cap == 0)
	{
		Node* i = a->sister->head;
		if (a->next) a->next->prev = a->prev;
		if (a->prev) a->prev->next = a->next;
		else         i->firstSaturated = a->next;
		a->next = i->firstNonsaturated;
		if (a->next) a->next->prev = a;
		a->prev = NULL;
		i->firstNonsaturated = a;
	}
	a->r_cap += delta;
}

template <typename FlowType, typename CostType> 
	inline FlowType MinCost<FlowType, CostType>::GetRCap(EdgeId e)
{
	Arc* a = &arcs[2*e];
	return a->r_cap;
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::SetRCap(Arc* a, FlowType new_rcap)
{
	assert(new_rcap >= 0);
#ifdef MINCOST_DEBUG
	a->cap_orig += new_rcap - a->r_cap;
#endif
	if (a->r_cap == 0)
	{
		Node* i = a->sister->head;
		if (a->next) a->next->prev = a->prev;
		if (a->prev) a->prev->next = a->next;
		else         i->firstSaturated = a->next;
		a->next = i->firstNonsaturated;
		if (a->next) a->next->prev = a;
		a->prev = NULL;
		i->firstNonsaturated = a;
	}
	a->r_cap = new_rcap;
	if (a->r_cap == 0)
	{
		Node* i = a->sister->head;
		if (a->next) a->next->prev = a->prev;
		if (a->prev) a->prev->next = a->next;
		else         i->firstNonsaturated = a->next;
		a->next = i->firstSaturated;
		if (a->next) a->next->prev = a;
		a->prev = NULL;
		i->firstSaturated = a;
	}
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::SetRCap(EdgeId e, FlowType new_rcap)
{
	SetRCap(&arcs[2*e], new_rcap);
}

template <typename FlowType, typename CostType> 
	inline FlowType MinCost<FlowType, CostType>::GetReverseRCap(EdgeId e)
{
	Arc* a = &arcs[2*e+1];
	return a->r_cap;
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::SetReverseRCap(EdgeId e, FlowType new_rcap)
{
	SetRCap(&arcs[2*e+1], new_rcap);
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PushFlow(Arc* a, FlowType delta)
{
	if (delta < 0) { a = a->sister; delta = -delta; }
	DecreaseRCap(a, delta);
	IncreaseRCap(a->sister, delta);
	a->head->excess += delta;
	a->sister->head->excess -= delta;
	cost += delta*a->cost;
	if (a->head->excess > 0 && !a->head->next)
	{
		a->head->next = firstActive;
		firstActive = a->head;
	}
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PushFlow(EdgeId e, FlowType delta)
{
	PushFlow(&arcs[2*e], delta);
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::UpdateCost(EdgeId e, FlowType cap_orig, CostType delta)
{
	Arc* a = &arcs[2*e];
	cost += delta*(cap_orig-a->r_cap);
	a->cost += delta;
	a->sister->cost = -a->cost;

	if (a->GetRCost() > 0) a = a->sister;
	if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
}

///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////

template <typename FlowType, typename CostType> 
	inline MinCost<FlowType, CostType>::PriorityQueue::PriorityQueue()
{
	N = 0;
	arraySize = 16;
	array = (Item*) malloc(arraySize*sizeof(Item));
}

template <typename FlowType, typename CostType> 
	inline MinCost<FlowType, CostType>::PriorityQueue::~PriorityQueue()
{
	free(array);
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PriorityQueue::Reset()
{
	N = 0;
}

template <typename FlowType, typename CostType> 
	inline CostType MinCost<FlowType, CostType>::PriorityQueue::GetKey(Node* i)
{
	return array[i->heap_ptr].key;
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PriorityQueue::Swap(int k1, int k2)
{
	Item* a = array+k1;
	Item* b = array+k2;
	a->i->heap_ptr = k2;
	b->i->heap_ptr = k1;
	Node* i = a->i;   a->i   = b->i;   b->i   = i;
	CostType key = a->key; a->key = b->key; b->key = key;
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PriorityQueue::Add(Node* i, CostType key)
{
	if (N == arraySize)
	{
		arraySize *= 2;
		array = (Item*) realloc(array, arraySize*sizeof(Item));
	}
	int k = i->heap_ptr = N ++;
	array[k].i = i;
	array[k].key = key;
	while (k > 0)
	{
		int k2 = (k-1)/2;
		if (array[k2].key <= array[k].key) break;
		Swap(k, k2);
		k = k2;
	}
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PriorityQueue::DecreaseKey(Node* i, CostType key)
{
	int k = i->heap_ptr;
	array[k].key = key;
	while (k > 0)
	{
		int k2 = (k-1)/2;
		if (array[k2].key <= array[k].key) break;
		Swap(k, k2);
		k = k2;
	}
}

template <typename FlowType, typename CostType> 
	inline typename MinCost<FlowType, CostType>::Node* MinCost<FlowType, CostType>::PriorityQueue::RemoveMin(CostType& key)
{
	if (N == 0) return NULL;

	Swap(0, N-1);
	N --;

	int k = 0;
	while ( 1 )
	{
		int k1 = 2*k + 1, k2 = k1 + 1;
		if (k1 >= N) break;
		int k_min = (k2 >= N || array[k1].key <= array[k2].key) ? k1 : k2;
		if (array[k].key <= array[k_min].key) break;
		Swap(k, k_min);
		k = k_min;
	}

	key = array[N].key;
	return array[N].i;
}


#endif
