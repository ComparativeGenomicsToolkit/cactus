#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MinCost.h"


template <typename FlowType, typename CostType> 
	MinCost<FlowType, CostType>::MinCost(int _nodeNum, int _edgeNumMax, void (*err_function)(const char *))
	: nodeNum(_nodeNum),
	  edgeNum(0),
	  edgeNumMax(_edgeNumMax),
	  counter(0),
	  cost(0),
	  error_function(err_function)
{
	nodes = (Node*) malloc(nodeNum*sizeof(Node));
	arcs = (Arc*) malloc(2*edgeNumMax*sizeof(Arc));
	if (!nodes || !arcs) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	memset(nodes, 0, nodeNum*sizeof(Node));
	memset(arcs, 0, 2*edgeNumMax*sizeof(Arc));
	firstActive = &nodes[nodeNum];
#ifdef MINCOST_DEBUG
	for (int i=0; i<nodeNum; i++) nodes[i].id = i;
#endif
}

template <typename FlowType, typename CostType> 
	MinCost<FlowType, CostType>::~MinCost()
{
	free(nodes);
	free(arcs);
}

template <typename FlowType, typename CostType> 
	void MinCost<FlowType, CostType>::Init()
{
	Node* i;
	Arc* a;

	for (a=arcs; a<arcs+2*edgeNum; a++)
	{
		if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
	}

	Node** lastActivePtr = &firstActive;
	for (i=nodes; i<nodes+nodeNum; i++)
	{
		if (i->excess > 0)
		{
			*lastActivePtr = i;
			lastActivePtr = &i->next;
		}
		else i->next = NULL;
	}
	*lastActivePtr = &nodes[nodeNum];
}


template <typename FlowType, typename CostType> 
	FlowType MinCost<FlowType, CostType>::Augment(Node* start, Node* end)
{
	FlowType delta = (start->excess < -end->excess) ? start->excess : -end->excess;
	Arc* a;

	for (a=end->parent; a; a=a->sister->head->parent)
	{
		if (delta > a->r_cap) delta = a->r_cap;
	}
	assert(delta > 0);

	end->excess += delta;
	for (a=end->parent; a; a=a->head->parent)
	{
		DecreaseRCap(a, delta);
		a = a->sister;
		IncreaseRCap(a, delta);
	}
	start->excess -= delta;

	return delta;
}

template <typename FlowType, typename CostType> 
	void MinCost<FlowType, CostType>::Dijkstra(Node* start)
{
	assert(start->excess > 0);

	Node* i;
	Node* j;
	Arc* a;
	CostType d;
	Node* permanentNodes;

	int FLAG0 = ++ counter; // permanently labeled nodes
	int FLAG1 = ++ counter; // temporarily labeled nodes

	start->parent = NULL;
	start->flag = FLAG1;
	queue.Reset();
	queue.Add(start, 0);

	permanentNodes = NULL;

	while ( (i=queue.RemoveMin(d)) )
	{
		if (i->excess < 0)
		{
			FlowType delta = Augment(start, i);
			cost += delta*(d - i->pi + start->pi);
			for (i=permanentNodes; i; i=i->next_permanent) i->pi += d;
			break;
		}

		i->pi -= d;
		i->flag = FLAG0;
		i->next_permanent = permanentNodes;
		permanentNodes = i;

		for (a=i->firstNonsaturated; a; a=a->next)
		{
			j = a->head;
			if (j->flag == FLAG0) continue;
			d = a->GetRCost();
			if (j->flag == FLAG1)
			{
				if (d >= queue.GetKey(j)) continue;
				queue.DecreaseKey(j, d);
			}
			else
			{
				queue.Add(j, d);
				j->flag = FLAG1;
			}
			j->parent = a;
		}

	}
}


template <typename FlowType, typename CostType> 
	CostType MinCost<FlowType, CostType>::Solve()
{
	Node* i;
	//Init();
	while ( 1 )
	{
		i = firstActive;
		if (i == &nodes[nodeNum]) break;
		firstActive = i->next;
		i->next = NULL;
		if (i->excess > 0)
		{
			Dijkstra(i);
			if (i->excess > 0 && !i->next) 
			{ 
				i->next = firstActive; 
				firstActive = i; 
			}
		}
	}
#ifdef MINCOST_DEBUG
	TestOptimality();
	TestCosts();
#endif

	return cost;
}


template <typename FlowType, typename CostType> 
	void MinCost<FlowType, CostType>::TestOptimality()
{
	Node* i;
	Arc* a;

	for (i=nodes; i<nodes+nodeNum; i++)
	{
		if (i->excess != 0)
		{
			assert(0);
		}
		for (a=i->firstSaturated; a; a=a->next)
		{
			if (a->r_cap != 0)
			{
				assert(0);
			}
		}
		for (a=i->firstNonsaturated; a; a=a->next)
		{
			CostType c = a->GetRCost();
			if (a->r_cap <= 0 || a->GetRCost() < -1e-5)
			{
				assert(0);
			}
		}
	}
}

#ifdef MINCOST_DEBUG

template <typename FlowType, typename CostType> 
	void MinCost<FlowType, CostType>::TestCosts()
{
	Arc* a;

	CostType _cost = 0;

	for (a=arcs; a<arcs+2*edgeNum; a+=2)
	{
		assert(a->r_cap + a->sister->r_cap == a->cap_orig + a->sister->cap_orig);
		_cost += a->cost*(a->cap_orig - a->r_cap);
	}

	CostType delta = cost - _cost;
	if (delta < 0) delta = -delta;
	if (delta >= 1e-5)
	{
		assert(0);
	}
}

#endif




///////////////////////////////////////////////////////////////////////////////////////

#define FLOW_INFTY ((int)0x00fffffff)

template <typename CostType> 
	DualMinCost<CostType>::DualMinCost(int _nodeNum, int _edgeNumMax)
	: MinCost<int,CostType>(_nodeNum+1, _edgeNumMax+2*_nodeNum)
{
	source = _nodeNum;
}

template <typename CostType> 
	DualMinCost<CostType>::~DualMinCost()
{
}

template <typename CostType> 
	void DualMinCost<CostType>::AddUnaryTerm(NodeId i, int objective_coef)
{
	MinCost<int, CostType>::AddNodeExcess(i, objective_coef);
	MinCost<int, CostType>::AddNodeExcess(source, -objective_coef);
}

template <typename CostType> 
	void DualMinCost<CostType>::SetLowerBound(NodeId i, CostType cmin)
{
	AddEdge(i, source, FLOW_INFTY, 0, -cmin);
}

template <typename CostType> 
	void DualMinCost<CostType>::SetUpperBound(NodeId i, CostType cmax)
{
	AddEdge(source, i, FLOW_INFTY, 0, cmax);
}

template <typename CostType> 
	void DualMinCost<CostType>::AddConstraint(NodeId i, NodeId j, CostType cmax)
{
	AddEdge(i, j, FLOW_INFTY, 0, cmax);
}

template <typename CostType> 
	void DualMinCost<CostType>::Solve()
{
	MinCost<int,CostType>::Solve();
}

template <typename CostType> 
	CostType DualMinCost<CostType>::GetSolution(NodeId i)
{
	return MinCost<int, CostType>::nodes[source].pi - MinCost<int, CostType>::nodes[i].pi;
}



#include "instances.inc"


