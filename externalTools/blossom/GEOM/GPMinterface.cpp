#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "GeomPerfectMatching.h"
#include "GPMkdtree.h"


GeomPerfectMatching::GeomPerfectMatching(int nodeNum, int _DIM)
	: DIM(_DIM),
	  node_num(0),
	  node_num_max(nodeNum),
	  edge_num(0)
{
	if (node_num_max < 1) { printf("too few nodes\n"); exit(1); }
	if (node_num_max & 1) { printf("# of points is odd: perfect matching cannot exist\n"); exit(1); }
	nodes = (Node*) malloc(node_num_max*sizeof(Node));
	memset(nodes, 0, node_num_max*sizeof(Node));
	edges = new Block<Edge>(512);
	coords = (REAL*) malloc((DIM+1)*node_num_max*sizeof(REAL));
	sums = coords + DIM*node_num_max;
	matching = (int*) malloc(node_num_max*sizeof(int));
	int i;
	for (i=0; i<node_num_max; i++) matching[i] = i;
}

GeomPerfectMatching::~GeomPerfectMatching()
{
	free(nodes);
	delete edges;
	free(coords);
	free(matching);
}

GeomPerfectMatching::PointId GeomPerfectMatching::AddPoint(REAL* coord)
{
	if (node_num >= node_num_max)
	{
		printf("Error: you are trying to add too many points!\n");
		exit(1);
	}
	memcpy(coords+DIM*node_num, coord, DIM*sizeof(REAL));
	return node_num ++;
}

void GeomPerfectMatching::AddInitialEdge(PointId _i, PointId _j)
{
	assert(_i>=0 && _i<node_num_max && _j>=0 && _j<node_num_max && _i!=_j);
	if (_j < _i) { int _k = _i; _i = _j; _j = _k; }
	Node* i = nodes + _i;
	Node* j = nodes + _j;
	Edge* e = edges->New();
	edge_num ++;

	e->head[1] = _i;
	e->head[0] = _j;
	e->next[0] = i->first[0];
	e->next[1] = j->first[1];
	i->first[0] = e;
	j->first[1] = e;
}


GeomPerfectMatching::REAL GeomPerfectMatching::ComputeCost(PointId* matching)
{
	if (node_num != node_num_max) { printf("ComputeCost() cannot be called before all points have been added!\n"); exit(1); }

	REAL cost = 0;
	int i;
	for (i=0; i<node_num; i++)
	{
		if (matching[i]==i || matching[i]<0 || matching[i]>=node_num || matching[matching[i]]!=i)
		{
			printf("ComputeCost(): not a valid matching!\n");
			exit(1);
		}
		if (matching[i] > i)
		{
			cost += Dist(i, matching[i]);
		}
	}
	return cost;
}


