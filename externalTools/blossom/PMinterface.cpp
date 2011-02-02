#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PMimplementation.h"


PerfectMatching::PerfectMatching(int nodeNum, int edgeNumMax)
	: node_num(nodeNum),
	  edge_num(0),
	  edge_num_max(edgeNumMax),
	  trees(NULL),
	  tree_num_max(0),
	  removed_first(NULL),
	  blossom_num(0),
	  removed_num(0),
	  first_solve(true)
{
	if (node_num & 1) { printf("# of nodes is odd: perfect matching cannot exist\n"); exit(1); }
	nodes = (Node*) malloc((node_num+1)*sizeof(Node));
	edges_orig = (char*) malloc(edge_num_max*sizeof(Edge)+1);
	edges = (Edge*) ( ( ((POINTER_TYPE)edges_orig) & 1 ) ? (edges_orig + 1) : edges_orig );
	memset(nodes, 0, (node_num+1)*sizeof(Node));

	blossoms   = new DBlock<Node>(256);
	tree_edges = new DBlock<TreeEdge>(256);
	expand_tmp_list = new Block<ExpandTmpItem>(256);
	pq_buf = PriorityQueue<REAL>::AllocateBuf();
}


void PerfectMatching::Save(char* filename, int format)
{
	if (!first_solve) { printf("Save() cannot be called after Solve()!\n"); exit(1); }
	int e;
	FILE* fp = fopen(filename, "w");
	if (!fp) { printf("Can't open %s\n", filename); exit(1); }
	if (format == 0)
	{
		fprintf(fp, "p edge %d %d\n", node_num, edge_num);
		for (e=0; e<edge_num; e++)
		{
			fprintf(fp, "e %d %d %d\n", 1+(int)(edges[e].head0[1]-nodes), 1+(int)(edges[e].head0[0]-nodes), (int)edges[e].slack/COST_FACTOR);
		}
	}
	else
	{
		fprintf(fp, "%d %d\n", node_num, edge_num);
		for (e=0; e<edge_num; e++)
		{
			fprintf(fp, "%d %d %d\n", (int)(edges[e].head0[1]-nodes), (int)(edges[e].head0[0]-nodes), (int)edges[e].slack/COST_FACTOR);
		}
	}
	fclose(fp);
}


PerfectMatching::~PerfectMatching()
{
	free(nodes);
	free(edges_orig);
	delete blossoms;
	delete tree_edges;
	delete expand_tmp_list;
	if (trees) free(trees);
	PriorityQueue<REAL>::DeallocateBuf(pq_buf);
}


PerfectMatching::EdgeId PerfectMatching::AddEdge(NodeId _i, NodeId _j, REAL cost)
{
	if (_i<0 || _i>=node_num || _j<0 || _j>node_num || _i==_j)
	{
		printf("wrong node id's! (%d,%d)\n", _i, _j); exit(1);
	}
	if (edge_num >= edge_num_max) ReallocateEdges();
	Node* i = nodes + _i;
	Node* j = nodes + _j;
	Edge* a = edges + edge_num;

	ADD_EDGE(i, a, 0);
	ADD_EDGE(j, a, 1);
	a->head0[0] = j;
	a->head0[1] = i;

	a->slack = cost*COST_FACTOR;
	PriorityQueue<REAL>::ResetItem(a);

	return edge_num ++;
}

int PerfectMatching::GetSolution(EdgeId e)
{
	assert(e>=0 && e<edge_num);
	Edge* a = edges + e;
	return (a->head0[1]->match == EDGE_DIR_TO_ARC(a, 0)) ? 1 : 0;
}

PerfectMatching::NodeId PerfectMatching::GetMatch(NodeId i)
{
	assert(i>=0 && i<node_num);
	return (int)(ARC_HEAD0(nodes[i].match)-nodes);
}

void PerfectMatching::GetRealEndpoints(Edge* a, Node*& tail, Node*& head)
{
	Node* i;
	Node* j;
	int delta = 0;

	for (i=a->head0[1]; !i->is_outer; i=i->blossom_parent, delta--) {}
	for (j=a->head0[0]; !j->is_outer; j=j->blossom_parent, delta++) {}
	if ( i == j )
	{
		i = a->head0[1];
		j = a->head0[0];
		while ( delta < 0 ) { i = i->blossom_parent; delta ++; }
		while ( delta > 0 ) { j = j->blossom_parent; delta --; }
		while ( i->blossom_parent != j->blossom_parent )
		{
			i = i->blossom_parent;
			j = j->blossom_parent; 
		}
	}
	tail = i;
	head = j;
	assert((i->is_outer && j->is_outer) || (i->blossom_parent==j->blossom_parent && !i->is_outer && !j->is_outer));
}

void PerfectMatching::ReallocateEdges()
{
	printf("Warning: reallocating edges. Increasing edge_num_max in the constructor may improve memory efficiency!\n");
	edge_num_max = edge_num_max*3/2 + 16;
	char* edges_orig_old = edges_orig;
	Edge* edges_old = edges;
	edges_orig = (char*) realloc(edges_orig_old, edge_num_max*sizeof(Edge)+1);
	edges = (Edge*) ( ( ((POINTER_TYPE)edges_orig_old) & 1 ) ? (edges_orig + 1) : edges_orig );
	if ( ((POINTER_TYPE)edges) & 1 )
	{
		char* edges_orig_old2 = edges_orig;
		Edge* edges_old2 = edges;

		edges_orig = (char*) malloc(edge_num_max*sizeof(Edge)+1);
		edges = (Edge*) ( ( ((POINTER_TYPE)edges_orig_old) & 1 ) ? (edges_orig + 1) : edges_orig );
		memcpy(edges, edges_old2, edge_num*sizeof(Edge));
		free(edges_orig_old2);
	}

#define UPDATE_EDGE_PTR(ptr) ptr = (Edge*)((char*)(ptr) + ((char*)edges - (char*)edges_old))
#define UPDATE_ARC_PTR(ptr) ptr = (Arc*)((char*)(ptr) + ((char*)edges - (char*)edges_old))

	Node* i;
	Edge* a;
	for (a=edges; a<edges+edge_num; a++)
	{
		if (a->next[0]) UPDATE_EDGE_PTR(a->next[0]);
		if (a->next[1]) UPDATE_EDGE_PTR(a->next[1]);
		if (a->prev[0]) UPDATE_EDGE_PTR(a->prev[0]);
		if (a->prev[1]) UPDATE_EDGE_PTR(a->prev[1]);
	}
	if (first_solve)
	{
		for (i=nodes; i<nodes+node_num; i++)
		{
			if (i->first[0]) UPDATE_EDGE_PTR(i->first[0]);
			if (i->first[1]) UPDATE_EDGE_PTR(i->first[1]);
		}
	}
	else
	{
		Node* i0;
		for (i0=nodes; i0<nodes+node_num; i0++)
		{
			i = i0;
			while ( 1 )
			{
				if (i->is_outer)
				{
					UPDATE_ARC_PTR(i->match);
					if (i->first[0]) UPDATE_EDGE_PTR(i->first[0]);
					if (i->first[1]) UPDATE_EDGE_PTR(i->first[1]);
					break;
				}
				UPDATE_ARC_PTR(i->blossom_sibling);
				if (i->first[0]) UPDATE_EDGE_PTR(i->first[0]);
				if (i->first[1]) UPDATE_EDGE_PTR(i->first[1]);

				i = i->blossom_parent;
				if (i->is_outer) { if ( i->is_marked) break; i->is_marked = 1; }
				else             { if (!i->is_marked) break; i->is_marked = 0; }
			}
		}
		for (i0=nodes; i0<nodes+node_num; i0++)
		{
			i = i0;
			while ( 1 )
			{
				if (i->is_outer) break;

				i = i->blossom_parent;
				if (i->is_outer) { if (!i->is_marked) break; i->is_marked = 0; }
				else             { if ( i->is_marked) break; i->is_marked = 1; }
			}
		}
	}
}

int PerfectMatching::GetBlossomNum()
{
	return blossom_num;
}

void PerfectMatching::GetDualSolution(int* blossom_parents, REAL* twice_y)
{
	int _i0, id = node_num;
	int* child_ptr;
	Node* i0;
	Node* i;
	int* tmp_array = new int[blossom_num];

	int* tmp_array_ptr = tmp_array;
	for (_i0=0, i0=nodes; _i0<node_num; _i0++, i0++)
	{
		twice_y[_i0] = i0->y;
		if (i0->is_outer)
		{
			blossom_parents[_i0] = -1;
			continue;
		}
		child_ptr = &blossom_parents[_i0];
		i = i0->blossom_parent;
		while ( 1 )
		{
			if (i->is_marked)
			{
				*child_ptr = i->lca_preorder;
				break;
			}
			i->is_marked = 1;
			*tmp_array_ptr ++ = i->lca_preorder;
			*child_ptr = i->lca_preorder = id ++;
			child_ptr = &blossom_parents[i->lca_preorder];
			twice_y[i->lca_preorder] = i->y;
			if (i->is_outer)
			{
				*child_ptr = -1;
				break;
			}
			i = i->blossom_parent;
		}
	}

	assert(id == node_num+blossom_num && tmp_array_ptr == tmp_array + blossom_num);

	tmp_array_ptr = tmp_array;
	for (_i0=0, i0=nodes; _i0<node_num; _i0++, i0++)
	{
		if (i0->is_outer) continue;
		i = i0->blossom_parent;
		while ( 1 )
		{
			if (!i->is_marked) break;
			i->is_marked = 0;
			i->lca_preorder = *tmp_array_ptr ++;
			if (i->is_outer) break;
			i = i->blossom_parent;
		}
	}

	delete [] tmp_array;
}
