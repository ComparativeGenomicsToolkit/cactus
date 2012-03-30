#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PMimplementation.h"


struct PerfectMatching::LCATreeX : LCATree
{
	LCATreeX(int size) : LCATree(size) { rev_mapping = new Node*[size]; }
	~LCATreeX() { delete [] rev_mapping; }
	Node**	rev_mapping;
};

void PerfectMatching::StartUpdate()
{
	Node* i0;
	Node* i;
	Node* j;
	Node* b;

	while ((i=removed_first))
	{
		removed_first = i->tree_sibling_next;
		blossoms->Delete(i);
		removed_num --;
	}

	Edge* a;
	Edge* selfloop_first = NULL;
	Edge* selfloop_last = NULL;

	for (i0=nodes; i0<nodes+node_num; i0++)
	{
		i0->is_processed = 0;
		if (i0->is_outer) continue;

		i0->is_tree_root = 0;
		i0->blossom_ptr = NULL;
		i = i0;
		while ( 1 )
		{
			j = i->blossom_parent;
			j->is_processed = 0;
			if (j->is_outer) { j->first_tree_child = i; break; }
			if (j->is_marked) break;
			if ((a=j->blossom_selfloops))
			{
				if (selfloop_last) selfloop_last->next[1] = a;
				else               selfloop_first         = a;
				selfloop_last = a;
				a->next[1] = NULL;
			}
			j->blossom_ptr = i;
			i = j;
		}
		b = (i->blossom_parent->is_outer) ? i->blossom_parent : i->blossom_parent->blossom_grandparent;
#ifdef LCA_REPAIRS
		if (!b->is_marked)
		{
			b->lca_size = 1;
			b->is_marked = 1;
		}
#endif
		while ( 1 )
		{
#ifdef LCA_REPAIRS
			b->lca_size ++;
#endif
			ARC_TO_EDGE_PTR(i->blossom_sibling)->y_saved = i->y;
			i->y += i->blossom_parent->y;
			if (!i->is_blossom) break;
			i->is_marked = 1;
			j = i;
			i = i->blossom_ptr;
			j->blossom_grandparent = b;
		}
		i->blossom_grandparent = b;
	}

#ifdef LCA_REPAIRS
	for (i0=nodes; i0<nodes+node_num; i0++)
	{
		if (i0->is_outer) continue;
		b = i0->blossom_grandparent;
		if (!b->is_marked) continue;
		b->is_marked = 0;
		LCATreeX* lca = new LCATreeX(b->lca_size);
		b->blossom_ptr = b->first_tree_child;
		i = b;
		while ( 1 )
		{
			if (i->blossom_ptr) i = i->blossom_ptr;
			else
			{
				while ( 1 )
				{
					if (i->is_outer) break;
					i->lca_preorder = lca->Add(i, i->blossom_parent);
					lca->rev_mapping[i->lca_preorder] = i;
					i = ARC_HEAD(i->blossom_sibling);
					if (i != i->blossom_parent->blossom_ptr) break;
					i = i->blossom_parent;
				}
				if (i->is_outer)
				{
					lca->AddRoot(i);
					break;
				}
			}
		}
		b->lca = lca;
	}
#endif

	while ((a=selfloop_first))
	{
		selfloop_first = a->next[1];
		do
		{
			Edge* a_next = a->next[0];

#ifdef LCA_REPAIRS
			int _i = a->head0[1]->lca_preorder;
			int _j = a->head0[0]->lca_preorder;
			Node* b = a->head0[1]->blossom_grandparent;
			b->lca->GetPenultimateNodes(_i, _j);
			i = b->lca->rev_mapping[_i];
			j = b->lca->rev_mapping[_j];
#else
			GetRealEndpoints(a, i, j);
#endif
			ADD_EDGE(i, a, 0);
			ADD_EDGE(j, a, 1);
			a->slack -= 2*i->blossom_eps;
			a = a_next;
		} while ( a );
	}

	/*
	for (i0=nodes; i0<nodes+node_num; i0++)
	{
		if (i0->is_outer) continue;
		b = i0->blossom_grandparent;
		if (b->lca)
		{
			delete b->lca;
			b->lca = NULL;
		}
	}
	*/

	nodes[node_num].first_tree_child = NULL;
}

void PerfectMatching::FinishUpdate()
{
	Node* i0;
	Node* i;
	Node* j;
	Edge* a;
	EdgeIterator I;
	int dir;
	Tree* t;

	for (i0=nodes; i0<nodes+node_num; i0++)
	{
		if (i0->is_outer) continue;

#ifdef LCA_REPAIRS
		if (i0->blossom_grandparent->lca)
		{
			delete i0->blossom_grandparent->lca;
			i0->blossom_grandparent->lca = NULL;
		}
#endif

		//////////////////////////////////////////////////////////////
		if (!i0->blossom_grandparent->is_removed)
		{
			i = i0;
			do
			{
				i->y = ARC_TO_EDGE_PTR(i->blossom_sibling)->y_saved;
				i->is_marked = 0;
				i->blossom_selfloops = NULL;
				i = i->blossom_parent;
			} while (i->is_marked);
			continue;
		}
		//////////////////////////////////////////////////////////////

		i = i0->blossom_parent;
		while ( 1 )
		{
			if (i->is_removed && !i->is_outer) break;
			REAL y_parent = (i->is_outer) ? 0 : i->blossom_parent->y;
			for (dir=0; dir<2; dir++)
			{
				if (!i->first[dir]) continue;
				i->first[dir]->prev[dir]->next[dir] = NULL;
				Edge* a_next;
				for (a=i->first[dir]; a; a=a_next)
				{
					a_next = a->next[dir];
					j = a->head0[1-dir];
					ADD_EDGE(j, a, dir);
					a->slack += j->blossom_parent->y - y_parent;
				}
				i->first[dir] = NULL;
			}
			if (i->is_removed) break;

			j = i->blossom_parent;
			i->is_removed = 1;
			i->tree_sibling_next = removed_first;
			removed_first = i;
			i = j;
		}
		i0->y = ARC_TO_EDGE_PTR(i0->blossom_sibling)->y_saved;
		i0->is_outer = 1;
		i0->flag = 2;
		i0->is_tree_root = 1;
	}

	Node* blossom_list = nodes[node_num].first_tree_child;



	for (i=nodes; i<nodes+node_num; i++)
	{
		if (!i->is_tree_root) continue;
		i->first_tree_child = nodes[node_num].first_tree_child;
		nodes[node_num].first_tree_child = i;
		REAL slack_min = PM_INFTY;
		FOR_ALL_EDGES(i, a, dir, I)
		{
			if (slack_min > a->slack) slack_min = a->slack;
		}
		i->y += slack_min;
		FOR_ALL_EDGES(i, a, dir, I) a->slack -= slack_min;
	}

	tree_num = 0;
	for (i=nodes[node_num].first_tree_child; i!=blossom_list; i=i->first_tree_child)
	{
		tree_num ++;
		if (!i->is_tree_root) continue;
		FOR_ALL_EDGES(i, a, dir, I)
		{
			j = a->head[dir];
			if (a->slack <= 0 && j->is_tree_root)
			{
				i->is_tree_root = j->is_tree_root = 0;
				i->match = EDGE_DIR_TO_ARC(a, dir);
				j->match = EDGE_DIR_TO_ARC(a, 1-dir);
				tree_num -= 2;
				break;
			}
		}
	}
	for ( ; i; i=i->first_tree_child)
	{
		if (i->is_removed) { i->is_tree_root = 0; continue; }
		tree_num ++;
	}

	if (tree_num > tree_num_max)
	{
		if (trees) free(trees);
		tree_num_max = tree_num;
		trees = (Tree*) malloc(tree_num_max*sizeof(Tree));
	}
	t = trees;

	Node* last_root = &nodes[node_num];
	Node* i_next;
	for (i=nodes; i; i=i_next)
	{
		if (!i->is_blossom) i_next = (i<nodes+node_num) ? (i + 1) : blossom_list;
		else                i_next = i->first_tree_child;
		if (!i->is_tree_root) continue;

		i->flag = 0;
		i->first_tree_child = NULL;
		i->tree_sibling_prev = last_root;
		last_root->tree_sibling_next = i;
		last_root = i;
		i->tree = t;
		t->root = i;
		t->eps = 0;
		t->first[0] = t->first[1] = NULL;
		t->pq_current = NULL;
		t->pq00.Reset();
		t->pq0.Reset();
		t->pq_blossoms.Reset();
		t ++;
	}

	assert(t == trees + tree_num);
	last_root->tree_sibling_next = NULL;

	while ((i=removed_first))
	{
		removed_first = i->tree_sibling_next;
		blossoms->Delete(i);
		blossom_num --;
	}
}

PerfectMatching::REAL PerfectMatching::GetTwiceSum(NodeId i)
{
	assert(i>=0 && i<node_num);
	return nodes[i].y;
}

inline void PerfectMatching::ProcessNegativeEdge(Edge* a)
{
	int dir;
	Node* i;
	for (dir=0; dir<2; dir++)
	{
		i = a->head0[dir];
		if (i->is_outer)
		{
			if (!i->is_tree_root)
			{
				i->is_tree_root = 1;
				i = ARC_HEAD(i->match);
				assert(!i->is_tree_root && i->is_outer);
				i->is_tree_root = 1;
				if (i->is_blossom)
				{
					i->first_tree_child = nodes[node_num].first_tree_child;
					nodes[node_num].first_tree_child = i;
				}
			}
            return; 
		}
		if (i->blossom_grandparent->is_removed) return;
	}

	Node* b = i->blossom_grandparent;
	assert(b->is_outer);

	if (!b->is_tree_root)
	{
		b->is_tree_root = 1;
		i = ARC_HEAD(b->match);
		assert(!i->is_tree_root && i->is_outer);
		i->is_tree_root = 1;
		if (i->is_blossom)
		{
			i->first_tree_child = nodes[node_num].first_tree_child;
			nodes[node_num].first_tree_child = i;
		}
	}

	b->is_removed = 1;
	b->tree_sibling_next = removed_first;
	removed_first = b;
}

PerfectMatching::EdgeId PerfectMatching::AddNewEdge(NodeId _i, NodeId _j, REAL cost, bool do_not_add_if_positive_slack)
{
	assert(_i>=0 && _i<node_num && _j>=0 && _j<node_num && _i!=_j);
	if (edge_num >= edge_num_max) ReallocateEdges();
	Node* i = nodes + _i;
	Node* j = nodes + _j;
	Edge* a = edges + edge_num;

	a->slack = cost*COST_FACTOR;
	a->head0[0] = j;
	a->head0[1] = i;
	Node* bi = (i->is_outer) ? i : i->blossom_grandparent;
	Node* bj = (j->is_outer) ? j : j->blossom_grandparent;
	if (bi == bj)
	{
#ifdef LCA_REPAIRS
		int _i = i->lca_preorder;
		int _j = j->lca_preorder;
		bi->lca->GetPenultimateNodes(_i, _j);
		i = bi->lca->rev_mapping[_i];
		j = bi->lca->rev_mapping[_j];
#else
		GetRealEndpoints(a, i, j);
#endif
		a->slack += i->blossom_parent->y + j->blossom_parent->y;
	}
	else
	{
		i = bi;
		j = bj;
	}
	a->slack -= a->head0[0]->y + a->head0[1]->y;

	if (a->slack >= 0 && do_not_add_if_positive_slack) return -1;

	ADD_EDGE(i, a, 0);
	ADD_EDGE(j, a, 1);
	PriorityQueue<REAL>::ResetItem(a);

	if (a->slack < 0)
	{
		ProcessNegativeEdge(a);
	}

	return edge_num ++;
}

void PerfectMatching::UpdateCost(EdgeId e, REAL delta_cost)
{
	assert(e>=0 && e<edge_num);
	Edge* a = edges + e;
	a->slack += delta_cost*COST_FACTOR;
	if (a->slack == 0) return;
	if (a->slack > 0)
	{
		Node* i = a->head[1];
		Node* j = a->head[0];
		if (i->is_outer)
		{
			if (ARC_TO_EDGE_PTR(i->match) != a && ARC_TO_EDGE_PTR(j->match) != a) return;
		}
		else
		{
			if (ARC_TO_EDGE_PTR(i->blossom_sibling) != a && ARC_TO_EDGE_PTR(j->blossom_sibling) != a) return;
		}
	}
	ProcessNegativeEdge(a);
}

