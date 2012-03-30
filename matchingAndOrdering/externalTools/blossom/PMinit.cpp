#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PMimplementation.h"


void PerfectMatching::InitGreedy(bool allocate_trees)
{
	Node* i;
	int dir;
	Edge* a;
	EdgeIterator I;
	Tree* t = NULL;
	Node* last_root = &nodes[node_num];
	REAL slack_min;

	for (i=nodes; i<nodes+node_num; i++) i->y = PM_INFTY;
	for (a=edges; a<edges+edge_num; a++)
	{
		if (a->head[0]->y > a->slack) a->head[0]->y = a->slack;
		if (a->head[1]->y > a->slack) a->head[1]->y = a->slack;
	}
	for (a=edges; a<edges+edge_num; a++)
	{
		i = a->head[0];
		if (!i->is_outer)
		{
			i->is_outer = 1;
			i->y /= 2;
		}
		a->slack -= i->y;
		i = a->head[1];
		if (!i->is_outer)
		{
			i->is_outer = 1;
			i->y /= 2;
		}
		a->slack -= i->y;
	}

	tree_num = node_num;
	for (i=nodes; i<nodes+node_num; i++)
	{
		if (i->flag == 2) continue;
		slack_min = PM_INFTY;
		FOR_ALL_EDGES(i, a, dir, I) if (slack_min > a->slack) slack_min = a->slack;
		i->y += slack_min;
		FOR_ALL_EDGES(i, a, dir, I)
		{
			if (a->slack <= slack_min && i->flag == 0 && a->head[dir]->flag == 0)
			{
				i->flag = 2;
				a->head[dir]->flag = 2;
				i->match = EDGE_DIR_TO_ARC(a, dir);
				a->head[dir]->match = EDGE_DIR_TO_ARC(a, 1-dir);
				tree_num -= 2;
			}
			a->slack -= slack_min;
		}
	}
	if (allocate_trees)
	{
		if (tree_num > tree_num_max)
		{
			if (trees) free(trees);
			tree_num_max = tree_num;
			trees = (Tree*) malloc(tree_num_max*sizeof(Tree));
		}
		t = trees;
	}
	for (i=nodes; i<nodes+node_num; i++)
	{
		if (i->flag != 0) continue;
		i->is_tree_root = 1;
		i->first_tree_child = NULL;
		i->tree_sibling_prev = last_root;
		last_root->tree_sibling_next = i;
		last_root = i;
		if (allocate_trees)
		{
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
	}
	last_root->tree_sibling_next = NULL;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

PerfectMatching::Node* PerfectMatching::FindBlossomRootInit(Edge* a0)
{
	Node* i;
	Node* j;
	Node* _i[2];
	Node* r;
	int branch;

	_i[0] = ARC_HEAD(a0);
	_i[1] = ARC_TAIL(a0);
	branch = 0;
	while ( 1 )
	{
		if (!_i[branch]->is_outer) 
		{
			r = _i[branch]; 
			j = _i[1-branch];
			break; 
		}
		_i[branch]->is_outer = 0;
		if (_i[branch]->is_tree_root)
		{
			j = _i[branch];
			i = _i[1-branch];
			while (i->is_outer)
			{
				i->is_outer = 0;
				i = ARC_HEAD(i->match);
				i->is_outer = 0;
				i = ARC_HEAD(i->tree_parent);
			}
			r = i;
			break;
		}
		i = ARC_HEAD(_i[branch]->match);
		i->is_outer = 0;
		_i[branch] = ARC_HEAD(i->tree_parent);
		branch = 1 - branch;
	}
	i = r;
	while ( i != j )
	{
		i = ARC_HEAD(i->match);
		i->is_outer = 1;
		i = ARC_HEAD(i->tree_parent);
		i->is_outer = 1;
	}
	return r;
}

void PerfectMatching::ShrinkInit(Edge* a0, Node* tree_root)
{
	int branch, flag;
	Node* i;
	Node* j;
	Node* r;
	Arc* a_prev;
	Arc* aa;

	tree_root->flag = 2;
	i = tree_root->first_tree_child;
	if ( i )
	while ( 1 )
	{
		ARC_HEAD(i->match)->flag = 2;
		i->flag = 2;

		MOVE_NODE_IN_TREE(i);
	}

	r = FindBlossomRootInit(a0);

	if ( !r->is_tree_root )
	{
		j = ARC_HEAD(r->match);
		j->match = aa = j->tree_parent;
		i = ARC_HEAD(aa);
		while ( !i->is_tree_root )
		{
			j = ARC_HEAD(i->match);
			i->match = ARC_REV(aa);
			j->match = aa = j->tree_parent;
			i = ARC_HEAD(aa);
		}
		i->match = ARC_REV(aa);
	}

	tree_root->is_tree_root = 0;

	branch = 0;
	flag = 0;
	a_prev = EDGE_DIR_TO_ARC(a0, 0);
	i = ARC_HEAD(a_prev);
	while ( 1 )
	{
		Arc* a_next = (flag == 0) ? i->match : i->tree_parent;
		flag = 1 - flag;
		i->flag = 0;
		i->match = NULL;
		if (branch == 0)
		{
			i->blossom_sibling = a_next;
			if (i == r)
			{
				branch = 1;
				flag = 0;
				a_prev = ARC_REV(a0);
				i = ARC_HEAD(a_prev);
				if (i == r) break;
			}
			else
			{
				a_prev = i->blossom_sibling;
				i = ARC_HEAD(a_prev);
			}
		}
		else
		{
			i->blossom_sibling = ARC_REV(a_prev);
			a_prev = a_next;
			i = ARC_HEAD(a_prev);
			if (i == r) break;
		}
	}
	i->blossom_sibling = ARC_REV(a_prev);
}

void PerfectMatching::ExpandInit(Node* k)
{
	Node* i = ARC_HEAD(k->blossom_sibling);
	Node* j;

	while ( 1 )
	{
		i->flag = 2; i->is_outer = 1;
		if (i == k) break;
		i->match = i->blossom_sibling;
		j = ARC_HEAD(i->match);
		j->flag = 2; j->is_outer = 1;
		j->match = ARC_REV(i->match);
		i = ARC_HEAD(j->blossom_sibling);
	}
}

void PerfectMatching::AugmentBranchInit(Node* i0, Node* r)
{
	Node* tree_root_prev = r->tree_sibling_prev;
	Node* i;
	Node* j;
	Arc* aa;

	r->flag = 2;
	i = r->first_tree_child;
	if ( i )
	while ( 1 )
	{
		ARC_HEAD(i->match)->flag = 2;
		i->flag = 2;

		MOVE_NODE_IN_TREE(i);
	}
	i = i0;
	if ( !i0->is_tree_root )
	{
		j = ARC_HEAD(i0->match);
		j->match = aa = j->tree_parent;
		i = ARC_HEAD(aa);
		while ( !i->is_tree_root )
		{
			j = ARC_HEAD(i->match);
			i->match = ARC_REV(aa);
			j->match = aa = j->tree_parent;
			i = ARC_HEAD(aa);
		}
		i->match = ARC_REV(aa);
	}
	r->is_tree_root = 0;
	tree_root_prev->tree_sibling_next = r->tree_sibling_next;
	if (r->tree_sibling_next) r->tree_sibling_next->tree_sibling_prev = tree_root_prev;
	tree_num --;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

// true_slack(a) = slack(a) + ...

// i->flag=0, i->is_processed=1:              true_slack -= eps
// i->flag=1, i->match->head->is_processed=1: true_slack += eps - slack(i->match)

void PerfectMatching::InitGlobal()
{
	Node* i;
	Node* j;
	Node* r;
	Node* r2;
	Node* r3 = NULL; // initialize to prevent compiler warning
	Edge* a;
	EdgeIterator I;
	int dir;
	Tree TREE;
	enum { NONE, AUGMENT, SHRINK } flag;

	InitGreedy();

	for (i=nodes; i<nodes+node_num; i++) i->best_edge = NULL;

	PriorityQueue<REAL> pq;

	for (r=nodes[node_num].tree_sibling_next; r; )
	{
		r2 = r->tree_sibling_next;
		if (r2) r3 = r2->tree_sibling_next;
		i = r;

		pq.Reset();

		r->tree = &TREE;

		REAL eps = 0;
		Arc* critical_arc = NULL;
		REAL critical_eps = PM_INFTY; 
		flag = NONE;
		Node* branch_root = i;
	
		while ( 1 )
		{
			i->is_processed = 1;
			i->y -= eps;
			if (!i->is_tree_root) ARC_HEAD(i->match)->y += eps;

			FOR_ALL_EDGES(i, a, dir, I)
			{
				a->slack += eps;
				j = a->head[dir];

				if (j->tree == &TREE)
				{
					// same tree
					if (j->flag == 0)
					{
						REAL slack = a->slack;
						if (!j->is_processed) slack += eps;
						if (2*critical_eps > slack || critical_arc == NULL)
						{
							flag = SHRINK;
							critical_eps = slack/2;
							critical_arc = EDGE_DIR_TO_ARC(a, dir);
							if (critical_eps <= eps) break;
							//pq.DecreaseUpperBound(critical_eps);
						}
					}
				}
				else if (j->flag == 0)
				{
					// different tree
					if (critical_eps >= a->slack || critical_arc == NULL)
					{
						flag = AUGMENT;
						critical_eps = a->slack;
						critical_arc = EDGE_DIR_TO_ARC(a, dir);
						if (critical_eps <= eps) break;
						//pq.DecreaseUpperBound(critical_eps);
					}
				}
				else
				{
					// free node
					if (a->slack > eps)
					{
						if (a->slack < critical_eps)
						{
							if (j->best_edge == NULL)
							{
								j->best_edge = a;
								pq.Add(a);
							}
							else
							{
								if (a->slack < j->best_edge->slack)
								{
									pq.Decrease(j->best_edge, a, pq_buf);
									j->best_edge = a;
								}
							}
						}
					}
					else
					{
						assert(j->flag == 2 && !j->is_blossom && !ARC_HEAD(j->match)->is_blossom);
						if (j->best_edge) pq.Remove(j->best_edge, pq_buf);
						j->flag = 1;
						j->tree = i->tree;
						j->tree_parent = EDGE_DIR_TO_ARC(a, 1-dir);
						j = ARC_HEAD(j->match);
						if (j->best_edge) pq.Remove(j->best_edge, pq_buf);
						ADD_TREE_CHILD(i, j);
					}
				}
			}

			if (dir < 2 && a)
			{
				Edge* atmp = a;
				int dirtmp = dir;
				CONTINUE_FOR_ALL_EDGES(i, atmp, dirtmp, I) atmp->slack += eps;
				break;
			}

			// move i
			if (i->first_tree_child) i = i->first_tree_child;
			else
			{
				while (i != branch_root && !i->tree_sibling_next) { i = ARC_HEAD(i->match); i = ARC_HEAD(i->tree_parent); }
				if (i == branch_root)
				{
					PriorityQueue<REAL>::Item* q = pq.GetMin();
					if (q == NULL || q->slack >= critical_eps)
					{
						eps = critical_eps;
						break;
					}
					pq.Remove(q, pq_buf);
					a = (Edge*)q;
					dir = (a->head[0]->flag == 2) ? 0 : 1;
					j = a->head[0]; 
					Arc* aa = EDGE_DIR_TO_ARC(a, dir);
					eps = a->slack;
					assert(eps < critical_eps);

					// continue growth
					i = ARC_TAIL(aa);
					j = ARC_HEAD(aa);

					assert(j->flag == 2 && !j->is_blossom && !ARC_HEAD(j->match)->is_blossom);
					j->flag = 1;
					j->tree = i->tree;
					j->tree_parent = ARC_REV(aa);
					j = ARC_HEAD(j->match);
					if (j->best_edge) pq.Remove(j->best_edge, pq_buf);
					ADD_TREE_CHILD(i, j);
					i = branch_root = j;
					continue;
				}
				i = i->tree_sibling_next;
			}
		}

		// update slacks
		i = r;
		while ( 1 )
		{
			if (i->is_processed)
			{
				i->y += eps;
				if (!i->is_tree_root)
				{
					j = ARC_HEAD(i->match);
					j->y -= eps;
					REAL delta = eps - ARC_TO_EDGE_PTR(i->match)->slack;
					FOR_ALL_EDGES(j, a, dir, I) a->slack += delta;
					j->best_edge = NULL;
				}
				FOR_ALL_EDGES(i, a, dir, I)
				{
					if (!PriorityQueue<REAL>::isReset(a))
					{
						assert(a->head[dir]->flag == 2 && a->head[dir]->best_edge == a);
						a->head[dir]->best_edge = NULL;
						PriorityQueue<REAL>::ResetItem(a);
					}
					a->slack -= eps;
				}

				i->is_processed = 0;
			}
			else
			{
				if (!i->is_tree_root) ARC_HEAD(i->match)->best_edge = NULL;
			}
			i->best_edge = NULL;

			MOVE_NODE_IN_TREE(i);
		}

		i = ARC_TAIL(critical_arc);
		j = ARC_HEAD(critical_arc);
		if (flag == SHRINK)
		{
			// shrink
			ShrinkInit(ARC_TO_EDGE_PTR(critical_arc), r);
		}
		else
		{
			// augment
			AugmentBranchInit(i, r);
			if (j->is_outer)
			{
				AugmentBranchInit(j, j);
			}
			else
			{
				ExpandInit(j);
				tree_num --;
			}
			i->match = critical_arc;
			j->match = ARC_REV(critical_arc);
		}

		r = r2;
		if (r && !r->is_tree_root) r = r3;
	}

	if (tree_num > tree_num_max)
	{
		if (trees) free(trees);
		tree_num_max = tree_num;
		trees = (Tree*) malloc(tree_num_max*sizeof(Tree));
	}
	Tree* t = trees;
	for (r=nodes; r<nodes+node_num; r++)
	{
		if (!r->is_outer)
		{
			ExpandInit(r);
			r->is_tree_root = 1;
			r->flag = 0;
			r->first_tree_child = NULL;
			if (t == trees) { nodes[node_num].tree_sibling_next = r; r->tree_sibling_prev = &nodes[node_num]; }
			else            { (t-1)->root->tree_sibling_next = r;    r->tree_sibling_prev = (t-1)->root; }
			r->tree = t;
			t->root = r;
			t->eps = 0;
			t->first[0] = t->first[1] = NULL;
			t->pq_current = NULL;
			t->pq00.Reset();
			t->pq0.Reset();
			t->pq_blossoms.Reset();
			t ++;
		}
	}
	assert(t == trees+tree_num);
	if (t == trees) nodes[node_num].tree_sibling_next = NULL;
	else            (t-1)->root->tree_sibling_next = NULL;
}
