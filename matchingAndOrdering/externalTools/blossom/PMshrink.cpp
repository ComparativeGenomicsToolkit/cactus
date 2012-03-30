#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PMimplementation.h"


PerfectMatching::Node* PerfectMatching::FindBlossomRoot(Edge* a0)
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
		if (_i[branch]->is_marked)
		{
			r = _i[branch]; 
			j = _i[1-branch];
			break; 
		}
		_i[branch]->is_marked = 1;
		if (_i[branch]->is_tree_root)
		{
			j = _i[branch];
			i = _i[1-branch];
			while (!i->is_marked)
			{
				i->is_marked = 1;
				i = ARC_HEAD(i->match);
				GET_TREE_PARENT(i, i);
			}
			r = i;
			break;
		}
		i = ARC_HEAD(_i[branch]->match);
		GET_TREE_PARENT(i, _i[branch]);
		branch = 1 - branch;
	}
	i = r;
	while ( i != j )
	{
		i = ARC_HEAD(i->match);
		i = ARC_HEAD(i->tree_parent);
		i->is_marked = 0;
	}
	// clear is_marked and is_outer
	i = ARC_HEAD(a0);
	while (i != r)
	{
		i->is_marked = 0;
		i->is_outer = 0;
		i = ARC_HEAD(i->match);
		i->is_outer = 0;
		i = ARC_HEAD(i->tree_parent);
	}
	i = ARC_TAIL(a0);
	while (i != r)
	{
		i->is_marked = 0;
		i->is_outer = 0;
		i = ARC_HEAD(i->match);
		i->is_outer = 0;
		i = ARC_HEAD(i->tree_parent);
	}
	r->is_marked = 0;
	r->is_outer = 0;

	return r;
}


void PerfectMatching::Shrink(Edge* a0)
{
	//assert(a0->head[0]->is_outer && a0->head[1]->is_outer);
	//assert(a0->head[0]->flag == 0 && a0->head[1]->flag == 0);

	double start_time = get_time();

	int branch, dir;
	Node* r;
	Node* i;
	Node* j;
	Edge* a;
	Edge** a_inner_ptr;
	Arc* a_prev;
	Node* b = blossoms->New();
	Edge* a_augment = NULL;
	Edge* b_match;

	b->first[0] = b->first[1] = NULL;

	// set is_outer=0 for all nodes in the blossom
	r = FindBlossomRoot(a0);
	Tree* t = r->tree;
	REAL eps = t->eps;

	b->first_tree_child = NULL;
	i = ARC_HEAD(a0);
	branch = 0;
	while ( 1 )
	{
		if (i == r && branch) break;
		i->is_marked = 1;
		if (i == r)
		{
			branch = 1;
			i = ARC_TAIL(a0);
			continue;
		}

		// remove i from the list of children
		REMOVE_FROM_TREE(i);

		// move children of i to the list of children of b
		if (i->first_tree_child)
		{
			j = i->first_tree_child;
			if (!b->first_tree_child) b->first_tree_child = j;
			else
			{
				Node* j_last = j->tree_sibling_prev;
				j->tree_sibling_prev = b->first_tree_child->tree_sibling_prev;
				b->first_tree_child->tree_sibling_prev->tree_sibling_next = j;
				b->first_tree_child->tree_sibling_prev = j_last;
			}
		}

		// go to parent
		i = ARC_HEAD(i->match);
		i->is_marked = 1;
		if (i->is_blossom)
		{
			a = ARC_TO_EDGE_PTR(i->match);
			t->pq_blossoms.Remove(a, pq_buf);
			REAL tmp = a->slack; a->slack = i->y; i->y = tmp;
		}
		i = ARC_HEAD(i->tree_parent); 
	}

	// move children of r to the list of children of b
	if (i->first_tree_child)
	{
		j = i->first_tree_child;
		if (!b->first_tree_child) b->first_tree_child = j;
		else
		{
			Node* j_last = j->tree_sibling_prev;
			j->tree_sibling_prev = b->first_tree_child->tree_sibling_prev;
			b->first_tree_child->tree_sibling_prev->tree_sibling_next = j;
			b->first_tree_child->tree_sibling_prev = j_last;
		}
	}

	// init b
	b->is_removed = 0;
	b->is_outer = 1;
	b->flag = 0;
	b->is_blossom = 1;
	b->is_tree_root = r->is_tree_root;
	b->is_processed = 1;
	b->tree = t;
	b->y = -eps;
	b->is_marked = 0;

	// replace r with b in the tree
	b->tree_sibling_prev = r->tree_sibling_prev;
	b->tree_sibling_next = r->tree_sibling_next;
	Node* b_parent = NULL;
	if (!b->is_tree_root)
	{
		b_parent = ARC_HEAD(r->match); GET_TREE_PARENT(b_parent, b_parent); 
	}
	if (b->tree_sibling_prev->tree_sibling_next) b->tree_sibling_prev->tree_sibling_next = b;
	else b_parent->first_tree_child = b;
	if (b->tree_sibling_next) b->tree_sibling_next->tree_sibling_prev = b;
	else if (b_parent) b_parent->first_tree_child->tree_sibling_prev = b;

	if (b->is_tree_root)
	{
		b->tree->root = b;
		b_match = NULL;
	}
	else
	{
		b->match = r->match;
		b_match = ARC_TO_EDGE_PTR(b->match);
	}
	REAL b_match_slack = 0; // initialize to prevent compiler warning
	if (b_match && ARC_HEAD(b->match)->is_blossom)
	{
		b_match_slack = b_match->slack;
		b_match->slack = ARC_HEAD(b->match)->y;
	}

	// second pass over nodes in the blossom
	branch = 0;
	a_prev = EDGE_DIR_TO_ARC(a0, 0);
	i = ARC_HEAD(a_prev);
	while ( 1 )
	{
		// update Arc::next and Arc::head pointers
		if (i->flag == 0) i->y += eps;
		else              i->y -= eps;
		i->is_processed = 0;

		if (i->flag == 1)
		{
			Edge* a_prev;
			for (dir=0; dir<2; dir++)
			if (i->first[dir])
			{
				for (a_inner_ptr=&i->first[dir], a=*a_inner_ptr, a_prev=a->prev[dir], a_prev->next[dir]=NULL; a; a=*a_inner_ptr)
				{
					Node* j0 = a->head[dir];
					for (j=j0; !j->is_outer && !j->is_marked; j = j->blossom_parent) {}
					if (j != j0) { /*assert(j->flag == 0);*/ int dir_rev = 1 - dir; MOVE_EDGE(j0, j, a, dir_rev); }
					if (j->is_marked) // "inner" arc
					{
						a_inner_ptr = &a->next[dir];
						a->prev[dir] = a_prev;
						a_prev = a;

						if (j->flag == 1) a->slack += eps;
					}
					else // "boundary" arc
					{
						*a_inner_ptr = a->next[dir];
						ADD_EDGE(b, a, dir);

						if (j->flag == 0 && j->tree != t) 
						{
							j->tree->pq_current->pq01[1-j->tree->dir_current].Remove(a, pq_buf);
							if (a->slack + eps <= j->tree->eps) a_augment = a;
						}
						a->slack += 2*eps;
						if (j->flag == 2) t->pq0.Add(a);
						else if (j->flag == 0)
						{
							if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
							j->tree->pq_current->pq00.Add(a);
						}
						else if (j->tree != t)
						{
							if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
							j->tree->pq_current->pq01[j->tree->dir_current].Add(a);
						}
					}
				}
				if (i->first[dir])
				{
					a_prev->next[dir] = i->first[dir];
					i->first[dir]->prev[dir] = a_prev;
				}
			}
		}

		Arc* a_next = (i->flag == 0) ? i->match : i->tree_parent;
		i->blossom_parent = b;
		i->match = NULL;
		i->blossom_grandparent = b;
		i->blossom_selfloops = NULL;
		if (branch == 0)
		{
			i->blossom_sibling = a_next;
			if (i == r)
			{
				branch = 1;
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
	r->is_tree_root = 0;

	for (i=ARC_HEAD(r->blossom_sibling); ; i = ARC_HEAD(i->blossom_sibling))
	{
		i->is_marked = 0;
		i->blossom_eps = eps;
		if (i == r) break;
	}

	if (b_match)
	{
		if (ARC_HEAD(b->match)->is_blossom)
		{
			b_match->slack = b_match_slack;
		}
		dir = ARC_TO_EDGE_DIR(b->match);
		//assert(b_match->head[1-dir] == r);
		MOVE_EDGE(r, b, b_match, dir);
	}

	stat.shrink_count ++;
	blossom_num ++;
	stat.shrink_time += get_time() - start_time;

	if (a_augment) Augment(a_augment);
}

