#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PMimplementation.h"




inline void PerfectMatching::ProcessSelfloop(Node* b, Edge* a)
{
	int dir;
	Node* j;
	Node* prev[2];
	for (dir=0; dir<2; dir++)
	{
		j = a->head[dir];
		GET_PENULTIMATE_BLOSSOM(j);
		prev[dir] = j;
	}
	if (prev[0] != prev[1])
	{
		ADD_EDGE(prev[0], a, 1);
		ADD_EDGE(prev[1], a, 0);
		a->slack -= 2*prev[0]->blossom_eps;
	}
	else
	{
		a->next[0] = prev[0]->blossom_selfloops;
		prev[0]->blossom_selfloops = a;
	}

}



void PerfectMatching::Expand(Node* b)
{
	assert(b->is_blossom);
	assert(b->is_outer);
	assert(b->flag == 1);

	double start_time = get_time();

	Node* i;
	Node* j;
	Node* k;
	Edge* a;
	EdgeIterator I;
	int dir;
	ExpandTmpItem* tmp_item;
	Tree* t = b->tree;
	REAL eps = t->eps;
	Edge* a_augment = NULL;

	GET_TREE_PARENT(b, i);
	a = ARC_TO_EDGE_PTR(b->tree_parent);
	dir = ARC_TO_EDGE_DIR(b->tree_parent);

	j = a->head0[1-dir];
	GET_PENULTIMATE_BLOSSOM(j);
	MOVE_EDGE(b, j, a, dir);

	a = ARC_TO_EDGE_PTR(b->match);
	dir = ARC_TO_EDGE_DIR(b->match);
	k = a->head0[1-dir];
	GET_PENULTIMATE_BLOSSOM(k);
	MOVE_EDGE(b, k, a, dir);

	i = ARC_HEAD(k->blossom_sibling);
	while ( 1 )
	{
		tmp_item = expand_tmp_list->New();
		tmp_item->i = i; tmp_item->blossom_parent = i->blossom_parent; tmp_item->blossom_grandparent = i->blossom_grandparent;
		i->flag = 2;

		// blossom_selfloops
		i->is_outer = 1;
		while ((a=i->blossom_selfloops))
		{
			i->blossom_selfloops = a->next[0];
			ProcessSelfloop(i, a);
		}
		i->is_outer = 0;

		if (i == k) break;
		i->match = i->blossom_sibling;
		j = ARC_HEAD(i->match);
		tmp_item = expand_tmp_list->New();
		tmp_item->i = j; tmp_item->blossom_parent = j->blossom_parent; tmp_item->blossom_grandparent = j->blossom_grandparent;
		j->flag = 2;

		// blossom_selfloops
		j->is_outer = 1;
		while ((a=j->blossom_selfloops))
		{
			j->blossom_selfloops = a->next[0];
			ProcessSelfloop(j, a);
		}
		j->is_outer = 0;

		j->match = ARC_REV(i->match);
		i = ARC_HEAD(j->blossom_sibling);
	}
	k->match = b->match;
	i = ARC_TAIL(b->tree_parent);
	Arc* aa = i->blossom_sibling;
	i->flag = 1; i->tree = b->tree; i->y += b->tree->eps;
	i->tree_parent = b->tree_parent;
	if (i != k)
	{
		Node** i_ptr;
		if (i->match == aa)
		{
			i = ARC_HEAD(i->match);
			i_ptr = &j;
			while ( 1 )
			{
				aa = i->blossom_sibling;
				i->flag = 0; i->tree = b->tree; i->y -= t->eps;
				*i_ptr = i;
				i_ptr = &i->first_tree_child;
				i->tree_sibling_prev = i;
				i->tree_sibling_next = NULL;
				i = ARC_HEAD(aa);
				i->flag = 1; i->tree = b->tree; i->y += t->eps;
				i->tree_parent = ARC_REV(aa);
				if (i == k) break;
				i = ARC_HEAD(i->match);
			}
			*i_ptr = ARC_HEAD(k->match);
		}
		else
		{
			i = k;
			j = ARC_HEAD(k->match);
			do
			{
				i->tree_parent = i->blossom_sibling;
				i->flag = 1; i->tree = b->tree; i->y += b->tree->eps;
				i = ARC_HEAD(i->tree_parent);
				i->flag = 0; i->tree = b->tree; i->y -= b->tree->eps;
				i->first_tree_child = j;
				j = i;
				i->tree_sibling_prev = i;
				i->tree_sibling_next = NULL;
				i = ARC_HEAD(i->match);
			} while ( i->flag != 1 );
		}
		i = ARC_HEAD(k->match);

		j->tree_sibling_prev = i->tree_sibling_prev;
		j->tree_sibling_next = i->tree_sibling_next;
		if (i->tree_sibling_prev->tree_sibling_next) i->tree_sibling_prev->tree_sibling_next = j;
		else ARC_HEAD(b->tree_parent)->first_tree_child = j;
		if (i->tree_sibling_next) i->tree_sibling_next->tree_sibling_prev = j;
		else ARC_HEAD(b->tree_parent)->first_tree_child->tree_sibling_prev = j;

		i->tree_sibling_prev = i;
		i->tree_sibling_next = NULL;
	}

	// go through inner arcs
	i = k;
	while ( 1 )
	{
		// "-" node
		if (i->is_blossom)
		{
			a = ARC_TO_EDGE_PTR(i->match);
			REAL tmp = a->slack; a->slack = i->y; i->y = tmp;
			t->pq_blossoms.Add(a);
		}
		FOR_ALL_EDGES(i, a, dir, I)
		{
			j = a->head[dir];
			if (j->flag != 0) a->slack -= eps;
		}
		i->is_processed = 1;
		if (i->tree_parent == b->tree_parent) break;
		i = ARC_HEAD(i->tree_parent);
		// "+" node
		FOR_ALL_EDGES(i, a, dir, I)
		{
			j = a->head[dir];
			if (j->flag == 2)
			{
				a->slack += eps;
				t->pq0.Add(a);
			}
			else if (j->flag == 0 && i < j)
			{
				a->slack += 2*eps;
				t->pq00.Add(a);
			}
		}
		i->is_processed = 1;
		i = ARC_HEAD(i->match);
	}

	// go through boundary arcs
	for (tmp_item=expand_tmp_list->ScanFirst(); tmp_item; tmp_item=expand_tmp_list->ScanNext())
	{
		i = tmp_item->i;
		j = tmp_item->blossom_parent; tmp_item->blossom_parent = i->blossom_parent; i->blossom_parent = j;
		j = tmp_item->blossom_grandparent; tmp_item->blossom_grandparent = i->blossom_grandparent; i->blossom_grandparent = j;
	}
	for (dir=0; dir<2; dir++)
	{
		if (!b->first[dir]) continue;
		b->first[dir]->prev[dir]->next[dir] = NULL;

		Edge* a_next;
		for (a=b->first[dir]; a; a=a_next)
		{
			a_next = a->next[dir];
			i = a->head0[1-dir];
			GET_PENULTIMATE_BLOSSOM2(i);
			ADD_EDGE(i, a, dir);
			GET_OUTER_HEAD(a, dir, j);

			if (i->flag == 1) continue;

			if (j->flag == 0 && j->tree != t) j->tree->pq_current->pq01[1-j->tree->dir_current].Remove(a, pq_buf);

			if (i->flag == 2)
			{
				a->slack += eps;
				if (j->flag == 0) j->tree->pq0.Add(a);
			}
			else
			{
				a->slack += 2*eps;
				if      (j->flag == 2) t->pq0.Add(a);
				else if (j->flag == 0)
				{
					if (j->tree != t)
					{
						if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
						if (a->slack <= j->tree->eps + eps) a_augment = a;
					}
					j->tree->pq_current->pq00.Add(a);
				}
				else if (j->tree != t)
				{
					if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
					j->tree->pq_current->pq01[j->tree->dir_current].Add(a); 
				}
					
			}
		}
	}
	for (tmp_item=expand_tmp_list->ScanFirst(); tmp_item; tmp_item=expand_tmp_list->ScanNext())
	{
		i = tmp_item->i;
		i->blossom_parent = tmp_item->blossom_parent;
		i->blossom_grandparent = tmp_item->blossom_grandparent;
		i->is_outer = 1;
	}
	expand_tmp_list->Reset();

	b->tree_sibling_next = removed_first;
	removed_first = b;
	removed_num ++;
	if (4*removed_num > node_num) FreeRemoved();

	blossom_num --;
	stat.expand_count ++;

	stat.expand_time += get_time() - start_time;

	if (a_augment) Augment(a_augment);
}

void PerfectMatching::FreeRemoved()
{
	Node* i0;
	Node* i;
	for (i0=nodes; i0<nodes+node_num; i0++)
	{
		for (i=i0; !i->is_outer && !i->is_marked; i=i->blossom_parent)
		{
			i->is_marked = 1;
			if (i->blossom_grandparent->is_removed) i->blossom_grandparent = i->blossom_parent;
		}
	}
	for (i0=nodes; i0<nodes+node_num; i0++)
	{
		for (i=i0; !i->is_outer && i->is_marked; i=i->blossom_parent)
		{
			i->is_marked = 0;
		}
	}

	while ((i=removed_first))
	{
		removed_first = i->tree_sibling_next;
		blossoms->Delete(i);
		removed_num --;
	}

	assert(removed_num == 0);
}

