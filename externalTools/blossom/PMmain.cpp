#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PMimplementation.h"


void PerfectMatching::Finish()
{

#define IS_VALID_MATCH(i) ((Edge*)(i->match) >= edges && (Edge*)(i->match) < edges + edge_num)

	Node* i0;
	Node* i;
	Node* j;
	Node* k;
	Node* b;
	Node* b_prev;
	Node* b_prev_prev;

	for (i0=nodes; i0<nodes+node_num; i0++)
	{
		if (IS_VALID_MATCH(i0)) continue;
		b_prev = NULL;
		b = i0;
		do
		{
			b->blossom_grandparent = b_prev;
			b_prev = b;
			b = b->blossom_parent;
		} while (!IS_VALID_MATCH(b));

		b_prev_prev = b_prev->blossom_grandparent;
		while ( 1 )
		{
			for (k=ARC_TAIL0(b->match); k->blossom_parent!=b; k=k->blossom_parent) {}
			k->match = b->match;
			i = ARC_HEAD(k->blossom_sibling);
			while ( i != k )
			{
				i->match = i->blossom_sibling;
				j = ARC_HEAD(i->match);
				j->match = ARC_REV(i->match);
				i = ARC_HEAD(j->blossom_sibling);
			}

			b = b_prev;
			if (!b->is_blossom) break;
			b_prev = b_prev_prev;
			b_prev_prev = b_prev->blossom_grandparent;
		}
	}
}


void PerfectMatching::AddTreeEdge(Tree* t0, Tree* t1)
{
	TreeEdge* e = tree_edges->New();
	e->head[0] = t1;
	e->head[1] = t0;
	e->next[0] = t0->first[0];
	t0->first[0] = e;
	e->next[1] = t1->first[1];
	t1->first[1] = e;

	e->pq00.Reset();
	e->pq01[0].Reset();
	e->pq01[1].Reset();

	t1->pq_current = e;
	t1->dir_current = 0;
}

bool PerfectMatching::ProcessEdge00(Edge* a, bool update_boundary_edge)
{
	int dir;
	Node* j;
	Node* prev[2];
	Node* last[2];
	for (dir=0; dir<2; dir++)
	{
		if (a->head[dir]->is_outer) { prev[dir] = NULL; last[dir] = a->head[dir]; }
		else
		{
			j = a->head[dir];
			GET_PENULTIMATE_BLOSSOM(j);
			prev[dir] = j;
			last[dir] = prev[dir]->blossom_parent;
			//assert(last[dir]->is_outer);
		}
	}

	if (last[0] != last[1])
	{
		for (dir=0; dir<2; dir++)
		{
			j = a->head[dir];
			if (j != last[dir]) { int dir_rev = 1 - dir; MOVE_EDGE(j, last[dir], a, dir_rev); }
		}
		if (update_boundary_edge) a->slack -= 2*a->head[0]->tree->eps;
		return true;
	}

	if (prev[0] != prev[1])
	{
		for (dir=0; dir<2; dir++)
		{
			j = a->head[dir];
			if (j != prev[dir]) { int dir_rev = 1 - dir; MOVE_EDGE(j, prev[dir], a, dir_rev); }
		}
		a->slack -= 2*prev[0]->blossom_eps;
		return false;
	}

	for (dir=0; dir<2; dir++)
	{
		j = a->head[1-dir];
		REMOVE_EDGE(j, a, dir);
	}
	a->next[0] = prev[0]->blossom_selfloops;
	prev[0]->blossom_selfloops = a;
	return false;
}


inline void PerfectMatching::AugmentBranch(Node* i0)
{
	int dir;
	Tree* t = i0->tree;
	Node* r = t->root;
	Node* tree_root_prev = r->tree_sibling_prev;
	Node* i;
	Node* j;
	Edge* a;
	EdgeIterator I;
	Arc* aa;
	REAL eps = t->eps;
	PriorityQueue<REAL>::Item* q;
	TreeEdge* e;
	TreeEdgeIterator T;
	Tree* t2;

	t = r->tree;
	t->pq_current = t;

	FOR_ALL_TREE_EDGES_X(t, e, dir, T)
	{
		t2 = e->head[dir];
		e->head[1-dir] = NULL; // mark it for deletion

		t2->pq_current = e;
		t2->dir_current = dir;
	}

	i = r->first_tree_child;
	if ( i )
	while ( 1 )
	{
		Node* i0 = i;
		i = ARC_HEAD(i->match);
		if (i->is_processed)
		{
			if (i->is_blossom)
			{
				a = ARC_TO_EDGE_PTR(i->match);
				REAL tmp = a->slack; a->slack = i->y; i->y = tmp;
				PriorityQueue<REAL>::ResetItem(a);
			}
			FOR_ALL_EDGES(i, a, dir, I)
			{
				GET_OUTER_HEAD(a, dir, j);

				if (j->flag == 0 && j->is_processed)
				{
					if (j->tree != t)
					{
						a->slack += eps;
						if (PriorityQueue<REAL>::isReset(a)) j->tree->pq0.Add(a);
					}
				}
				else a->slack += eps;
			}
		}

		i = i0;
		MOVE_NODE_IN_TREE(i);
	}

	///////////////////////////////////////////////////////////////////

	FOR_ALL_TREE_EDGES(t, e, dir)
	{
		t2 = e->head[dir];
		t2->pq_current = NULL;

		e->pq01[1-dir].Merge(t2->pq0);
		for (q=e->pq00.GetFirst(); q; q=e->pq00.GetNext(q))
		{
			q->slack -= eps;
			int dir2;
			for (dir2=0; dir2<2; dir2++) GET_OUTER_HEAD((Edge*)q, dir2, j);
		}
		e->pq00.Merge(t2->pq0);
		for (q=e->pq01[dir].GetAndResetFirst(); q; q=e->pq01[dir].GetAndResetNext())
		{
			q->slack -= eps;
			int dir2;
			for (dir2=0; dir2<2; dir2++) GET_OUTER_HEAD((Edge*)q, dir2, j);
		}
	}
	for (q=t->pq0.GetAndResetFirst(); q; q=t->pq0.GetAndResetNext())
	{
		q->slack -= eps;
		int dir2;
		for (dir2=0; dir2<2; dir2++) GET_OUTER_HEAD((Edge*)q, dir2, j);
	}
	for (q=t->pq00.GetAndResetFirst(); q; q=t->pq00.GetAndResetNext())
	{
		ProcessEdge00((Edge*)q);
	}

	///////////////////////////////////////////////////////////////////

	r->flag = 2;
	r->is_processed = 0;
	i = r->first_tree_child;
	r->y += eps;
	if ( i )
	while ( 1 )
	{
		j = ARC_HEAD(i->match);
		j->flag = 2;
		i->flag = 2;
		j->is_processed = 0;
		i->is_processed = 0;
		j->y -= eps;
		i->y += eps;

		MOVE_NODE_IN_TREE(i);
	}

	///////////////////////////////////////////////////////////////////

	i = i0;
	if ( !i0->is_tree_root )
	{
		j = ARC_HEAD(i0->match);
		GET_TREE_PARENT(j, i);
		j->match = aa = j->tree_parent;
		while ( !i->is_tree_root )
		{
			j = ARC_HEAD(i->match);
			i->match = ARC_REV(aa);
			GET_TREE_PARENT(j, i);
			j->match = aa = j->tree_parent;
		}
		i->match = ARC_REV(aa);
	}
	r->is_tree_root = 0;
	tree_root_prev->tree_sibling_next = r->tree_sibling_next;
	if (r->tree_sibling_next) r->tree_sibling_next->tree_sibling_prev = tree_root_prev;
	tree_num --;
}


void PerfectMatching::Augment(Edge* a)
{
	Node* j;
	int dir;

	for (dir=0; dir<2; dir++)
	{
		GET_OUTER_HEAD(a, dir, j);
		AugmentBranch(j);
		j->match = EDGE_DIR_TO_ARC(a, 1-dir);
	}
	if (options.verbose)
	{
		int k = 1;
		while (k < tree_num) k *= 2;
		if (k == tree_num || tree_num<=8 || (tree_num<=64 && (tree_num%8)==0)) { printf("%d.", tree_num); fflush(stdout); }
	}
}

inline void PerfectMatching::GrowNode(Node* i)
{
	//assert(i->is_outer);
	//assert(i->flag == 0);

	Edge* a;
	EdgeIterator I;
	int dir;
	Node* j;
	Tree* t = i->tree;
	REAL eps = t->eps;
	Edge* a_augment = NULL;

	FOR_ALL_EDGES(i, a, dir, I)
	{
		GET_OUTER_HEAD(a, dir, j);

		if (j->flag == 2)
		{
			a->slack += eps; 
			if (a->slack > 0) 
			{ 
				t->pq0.Add(a); 
			}
			else
			{
				j->flag = 1;
				j->tree = i->tree;
				j->tree_parent = EDGE_DIR_TO_ARC(a, 1-dir);
				j->y += eps;
				j = ARC_HEAD(j->match);
				j->y -= eps;
				ADD_TREE_CHILD(i, j);
			}
		}
		else
		{
			if (j->flag == 0 && j->is_processed)
			{
				if (!PriorityQueue<REAL>::isReset(a)) j->tree->pq0.Remove(a, pq_buf);
				if (a->slack <= j->tree->eps && j->tree != t) a_augment = a;
				a->slack += eps;
				if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
				j->tree->pq_current->pq00.Add(a);
			}
			else
			{
				a->slack += eps;
				if (j->flag == 1 && j->tree != t)
				{
					if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
					j->tree->pq_current->pq01[j->tree->dir_current].Add(a);
				}
			}
		}
	}

	//assert(!i->is_processed);
	i->is_processed = 1;

	if (!i->is_tree_root)
	{
		j = ARC_HEAD(i->match);
		//assert(!j->is_processed);
		j->is_processed = 1;
		if (j->is_blossom)
		{
			a = ARC_TO_EDGE_PTR(i->match);
			REAL tmp = a->slack; a->slack = j->y; j->y = tmp;
			t->pq_blossoms.Add(a);
		}
	}

	if (a_augment) Augment(a_augment);

	stat.grow_count ++;
}



void PerfectMatching::GrowTree(Node* r, bool new_subtree)
{
	//assert(r->flag == 0);

	Node* i = r;
	Node* j;
	Node* stop = r->tree_sibling_next;
	if (new_subtree && r->first_tree_child) stop = r->first_tree_child;
	Edge* a;
	EdgeIterator I;
	int dir;
	Tree* t = r->tree;
	REAL eps = t->eps;
	int tree_num0 = tree_num;

	while ( 1 )
	{
		if (!i->is_tree_root)
		{
			// process "-" node
			i = ARC_HEAD(i->match);
			FOR_ALL_EDGES(i, a, dir, I)
			{
				GET_OUTER_HEAD(a, dir, j);

				if (j->flag == 2) a->slack -= eps;
				else
				{
					if (j->flag == 0 && j->is_processed)
					{
						if (!PriorityQueue<REAL>::isReset(a)) j->tree->pq0.Remove(a, pq_buf);
						a->slack -= eps;
						if (j->tree != t)
						{
							if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
							j->tree->pq_current->pq01[1-j->tree->dir_current].Add(a);
						}
					}
					else a->slack -= eps;
				}
			}
			i = ARC_HEAD(i->match);
		}
		// process "+" node
		GrowNode(i);
		if (tree_num != tree_num0) break;

		if (i->first_tree_child) i = i->first_tree_child;
		else
		{
			while (i != r && !i->tree_sibling_next) { i = ARC_HEAD(i->match); GET_TREE_PARENT(i, i); }
			i = i->tree_sibling_next;
		}
		if (i == stop) break;
	}
}

void PerfectMatching::Solve(bool finish)
{
	Node* i;
	Node* j;
	Node* r;
	Node* r2;
	Node* r3 = NULL; // initialize to prevent compiler warning
	PriorityQueue<REAL>::Item* q;
	Edge* a;
	Tree* t;
	Tree* t2;
	TreeEdge* e;
	TreeEdgeIterator T;
	int dir;
	REAL eps;

	double start_time = get_time();

	if (IS_INT)
	{
		if (options.dual_greedy_update_option == 2)
		{
			printf("Fixed eps approach can only be used with floating point REAL!\n");
			printf("Change REAL to double in PerfectMatching.h and recompile\n");
			exit(1);
		}
		if (options.dual_LP_threshold > 0)
		{
			printf("LP approach can only be used with floating point REAL!\n");
			printf("Change REAL to double in PerfectMatching.h and recompile\n");
			exit(1);
		}
	}
	if (options.verbose) { printf("perfect matching with %d nodes and %d edges\n", node_num, edge_num); fflush(stdout); }

	if (first_solve)
	{
		if (options.verbose) { printf("    starting init..."); fflush(stdout); }
		if (options.fractional_jumpstart) InitGlobal();
		else                              InitGreedy();
		if (options.verbose) printf("done [%.3f secs]. ", get_time() - start_time);
		first_solve = false;
	}
	else if (options.verbose) printf("    solving updated problem. ");

	if (options.verbose) { printf("%d trees\n    .", tree_num); fflush(stdout); }

	memset(&stat, 0, sizeof(Stat));

	///////////////////////////////////////////////////////
	//       first pass - initialize auxiliary graph     //
	///////////////////////////////////////////////////////

	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		//assert(!r->is_processed);
		t = r->tree;
		//assert(!t->first[0] && !t->first[1]);

		EdgeIterator I;
		FOR_ALL_EDGES(r, a, dir, I)
		{
			j = a->head[dir];
			if (j->flag == 2) t->pq0.Add(a);
			else if (j->is_processed)
			{
				//assert(j->flag == 0);
				if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
				j->tree->pq_current->pq00.Add(a);
			}
		}
		r->is_processed = 1;
		FOR_ALL_TREE_EDGES(t, e, dir) e->head[dir]->pq_current = NULL;
	}

	///////////////////////////////////////////////////////
	//                  main loop                        //
	///////////////////////////////////////////////////////

	while ( 1 )
	{
		int tree_num0 = tree_num;
		Stat stat0 = stat;
		REAL delta = 0;

		for (r=nodes[node_num].tree_sibling_next; r; )
		{
			r2 = r->tree_sibling_next;
			if (r2) r3 = r2->tree_sibling_next;
			t = r->tree;

			int tree_num1 = tree_num;

			//////////////////////////////////////////////////////////////////////
			// step 1 - traversing auxiliary graph, setting pq_current pointers //
			//////////////////////////////////////////////////////////////////////
			t->pq_current = t;
			if (options.update_duals_before)
			{
				eps = PM_INFTY;
				Edge* a_augment = NULL;
				REAL eps_augment = PM_INFTY;
				if ((q=t->pq0.GetMin())) eps = q->slack;
				if ((q=t->pq_blossoms.GetMin()) && eps > q->slack) eps = q->slack;
				while ((q=t->pq00.GetMin()))
				{
					if (ProcessEdge00((Edge*)q, false)) break;
					t->pq00.Remove(q, pq_buf);
				}
				if (q && 2*eps > q->slack) eps = q->slack/2;
				FOR_ALL_TREE_EDGES_X(t, e, dir, T)
				{
					t2 = e->head[dir];
					t2->pq_current = e;
					t2->dir_current = dir;
					if ((q=e->pq00.GetMin()) && (!a_augment || eps_augment > q->slack-t2->eps)) { a_augment = (Edge*)q; eps_augment = q->slack-t2->eps; }
					if ((q=e->pq01[dir].GetMin()) && eps > q->slack+t2->eps) eps = q->slack+t2->eps;
				}
				if (eps > eps_augment) eps = eps_augment;
				if (eps > t->eps)
				{
					delta += eps - t->eps;
					t->eps = eps;
				}
				if (a_augment && eps_augment <= t->eps) Augment(a_augment);
			}
			else
			{
				FOR_ALL_TREE_EDGES_X(t, e, dir, T)
				{
					t2 = e->head[dir];
					t2->pq_current = e;
					t2->dir_current = dir;

					if ((q=e->pq00.GetMin()) && (q->slack - t->eps <= t2->eps))
					{
						Augment((Edge*)q);
						break;
					}
				}
			}

			/////////////////////////////////
			//   step 2 - growing tree     //
			/////////////////////////////////
			eps = t->eps;
			REAL twice_eps = 2*eps;

			while ( tree_num1 == tree_num )
			{
				if ((q=t->pq0.GetMin()) && q->slack <= t->eps)
				{
					a = (Edge*)q;
					dir = (a->head[1]->flag == 2 && a->head[1]->is_outer) ? 1 : 0;
					GET_OUTER_HEAD(a, 1-dir, i);
					j = a->head[dir];
					//assert(i->flag==0 && j->flag==2 && i->is_outer && j->is_outer && i->tree==t);

					j->flag = 1;
					j->tree = i->tree;
					j->tree_parent = EDGE_DIR_TO_ARC(a, 1-dir);
					j->y += eps;
					j = ARC_HEAD(j->match);
					j->y -= eps;
					ADD_TREE_CHILD(i, j);

					GrowTree(j, true);
				}
				else if ((q=t->pq00.GetMin()) && q->slack <= twice_eps)
				{
					t->pq00.Remove(q, pq_buf);
					a = (Edge*)q;
					if (ProcessEdge00(a)) Shrink(a);
				}
				else if ((q=t->pq_blossoms.GetMin()) && q->slack <= eps)
				{
					t->pq_blossoms.Remove(q, pq_buf);
					a = (Edge*)q;
					j = (a->head[0]->flag == 1) ? a->head[0] : a->head[1];
					REAL tmp = a->slack; a->slack = j->y; j->y = tmp;
					Expand(j);
				}
				else break;
			}

			///////////////////////////////////////////////////////////////////////
			// step 3 - traversing auxiliary graph, clearing pq_current pointers //
			///////////////////////////////////////////////////////////////////////
			if ( tree_num1 == tree_num )
			{
				t->pq_current = NULL;
				if (options.update_duals_after)
				{
					eps = PM_INFTY;
					Edge* a_augment = NULL;
					REAL eps_augment = PM_INFTY;
					if ((q=t->pq0.GetMin())) eps = q->slack;
					if ((q=t->pq_blossoms.GetMin()) && eps > q->slack) eps = q->slack;
					while ((q=t->pq00.GetMin()))
					{
						if (ProcessEdge00((Edge*)q, false)) break;
						t->pq00.Remove(q, pq_buf);
					}
					if (q && 2*eps > q->slack) eps = q->slack/2;
					FOR_ALL_TREE_EDGES(t, e, dir)
					{
						t2 = e->head[dir];
						e->head[dir]->pq_current = NULL;
						if ((q=e->pq00.GetMin()) && (!a_augment || eps_augment > q->slack-t2->eps)) { a_augment = (Edge*)q; eps_augment = q->slack-t2->eps; }
						if ((q=e->pq01[dir].GetMin()) && eps > q->slack+t2->eps) eps = q->slack+t2->eps;
					}
					if (eps > eps_augment) eps = eps_augment;
					bool progress = false;
					if (eps > t->eps)
					{
						delta += eps - t->eps;
						t->eps = eps;
						progress = true;
					}
					if (a_augment && eps_augment <= t->eps) Augment(a_augment);
					else if (progress && tree_num >= options.single_tree_threshold*node_num)
					{
						// continue with the same tree
						r = t->root;
						continue;
					}
				}
				else
				{
					FOR_ALL_TREE_EDGES(t, e, dir) e->head[dir]->pq_current = NULL;
				}
			}

			///////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////

			r = r2;
			if (r && !r->is_tree_root) r = r3;
		}

		if (tree_num == 0) break;

		if ( tree_num == tree_num0 )
		  //&& stat.grow_count == stat0.grow_count 
		  //&& stat.shrink_count == stat0.shrink_count 
		  //&& stat.expand_count == stat0.expand_count )
		{
			if (!UpdateDuals())
			{
				if (!IS_INT && delta <= PM_THRESHOLD) // for numerical stability
				{
					//CommitEps();
					int dual_greedy_update_option = options.dual_greedy_update_option;
					options.dual_greedy_update_option = 2;
					UpdateDuals();
					options.dual_greedy_update_option = dual_greedy_update_option;
				}
			}
		}
	}

	if (finish) Finish();

	if (options.verbose)
	{
		printf("\ndone [%.3f secs]. %d grows, %d expands, %d shrinks\n", get_time()-start_time, stat.grow_count, stat.expand_count, stat.shrink_count); 
		printf("    expands: [%.3f secs], shrinks: [%.3f secs], dual updates: [%.3f secs]\n", stat.expand_time, stat.shrink_time, stat.dual_time); 
		fflush(stdout); 
	}
}



