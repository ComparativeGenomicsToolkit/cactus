#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PMimplementation.h"
#include "MinCost/MinCost.h"


void PerfectMatching::ComputeEpsGlobal()
{
	Node* r;
	PriorityQueue<REAL>::Item* q;
	Tree* t;
	Tree* t2;
	TreeEdge* e;
	int i, j, k, N = 0, E = 0;

	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		t = r->tree;
		t->id = N;
		N += 2;
		for (k=0; k<2; k++)
		for (e=t->first[k]; e; e=e->next[k]) E += 6;
	}
	DualMinCost<REAL>* m = new DualMinCost<REAL>(N, E);
	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		t = r->tree;
		i = t->id;
		m->AddUnaryTerm(i, -1);
		m->SetLowerBound(i, 0);
		m->AddUnaryTerm(i+1, 1);
		m->SetUpperBound(i+1, 0);

		if (t->eps_delta < PM_INFTY)
		{
			m->SetUpperBound(i, t->eps_delta);
			m->SetLowerBound(i+1, -t->eps_delta);
		}
		for (e=t->first[0]; e; e=e->next[0])
		{
			t2 = e->head[0];
			if (t2 == NULL) continue;
			j = e->head[0]->id;
			if ((q=e->pq01[0].GetMin()))
			{
				m->AddConstraint(j, i, q->slack - t->eps + t2->eps);
				m->AddConstraint(i+1, j+1, q->slack - t->eps + t2->eps);
			}
			if ((q=e->pq01[1].GetMin()))
			{
				m->AddConstraint(i, j, q->slack - t2->eps + t->eps);
				m->AddConstraint(j+1, i+1, q->slack - t2->eps + t->eps);
			}
			if ((q=e->pq00.GetMin()))
			{
				m->AddConstraint(i+1, j, q->slack - t->eps - t2->eps);
				m->AddConstraint(j+1, i, q->slack - t->eps - t2->eps);
			}
		}
	}
	m->Solve();
	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		t = r->tree;
		i = t->id;
		t->eps_delta = (m->GetSolution(i) - m->GetSolution(i+1))/2;
	}
	delete m;
}

void PerfectMatching::ComputeEpsSingle()
{
	Node* r;
	PriorityQueue<REAL>::Item* q;
	Tree* t;
	Tree* t2;
	TreeEdge* e;
	REAL eps = PM_INFTY;

	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		t = r->tree;
		if (eps > t->eps_delta) eps = t->eps_delta;
		for (e=t->first[0]; e; e=e->next[0])
		{
			t2 = e->head[0];
			if ((q=e->pq00.GetMin()) && 2*eps > q->slack-t->eps-t2->eps)
			{
				eps = (q->slack-t->eps-t2->eps)/2;
			}
		}
	}
	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		r->tree->eps_delta = eps;
	}
}


void PerfectMatching::ComputeEpsCC()
{
	Node* r;
	PriorityQueue<REAL>::Item* q;
	Tree* t0;
	Tree* t;
	Tree* t2;
	Tree* t_next;
	TreeEdge* e;
	REAL eps, eps2;
	Tree* queue_last;
	int dir;
	Tree* FIXED_TREE = trees-1;
	int component_num = 0;
	TreeEdge** e_ptr;

	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		t0 = r->tree;
		t0->next = NULL;
	}
	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		t0 = r->tree;
		if (t0->next) continue;
		eps = t0->eps_delta;

		t0->next = queue_last = t = t0;
		while ( 1 )
		{
			for (dir=0; dir<2; dir++)
			for (e_ptr=&t->first[dir], e=*e_ptr; e; e=*e_ptr)
			{
				t2 = e->head[dir];
				if (t2 == NULL) { *e_ptr = e->next[dir]; tree_edges->Delete(e); continue; }
				e_ptr = &e->next[dir];

				REAL eps00 = ((q=e->pq00.GetMin())) ? (q->slack - t->eps - t2->eps) : PM_INFTY;
				if (t2->next && t2->next != FIXED_TREE)
				{
					if (2*eps > eps00) eps = eps00/2;
					continue;
				}

				REAL eps01[2];
				eps01[dir]   = ((q=e->pq01[dir].GetMin()))   ? (q->slack - t->eps + t2->eps) : PM_INFTY;
				eps01[1-dir] = ((q=e->pq01[1-dir].GetMin())) ? (q->slack - t2->eps + t->eps) : PM_INFTY;

				if (t2->next == FIXED_TREE) eps2 = t2->eps_delta;
				else if (eps01[0] > 0 && eps01[1] > 0) eps2 = 0;
				else
				{
					queue_last->next = t2;
					queue_last = t2;
					t2->next = t2;
					if (eps > eps00) eps = eps00;
					if (eps > t2->eps_delta) eps = t2->eps_delta;
					continue;
				}
				if (eps > eps00 - eps2)      eps = eps00 - eps2;
				if (eps > eps2 + eps01[dir]) eps = eps2 + eps01[dir];
			}

			if (t->next == t) break;
			t = t->next;
		}
		for (t=t0; ; t=t_next)
		{
			t->eps_delta = eps;
			t_next = t->next;
			t->next = FIXED_TREE;
			if (t_next == t) break;
		}
		component_num ++;
	}
	//printf("%d CCs ", component_num);
}


void PerfectMatching::ComputeEpsSCC()
{
	PriorityQueue<REAL>::Item* q;
	Node* r;
	Tree* t0;
	Tree* t;
	Tree* t2;
	TreeEdge* e;
	TreeEdge** e_ptr;
	REAL eps;
	int c, dir;

	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		t0 = r->tree;
		t0->dfs_parent = NULL;

		for (dir=0; dir<2; dir++)
		for (e_ptr=&t0->first[dir], e=*e_ptr; e; e=*e_ptr)
		{
			t2 = e->head[dir];
			if (t2 == NULL) { *e_ptr = e->next[dir]; tree_edges->Delete(e); continue; }
			e_ptr = &e->next[dir];
		}
	}
	Tree* stack = NULL;

	// first DFS
	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		t0 = r->tree;
		if (t0->dfs_parent) continue;
		t = t0;
		e = (t->first[0]) ? t->first[0] : t->first[1];
		t->dfs_parent = (TreeEdge*)trees;
		while ( 1 )
		{
			if (e == NULL)
			{
				t->next = stack;
				stack = t;

				if (t == t0) break;

				e = t->dfs_parent;
				if (t == e->head[0]) { t = e->head[1]; e = (e->next[0]) ? e->next[0] : t->first[1]; }
				else                 { t = e->head[0]; e = e->next[1]; }
				continue;
			}

			if (e->head[1] == t)
			{
				if (e->head[0]->dfs_parent || !(q=e->pq01[0].GetMin()) || q->slack - t->eps + e->head[0]->eps > 0) { e = (e->next[0]) ? e->next[0] : t->first[1]; continue; }
				t = e->head[0];
			}
			else
			{
				if (e->head[1]->dfs_parent || !(q=e->pq01[1].GetMin()) || q->slack - t->eps + e->head[1]->eps > 0) { e = e->next[1]; continue; }
				t = e->head[1];
			}
			t->dfs_parent = e;
			e = (t->first[0]) ? t->first[0] : t->first[1];
		}
	}

	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next) r->tree->dfs_parent = NULL;

	int component_num = 0;
	while (stack)
	{
		t0 = stack;
		stack = t0->next;
		if (t0->dfs_parent) continue;
		t = t0;
		e = (t->first[0]) ? t->first[0] : t->first[1];
		t->dfs_parent = (TreeEdge*)trees;
		while ( 1 )
		{
			if (e == NULL)
			{
				e = t->dfs_parent;
				t->dfs_parent = (TreeEdge*)((char*)trees + component_num);
				if (t == t0) break;
				if (t == e->head[0]) { t = e->head[1]; e = (e->next[0]) ? e->next[0] : t->first[1]; }
				else                 { t = e->head[0]; e = e->next[1]; }
				continue;
			}

			if (e->head[1] == t)
			{
				if (e->head[0]->dfs_parent || !(q=e->pq01[1].GetMin()) || q->slack - e->head[0]->eps + t->eps > 0) { e = (e->next[0]) ? e->next[0] : t->first[1]; continue; }
				t = e->head[0];
			}
			else
			{
				if (e->head[1]->dfs_parent || !(q=e->pq01[0].GetMin()) || q->slack - e->head[1]->eps + t->eps > 0) { e = e->next[1]; continue; }
				t = e->head[1];
			}
			t->dfs_parent = e;
			e = (t->first[0]) ? t->first[0] : t->first[1];
		}
		component_num ++;
	}

	Tree** array = new Tree*[component_num];
	memset(array, 0, component_num*sizeof(Tree*));

	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		t = r->tree;
		t->id = (int)((char*)t->dfs_parent - (char*)trees);
		t->next = array[t->id];
		array[t->id] = t;
	}

	for (c=component_num-1; c>=0; c--)
	{
		eps = PM_INFTY;
		for (t=array[c]; t; t=t->next)
		{
			if (eps > t->eps_delta) eps = t->eps_delta;
			FOR_ALL_TREE_EDGES(t, e, dir)
			{
				t2 = e->head[dir];
				REAL eps00 = (q=e->pq00.GetMin()) ? (q->slack-t->eps-t2->eps) : PM_INFTY;
				REAL eps01[2];
				eps01[dir]   = ((q=e->pq01[dir].GetMin()))   ? (q->slack - t->eps + t2->eps) : PM_INFTY;
				eps01[1-dir] = ((q=e->pq01[1-dir].GetMin())) ? (q->slack - t2->eps + t->eps) : PM_INFTY;
				if (t2->id < c)
				{
					if (eps > eps01[dir]) eps = eps01[dir];
					if (eps > eps00) eps = eps00;
				}
				else if (t2->id == c)
				{
					if (2*eps > eps00) eps = eps00 / 2;
				}
				else
				{
					if (eps > eps01[dir] + t2->eps_delta) eps = eps01[dir] + t2->eps_delta;
					if (eps > eps00 - t2->eps_delta)  eps = eps00 - t2->eps_delta;
				}
			}
		}
		for (t=array[c]; t; t=t->next) t->eps_delta = eps;
	}

	delete [] array;
	//printf("%d SCCs ", component_num);
}

void PerfectMatching::CommitEps()
{
	printf("CommitEps()\n");
	Node* i;
	Node* j;
	Node* r;
	int dir;
	Edge* a;
	EdgeIterator I;
	Tree* t;
	TreeEdge* e;
	TreeEdge** e_ptr;
	REAL eps, eps2;
	PriorityQueue<REAL>::Item* q;

	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		t = r->tree;
		eps = t->eps;

		i = r;
		while ( 1 )
		{
			i->y += eps;
			if (!i->is_tree_root)
			{
				Node* i0 = i;
				i = ARC_HEAD(i0->match);
				if (i->is_blossom) ARC_TO_EDGE_PTR(i0->match)->slack -= eps;
				else               i->y                              -= eps;
				FOR_ALL_EDGES(i, a, dir, I)
				{
					GET_OUTER_HEAD(a, dir, j);

					a->slack += eps;
					if (j->flag == 0) a->slack -= j->tree->eps;
				}
				i = i0;
			}

			MOVE_NODE_IN_TREE(i);
		}

		t->pq0.Update(-eps);

		PriorityQueue<REAL> pq00 = t->pq00;
		t->pq00.Reset();
		for (q=pq00.GetAndResetFirst(); q; q=pq00.GetAndResetNext())
		{
			a = (Edge*)q;
			if (ProcessEdge00(a)) t->pq00.Add(a);
		}

		for (e_ptr=&t->first[0], e=*e_ptr; e; e=*e_ptr)
		{
			if (e->head[0] == NULL) { *e_ptr = e->next[0]; tree_edges->Delete(e); continue; }
			e_ptr = &e->next[0];

			eps2 = e->head[0]->eps;
			e->pq00.Update( - eps - eps2 );
		}
	}

	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next) r->tree->eps = 0;
}

bool PerfectMatching::UpdateDuals()
{
	Node* r;

	double start_time = get_time();

	////////////////////////////////////////////////////////////////////////////////////
	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		Tree* t = r->tree;
		PriorityQueue<REAL>::Item* q;
		REAL eps = PM_INFTY;
		if ((q=t->pq0.GetMin())) eps = q->slack;
		if ((q=t->pq_blossoms.GetMin()) && eps > q->slack) eps = q->slack;
		while ((q=t->pq00.GetMin()))
		{
			if (ProcessEdge00((Edge*)q, false)) break;
			t->pq00.Remove(q, pq_buf);
		}
		if (q && 2*eps > q->slack) eps = q->slack/2;
		t->eps_delta = eps - t->eps;
	}

	if (tree_num >= options.dual_LP_threshold*node_num)
	{
		if      (options.dual_greedy_update_option == 0) ComputeEpsCC();
		else if (options.dual_greedy_update_option == 1) ComputeEpsSCC();
		else                                             ComputeEpsSingle();
	}
	else ComputeEpsGlobal();

	REAL delta = 0;
	for (r=nodes[node_num].tree_sibling_next; r; r=r->tree_sibling_next)
	{
		if (r->tree->eps_delta > 0)
		{
			delta += r->tree->eps_delta;
			r->tree->eps += r->tree->eps_delta;
		}
	}

	stat.dual_time += get_time() - start_time;

	return (delta > PM_THRESHOLD);
}
