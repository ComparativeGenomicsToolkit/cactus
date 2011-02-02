#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "GeomPerfectMatching.h"
#include "GPMkdtree.h"
#include "../timer.h"



GeomPerfectMatching::REAL GeomPerfectMatching::SolveComplete()
{
	if (node_num != node_num_max) { printf("ComputeCost() cannot be called before all points have been added!\n"); exit(1); }

	PointId p, q;
	int e = 0, E = node_num*(node_num-1)/2;
	PerfectMatching* pm = new PerfectMatching(node_num, E);
	for (p=0; p<node_num; p++)
	{
		for (q=p+1; q<node_num; q++)
		{
			pm->AddEdge(p, q, Dist(p, q));
		}
	}
	pm->options = options;
	pm->Solve();
	for (p=0; p<node_num; p++)
	{
		for (q=p+1; q<node_num; q++)
		{
			if (pm->GetSolution(e++))
			{
				matching[p] = q;
				matching[q] = p;
			}
		}
	}
	delete pm;
	return ComputeCost(matching);
}

GeomPerfectMatching::REAL GeomPerfectMatching::Solve()
{
	double start_time = get_time();
	double perfect_matching_time = 0;
	double negative_edges_time = 0;
	if (options.verbose) { printf("starting geometric matching with %d points\n", node_num); fflush(stdout); }
	PointId p, q;
	Edge* e;
	int _e;
	int iter;
	bool success = false;
	PerfectMatching* pm = NULL;
	GPMKDTree* kd_tree;

	double init_matching_time = get_time();

	if (gpm_options.init_Delaunay) InitDelaunay();
	if (gpm_options.init_KNN > 0) InitKNN(gpm_options.init_KNN);
	if (gpm_options.init_greedy) CompleteInitialMatching();

	init_matching_time = get_time() - init_matching_time;

	graph_update_time = 0;

	int iter_max = gpm_options.iter_max;
	for (iter=0; iter_max<=0 || iter<iter_max; iter++)
	{
		if (pm)
		{
			double negative_edges_start_time = get_time();
			int edge_num0 = edge_num;
			pm->StartUpdate();
			for (p=0; p<node_num; p++)
			{
				PerfectMatching::REAL s = pm->GetTwiceSum(p);
				if ( ((REAL)1 / 2) == 0 && ((PerfectMatching::REAL)1 / 2) != 0 ) sums[p] = (REAL)ceil((double)s);
				else                                                             sums[p] = (REAL)s;
			}
			if (options.verbose) { printf("building kd_tree..."); fflush(stdout); }
			{
				kd_tree = new GPMKDTree(DIM+1, node_num, coords, this);
			}
			if (options.verbose) { printf(" done. Now adding negative edges:\n    "); fflush(stdout); }

			for (p=0; p<node_num; p++)
			{
				if (options.verbose && (p%(node_num/72)==0)) { printf("+"); fflush(stdout); }
				for (e=nodes[p].first[0]; e; e=e->next[0]) nodes[e->head[0]].is_marked = 1;
				for (e=nodes[p].first[1]; e; e=e->next[1]) nodes[e->head[1]].is_marked = 1;
				kd_tree->AddNegativeEdges(p, pm);
				for (e=nodes[p].first[0]; e; e=e->next[0]) nodes[e->head[0]].is_marked = 0;
				for (e=nodes[p].first[1]; e; e=e->next[1]) nodes[e->head[1]].is_marked = 0;
			}
			delete kd_tree;
			//if (edge_num - edge_num0 > node_num / 32)
			if ( 0 ) // always reuse previous computation
			{
				delete pm;
				pm = NULL;
			}
			else
			{
				pm->FinishUpdate();
				if (edge_num0 == edge_num) success = true;
			}
			if (options.verbose) { printf("\ndone (%d edges added)\n", edge_num-edge_num0); fflush(stdout); }
			negative_edges_time += get_time() - negative_edges_start_time;
		}
		if (!pm)
		{
			int E = 5*node_num;
			if (E < 5*edge_num/4) E = 5*edge_num/4;
			pm = new PerfectMatching(node_num, E);
			for (e=edges->ScanFirst(); e; e=edges->ScanNext())
			{
				p = e->head[1]; q = e->head[0];
				pm->AddEdge(p, q, Dist(p, q));
			}
		}
		if (options.verbose) printf("iter %d: ", iter+1);
		pm->options = options;
		double perfect_matching_start = get_time();
		pm->Solve();
		perfect_matching_time += get_time() - perfect_matching_start;
		if (success) break;
	}

	for (_e=0, e=edges->ScanFirst(); e; _e++, e=edges->ScanNext())
	{
		if (pm->GetSolution(_e))
		{
			p = e->head[1]; q = e->head[0];
			matching[p] = q;
			matching[q] = p;
		}
	}
	delete pm;
	REAL cost = ComputeCost(matching);
	if (options.verbose)
	{
		printf("geometric matching finished [%.3f secs]. cost=%.1f \n", get_time()-start_time, (double)cost); 
		printf("    selecting initial edges: [%.3f secs], perfect matching: [%.3f secs]\n", init_matching_time, perfect_matching_time);
		printf("    pricing: [%.3f secs] including graph updates: [%.3f secs]\n", negative_edges_time, graph_update_time); 
		fflush(stdout);
	}
	return cost;
}

