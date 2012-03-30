#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "GeomPerfectMatching.h"
#include "GPMkdtree.h"




// greedy procedure to make sure that a perfect matching exists:
// 1. construct a matching among existing edges
//    (greedy procedure: pick a node, check whether there are edges leading
//     to unmatched nodes, if there are pick the edge with the smallest length).
// 2. take remaining unmatched nodes, construct kd-tree for them,
//    assign an ordering to nodes (last visited time during left-most depth-first search),
//    add edges between consecutive nodes (2*i,2*i+1)
void GeomPerfectMatching::CompleteInitialMatching()
{
	if (options.verbose) printf("adding edges to make sure that a perfect matching exists...");
	PointId p, q;
	Edge* e;
	double len, len_min;
	int unmatched_num = 0, edge_num0 = edge_num;

	// construct greedy matching
	for (p=0; p<node_num; p++)
	{
		if (nodes[p].is_marked) continue;
		q = -1;
		for (e=nodes[p].first[0]; e; e=e->next[0])
		{
			if (nodes[e->head[0]].is_marked) continue;
			len = Dist2(p, e->head[0]);
			if (q < 0 || len_min > len)
			{
				q = e->head[0];
				len_min = len;
			}
		}
		if (q >= 0)
		{
			nodes[p].is_marked = nodes[q].is_marked = 1;
		}
		else unmatched_num ++;
	}

	if (unmatched_num == 0)
	{
		for (p=0; p<node_num; p++) nodes[p].is_marked = 0;
		return;
	}

	//printf("%d unmatched\n", unmatched_num);

	REAL* unmatched_coords = new REAL[unmatched_num*DIM];
	int* rev_mapping = new int[unmatched_num];
	unmatched_num = 0;
	for (p=0; p<node_num; p++)
	{
		if (nodes[p].is_marked) nodes[p].is_marked = 0;
		else
		{
			memcpy(unmatched_coords+unmatched_num*DIM, coords+p*DIM, DIM*sizeof(REAL));
			rev_mapping[unmatched_num ++] = p;
		}
	}

	GPMKDTree* kd_tree = new GPMKDTree(DIM, unmatched_num, unmatched_coords, this);
	kd_tree->AddPerfectMatching(rev_mapping);

	delete kd_tree;
	delete [] unmatched_coords;
	delete [] rev_mapping;

	if (options.verbose) printf("done (%d edges)\n", edge_num-edge_num0);
}

void GeomPerfectMatching::InitKNN(int K)
{
	if (node_num != node_num_max) { printf("InitKNN() cannot be called before all points have been added!\n"); exit(1); }
	if (options.verbose) printf("adding K nearest neighbors (K=%d)\n", K);

	int dir, k;
	PointId p;
	Edge* e;

	if (K > node_num - 1) K = node_num - 1;

	GPMKDTree* kd_tree = new GPMKDTree(DIM, node_num, coords, this);
	PointId* neighbors = new PointId[K];

	for (p=0; p<node_num; p++)
	{
		for (dir=0; dir<2; dir++)
		for (e=nodes[p].first[dir]; e; e=e->next[dir])
		{
			nodes[e->head[dir]].is_marked = 1;
		}

		kd_tree->ComputeKNN(p, K, neighbors);
		for (k=0; k<K; k++)
		{
			if (nodes[neighbors[k]].is_marked) continue;
			AddInitialEdge(p, neighbors[k]);
			nodes[neighbors[k]].is_marked = 1;
		}

		for (dir=0; dir<2; dir++)
		for (e=nodes[p].first[dir]; e; e=e->next[dir])
		{
			nodes[e->head[dir]].is_marked = 0;
		}
	}

	delete kd_tree;
	delete [] neighbors;
}

#ifdef DELAUNAY_TRIANGLE

#ifdef _MSC_VER
#pragma warning(disable: 4311)
#pragma warning(disable: 4312)
#endif

extern "C" {
#define ANSI_DECLARATORS
#define TRILIBRARY
#define NO_TIMER
#define main NO_MAIN_FUNCTION
#include "../triangle/triangle.c"
}

void GeomPerfectMatching::InitDelaunay()
{
	if (node_num < 16) return;
	if (options.verbose) printf("adding edges in Delaunay triangulation\n");

	int k;

	struct triangulateio in, out, vorout;
	in.numberofpoints = node_num;
	in.numberofpointattributes = 0;
	in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
	for (k=0; k<2*node_num; k++) in.pointlist[k] = coords[k];
	in.pointattributelist = NULL;
	in.pointmarkerlist = NULL;
	in.numberofsegments = 0;
	in.numberofholes = 0;
	in.numberofregions = 0;
	in.regionlist = 0;

	out.pointlist = (REAL *) NULL;
	out.pointattributelist = (REAL *) NULL;
	out.pointmarkerlist = (int *) NULL;
	out.trianglelist = (int *) NULL;
	out.triangleattributelist = (REAL *) NULL;
	out.neighborlist = (int *) NULL;
	out.segmentlist = (int *) NULL;
	out.segmentmarkerlist = (int *) NULL;
	out.edgelist = (int *) NULL;
	out.edgemarkerlist = (int *) NULL;

	vorout.pointlist = (REAL *) NULL;
	vorout.pointattributelist = (REAL *) NULL;
	vorout.edgelist = (int *) NULL;
	vorout.normlist = (REAL *) NULL;

	triangulate("pczAevn", &in, &out, &vorout);

	free(in.pointlist);
	free(out.pointlist);
	free(out.pointmarkerlist);
	free(out.trianglelist);
	free(out.neighborlist);
	free(out.segmentlist);
	free(out.segmentmarkerlist);
	free(out.edgemarkerlist);
	free(vorout.pointlist);
	free(vorout.pointattributelist);
	free(vorout.edgelist);
	free(vorout.normlist);

	for (k=0; k<out.numberofedges; k++) AddInitialEdge(out.edgelist[2*k], out.edgelist[2*k+1]);

	free(out.edgelist);
}



#else

void GeomPerfectMatching::InitDelaunay()
{
	printf("You need to download the 'Triangle' software from \n\thttp://www.cs.cmu.edu/~quake/triangle.html ,\nextract it to the directory GeomPerfectMatching and define DELAUNARY_TRIANG in GeomPerfectMatching.h\n");
	exit(1);
}

#endif

