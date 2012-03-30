#include <stdio.h>
#include "PerfectMatching.h"
#include "GEOM/GeomPerfectMatching.h"

void LoadFile(int& node_num, int& edge_num, int*& edges, int*& weights, char* filename)
{
	int e = 0;
	char LINE[1000];
	FILE* fp = fopen(filename, "r");
	if (!fp) { printf("Can't open %s\n", filename); exit(1); }

	int format = -1; // 0: DIMACS format. node id's start from 1
	                 // 1: simpler format (without "p" and "e"). node id's start from 0

	edge_num = -1;
	while (fgets(LINE, sizeof(LINE)-1, fp))
	{
		if (LINE[0] == 'c') continue;
		if (format < 0)
		{
			if (LINE[0] == 'p')
			{
				format = 0;
				if (sscanf(LINE, "p edge %d %d\n", &node_num, &edge_num) != 2) { printf("%s: wrong format #1\n", filename); exit(1); }
			}
			else
			{
				format = 1;
				if (sscanf(LINE, "%d %d\n", &node_num, &edge_num) != 2) { printf("%s: wrong format #1\n", filename); exit(1); }
			}

			//////////////////////////////////////////////////////////////////////////////////
			if (node_num <= 0 || edge_num < 0) { printf("# of nodes and edges should be positive\n"); exit(1); }
			if (node_num & 1) { printf("# of nodes is odd: perfect matching cannot exist\n"); exit(1); }
			edges = new int[2*edge_num];
			weights = new int[edge_num];
			//////////////////////////////////////////////////////////////////////////////////
		}
		else
		{
			int i, j;
			char* ptr = LINE;
			if (format == 0) { if (LINE[0] != 'e') continue; ptr = &LINE[1]; }
			else             ptr = &LINE[0];

			int len;
			if (sscanf(ptr, "%d %d %d\n", &i, &j, &len) != 3) continue;
			if (format == 0) { i --; j --; }
			edges[2*e] = i;
			edges[2*e+1] = j;
			weights[e] = len;
			e ++;
		}
	}

	if (e != edge_num) { printf("%s: wrong format #3\n", filename); exit(1); }
	fclose(fp);
}

void SaveMatching(int node_num, PerfectMatching* pm, char* filename)
{
	FILE* fp = fopen(filename, "w");
	if (!fp) { printf("Can't open %s\n", filename); exit(1); }
	fprintf(fp, "%d %d\n", node_num, node_num/2);
	int i, j;
	for (i=0; i<node_num; i++)
	{
		j = pm->GetMatch(i);
		if (i < j) fprintf(fp, "%d %d\n", i, j);
	}
	fclose(fp);
}

void LoadGeomFile(int& node_num, int*& x_array, int*& y_array, char* filename)
{
	int i = 0, i_tmp, x, y, DIM = 0;
	char LINE[1000];

	x_array = y_array = NULL;

	FILE* fp = fopen(filename, "r");
	if (!fp) { printf("Can't open %s\n", filename); exit(1); }

	while (fgets(LINE, sizeof(LINE)-1, fp))
	{
		if (sscanf(LINE, "DIMENSION : %d", &node_num) == 1)
		{
			if (node_num < 1) { printf("too few nodes\n"); exit(1); }
			if (node_num & 1) { printf("# of points is odd: perfect matching cannot exist\n"); exit(1); }
			if (x_array) { printf("wrong format\n"); exit(1); }
			x_array = new int[node_num];
			y_array = new int[node_num];
			continue;
		}
		if (sscanf(LINE, "EDGE_WEIGHT_TYPE : EUC_%dD", &DIM) == 1)
		{
			if (DIM != 2) { printf("only EUC_2D is supported"); exit(1); }
			continue;
		}
		if (sscanf(LINE, "%d %d %d", &i_tmp, &x, &y) == 3)
		{
			i_tmp --;
			if (i_tmp != i ++ || i > node_num) { printf("wrong number of points\n"); exit(1); }
			x_array[i_tmp] = x;
			y_array[i_tmp] = y;
			continue;
		}
		printf("%s", LINE);
	}
	fclose(fp);
	if (i != node_num || !x_array || DIM != 2) { printf("wrong format\n"); exit(1); }
}


void SaveMatching(int node_num, GeomPerfectMatching* gpm, char* filename)
{
	FILE* fp = fopen(filename, "w");
	if (!fp) { printf("Can't open %s\n", filename); exit(1); }
	fprintf(fp, "%d %d\n", node_num, node_num/2);
	int i, j;
	for (i=0; i<node_num; i++)
	{
		j = gpm->GetMatch(i);
		if (i < j)
		{
			GeomPerfectMatching::REAL len = gpm->Dist(i, j);
			if ( ((GeomPerfectMatching::REAL)1 / 2) == 0 ) fprintf(fp, "%d %d %d\n", i, j, len);
			else                                           fprintf(fp, "%d %d %f\n", i, j, (double)len);
		}
	}
	fclose(fp);
}


void ShowUsage()
{
	printf("Usage: see USAGE.TXT\n");
	exit(1);
}


int main(int argc, char* argv[])
{
	struct PerfectMatching::Options options;
	struct GeomPerfectMatching::GPMOptions gpm_options;
	char* filename = NULL;
	char* geom_filename = NULL;
	char* save_filename = NULL;
	bool check_perfect_matching = false;
	int i, e, node_num, edge_num;
	int* edges;
	int* weights;

	for (i=1; i<argc; i++)
	{
		if (argv[i][0] != '-') { printf("Unknown option: %s\n", argv[i]); ShowUsage(); }
		switch (argv[i][1])
		{
			case 'e':
				if (filename || argv[i][2] || ++i == argc) ShowUsage();
				filename = argv[i];
				break;
			case 'g':
				if (geom_filename || argv[i][2] || ++i == argc) ShowUsage();
				geom_filename = argv[i];
				break;
			case 'j':
				if (argv[i][2]) ShowUsage();
				options.fractional_jumpstart = false;
				break;
			case 'm':
				options.dual_LP_threshold = atof(&argv[i][2]);
				if (options.dual_LP_threshold<0 || options.dual_LP_threshold>1) ShowUsage();
				break;
			case 'd':
				options.dual_greedy_update_option = atoi(&argv[i][2]);
				if (options.dual_greedy_update_option<1 || options.dual_greedy_update_option>2) ShowUsage();
				break;
			case 'b':
				if (argv[i][2]) ShowUsage();
				options.update_duals_before = true;
				break;
			case 'a':
				if (argv[i][2]) ShowUsage();
				options.update_duals_after = true;
				break;
			case 'D':
				if (argv[i][2]) ShowUsage();
				gpm_options.init_Delaunay = false;
				break;
			case 'K':
				gpm_options.init_KNN = atoi(&argv[i][2]);
				if (gpm_options.init_KNN<0) ShowUsage();
				break;
			case 'I':
				if (argv[i][2]) ShowUsage();
				gpm_options.init_greedy = false;
				break;
			case 'T':
				gpm_options.iter_max = atoi(&argv[i][2]);
				if (gpm_options.iter_max<0) ShowUsage();
				break;
			case 'w':
				if (save_filename || argv[i][2] || ++i == argc) ShowUsage();
				save_filename = argv[i];
				break;
			case 'c':
				if (argv[i][2]) ShowUsage();
				check_perfect_matching = true;
				break;
			case 'V':
				if (argv[i][2]) ShowUsage();
				options.verbose = false;
				break;
			default:
				printf("Unknown option: %s\n", argv[i]);
				ShowUsage();
				break;
		}
		
	}

	if (!filename && !geom_filename) ShowUsage();

	if (filename) LoadFile(node_num, edge_num, edges, weights, filename);

	if (!geom_filename)
	{
		PerfectMatching *pm = new PerfectMatching(node_num, edge_num);
		for (e=0; e<edge_num; e++) pm->AddEdge(edges[2*e], edges[2*e+1], weights[e]);
		pm->options = options;
		pm->Solve();
		if (check_perfect_matching)
		{
			int res = CheckPerfectMatchingOptimality(node_num, edge_num, edges, weights, pm);
			printf("check optimality: res=%d (%s)\n", res, (res==0) ? "ok" : ((res==1) ? "error" : "fatal error"));
		}
		double cost = ComputePerfectMatchingCost(node_num, edge_num, edges, weights, pm);
		printf("cost = %.1f\n", cost);
		if (save_filename) SaveMatching(node_num, pm, save_filename);
		delete pm;
	}
	else
	{
		int geom_node_num;
		int* x_array;
		int* y_array;
		LoadGeomFile(geom_node_num, x_array, y_array, geom_filename);
		GeomPerfectMatching *gpm = new GeomPerfectMatching(geom_node_num, 2);
		for (i=0; i<geom_node_num; i++)
		{
			GeomPerfectMatching::REAL coord[2];
			coord[0] = x_array[i];
			coord[1] = y_array[i];
			gpm->AddPoint(coord);
		}
		delete [] x_array;
		delete [] y_array;
		gpm->options = options;
		gpm->gpm_options = gpm_options;
		if (filename)
		{
			if (node_num != geom_node_num) { printf("%s and %s don't match!\n", geom_filename, filename); exit(1); }
			for (e=0; e<edge_num; e++)
			{
				if (weights[e] != gpm->Dist(edges[2*e], edges[2*e+1]))
				{ printf("edge lengths in %s and %s don't match!\n", geom_filename, filename); exit(1); }
				gpm->AddInitialEdge(edges[2*e], edges[2*e+1]);
			}
		}
		gpm->Solve();
		if (save_filename) SaveMatching(geom_node_num, gpm, save_filename);
		delete gpm;
	}

	if (filename)
	{
		delete [] edges;
		delete [] weights;
	}

	return 0;
}


