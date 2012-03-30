#include <stdio.h>
#include "PerfectMatching.h"
#include "LCA.h"

struct Node
{
	PerfectMatching::REAL sum; // = twice_y[i] + twice_y[i->parent] + twice_y[i->parent->parent] + ...
	Node* match;
	Node* parent;
	Node* child;
	Node* sibling;
	int lca_preorder;
};


int CheckPerfectMatchingOptimality(int node_num, int edge_num, int* edges, int* weights, PerfectMatching* pm, PerfectMatching::REAL threshold)
{
	int _i, _j, _e;
	Node* i;
	Node* j;
	int blossom_num = pm->GetBlossomNum();
	int* blossom_parents = new int[node_num+blossom_num];
	PerfectMatching::REAL* twice_y = new PerfectMatching::REAL[node_num+blossom_num];

	PerfectMatching::REAL y_blossom_min = 0;
	PerfectMatching::REAL slack_min = 0;
	PerfectMatching::REAL active_slack_max = 0;

	// step 1 - read dual solution and construct tree
	pm->GetDualSolution(blossom_parents, twice_y);
	Node* nodes = new Node[node_num+blossom_num+1];
	memset(nodes, 0, (node_num+blossom_num+1)*sizeof(Node));
	Node* ROOT = nodes+node_num+blossom_num;
	for (_i=0, i=nodes; _i<node_num+blossom_num; _i++, i++)
	{
		i->sum = twice_y[_i];
		if (_i >= node_num && y_blossom_min > i->sum) y_blossom_min = i->sum;
		if (blossom_parents[_i] >= 0)
		{
			if (blossom_parents[_i]<node_num || blossom_parents[_i]>=node_num+blossom_num)
			{
				delete [] nodes;
				delete [] blossom_parents;
				delete [] twice_y;
				return 2;
			}
			i->parent = nodes + blossom_parents[_i];
			i->sibling = i->parent->child;
			i->parent->child = i;
		}
	}
	delete [] blossom_parents;
	delete [] twice_y;

	for (i=nodes; i<nodes+node_num+blossom_num; i++)
	{
		if (!i->parent)
		{
			i->parent = ROOT;
			i->sibling = ROOT->child;
			ROOT->child = i;
		}
	}

	LCATree* lca_tree = new LCATree(node_num+blossom_num+1);
	Node** rev_mapping = new Node*[node_num+blossom_num];

	i = ROOT;
	while ( 1 )
	{
		if (i->child)
		{
			if (i < nodes+node_num) { delete [] nodes; delete lca_tree; delete [] rev_mapping; return 2; }
			i->child->sum += i->sum;
			i = i->child;
		}
		else
		{
			if (i >= nodes+node_num) { delete [] nodes; delete lca_tree; delete [] rev_mapping; return 2; }
			while ( 1 )
			{
				i->lca_preorder = lca_tree->Add(i, i->parent);
				rev_mapping[i->lca_preorder] = i;
				if (i->sibling) break;
				i = i->parent;
				if (i == ROOT)
				{
					i->lca_preorder = lca_tree->AddRoot(i);
					break;
				}
			}
			if (i == ROOT) break;
			i = i->sibling;
			i->sum += i->parent->sum;
		}
	}

	int matched_num = 0;
	for (_e=0; _e<edge_num; _e++)
	{
		_i = edges[2*_e];
		_j = edges[2*_e+1];
		if (_i<0 || _j<0 || _i>=node_num || _j>=node_num || _i==_j) { delete [] nodes; delete lca_tree; delete [] rev_mapping; return 2; }

		int lca_i = nodes[_i].lca_preorder;
		int lca_j = nodes[_j].lca_preorder;
		lca_tree->GetPenultimateNodes(lca_i, lca_j);
		i = rev_mapping[lca_i];
		j = rev_mapping[lca_j];
		PerfectMatching::REAL twice_slack = 2*weights[_e] - (nodes[_i].sum - i->parent->sum) - (nodes[_j].sum - j->parent->sum);
		if (slack_min > twice_slack) slack_min = twice_slack;
		if (pm->GetSolution(_e))
		{
			if (pm->GetMatch(_i)!=_j || pm->GetMatch(_j)!=_i || i->match || j->match) { delete [] nodes; delete lca_tree; delete [] rev_mapping; return 2; }
			i->match = j;
			j->match = i;
			if (active_slack_max < twice_slack) active_slack_max = twice_slack;
			matched_num += 2;
		}
	}

	delete [] nodes;
	delete lca_tree;
	delete [] rev_mapping;

	if (matched_num != node_num) return 2;

	if (y_blossom_min < -threshold || slack_min < -threshold || active_slack_max > threshold)
	{
		printf("ERROR in CheckPerfectMatchingOptimality():\n");
		if ( ((PerfectMatching::REAL)1 / 2) == 0 )
			printf("\ty_blossom_min=%d\n\tslack_min=%d\n\tactive_slack_max=%d\n", (int)y_blossom_min, (int)slack_min, (int)active_slack_max);
		else
			printf("\ty_blossom_min=%.15f\n\tslack_min=%.15f\n\tactive_slack_max=%.15f\n", (double)y_blossom_min, (double)slack_min, (double)active_slack_max);
		return 1;
	}

	return 0;
}

double ComputePerfectMatchingCost(int node_num, int edge_num, int* edges, int* weights, PerfectMatching* pm)
{
	int i;
	int j;
	int e;
	double cost = 0;

	int* nodes = new int[node_num];
	memset(nodes, 0, node_num*sizeof(int));
	for (e=0; e<edge_num; e++)
	{
		if (pm->GetSolution(e))
		{
			i = edges[2*e];
			j = edges[2*e+1];
			nodes[i] ++;
			nodes[j] ++;
			cost += weights[e];
		}
	}
	for (i=0; i<node_num; i++)
	{
		if (nodes[i] != 1)
		{
			printf("ComputeCost(): degree = %d!\n", nodes[i]);
			exit(1);
		}
	}
	delete [] nodes;
	return cost;
}

