#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "GPMkdtree.h"
#include "../timer.h"


// 'mapping' must be of size 2*N. (array[N], ... array[2*N-1] is used as a temporary buffer).
// After the call  array[mapping[0]] <= array[mapping[1]] <= ... <= array[mapping[N-1]].
// array is not modified.
template <typename Type> inline void sort(Type* array, int array_skip, int* mapping, int N)
{
	// mergesort
	int i;
	int* mappingSrc = mapping;
	int* mappingDst = mapping + N;
	int* pSrc1;
	int* pSrc2;
	int* pSrc1End;
	int* pSrc2End;
	int* pDst;

	for (i=0; i<(N&(~1)); i+=2)
	{
		if (array[array_skip*i] < array[array_skip*(i+1)])
		{
			mappingSrc[i]   = i;
			mappingSrc[i+1] = i+1;
		}
		else
		{
			mappingSrc[i]   = i+1;
			mappingSrc[i+1] = i;
		}
	}
	if (i != N) mappingSrc[i] = i;

	int step;
	for (step=2; step<N; step*=2)
	{
		pSrc2End = mappingSrc;
		pDst = mappingDst;
		while ( 1 )
		{
			pSrc1 = pSrc2End;
			pSrc1End = pSrc1 + step;
			if (pSrc1End >= mappingSrc + N)
			{
				memcpy(pDst, pSrc1, (int)((char*)(mappingSrc + N) - (char*)pSrc1));
				break;
			}
			pSrc2 = pSrc1End;
			pSrc2End = pSrc2 + step;
			if (pSrc2End > mappingSrc + N) pSrc2End = mappingSrc + N;
			while ( 1 )
			{
				if (array[(array_skip)*(*pSrc1)] < array[array_skip*(*pSrc2)])
				{
					*pDst ++ = *pSrc1 ++;
					if (pSrc1 == pSrc1End)
					{
						memcpy(pDst, pSrc2, (int)((char*)pSrc2End - (char*)pSrc2));
						pDst = (int*) ((char*)pDst + (int)((char*)pSrc2End - (char*)pSrc2));
						break;
					}
				}
				else
				{
					*pDst ++ = *pSrc2 ++;
					if (pSrc2 == pSrc2End)
					{
						memcpy(pDst, pSrc1, (int)((char*)pSrc1End - (char*)pSrc1));
						pDst = (int*) ((char*)pDst + (int)((char*)pSrc1End - (char*)pSrc1));
						break;
					}
				}
			}
		}
		pDst = mappingDst;
		mappingDst = mappingSrc;
		mappingSrc = pDst;
	}
	if (mappingSrc != mapping) memcpy(mapping, mappingSrc, N*sizeof(int));
}

//////////////////////////////////////////////////////////////////////////////////////////

#define NEIGHBOR_PARENT(k) (((k)-1)>>1)
#define NEIGHBOR_FIRST_CHILD(k) (((k)<<1)+1)

Neighbors::Neighbors()
{
	K_max = 0;
	dist_array = NULL;
}

Neighbors::~Neighbors()
{
	if (dist_array) delete [] dist_array;
}

void Neighbors::Init(int _K, PointId* _array)
{
	K = _K;
	array = _array;
	num = 0;
	if (K > K_max)
	{
		if (dist_array) delete [] dist_array;
		K_max = K;
		dist_array = new double[K_max];
	}
}

inline void Neighbors::Swap(int k1, int k2)
{
	PointId p = array[k1]; array[k1] = array[k2]; array[k2] = p;
	double d = dist_array[k1]; dist_array[k1] = dist_array[k2]; dist_array[k2] = d;
}

inline void Neighbors::Add(PointId p, double dist)
{
	int k;
	if (num < K)
	{
		k = num ++;
		array[k] = p;
		dist_array[k] = dist;
		while ( k > 0 )
		{
			int k_parent = NEIGHBOR_PARENT(k);
			if (dist_array[k] <= dist_array[k_parent]) break;
			Swap(k, k_parent);
			k = k_parent;
		}
	}
	else
	{
		if (dist_array[0] <= dist) return;
		array[0] = p;
		dist_array[0] = dist;
		k = 0;
		while ( 1 )
		{
			int k_child = NEIGHBOR_FIRST_CHILD(k);
			if (k_child >= K) break;
			if (k_child+1 < K && dist_array[k_child+1] > dist_array[k_child]) k_child ++;
			if (dist_array[k] >= dist_array[k_child]) break;
			Swap(k, k_child);
			k = k_child;
		}
	}
	//for (k=1; k<num; k++) assert(dist_array[k] <= dist_array[NEIGHBOR_PARENT(k)]);
}

inline double Neighbors::GetMax()
{
	assert(num > 0);
	return dist_array[0];
}

//////////////////////////////////////////////////////////////////////////////////////////

GPMKDTree::GPMKDTree(int _D, int _point_num, REAL* coords, GeomPerfectMatching* _GPM)
	: D(_D), DIM(_GPM->DIM), point_num(_point_num), GPM(_GPM)
{
	Node* i;
	Node* j;
	int d, d0, k;
	int* mapping = new int[(D+2)*point_num];
	int* buf = mapping + D*point_num;
	int* marking = buf + point_num;
	memset(marking, 0, point_num*sizeof(int));
	int* ptr = mapping;

	int node_num_max = 4*point_num/3+2;
	nodes = (Node*)malloc(node_num_max*sizeof(Node));
	rev_mapping = (Node**)malloc(point_num*sizeof(Node*));
	memset(rev_mapping, 0, point_num*sizeof(Node*));

	REAL** coords_array = new REAL*[D];
	int* skip_array = new int[D];
	for (d=0; d<DIM; d++) { coords_array[d] = coords + d; skip_array[d] = DIM; }
	if (d < D)            { coords_array[d] = GPM->sums; skip_array[d] = 1; }

	for (d=0; d<D; d++)
	{
		sort<REAL>(coords_array[d], skip_array[d], ptr, point_num);
		if (d == DIM) sum_max = GPM->sums[ptr[point_num-1]];
		ptr += point_num;
	}

	nodes[0].parent = NULL;
	nodes[0].order = 0;
	nodes[0].d = point_num;
	nodes[0].first_child = NULL;
	node_num = 1;
	Node* first_unprocessed = &nodes[0];
	while ( (i=first_unprocessed) )
	{
		first_unprocessed = i->first_child;

		int start = i->order;
		int num0 = i->d, num;
		if ((DIM==D && num0<=2) || (DIM<D && num0 <= CHILDREN_MAX))
		{
			// leaf
			i->d = -num0;
			for (k=0; k<num0; k++)
			{
				i->points[k] = mapping[start+k];
				rev_mapping[mapping[start+k]] = i;
			}
			continue;
		}

		// not a leaf.
		if (node_num + 2 > node_num_max)
		{
			node_num_max = 3*node_num_max/2 + 16;
			Node* nodes_old = nodes;
			nodes = (Node*)realloc(nodes, node_num_max*sizeof(Node));
#define UPDATE_NODE_PTR(ptr) ptr = (Node*)((char*)ptr + ((char*)nodes-(char*)nodes_old))
			UPDATE_NODE_PTR(i);
			if (first_unprocessed) UPDATE_NODE_PTR(first_unprocessed);
			for (k=0; k<node_num; k++)
			{
				if (nodes[k].parent)                       UPDATE_NODE_PTR(nodes[k].parent);
				if (nodes[k].d>=0 && nodes[k].first_child) UPDATE_NODE_PTR(nodes[k].first_child);
			}
			for (k=0; k<point_num; k++) if (rev_mapping[k]) UPDATE_NODE_PTR(rev_mapping[k]);
		}

		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		// choose split dimension d0 and number of items 'num' in the first group.

		//d0 = (i->parent) ? ((i->parent->d + 1) % D) : 0;
		//num = num0/2;

		const int FRACTION = 20;
		int num_min = 1;      if (num_min < num0/FRACTION) num_min = num0/FRACTION;
		int num_max = num0-1; if (num_max > (FRACTION-1)*num0/FRACTION) num_max = (FRACTION-1)*num0/FRACTION;
		int num_max_DIM = num0-1;

		if (D>DIM && (!i->parent || !i->parent->parent || !i->parent->parent->parent)) d0 = D-1;
		else
		{
			d0 = -1;
			REAL diff_max = 0;
			for (d=0; d<D; d++)
			{
				ptr = mapping + d*point_num;

				REAL c_min, c_max;
				c_min = coords_array[d][ptr[start +          num_min             ]*skip_array[d]];
				c_max = coords_array[d][ptr[start + ((d<DIM)?num_max:num_max_DIM)]*skip_array[d]];
				REAL diff = c_max - c_min;
				//if (d == DIM) diff *= 2;
				if (d0<0 || diff_max < diff) { d0 = d; diff_max = diff; }
			}
		}

		ptr = mapping + d0*point_num;
		REAL c_min, c_max;
		if (d < DIM)
		{
			c_min = coords_array[d0][ptr[start+num_min]*skip_array[d0]];
			c_max = coords_array[d0][ptr[start+num_max]*skip_array[d0]];
		}
		else
		{
			c_min = coords_array[d0][ptr[start+1]*skip_array[d0]];
			c_max = coords_array[d0][ptr[start+num0-1]*skip_array[d0]]; // for this dimension upper tail is important, so don't discard it!
		}
		REAL split = (c_min + c_max) / 2;
		for (num=num_min; num<num_max; num++)
		{
			if (coords_array[d0][ptr[start+num]*skip_array[d0]] > split) break;
		}

		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////

		ptr = mapping + d0*point_num;
		for (k=start; k<start+num; k++) marking[ptr[k]] = 1;
		for (d=0; d<D; d++)
		{
			if (d == d0) continue;
			ptr = mapping + d*point_num;
			int* bufa = buf;
			int* bufb = buf + num;
			for (k=start; k<start+num0; k++)
			{
				if (marking[ptr[k]]) *bufa ++ = ptr[k];
				else                 *bufb ++ = ptr[k];
			}
			memcpy(ptr+start, buf, num0*sizeof(int));
		}
		ptr = mapping + d0*point_num;
		for (k=start; k<start+num; k++) marking[ptr[k]] = 0;

		i->d = d0;
		PointId p = ptr[start + num - ((num >= num0-num) ? 1 : 0)];
		i->coord = coords_array[d0][p*skip_array[d0]];
		i->first_child = j = &nodes[node_num];
		node_num += 2;
		j->parent = (j+1)->parent = i;
		j->order = start;
		j->d = num;
		(j+1)->order = start+num;
		(j+1)->d = num0-num;
		(j+1)->first_child = first_unprocessed;
		j->first_child = j+1;
		first_unprocessed = j;
	}

	delete [] coords_array;
	delete [] skip_array;
	delete [] mapping;

	////////////////////////////////////////////////////////////////////////////////////
	// set ordering and depth_max
	int depth = 0, depth_max = 0;

	// set ordering
	i = &nodes[0];
	k = 0;
	if (D > DIM)
	{
		// ordering for AddNegativeEdges() - tree preorder
		while ( 1 )
		{
			if (!IS_LEAF(i))
			{
				i = i->first_child;
				depth ++;
				if (depth_max < depth) depth_max = depth;
			}
			else
			{
				while ( 1 )
				{
					i->order = k ++;

					if (!i->parent) break;
					if (i->parent->first_child == i) { i ++; break; }
					i = i->parent;
					depth --;
				}
				if (!i->parent) break;
			}
		}
	}
	else
	{
		// compute tree inorder - useful for nearest neighbor search so that that branches close to the input node are explored first
		while ( 1 )
		{
			if (!IS_LEAF(i))
			{
				i = i->first_child;
				depth ++;
				if (depth_max < depth) depth_max = depth;
			}
			else
			{
				i->order = k ++;
				while ( 1 )
				{
					if (!i->parent) break;
					if (i->parent->first_child == i) { i->parent->order = k ++; i ++; break; }
					i = i->parent;
					depth --;
				}
				if (!i->parent) break;
			}
		}
	}

	traversing_buf = (REAL*) malloc((D + depth_max + 2)*sizeof(REAL));
}

GPMKDTree::~GPMKDTree()
{
	free(nodes);
	free(rev_mapping);
	free(traversing_buf);
}

void GPMKDTree::AddPerfectMatching(PointId* rev_mapping)
{
	Node* i;
	int k;
	PointId p, q = -1;
	i = &nodes[0];
	do
	{
		if (IS_LEAF(i))
		{
			for (k=0; k<-i->d; k++)
			{
				p = i->points[k];
				if (q < 0) q = p;
				else { GPM->AddInitialEdge(rev_mapping[p], rev_mapping[q]); q = -1; }
			}
		}
		else
		{
			i = i->first_child;
			continue;
		}

		while ( i->parent )
		{
			if (i->parent->first_child == i) { i ++; break; }
			i = i->parent;
		}
	} while (i->parent);
}

//////////////////////////////////////////////////////////////////////////////////////////

#define MOVE_DOWN_LEFT(i)\
	{\
		*stack ++ = current_diff[i->d];\
		if (current_diff[i->d] <= 0 && (diff = i->coord - coord0[i->d]) < 0) current_diff[i->d] = diff;\
		i = i->first_child;\
	}

#define MOVE_DOWN_RIGHT(i)\
	{\
		*stack ++ = current_diff[i->d];\
		if (current_diff[i->d] >= 0 && (diff = i->coord - coord0[i->d]) > 0) current_diff[i->d] = diff;\
		i = i->first_child+1;\
	}

#define MOVE_LEFT(i)\
	{\
		int d_prev = i->parent->d;\
		current_diff[d_prev] = stack[-1];\
		if (current_diff[d_prev] <= 0 && (diff = i->parent->coord - coord0[d_prev]) < 0) current_diff[d_prev] = diff;\
		i --;\
	}

#define MOVE_RIGHT(i)\
	{\
		int d_prev = i->parent->d;\
		current_diff[d_prev] = stack[-1];\
		if (current_diff[d_prev] >= 0 && (diff = i->parent->coord - coord0[d_prev]) > 0) current_diff[d_prev] = diff;\
		i ++;\
	}

#define MOVE_UP(i)\
	{\
		i = i->parent;\
		current_diff[i->d] = *(-- stack);\
	}

//////////////////////////////////////////////////////////////////////////////////////////

#define MOVE_DOWN_LEFT_X(i)\
	{\
		*stack ++ = current_diff[i->d];\
		if (i->d < D-1) { if (current_diff[i->d] <= 0 && (diff = i->coord - coord0[i->d]) < 0) current_diff[i->d] = diff; }\
		else current_diff[i->d] = i->coord;\
		i = i->first_child;\
	}

#define MOVE_DOWN_RIGHT_X(i)\
	{\
		*stack ++ = current_diff[i->d];\
		if (i->d < D-1) { if (current_diff[i->d] >= 0 && (diff = i->coord - coord0[i->d]) > 0) current_diff[i->d] = diff; }\
		i = i->first_child+1;\
	}

#define MOVE_LEFT_X(i)\
	{\
		int d_prev = i->parent->d;\
		current_diff[d_prev] = stack[-1];\
		if (d_prev < D-1) { if (current_diff[d_prev] <= 0 && (diff = i->parent->coord - coord0[d_prev]) < 0) current_diff[d_prev] = diff; }\
		else current_diff[d_prev] = i->parent->coord;\
		i --;\
	}

#define MOVE_RIGHT_X(i)\
	{\
		int d_prev = i->parent->d;\
		current_diff[d_prev] = stack[-1];\
		if (d_prev < D-1) { if (current_diff[d_prev] >= 0 && (diff = i->parent->coord - coord0[d_prev]) > 0) current_diff[d_prev] = diff; }\
		i ++;\
	}

#define MOVE_UP_X(i)\
	{\
		i = i->parent;\
		current_diff[i->d] = *(-- stack);\
	}

//////////////////////////////////////////////////////////////////////////////////////////

/*
void GPMKDTree::ComputeKNN(PointId p0, int K, PointId* neighbors_array)
{
	int p;
	REAL* coord0 = GPM->coords + p0*DIM;

	neighbors.Init(K, neighbors_array);

	for (p=0; p<point_num; p++)
	if (p != p0)
	{
		neighbors.Add(p, GPM->Dist(coord0, GPM->coords+p*DIM));
	}
}
*/

void GPMKDTree::ComputeKNN(PointId p0, int K, PointId* neighbors_array)
{
	int neighbor_num = 0;
	Node* i = rev_mapping[p0];
	int k, order0 = i->order;
	REAL diff;
	REAL* coords = GPM->coords;
	REAL* coord0 = GPM->coords + p0*DIM;
	REAL* current_diff = traversing_buf;
	REAL* stack = traversing_buf + D;

	neighbors.Init(K, neighbors_array);

	for (k=0; k<D; k++) current_diff[k] = 0;
	i = &nodes[0];
	do
	{
		if (IS_LEAF(i))
		{
			for (k=0; k<-i->d; k++)
			{
				PointId p = i->points[k];
				if (p == p0) continue;
				double dist2;
				REAL* coord2 = coords+p*DIM;
				GPM_GET_DIST2(dist2, coord0, coord2);
				neighbors.Add(p, dist2);
			}
		}
		else
		{
			if (neighbors.GetNum() < K || GPM->Norm2(current_diff) < neighbors.GetMax())
			{
				if (i->order > order0)
				{
					MOVE_DOWN_LEFT(i);
				}
				else
				{
					MOVE_DOWN_RIGHT(i);
				}
				continue;
			}
		}

		while ( i->parent )
		{
			if (i->parent->order > order0)
			{
				if (i->parent->first_child == i)
				{
					MOVE_RIGHT(i);
					break;
				}
			}
			else
			{
				if (i->parent->first_child != i)
				{
					MOVE_LEFT(i);
					break;
				}
			}
			MOVE_UP(i);
		}
	} while ( i->parent );
}


//////////////////////////////////////////////////////////////////////////////////////////

/*
void GPMKDTree::AddNegativeEdges(PointId p, PerfectMatching* pm)
{
	PointId q;
	for (q=p+1; q<node_num; q++)
	{
		if (GPM->nodes[q].is_marked) continue;
		REAL len = GPM->Dist(p, q);
		if (2*len - GPM->sums[p] - GPM->sums[q] < 0)
		{
			if (pm->AddNewEdge(p, q, len, true)>=0)
			{
				GPM->AddInitialEdge(p, q);
			}
		}
	}
}
*/


void GPMKDTree::AddNegativeEdges(PointId p0, PerfectMatching* pm)
{
	Node* i = rev_mapping[p0];
	int k, order0 = i->order;
	bool check;
	REAL diff;
	REAL* coords = GPM->coords;
	REAL* coord0 = coords + p0*DIM;
	REAL* sums = GPM->sums;
	REAL sum0 = sums[p0];
	REAL* current_diff = traversing_buf;
	REAL* stack = traversing_buf + D;

	for (k=0; i->points[k]!=p0; k++) {}
	for (k++; k<-i->d; k++)
	{
		PointId p = i->points[k];
		
		REAL len = GPM->Dist(coord0, GPM->coords+p*DIM);
		if (2*len - GPM->sums[p] < GPM->sums[p0])
		{
			double start_time = get_time();
			if (pm->AddNewEdge(p0, p, len, true)>=0) GPM->AddInitialEdge(p0, p);
			GPM->graph_update_time += get_time() - start_time;
		}
	}

	for (k=0; k<DIM; k++) current_diff[k] = 0;
	current_diff[k] = sum_max;
	i = &nodes[0];
	do
	{
		if (i->order > order0)
		{
			if (IS_LEAF(i))
			{
				for (k=0; k<-i->d; k++)
				{
					PointId p = i->points[k];
					if (!GPM->nodes[p].is_marked)
					{
						//REAL* coord2 = coords+p*DIM;
						//REAL threshold = sums[p0]+sums[p];
						//GPM_CHECK_DIST(check, coord0, coord2, threshold);
						//if (check)
						{
							REAL len = GPM->Dist(coord0, GPM->coords+p*DIM);
							if (2*len - GPM->sums[p] < GPM->sums[p0])
							{
								double start_time = get_time();
								if (pm->AddNewEdge(p0, p, len, true)>=0) GPM->AddInitialEdge(p0, p);
								GPM->graph_update_time += get_time() - start_time;
							}
						}
					}
				}
			}
			else
			{
				REAL threshold = current_diff[D-1] + sum0;
				GPM_CHECK_NORM(check, current_diff, threshold);
				if (check)
				{
					MOVE_DOWN_LEFT_X(i);
					continue;
				}
			}
		}

		while ( i->parent )
		{
			if (i->parent->first_child == i)
			{
				MOVE_RIGHT_X(i);
				break;
			}
			MOVE_UP_X(i);
		}
	} while (i->parent);
}
