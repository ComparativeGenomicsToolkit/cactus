/*
    GPMkdtree.h - kd-tree data structure for pricing in complete geometric instances

    Copyright 2008 Vladimir Kolmogorov (vnk@adastral.ucl.ac.uk)

    This software can be used for research purposes only. Commercial use is prohibited.
    Public redistribution of the code or its derivatives is prohibited.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#ifndef NJKASKJTASJNBAJSNRBAJS
#define NJKASKJTASJNBAJSNRBAJS

#include "GeomPerfectMatching.h"

struct Neighbors
{
	typedef GeomPerfectMatching::REAL REAL;
	typedef GeomPerfectMatching::PointId PointId;

	Neighbors();
	~Neighbors();

	int GetNum() { return num; }

	void Init(int K, PointId* array);
	void Add(PointId p, double dist);
	double GetMax();

private:
	void Swap(int k1, int k2);
	PointId* array;
	double* dist_array;
	int num, K, K_max;
};



struct GPMKDTree
{
	typedef GeomPerfectMatching::REAL REAL;
	typedef GeomPerfectMatching::PointId PointId;

	GPMKDTree(int D, int point_num, REAL* coords, GeomPerfectMatching* GPM);
	~GPMKDTree();

	// if D == DIM
	void AddPerfectMatching(PointId* rev_mapping);
	void ComputeKNN(PointId p, int K, PointId* neighbors);

	// if D == DIM+1
	void AddNegativeEdges(PointId p, PerfectMatching* pm);

	//////////////////////////////////////////////////////////////////////////

public:
#define CHILDREN_MAX 2
	struct Node
	{
		Node* parent;
		int d; // split dimension. d<0 indicates a leaf.
#define IS_LEAF(i) ((i)->d < 0)
		union
		{
			struct // for non-leaves
			{
				REAL coord;
				Node* first_child; // the second child is first_child+1
			};
			struct // for leaves
			{
				PointId points[CHILDREN_MAX]; // the number of points is -d
			};
		};
		int order;
	}* nodes;


	//////////////////////////////////////////////////////////////////////////
	Node** rev_mapping;

	int D, DIM, point_num, node_num;
	REAL sum_max;
	REAL* traversing_buf;
	GeomPerfectMatching* GPM;

	Neighbors neighbors;
};




#endif
