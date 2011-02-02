/*
    GeomPerfectMatching.h - computing min cost perfect matching in complete geometric instances (with edge weights equal to Euclidean distances)

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

#ifndef ASKHALKSJBRASMNABFAJSTAS
#define ASKHALKSJBRASMNABFAJSTAS

#include <assert.h>
#include <math.h>
#include "../PerfectMatching.h"

//#define DELAUNAY_TRIANGLE

struct GPMKDTree;

class GeomPerfectMatching
{
public:
	typedef int REAL; // if you change it to double, you should also change PerfectMatching::REAL to double!
	typedef int PointId;

	// pointNum must be a positive even number.
	GeomPerfectMatching(int pointNum, int DIM);
	~GeomPerfectMatching();

	// Must be called exactly point_num times (before anything else can be called).
	// First call returns 0, second returns 1, and so on.
	// coord must be an array of size DIM. (This array is read into internal memory.)
	PointId AddPoint(REAL* coord);


	////////////// SOLVING //////////////

	// solves perfect matching in the complete graph. Inefficient, just for testing.
	REAL SolveComplete(); 


	// options for Solve()
	struct GPMOptions
	{
		GPMOptions() : init_Delaunay(true), init_KNN(0), init_greedy(true), iter_max(0) {}

		// three variables below determine the initial subset of edges. 
		// Delaunay initialization seems to be more robust than K nearest neighbors
		// (but you need to download the "Triangle" package of Shewchuk
		// from http://www.cs.cmu.edu/~quake/triangle.html , extract it to the directory 'triangle'
		// and define DELAUNAY_TRIANGLE above). 
		bool	init_Delaunay;  // add Delaunay triangulation edges.
		int		init_KNN;       // use init_KNN nearest neighbors for each point.
		bool	init_greedy;    // add edges greedily to make sure that a perfect matching exists (see comments before CompleteInitialMatching() in GPMinit.cpp for details)

		int		iter_max;   // If iter_max <= 0 then adds subsets of edges until an optimal solution is found.
				            // Otherwise runs at most iter_max iterations, so the solution may be suboptimal. 
				            // (iter_max=1 runs perfect matching just for the initial subset).
	};
	struct PerfectMatching::Options options;
	struct GPMOptions gpm_options;

	REAL Solve(); 

	// You can also specify the initial subset manually
	void AddInitialEdge(PointId i, PointId j);

	//////////// READING RESULTS //////////
	// Can be called after Solve() or SolveComplete(). Returns pointer to an array 'matching' 
	// of size node_num. (matching[i] is the point corresponding to point i).
	// User must not modify this array.
	PointId GetMatch(PointId p) { return matching[p]; }



	REAL Dist(REAL* coord1, REAL* coord2);
	REAL Dist(PointId p, PointId q);


	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
private:

friend struct GPMKDTree;

	struct Edge
	{
		PointId		head[2];
		Edge*		next[2];
	};
	struct Node
	{
		Edge*		first[2];
		int			is_marked;
	};

	Node*			nodes;
	Block<Edge>*	edges;
	REAL*			coords; // array of size DIM*node_num_max
	REAL*			sums; // array of size node_num_max
	PointId*		matching; // array of size node_num_max

	int DIM;
	int node_num, node_num_max;
	int edge_num;

	double graph_update_time;

// below x denotes REAL, X denotes double, c,c1,c2 are coords (REAL*)
#define GPM_ROUND(X) (REAL)( ( ((REAL)1 / 2) == 0 ) ? ((X)+0.5) : (X) )
#define GPM_GET_NORM2(X, c) { int d; X = 0; for (d=0; d<DIM; d++) X += ((double)(c)[d])*(c)[d]; }
#define GPM_GET_NORM(x, c) { double X; GPM_GET_NORM2(X, c); X = sqrt(X); x = GPM_ROUND(X); }
#define GPM_GET_DIST2(X, c1, c2) { int d; X = 0; for (d=0; d<DIM; d++) X += ((double)((c1)[d]-(c2)[d]))*((c1)[d]-(c2)[d]); }
#define GPM_GET_DIST(x, c1, c2) { double X; GPM_GET_DIST2(X, c1, c2); X = sqrt(X); x = GPM_ROUND(X); }
// if returns false, then 2*||coord||-treshold >= 0
#define GPM_CHECK_NORM(result, c, x)\
	{\
		if (threshold <= 0) result = false;\
		else\
		{\
			double X;\
			GPM_GET_NORM2(X, c);\
			result = (4*X < (double)x*threshold);\
		}\
	}
#define GPM_CHECK_DIST(result, c1, c2, x)\
	{\
		if (threshold <= 0) result = false;\
		else\
		{\
			double X;\
			GPM_GET_DIST2(X, c1, c2);\
			result = (4*X < (double)x*threshold);\
		}\
	}

	double Norm2(REAL* coord)                { double norm2; GPM_GET_NORM2(norm2, coord); return norm2; }
	REAL Norm(REAL* coord)                   { REAL norm;    GPM_GET_NORM (norm,  coord); return norm;  }
	double Dist2(REAL* coord1, REAL* coord2) { double dist2; GPM_GET_DIST2(dist2, coord1, coord2); return dist2; }
	double Dist2(PointId p, PointId q) { return Dist2(coords+DIM*p, coords+DIM*q); }

	void CompleteInitialMatching(); // add edges so that a perfect matching is guaranteed to exist
	void InitKNN(int K);
	void InitDelaunay();

	REAL ComputeCost(PointId* matching);
};

inline GeomPerfectMatching::REAL GeomPerfectMatching::Dist(REAL* coord1, REAL* coord2)
{
	REAL dist;
	GPM_GET_DIST (dist,  coord1, coord2);
	return dist; 
}

inline GeomPerfectMatching::REAL GeomPerfectMatching::Dist(PointId p, PointId q)
{
	return Dist(coords+DIM*p,  coords+DIM*q);
}

#endif
