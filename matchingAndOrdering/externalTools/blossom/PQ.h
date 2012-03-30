/*
    PQ.h - implements pairing heaps priority queue (with multipass 'delete min')

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

#ifndef HFKSJHFKJHARBABDAKFAF
#define HFKSJHFKJHARBABDAKFAF

// exactly one flag must be defined
//#define PQ_MULTIPASS
#define PQ_INTERLEAVED_MULTIPASS

#include <string.h>

template <typename REAL> class PriorityQueue
{
public:
	struct Item
	{
		REAL	slack;

		Item*	parentPQ;
		union
		{
			struct
			{
				Item*	leftPQ;
				Item*	rightPQ;
			};
			REAL	y_saved; // used in repairs
		};
	};
	static void* AllocateBuf();
	static void DeallocateBuf(void* buf);

	static void ResetItem(Item* i);
	static bool isReset(Item* i);

	//////////////////////////////////////////////////////////

	void Reset();
	void Add(Item* i);
#define Remove(i, buf) _Remove(i)
	void _Remove(Item* i);
	void Decrease(Item* i_old, Item* i_new, void* buf);
	Item* GetMin();

	//////////////////////////////////////////////////////////

	void Update(REAL delta);
	void Merge(PriorityQueue<REAL>& dest);

	// traversing items in the order they are stored (irrespective of slack).
	// The caller must go through all items, no other member functions can be called during the scan.
	Item* GetAndResetFirst();
	Item* GetAndResetNext();

	Item* GetFirst();
	Item* GetNext(Item* i);

	//////////////////////////////////////////////////////////

private:
	struct Buf
	{
	};
	Item*	rootPQ;
	void RemoveRoot();
};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

template <typename REAL> inline void* PriorityQueue<REAL>::AllocateBuf()
{
	return NULL;
}

template <typename REAL> inline void PriorityQueue<REAL>::DeallocateBuf(void* _buf)
{
}

template <typename REAL> inline void PriorityQueue<REAL>::ResetItem(Item* i) 
{ 
	i->parentPQ = NULL;
}

template <typename REAL> inline bool PriorityQueue<REAL>::isReset(Item* i) 
{ 
	return (i->parentPQ == NULL);
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

template <typename REAL> inline void PriorityQueue<REAL>::Reset() 
{ 
	rootPQ = NULL; 
}

/*
template <typename REAL> inline void PriorityQueue<REAL>::RemoveRoot()
{
	Item* r = rootPQ;
	PriorityQueue<REAL> pq;
	pq.rootPQ = rootPQ;
	rootPQ = NULL;
	Item* i;
	for (i=pq.GetAndResetFirst(); i; i=pq.GetAndResetNext())
	{
		if (i != r) Add(i);
	}
	r->parentPQ = NULL;
}
*/

// sets i = merge(i, j). Ignores parentPQ and rightPQ for i and j.
#define MERGE_PQ(i, j)\
	{\
		if (i->slack <= j->slack)\
		{\
			j->rightPQ = i->leftPQ;\
			if (j->rightPQ) j->rightPQ->parentPQ = j;\
			j->parentPQ = i;\
			i->leftPQ = j;\
		}\
		else\
		{\
			i->rightPQ = j->leftPQ;\
			if (i->rightPQ) i->rightPQ->parentPQ = i;\
			i->parentPQ = j;\
			j->leftPQ = i;\
			i = j;\
		}\
	}

template <typename REAL> inline void PriorityQueue<REAL>::RemoveRoot()
{
	Item* i = rootPQ->leftPQ;
	rootPQ->parentPQ = NULL;
	if (i)
	{
#ifdef PQ_MULTIPASS
		while ( i->rightPQ )
		{
			Item** prev_ptr = &rootPQ;
			while ( 1 )
			{
				if (i->rightPQ)
				{
					Item* j = i->rightPQ;
					Item* next = j->rightPQ;
					MERGE_PQ(i, j);
					*prev_ptr = i;
					if (!next) { i->rightPQ = NULL; break; }
					prev_ptr = &i->rightPQ;
					i = next;
				}
				else
				{
					*prev_ptr = i;
					i->rightPQ = NULL;
					break;
				}
			}
			i = rootPQ;
		}
#endif

#ifdef PQ_INTERLEAVED_MULTIPASS
		while ( i->rightPQ )
		{
			Item* prev = NULL;
			while ( i )
			{
				Item* next;
				if (i->rightPQ)
				{
					Item* j = i->rightPQ;
					next = j->rightPQ;
					MERGE_PQ(i, j);
				}
				else next = NULL;
				i->rightPQ = prev;
				prev = i;
				i = next;
			}
			i = prev;
		}
#endif
		i->parentPQ = i;
	}
	rootPQ = i;
}

template <typename REAL> inline void PriorityQueue<REAL>::Add(Item* i)
{
	if (!rootPQ)
	{
		rootPQ = i;
		i->parentPQ = i;
		i->leftPQ = i->rightPQ = NULL;
	}
	else if (i->slack <= rootPQ->slack)
	{
		rootPQ->parentPQ = i;
		i->leftPQ = rootPQ;
		i->rightPQ = NULL;
		rootPQ = i;
		i->parentPQ = i;
	}
	else
	{
		i->leftPQ = NULL;
		i->rightPQ = rootPQ->leftPQ;
		if (i->rightPQ) i->rightPQ->parentPQ = i;
		rootPQ->leftPQ = i;
		i->parentPQ = rootPQ;
	}
}


template <typename REAL> inline void PriorityQueue<REAL>::_Remove(Item* i)
{
	Item* p = i->parentPQ;
	if (p == i) RemoveRoot();
	else
	{
		if (i->rightPQ) i->rightPQ->parentPQ = p;
		if (p->leftPQ == i) p->leftPQ  = i->rightPQ;
		else                p->rightPQ = i->rightPQ;
		if (i->leftPQ)
		{
			i->parentPQ = i;
			i->rightPQ = NULL;
			PriorityQueue<REAL> pq;
			pq.rootPQ = i;
			pq.RemoveRoot();
			pq.Merge(*this);
		}
		else i->parentPQ = NULL;
	}
}

template <typename REAL> inline void PriorityQueue<REAL>::Decrease(Item* i_old, Item* i_new, void* _buf)
{
	if (i_old->parentPQ == i_old)
	{
		if (i_old != i_new)
		{
			rootPQ = i_new;
			i_new->parentPQ = i_new;
			i_new->leftPQ = i_old->leftPQ;
			i_new->rightPQ = NULL;
			if (i_new->leftPQ) i_new->leftPQ->parentPQ = i_new;
			i_old->parentPQ = NULL;
		}
	}
	else
	{
		Remove(i_old, _buf);
		Add(i_new);
	}
}

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetMin()
{
	return rootPQ;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



template <typename REAL> inline void PriorityQueue<REAL>::Merge(PriorityQueue<REAL>& dest)
{
	if (!rootPQ) return;
	if (!dest.rootPQ) dest.rootPQ = rootPQ;
	else
	{
		if (rootPQ->slack < dest.rootPQ->slack)
		{
			Item* j = rootPQ; rootPQ = dest.rootPQ; dest.rootPQ = j;
		}
		rootPQ->rightPQ = dest.rootPQ->leftPQ;
		if (rootPQ->rightPQ) rootPQ->rightPQ->parentPQ = rootPQ;
		rootPQ->parentPQ = dest.rootPQ;
		dest.rootPQ->leftPQ = rootPQ;
	}
	rootPQ = NULL;
}



template <typename REAL> inline void PriorityQueue<REAL>::Update(REAL delta)
{
	if (!rootPQ) return;

	Item* i = rootPQ;
	while (i->leftPQ) i = i->leftPQ;

	while ( 1 )
	{
		i->slack += delta;

		if (i->rightPQ)
		{
			i = i->rightPQ;
			while (i->leftPQ) i = i->leftPQ;
		}
		else
		{
			while ( 1 )
			{
				Item* j = i;
				i = i->parentPQ;
				if (i == j) return;
				if (i->leftPQ == j) break;
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetAndResetFirst()
{
	if (!rootPQ) return NULL;
	return GetAndResetNext();
}

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetAndResetNext()
{
	if (!rootPQ) return NULL;
	Item* result = rootPQ;
	result->parentPQ = NULL;
	Item* i = rootPQ->leftPQ;
	if (!i) rootPQ = result->rightPQ;
	else
	{
		rootPQ = i;
		while (i->rightPQ) i = i->rightPQ;
		i->rightPQ = result->rightPQ;
	}
	return result;
}

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetFirst()
{
	if (!rootPQ) return NULL;
	Item* i = rootPQ;
	while (i->leftPQ) i = i->leftPQ;
	return i;
}

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetNext(Item* i)
{
	if (i->rightPQ)
	{
		i = i->rightPQ;
		while (i->leftPQ) i = i->leftPQ;
		return i;
	}
	while ( 1 )
	{
		Item* j = i;
		i = i->parentPQ;
		if (i == j) return NULL;
		if (i->leftPQ == j) return i;
	}
}

#endif
