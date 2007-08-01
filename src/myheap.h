//Copyright 1998 President and Fellows of Harvard University
//heap.h

#if !defined(HEAP)
#define HEAP
/*This heap maintains the following data structures:

	kvec[i] holds the values of the keys i=1..num <=max)
	h[i] maintains the keys in a order from h[1] to h[num] such that
		h[i]>h[parent(i)] and h[1] is the lowest
	pos[i] is the inverse of h, that is, pos[h[i]]=i

By only asserting that the heap must maintain this property with each operation,
a sequence of keys associated with values in ascending order may be returned
without doing a full sort on the sequence.  The parent function is controlled
by the base value.  If set to 1, then a full sort is done each time.  To fully
benefit from this algorithm a value of at least two should be used.*/
#include "standard.h"

typedef double keytyp;

class Heap{
	int	max;/* max number of items in heap */
	int	num;/* number of items in heap */
	int	base;			/* base of heap */
	int	*h;			/* {h[1],...,h[n]} is set of items */
	int	*pos;			/* pos[i] gives position of i in h */
	keytyp	*kvec;			/* kvec[i] is key of item i */
	int	minchild(int i);	/* returm smallest child of item */
	void siftup(int i ,int x);/* move item up to restore heap order */
	void siftdown(int i,int x);	/* move item down to restore heap order */
	void heap_error(char *s);
	/* parent of item, leftmost and rightmost children */
	int parent(int x) {return (x+base-2)/base;}
	//guarantees that parent(2)=1
	int left(int x) {return base*(x-1)+2;}//inverse of parent
	int right(int x) {return base*x+1;}
public:	
	Heap(int mx,int bs);
	~Heap();
	void insert(int i,keytyp k);	/* insert item with specified key */
	int	remove(int i);	/* remove item from heap */
	int	del_min();		/* delete and return smallest item */
	void renew_heap();
	keytyp value(int i){return kvec[i];}
};

#endif


