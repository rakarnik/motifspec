//Copyright 1998 President and Fellows of Harvard University
//heap.cpp

#include "myheap.h"

Heap::Heap(int N,int D){
	/* Initialize a heap to store items in {1,...,N}.*/
	int i;
	
	max = N;
	base = D;
	num = 0;
	h=new int[N+1];
	pos=new int[N+1];
	kvec=new keytyp[N+1];
	for (i=1; i<= N; i++) pos[i] = -1;
}

Heap::~Heap(){
	delete[] h;
	delete[] pos;
	delete[] kvec;
}

void Heap::renew_heap(){
	int i;

	num=0;
	for(i=1;i<=max;i++) pos[i]=-1;
}

void Heap::insert(int i,keytyp k)
/* insert item i with specified key */
{
	if(i<1) heap_error("Heap::insert: fatal! item to insert is < 1");
	if(pos[i]!=-1) remove(i);
	kvec[i] = k; num++; siftup(i,num);
}

int Heap::remove(int i)
/* Remove item i from heap. */
{
	int j;
	if(pos[i]==-1) return -1;
	j = h[num--];
	if (i != j && kvec[j] <= kvec[i]) siftup(j,pos[i]);
	else if (i != j && kvec[j]>kvec[i]) siftdown(j,pos[i]);
	pos[i] = -1;
	return i;
}

int	Heap::del_min()
/* delete and return item with smallest key */
{
	int i;
	if (num == 0) return -1;
	i = h[1];
	remove(h[1]);
	return i;
}


void Heap::siftup(int i ,int x)
/* Shift i up from position x to restore heap order.*/
{
	int px = parent(x);
	while (x > 1 && kvec[h[px]] > kvec[i]) {
		h[x] = h[px]; pos[h[x]] = x;
		x = px; px = parent(x);
	}
	h[x] = i; pos[i] = x;
}

void Heap::siftdown(int i,int x)
/* Shift i down from position x to restore heap order.*/
{
	int cx = minchild(x);
	while (cx != -1 && kvec[h[cx]] < kvec[i]) {
		h[x] = h[cx]; pos[h[x]] = x;
		x = cx; cx = minchild(x);
	}
	h[x] = i; pos[i] = x;
}

int Heap::minchild(int x)
/* Return the position of the child of the item at position x
   having minimum key. */
//min(left...right)
{
	int y, minc;
	if ((minc = left(x)) > num) return -1;
	for (y = minc + 1; y <= right(x) && y <= num; y++) {
		if (kvec[h[y]] < kvec[h[minc]]) minc = y;
	}
	return minc;
}

void Heap::heap_error(char *s){cerr<<"dheap: "<<s<<'\n';exit(3);}

