//Copyright 1998 President and Fellows of Harvard University
//genome.h

#ifndef _Random
#define _Random
#include "standard.h"

template<class C> class Random{
  int r_seed;
  C r_min;
  C r_max;
  double r_range, r_denom;
  void set_type(){
    r_denom=RAND_MAX;
    if((C)(RAND_MAX/(RAND_MAX+1.0))==0){
      //integer math, increase range (in set_range), denom by one
      r_denom+=1.0;
    }
  }

 public:

  Random(int i=-1){
    set_type();
    set_seed(i);
  }

  Random(C min, C max, int i=-1){
    set_type();
    set_seed(i);
    set_range(min, max);
    cerr<<r_min<<'\t'<<r_max<<'\t'<<r_range<<'\t'<<r_seed<<'\n';
  }

  void set_range(C min, C max){
    r_min=min;
    r_max=max;
    r_range=r_max-r_min;
    if(r_range>RAND_MAX){
      cerr<<"Range is greater than RAND_MAX.  A specialized implementation is needed.\n";
      cerr<<" Continuing anyway\n\n";
    }
    if(r_denom>RAND_MAX){
      r_range+=1.0;//integer math
    }
  }

  void set_seed(int i=-1){
    r_seed=i;
    if(r_seed==-1){r_seed=(unsigned)time(NULL);}
    while(r_seed<0||r_seed>RAND_MAX){
      cerr<<"Improper seed to random number generator "<<r_seed<<'\n';
      r_seed=(unsigned)time(NULL);
      cerr<<" Continuing with seed "<<r_seed<<"\n\n";
    }
    srand(r_seed);
    rand();rand();rand();rand();rand();rand();rand();
    //some implementations are not very random for the first few calls
    //so the first seven are thrown out.
    //for floats, 0.0 is possible, 1.0 is not.  This is easier, and 
    //probably safer wrt division errors leading to >1.0 results.
  }

  int seed() const {return r_seed;}

  C rnum(){
    return (C)( r_range*rand()/r_denom )+ r_min ;
  }
  //integer translation (1-5): (int)(5.0*rand()/(RAND_MAX+1)) + 1
  //double translation (1.0-5.0): (double)(4.0*rand()/RAND_MAX) + 1.0

};

#endif
