//Copyright 1999 President and Fellows of Harvard University
//seqset.h

#ifndef _seqset
#define _seqset
#include "standard.h"

class Seqset{
  int ss_num_seqs;
  int* ss_len_seq;
  char** ss_seq;//[ace_num_seqs][ace_len_seq[i]] ~ NACGTN=012345
	float gc_genome;
	float* gc;

 public:
 	Seqset();
	Seqset(const vector<string>& v);
  ~Seqset();
  int num_seqs() const { return ss_num_seqs; }
  int len_seq(const int i) const { return ss_len_seq[i]; }
  char** seq_ptr() const { return ss_seq; }
	float gcgenome() const { return gc_genome; }
	float gccontent(const int i) const { return gc[i]; }
};

#endif
