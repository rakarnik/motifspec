//Copyright 1999 President and Fellows of Harvard University
//seqset.h

#ifndef _seqset
#define _seqset
#include "standard.h"

class Seqset{
	int ss_num_seqs;
  vector<vector <int> > ss_seq;
	
 public:
 	Seqset();
	Seqset(const vector<string>& v);
  int num_seqs() const { return ss_num_seqs; }                  // Return number of sequences in this set
  int len_seq(const int i) const { return ss_seq[i].size(); }   // Return length of a specified sequence
  vector<vector<int> > const& seq() const { return ss_seq; }    // Return reference to sequence data
};

#endif
