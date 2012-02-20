//Copyright 1999 President and Fellows of Harvard University
//seqset.h

#ifndef _seqset
#define _seqset
#include "standard.h"

class Seqset{
	int nseqs;
	vector<vector <char> > seqs;
	vector<int> seq_lens;
	
public:
	Seqset();
	Seqset(const vector<string>& v);
	int num_seqs() const { return nseqs; }                     // Return number of sequences in this set
	int len_seq(const int i) const { return seq_lens[i]; }     // Return length of a specified sequence
	vector<vector<char> > const& seq() const { return seqs; }   // Return reference to sequence data
};

#endif
