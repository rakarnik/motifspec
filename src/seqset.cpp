//Copyright 1999 President and Fellows of Harvard University
//seqset.cpp

#include "seqset.h"

Seqset::Seqset() : 
ss_num_seqs(0),
ss_seq(ss_num_seqs) {
}

Seqset::Seqset(const vector<string>& v) :
ss_num_seqs(v.size()),
ss_seq(ss_num_seqs) {
	map<char, int> nt;
	nt['A'] = nt['a'] = 0;
	nt['C'] = nt['c'] = 1;
	nt['G'] = nt['g'] = 2;
	nt['T'] = nt['t'] = 3;
  int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = v[i].length();
		ss_seq[i].reserve(len);
    for(int j = 0; j < len; j++) {
			ss_seq[i].push_back(nt[v[i][j]]);
    }
	}
}
