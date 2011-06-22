//Copyright 1999 President and Fellows of Harvard University
//seqset.cpp

#include "seqset.h"

Seqset::Seqset() : 
nseqs(0),
seqs(nseqs),
seq_lens(nseqs) {
}

Seqset::Seqset(const vector<string>& v) :
nseqs(v.size()),
seqs(nseqs),
seq_lens(nseqs) {
	map<char, int> nt;
	nt['A'] = nt['a'] = 0;
	nt['C'] = nt['c'] = 1;
	nt['G'] = nt['g'] = 2;
	nt['T'] = nt['t'] = 3;
	for(int i = 0; i < nseqs; i++) {
		seq_lens[i] = v[i].length();
		seqs[i].reserve(seq_lens[i]);
		for(int j = 0; j < seq_lens[i]; j++) {
			seqs[i].push_back(nt[v[i][j]]);
		}
	}
}
