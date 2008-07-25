//Copyright 1999 President and Fellows of Harvard University
//seqset.cpp

#include "seqset.h"

Seqset::Seqset() {
  map<char,int> code;
  code['n'] = code['N'] = 0;
  code['a'] = code['A'] = 1;
  code['c'] = code['C'] = 2;
  code['g'] = code['G'] = 3;
  code['t'] = code['T'] = 4;
	ss_num_seqs = 0;
  ss_len_seq = new int[ss_num_seqs];
  ss_seq = new char*[ss_num_seqs];
  gc = new float[ss_num_seqs];
}

Seqset::Seqset(const vector<string>& v) {
  map<char,int> code;
  code['n'] = code['N'] = 0;
  code['a'] = code['A'] = 1;
  code['c'] = code['C'] = 2;
  code['g'] = code['G'] = 3;
  code['t'] = code['T'] = 4;
	ss_num_seqs = v.size();
  ss_len_seq = new int[ss_num_seqs];
  ss_seq = new char*[ss_num_seqs];
  gc = new float[ss_num_seqs];
	for(int i = 0; i < ss_num_seqs; i++) {
    ss_len_seq[i] = v[i].length();
    ss_seq[i] = new char[ss_len_seq[i]];
		gc[i] = 0.0;
    for(int j = 0; j < ss_len_seq[i]; j++){
      ss_seq[i][j] = code[v[i][j]];
			if(code[v[i][j]] == 2 || code[v[i][j]] == 3) {
				gc[i]++;
			}
    }
		gc[i] /= ss_len_seq[i];
  }
}

Seqset::~Seqset(){
	for(int i=0;i<ss_num_seqs;i++){
    delete [] ss_seq[i];
  }
  delete [] ss_seq;
  delete [] ss_len_seq;
	delete [] gc;
}

