//Copyright 1999 President and Fellows of Harvard University
//seqset.h

#ifndef _seqset
#define _seqset
#include "standard.h"

class Seqset{
	int ss_num_seqs;
  int total_seq_len;
  vector<vector <int> > ss_seq;
	float gc_genome;
	vector<float> gc;
	vector<float> bgmodel0;
	vector<float> bgmodel1;
	vector<float> bgmodel2;
	vector<float> bgmodel3;
	vector<vector <float> > wbgscores;
	vector<vector <float> > cbgscores;
	
	void train_background3();                                     // Train third order background model
	void train_background2();                                     // Train second order background model
	void train_background1();                                     // Train first order background model
	void train_background0();                                     // Train zeroeth order background model
	void calc_bg_scores();                                        // Calculate scores according to background model
	                                                              // for all sequences
	
 public:
 	Seqset();
	Seqset(const vector<string>& v);
  ~Seqset();
	int num_seqs() const { return ss_num_seqs; }                  // Return number of sequences in this set
  int len_seq(const int i) const { return ss_seq[i].size(); }   // Return length of a specified sequence
  vector<vector<int> > const& seq() const { return ss_seq; }    // Return reference to sequence data
	float gcgenome() const { return gc_genome; }                  // Return GC content of the set
	float gccontent(const int i) const { return gc[i]; }          // Return GC content of a specified sequence
	vector<vector <float> > const& get_wbgscores() const { return wbgscores; }
	vector<vector <float> > const& get_cbgscores() const { return cbgscores; }
};

#endif
