#ifndef _motifsearchsubset
#define _motifsearchsubset

#include "motifsearch.h"

class MotifSearchSubset : public MotifSearch {
	const vector<string>& subset;
	
public:
	MotifSearchSubset(const vector<string>& names, const vector<string>& seqs,
									const int nc, const int order, const double sim_cut, const int maxm,
									const vector<string>& sub);
	void set_final_params();
	void reset_search_space();
	void adjust_search_space();
	void set_search_space_cutoff(const int phase);
	double score();
	int search_for_motif(const int worker, const int iter, const string outfile);
};

#endif

