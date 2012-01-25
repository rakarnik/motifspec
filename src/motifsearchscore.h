#ifndef _motifsearchscore
#define _motifsearchscore

#include "motifsearch.h"

class MotifSearchScore : public MotifSearch {
	vector<float>& scores;
	float mean_sc;
	float stdev_sc;
	vector<float> cumul_scores;
	vector<struct idscore> scranks;

public:
	MotifSearchScore(const vector<string>& names, const vector<string>& seqs,
			const int nc, const int order, const double sim_cut, const int maxm,
			vector<float>& sctab);
	void reset_search_space();
	void adjust_search_space();
	void set_search_space_cutoff(const int phase);
	void calc_matrix(double* score_matrix);
	int search_for_motif(const int worker, const int iter, const string outfile);
	void print_status(ostream& out, const int i, const int phase);
};

#endif

