#include "motifsearch.h"

class MotifSearchExpr : public MotifSearch {
	vector<vector <float> >& expr;                                // Expression data
	int npoints;                                                  // Number of experiments
	vector<float> mean;                                           // Mean expression pattern
	vector<float> expscores;                                      // Expression scores
	vector<struct idscore> expranks;                              // Ranked list of sequences by expression score
	
	void calc_mean();                                             // Calculate the mean expression pattern
	void compute_expr_scores();                  // Calculate expression scores for sequences
	
public:
	MotifSearchExpr(const vector<string>& names, const vector<string>& seqs,
                  const int nc, const int order, const double sim_cut,
									vector<vector <float> >& exprtab, const int npts);
	void reset_search_space();
	void adjust_search_space();
	void set_search_space_cutoff(const int phase);
	int search_for_motif(const int worker, const int iter, const string outfile);
  void print_status(ostream& out, const int i, const int phase);
};