#ifndef _semodel
#define _semodel

#include "standard.h"
#include "Random.h"
#include "seqset.h"
#include "bgmodel.h"
#include "archivesites.h"

// Search types
#define UNDEFINED 0
#define EXPRESSION 1
#define SUBSET 2
#define SCORE 3

struct SEParams{
  int expect;						           // number of expected sites
  double weight;				           // fractional weight on priors
  double psfact;				           // psfact * numsites=npseudo
  double npseudo;				           // number of pseudo counts
  double backfreq[4];		           // array for gc content
  double pseudo[4];                // pseudocounts for any frequency calculations
  int maxlen;						           // maximum length of sites
  int npass;
  int minpass;
  int nruns;
	float minprob[3];
  bool fragment;
  int seed;
  double select;
  int flanking;
  int undersample;
  int oversample;
	int minsize;
	float mincorr;
};

/* Integrated model */
class SEModel {
	/* Common */
	int search_type;
	vector<string> nameset;
	int ngenes;
	vector<bool> possible;
	
	/* Sequence model */
	SEParams separams;
  bool verbose;
  Random<int> ran_int;
  Random<double> ran_dbl;
  Seqset seqset;
	BGModel bgmodel;
  Motif motif;
	Motif select_sites;
  ArchiveSites archive;
	int members;
  
  /* Expression model */
	vector<vector <float> >& expr;
	int npoints;
	vector<float> mean;
	
	/* Subset */
	const vector<string>& subset;

	/* Score */
	vector<float>& scores;
  vector<float> cumul_scores;
	
	struct idscore {
		int id;
		double score;
	};
	
	struct iscomp {
		bool operator() (struct idscore is1, struct idscore is2) { return (is1.score > is2.score); }
	} isc;
	
	vector<double> seqscores;
	vector<int> bestpos;
	vector<struct idscore> seqranks;
	vector<double> expscores;
	vector<struct idscore> expranks;
	vector<struct idscore> scranks;

	double score_site(double* score_matrix, const int c, const int p, const bool s);
	void set_cutoffs();
	void set_seq_cutoff(const int phase);
	void set_seq_cutoff_expr(const int phase);
	void set_seq_cutoff_subset(const int phase);
	void set_seq_cutoff_score(const int phase);
	void set_expr_cutoff();
	void set_score_cutoff();
	void print_possible(ostream& out);
	
 public:
	/* Return codes for search */
	static const int BAD_SEARCH_SPACE = 1;
	static const int TOO_SIMILAR = 2;
	static const int BAD_SEED = 3;
	static const int TOO_FEW_SITES = 4;
	static const int TOO_MANY_SITES = 5;

	/* General */
	SEModel(const int s_type, const vector<string>& names, const vector<string>& seqs,
					vector<vector <float> >& exprtab, const vector<string>& sub, vector<float>& sctab,
					const int npts, const int nc, const int order, const double sim_cut);
	void modify_params(int argc, char *argv[]);
  double get_best_motif(int i=0);
  void output_params(ostream &fout);
  void set_default_params();//modifiable
  void set_final_params();//derived from others
	SEParams get_params() const { return separams; };
  void ace_initialize();
	vector<string> names() const { return nameset; };
	ArchiveSites& get_archive() { return archive; };
	
	/* Manage membership */
	void add_possible(const int gene);                            // Add gene as a potential member
	void remove_possible(const int gene);                         // Remove gene as a potential member
	void clear_all_possible();                                    // Remove all genes as potential members
	void reset_possible();                                        // Adjust potential membership back to initial conditions
	bool is_possible(const int gene) const;									      // Return whether this gene is a potential member
	int possible_size() const;															      // Return number of potential members
	int total_positions() const;														      // Return total number of possible positions
	int possible_positions() const;													      // Return number of potential positions
	void clear_sites();																			      // Remove all sites currently in model
	bool is_member(const int gene) const;										      // Return whether the gene is a member
	int size() const;                                             // Return number of genes
	int motif_size() const;																	      // Return number of sites
	void genes(int* genes) const;                                 // Return the genes that are assigned to this model
	
	/* Sequence model*/
	Seqset& get_seqset() { return seqset; }                       // Return the set of sequences
	void calc_matrix(double* score_matrix);                       // Calculate the PWM for the current set of sites
	double score();                                               // Calculate the score of the current model
	double matrix_score();                                        // Calculate the entropy score for the current matrix
	double over_score();                                          // Calculate the overrepresentation score
	double spec_score();                                          // Calculate the specificity score
	
	/* Expression model */
	void calc_mean();                                             // Calculate the mean for this model
	vector<float>& get_mean() { return mean; };                   // Get the mean for this model
	float get_corr_with_mean(const vector<float>& pattern) const; // Calculate the correlation between the model mean and 'pattern
	
	/* Algorithm steps */
	void seed_random_site();
  void seed_high_scoring_site();
  void single_pass(const double seqcut = 0.0, bool greedy = false);
	void single_pass_select(const double seqcut = 0.0, bool greedy = false);
	void compute_seq_scores();
	void compute_seq_scores_minimal();
	void compute_expr_scores();
	bool column_sample(const bool add, const bool remove);
	void adjust_search_space();
	int search_for_motif(const int worker, const int iter, const string outfile);
  int search_for_motif_expr(const int worker, const int iter, const string outfile);
  int search_for_motif_subset(const int worker, const int iter, const string outfile);
  int search_for_motif_score(const int worker, const int iter, const string outfile);
	bool consider_motif(const char* filename);
	
	/* Output */
	void full_output(ostream &fout);                        // Output all motifs stored in archive_sites
  void full_output(char *name);
	void print_status(ostream& out, const int i, const int phase);
};

#endif
