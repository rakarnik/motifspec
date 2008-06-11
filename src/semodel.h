#ifndef _semodel
#define _semodel

#include "standard.h"
#include "Random.h"
#include "seqset.h"
#include "sites.h"
#include "archivesites.h"
#include "myheap.h"

struct SEParams{
  int expect;						//number of expected sites
  double weight;				//fractional weight on priors
  double psfact;				//psfact*numsites=npseudo
  double npseudo;				//number of pseudo counts
  double gcback;				//genome gc content
  double backfreq[6];		//array for gc content
  double pseudo[6];
  int maxlen;						//maximum length of sites
  int npass;
  int minpass[3];
  int nruns;
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
	vector<string> nameset;
	int ngenes;
	bool* possible;
	int npossible;
	
	/* Sequence model */
	SEParams separams;
  int max_motifs;
  double map_cutoff;
  double sim_cutoff;
  bool verbose;
  Random<int> ran_int;
  Random<double> ran_dbl;
  Seqset seqset;
  Sites sites;
  Sites select_sites;
	Sites print_sites;
  ArchiveSites archive;
	int members;
  int* freq_matrix;
  double* score_matrix;
  
  /* Expression model */
	float** expr;
	int npoints;
	float* mean;
	float* stdev;
	float** pcorr;
	
	struct hitscore {
		int seq;
		double score;
	};
	
	struct hscomp {
		bool operator() (struct hitscore hp1, struct hitscore hp2) { return (hp1.score > hp2.score); }
	} hsc;
	
	double *seqscores;
	struct hitscore *seqranks;
	double *expscores;
	struct hitscore *expranks;
	void set_seq_cutoffs(const int i);
	void set_expr_cutoffs(const int i);
	void print_possible(ostream& out);
	
 public:
	/* Return codes for search */
	static const int BAD_SEARCH_SPACE = 1;
	static const int TOO_SIMILAR = 2;
	static const int BAD_SEED = 3;
	static const int TOO_FEW_SITES = 4;
		
	/* General */
	SEModel(const vector<string>& seqs, float** exprtab, const int numexpr, const vector<string>& names, const int nc = 10, const int bf = 1000, const double map_cut = -20.0, const double sim_cut = 0.8);
  ~SEModel();
	void modify_params(int argc, char *argv[]);
  double get_best_motif(int i=0);
  void output_params(ostream &fout);
  void set_default_params();//modifiable
  void set_final_params();//derived from others
	SEParams get_params() const { return separams; };
  void ace_initialize();
	vector<string> names() const { return nameset; };
	ArchiveSites* get_archive() { return &archive; };
	
	/* Manage membership */
	void add_possible(const int gene);                      // Add gene as a potential member
	void remove_possible(const int gene);                   // Remove gene as a potential member
	void clear_all_possible();                              // Remove all genes as potential members
	bool is_possible(const int gene) const;									// Return whether this gene is a potential member
	int possible_size() const;															// Return number of potential members
	int possible_positions() const;													// Return number of potential positions
	void clear_sites();																			// Remove all sites currently in model
	bool is_member(const int gene) const;										// Return whether the gene is a member
	int size() const;                                       // Return number of genes
	int sites_size() const;																	// Return number of sites
	void genes(int* genes) const;                           // Return the genes that are assigned to this model
	
	/* Sequence model*/
	Seqset get_seqset() { return seqset; }                  // Return the set of sequences
	void calc_matrix();                                     // Calculate the PWM for the current set of sites
	string consensus() const;                               // Return the consensus sequence for the current set of sites
	double map_score();                                     // Calculate the MAP score
	double spec_score();                                    // Calculate the specificity score
	void orient_motif();                                    // Orient sites according to (A = 1.0, C = -1.5, G = 1.5, T = 1.0) scoring
	void orient_print_motif();                              // Orient print sites according to (A = 1.0, C = -1.5, G = 1.5, T = 1.0) scoring
	
	/* Expression model */
	void calc_mean();                                       // Calculate the mean for this model
	float* get_mean() { return mean; };                     // Get the mean for this model
	void calc_stdev();																			// Calculate the standard deviation for this model
	float get_avg_pcorr();																	// Calculate the average pairwise correlaion of current model
	float get_pcorr(const int g1, const int g2);            // Calculate the pairwise correlation for this pair of genes
	float get_corr_with_mean(const float* pattern) const;   // Calculate the correlation between the model mean and 'pattern
	
	/* Algorithm steps */
	void seed_random_site();
  void seed_biased_site();
  void single_pass(const double minprob = 0.0, bool greedy = false);
	void single_pass_select(const double minprob = 0.0);
	void compute_seq_scores(const bool sample);
	void compute_seq_scores() { compute_seq_scores(true); };
	void compute_expr_scores();
	void compute_scores() { compute_seq_scores(); compute_expr_scores(); };
	bool column_sample(const int c, const bool sample);
  bool column_sample(const int c) { return column_sample(c, true); }
  bool column_sample(const bool sample) { return column_sample(1000, sample); }
  bool column_sample() { return column_sample(1000,true); }
  void optimize_columns();
  void optimize_sites();
	void expand_search_around_mean(const double corr_cutoff);
	void expand_search_min_pcorr(const double corr_cutoff);
	void expand_search_avg_pcorr();
	void search_for_motif(const int worker, const int iter);
	bool consider_motif(const char* filename);
	
	/* Output */
	void full_output(ostream &fout);                        // Output all motifs stored in archive_sites
  void full_output(char *name);
	void print_status(ostream& out, const int i, const int phase);
};

#endif

