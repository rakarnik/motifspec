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
  double sitecut[4];
  int nruns;
  bool fragment;
  int seed;
  double select;
  int flanking;
  int undersample;
  int oversample;
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
  int *freq_matrix;
  double *score_matrix;
  double **site_bias;

  /* Expression model */
	float** expr;
	int npoints;
	float* mean;
	float* stdev;
	float** pcorr;
	double corr_cutoff;
	
	void set_cutoffs();
	void print_possible(ostream& out);

 public:
	/* Return codes for search */
	static const int BAD_SEARCH_SPACE = 1;
	static const int TOO_SIMILAR = 2;
	static const int BAD_SEED = 3;
	static const int TOO_FEW_SITES = 4;
		
	/* General */
	SEModel() {};
  ~SEModel();
	void init(const vector<string>& seqs, float** exprtab, const int numexpr, const vector<string>& names, const int nc = 10, const int bf = 1000, const double map_cut = -20.0, const double sim_cut = 0.8);
	void modify_params(int argc, char *argv[]);
  double get_best_motif(int i=0);
  void output_params(ostream &fout);
  void set_default_params();//modifiable
  void set_final_params();//derived from others
	SEParams get_params() const { return separams; };
  void ace_initialize();
	vector<string> names() const { return nameset; };
	
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
	void calc_matrix();
	string consensus();
	double map_score();
	void orient_motif();
	void orient_print_motif();
	
	/* Expression model */
	void calc_mean();                                       // Calculate the mean for this model
	float* get_mean();                                      // Get the mean for this model
	void calc_stdev();																			// Calculate the standard deviation for this model
	float get_pcorr(const int g1, const int g2);            // Calculate the pairwise correlation for this pair of genes
	float get_corr_with_mean(const float* pattern) const;   // Calculate the correlation between the model mean and 'pattern
	float prob_gene_given_model(int g) const;               // Calculate probability of gene belonging to this model
	
	/* Algorithm steps */
	void seed_random_site();
  void seed_biased_site();
  void single_pass(const double minprob = 0.0);
	void single_pass_select(const double minprob = 0.0);
	bool column_sample(const int c, const bool sample);
  bool column_sample(const int c){return column_sample(c,true);}
  bool column_sample(const bool sample) {return column_sample(1000,sample);}
  bool column_sample() {return column_sample(1000,true);}
  void optimize_columns();
  void optimize_sites();
	void expand_search_around_mean(const double corr_cutoff);
	void expand_search_min_pcorr(const double corr_cutoff);
	void expand_search_avg_pcorr(const double corr_cutoff);
	void search_for_motif(const double minsize, const double mincorr);
	int search_for_motif_near_seed(const double minsize, const double mincorr);
	
	/* Output */
	void output(ostream &fout);
	void full_output(ostream &fout);
  void full_output(char *name);
	void print_status(ostream& out, const int i, const int phase, const double cutoff, const double sc);
};

#endif

