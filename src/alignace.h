//Copyright 1998 President and Fellows of Harvard University
//alignace.h

#ifndef _alignace
#define _alignace
#include "standard.h"
#include "Random.h"
#include "seqset.h"
#include "sites.h"
#include "archivesites.h"
#include "myheap.h"

struct AlignACEParams{
  int ap_expect;//number of expected sites
  double ap_weight;//fractional weight on priors
  double ap_psfact;//psfact*numsites=npseudo
  double ap_npseudo;//number of pseudo counts
  double ap_gcback;//genome gc content
  double ap_backfreq[6];//array for gc content
  double ap_pseudo[6];
  int ap_maxlen;//maximum length of sites
  int ap_npass;
  int ap_minpass[3];
  double ap_sitecut[3];
  int ap_nruns;
  bool ap_fragment;
  int ap_seed;
  double ap_select;
  int ap_flanking;
  int ap_undersample;
  int ap_oversample;
};

class AlignACE{
 public:
  AlignACEParams ace_params;
  int ace_max_motifs;
  double ace_map_cutoff;
  double ace_sim_cutoff;
  bool ace_verbose;

  Random<int> ace_ran_int;
  Random<double> ace_ran_dbl;

  Seqset ace_seqset;
  Sites ace_sites;
  Sites ace_select_sites;
  ArchiveSites ace_archive;
	
	int* ace_membership;
	int ace_members;

  int *ace_freq_matrix;
  double *ace_score_matrix;
  double **ace_site_bias;

	AlignACE();
  ~AlignACE();
	void init(const vector<string>& v, const int nc=10, const int bf=100, const double map_cut=-20.0, const double sim_cut=0.8);
  void modify_params(int argc, char *argv[]);
  void doit();
  void output(ostream &fout);
  static void print_usage(ostream &fout);
  static void print_version(ostream &fout);
  double get_best_motif(int i=0);
  void full_output(ostream &fout);
  void full_output(char *name);
  void output_params(ostream &fout);
  void set_default_params();//modifiable
  void set_final_params();//derived from others
  void ace_initialize();
	void add_possible(int poss);
	void remove_possible(int poss);
	void clear_possible();
  void seed_random_sites(const int num);
	void seed_random_sites_restricted(const int num);
  void seed_biased_site();
  void calc_matrix();
	string consensus();
	bool consider_site(int gene, float corr, const double minprob = 0.0);
  void single_pass(const double minprob = 0.0);
	void single_pass_restricted(const double minprob = 0.0);
	void single_pass_select(const double minprob = 0.0);
	bool column_sample(const int c, const bool sample);
  bool column_sample(const int c){return column_sample(c,true);}
  bool column_sample(const bool sample) {return column_sample(-1,sample);}
  bool column_sample() {return column_sample(-1,true);}
  double map_score();
	double map_score_restricted();
  void optimize_columns();
  void optimize_sites();
  void orient_motif();
	void debug_check_columns();
};

void cAlignACE(int argc, char *argv[]);


#endif

