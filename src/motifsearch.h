#ifndef _motifsearch
#define _motifsearch

#include "standard.h"
#include "Random.h"
#include "fastmath.h"
#include "seqset.h"
#include "bgmodel.h"
#include "archivesites.h"
#include "searchparams.h"

class MotifSearch {
protected:
	/* Common */
	int search_type;
	vector<string> nameset;
	int ngenes;
	Random<int> ran_int;
	Random<double> ran_dbl;

	struct idscore {
		int id;
		double score;
	};
	
	struct iscomp {
		bool operator() (struct idscore is1, struct idscore is2) { return (is1.score > is2.score); }
	} isc;
	
	/* Sequence model */
	SearchParams params;
	Seqset seqset;
	BGModel bgmodel;
	Motif motif;
	Motif select_sites;
	ArchiveSites archive;
	vector<double> seqscores;
	vector<struct idscore> seqranks;
	vector<int> bestpos;
	vector<bool> beststrand;
	int members;
	
	double score_site(double* score_matrix, const int c, const int p, const bool s);
	void set_cutoffs();
	void set_seq_cutoff(const int phase);
	virtual void set_search_space_cutoff(const int phase) = 0;
	
public:
	/* Return codes for search */
	static const int BAD_SEARCH_SPACE = 1;
	static const int TOO_SIMILAR = 2;
	static const int BAD_SEED = 3;
	static const int TOO_FEW_SITES = 4;
	static const int TOO_MANY_SITES = 5;
	
	/* General */
	MotifSearch(const vector<string>& names, const vector<string>& seqs,
			const int nc, const int order, const double sim_cut, const int maxm);
	virtual ~MotifSearch() {};
	void modify_params(int argc, char *argv[]);
	double get_best_motif(int i=0);
	void output_params(ostream &fout);
	void set_default_params();
	virtual void set_final_params();
	SearchParams get_params() const { return params; };
	void ace_initialize();
	const vector<string>& names() const { return nameset; };
	ArchiveSites& get_archive() { return archive; };
	
	/* Manage search space */
	virtual void reset_search_space() = 0;                        // Set search space back to initial conditions
	virtual void adjust_search_space() = 0;                       // Set search space according to current cutoffs
	int total_positions() const;														      // Return total number of possible positions
	int positions_in_search_space() const;											  // Return number of potential positions
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
	
	/* Algorithm steps */
	void update_seq_count();
	void seed_random_site();
	void single_pass(bool greedy = false);
	void single_pass_select(bool greedy = false);
	void compute_seq_scores();
	void compute_seq_scores_minimal();
	virtual int search_for_motif(const int worker, const int iter, const string outfile) = 0;
	bool consider_motif(const char* filename);
	
	/* Output */
	void full_output(ostream &fout);                        // Output all motifs stored in archive_sites
	void full_output(char *name);
	virtual void print_status(ostream& out, const int i, const int phase);
};

#endif
