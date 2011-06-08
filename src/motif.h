#ifndef _motif
#define _motif
#include "standard.h"
#include "seqset.h"
#include "site.h"

class Motif {
	const Seqset& seqset;                    // set of sequences that this motif refers to
	int init_nc;                             // initial number of columns in this motif
	double* pseudo;                          // pseudocounts to use for score calculations
	int width;															 // width of the motif (including non-informative columns)
  int num_seqs;                            // total number of sequences in this set
  int max_width;                           // maximum width of this motif
  vector<Site> sitelist;                   // list of sites that comprise this motif
  vector<int> columns;                     // columns in this motif
	int num_seqs_with_sites;                 // number of sequences with sites
	vector<int> has_sites;                   // number of sites in each sequence
	vector<bool> possible;                   // whether or not each sequence is in search space
	
	double motif_score;                      // score of this motif
	int above_seqc;													 // number of sequences above sequence threshold
	int ssp_size;                            // number of sequences in search space
  int above_cutoffs;                       // number of sequences above both cutoffs
	double seq_cutoff;                       // sequence cutoff for this motif
	double expr_cutoff;                      // expression cutoff for this motif
	double score_cutoff;                     // score cutoff for this motif
	string iter;                             // iteration in which this motif was found
	int dejavu;                              // number of times this motif was seen

public:
	Motif();
  Motif(const Seqset& v, int nc, double* p);
  Motif(const Motif& m);
	Motif& operator= (const Motif& m);
  int number() const { return sitelist.size(); }
  int get_width() const { return width; }
  int ncols() const { return columns.size(); }
  int chrom(int i) const { return sitelist[i].chrom(); }
  int posit(int i) const { return sitelist[i].posit(); }
  bool strand(int i) const { return sitelist[i].strand(); }
  int get_max_width() const { return max_width; }
	bool is_open_site(const int c, const int p);
	int seqs_with_sites() const { return num_seqs_with_sites; }
	bool seq_has_site(const int c) const { return (has_sites[c] != 0); }
	double get_motif_score() const { return motif_score; }
	void set_motif_score(const double sc) { motif_score = sc; }
	void set_above_seqc(const int sc) { above_seqc = sc; }
	int get_above_seqc() const { return above_seqc; }
	int get_search_space_size() const { return ssp_size; }
	bool in_search_space(const int g) { return possible[g]; }
	void clear_search_space();                                 // Remove all genes from search space
	void add_to_search_space(const int g);                     // Add gene to search space
	void remove_from_search_space(const int g);                // Remove gene from search space
  void set_above_cutoffs(const int ac) { above_cutoffs = ac; }
  int get_above_cutoffs() const { return above_cutoffs; }
 	double get_seq_cutoff() const { return seq_cutoff; }
	void set_seq_cutoff(const double cutoff) { seq_cutoff = cutoff; }
	double get_expr_cutoff() const { return expr_cutoff; }
	void set_expr_cutoff(const double cutoff) { expr_cutoff = cutoff; }
	double get_score_cutoff() const { return score_cutoff; }
	void set_score_cutoff(const double cutoff) { score_cutoff = cutoff; }
	string get_iter() const { return iter; }
	void set_iter(const string it) { iter = it; }
	int get_dejavu() { return dejavu; }
	void set_dejavu(const int d) { dejavu = d; }
	void inc_dejavu() { dejavu++; }
	void add_site(const int c, const int p, const bool s);
  void clear_sites();
	void remove_all_sites();
  void calc_freq_matrix(int* fm) const;
  void freq_matrix_extended(vector<float>& fm) const;
	void calc_score_matrix(double* sm, double* pseudo) const;
	double score_site(double* score_matrix, const int c, const int p, const bool s) const;
	double compare(const Motif& other);
  int column(const int i) const { return columns[i]; };
	vector<int>::iterator first_column() { return columns.begin(); };
	vector<int>::iterator last_column() { return columns.end(); };
  bool column_freq(const int col, int *ret);
  int remove_col(const int c);
	void add_col(const int c);
	bool has_col(const int c);
  void shift_sites(const int shift);
  void flip_sites();
	void orient();
  int total_positions() const;
	int positions_in_search_space() const;
  void columns_open(int &l, int &r);
	string consensus() const;                                                // Return the consensus sequence for the current set of sites
	void read(istream& motin);                                               // Read list of sites from a stream
	void write(ostream& motout) const;                                       // Write list of sites to a stream
	void destroy();
	void print_columns(ostream& out);
	bool check_sites();
	void check_possible();
};

#endif

