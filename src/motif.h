#ifndef _motif
#define _motif
#include "standard.h"
#include "seqset.h"
#include "site.h"

class Motif {
	const Seqset& seqset;                    // set of sequences that this motif refers to
	int width;															 // width of the motif (including non-informative columns)
  int npseudo;                             // pseudocount
  int num_seqs;                            // total number of sequences in this set
  int max_width;                           // maximum width of this motif
  vector<Site> sitelist;                   // list of sites that comprise this motif
  vector<int> columns;                     // columns in this motif
	int num_seqs_with_sites;                 // number of sequences with sites
	vector<int> has_sites;                   // number of sites in each sequence
  
	double mapsc;
	double spec;
	int iter;
	int dejavu;
	double seq_cutoff;
	double expr_cutoff;

public:
	Motif();
  Motif(const Seqset& v, int nc = 10, int np = 1);
  Motif(const Motif& m);
	Motif& operator= (const Motif& m);
  int number() const { return sitelist.size(); }
  int get_width() const { return width; }
  int ncols() const { return columns.size(); }
  int chrom(int i) const { return sitelist[i].chrom(); }
  int posit(int i) const { return sitelist[i].posit(); }
  bool strand(int i) const { return sitelist[i].strand(); }
  int get_max_width() const { return max_width; }
	int get_iter() const { return iter; }
	void set_iter(const int it) { iter = it; }
	int get_dejavu() { return dejavu; }
	void set_dejavu(const int d) { dejavu = d; }
	void inc_dejavu() { dejavu++; }
	double get_seq_cutoff() const { return seq_cutoff; }
	void set_seq_cutoff(const double cutoff) { seq_cutoff = cutoff; }
	double get_expr_cutoff() const { return expr_cutoff; }
	void set_expr_cutoff(const double cutoff) { expr_cutoff = cutoff; }
	double get_map() const { return mapsc; }
	void set_map(const double sc) { mapsc = sc; }
  double get_spec() const { return spec; }
	void set_spec(const double s) { spec = s; }
	bool is_open_site(const int c, const int p);
	int seqs_with_sites() const { return num_seqs_with_sites; }
	bool seq_has_site(const int c) const { return (has_sites[c] != 0); }
  void add_site(const int c, const int p, const bool s);
  void clear_sites();
	void remove_all_sites();
  void calc_freq_matrix(int *fm);
  void freq_matrix_extended(double *fm) const;
	void calc_score_matrix(double *sm, double* pseudo);
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
  int positions_available() const;
	int positions_available(const vector<bool>& possible) const;
  void columns_open(int &l, int &r);
	string consensus() const;                                                // Return the consensus sequence for the current set of sites
	void read(istream& motin);                                               // Read list of sites from a stream
	void write(ostream& motout) const;                                       // Write list of sites to a stream
	void destroy();
	void print_columns(ostream& out);
	bool check_sites();
};

#endif

