//Copyright 1999 President and Fellows of Harvard University
//sites.h

#ifndef _sites
#define _sites
#include "standard.h"
#include "seqset.h"

class Seqset;

class Site {
	int chr;
	int pos;
	bool str;
public:
	Site() {};
	Site(const int chr, const int pos, const bool str) { this->chr = chr; this->pos = pos; this->str = str; };
	int chrom() const { return chr; };
	int posit() const { return pos; };
	bool strand() const { return str; };
	void shift(const int shift) { pos += shift; };
	void flip() { str = ! str; };
};

class Sites{
  int sites_depth;//ie 4 for acgt, 6 for nacgtn

	const Seqset& seqset;                          // set of sequences that this motif refers to
  int sites_num_seqs;                            // total number of sequences in this set
  vector<int> sites_len_seq;                     // length of the sequences in this motif
  int sites_max_width;                           // maximum width of this motif
  vector<Site> sitelist;                         // list of sites that comprise this motif
	int sites_num_seqs_with_sites;                 // number of sequences with sites
	vector<int> sites_has_sites;                   // number of sites in each sequence
  vector<int> columns;                           // columns in this motif
  
	double mapsc;
	double spec;
	int iter;
	int dejavu;
	double seq_cutoff;
	double expr_cutoff;

public:
  Sites(const Seqset& v, int nc = 10, int memx = 1, int dp = 6);
  Sites(const Sites& s);
	Sites& operator= (const Sites& s);
  int number() const { return sitelist.size(); }
  int width() const { return columns.back() + 1; }
  int depth() const { return sites_depth; }
  int ncols() const { return columns.size(); }
  int chrom(int i) const { return sitelist[i].chrom(); }
  int posit(int i) const { return sitelist[i].posit(); }
  bool strand(int i) const { return sitelist[i].strand(); }
  int max_width() const { return sites_max_width; }
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
	int seqs_with_sites() const { return sites_num_seqs_with_sites; }
	bool seq_has_site(const int c) const { return (sites_has_sites[c] != 0); }
  void add_site(const int c, const int p, const bool s);
  void clear_sites();
	void remove_all_sites();
  void calc_freq_matrix(int *fm);
  void freq_matrix_extended(double *fm) const;
  int column(const int i) const {return columns.at(i);}
  bool column_freq(const int col, int *ret);
  int remove_col(const int c);
	void add_col(const int c);
	bool has_col(const int c);
  int positions_available() const;
	int positions_available(const bool* possible) const;
  void shift_sites(const int l, const int r);
  void flip_sites();
	void orient();
  void columns_open(int &l, int &r);
	void write(ostream& motout) const;
	void read(istream& motin);
  void destroy();
	void print_columns(ostream& out);
};

#endif
