//Copyright 1999 President and Fellows of Harvard University
//sites.h

#ifndef _sites
#define _sites
#include "standard.h"
#include "seqset.h"

class Seqset;
class Sites{
	int sites_num;
  int sites_width;
  int sites_num_cols;
  int sites_depth;//ie 4 for acgt, 6 for nacgtn

  int sites_num_seqs; //need to deallocate memory
  int *sites_len_seq; //need for copy constructor
  int sites_max_num_sites;//used for copy assignment
  int sites_max_width;

  int *sites_chrom;
  int *sites_posit;
  bool *sites_strand;
  //these three are an unsorted compact list, just need to sort for output

	int sites_num_seqs_with_sites;
	int* sites_has_sites;

  int *sites_active_fwd;
  //columns 0..wide-1;fwd(0)=2nd column
  bool sites_alloc;

public:
  Sites(){sites_alloc=false;destroy();}
  void sites_init(const Sites& s);
  Sites(const vector<string>& v, int nc=10, int memx=1, int dp=6);
  Sites(const Sites& s);
	void init(const vector<string>& v, int nc=10, int memx=1, int dp=6);
  Sites& operator= (const Sites& s);
  void allocate_mem();
  ~Sites();
  int number() const {return sites_num;}
  int width() const {return sites_width;}
  int depth() const {return sites_depth;}
  int ncols() const {return sites_num_cols;}
  int chrom(int i) const {return sites_chrom[i];}
  int posit(int i) const {return sites_posit[i];}
  bool strand(int i) const {return sites_strand[i];}
  int max_width() const {return sites_max_width;}
  bool is_open_site(const int c, const int p);
	int seqs_with_sites() const { return sites_num_seqs_with_sites; }
	bool seq_has_site(const int c) const { return (sites_has_sites[c] > 0); }; 
  void add_site(const int c, const int p, const bool s);
  void remove_site(const int c, const int p);
  void remove_all_sites();
  void clear_sites();
  void calc_freq_matrix(const Seqset& b, int *fm);
  void freq_matrix_extended(const Seqset& b, double *fm) const;
  int next_column(const int i){return sites_active_fwd[i];}
  bool column_freq(const int col, const Seqset& s, int *ret);
  int remove_col(const int c);
  void add_col(const int c);
  int positions_available() const;
	int positions_available(const bool* possible) const;
  void shift_sites(const int l, const int r);
  void flip_sites();
  void columns_open(int &l, int &r);
  void destroy();
};

double CompareACE(const Sites &s1, const Sites &s2, const Seqset &t1, const Seqset &t2);

#endif
