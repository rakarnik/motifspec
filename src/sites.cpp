//Copyright 1999 President and Fellows of Harvard University
//sites.cpp

#include "sites.h"

Sites::Sites(const Seqset& s, int nc, int mx, int dp) :
sites_depth(dp),
sites_num_seqs(s.num_seqs()),
sites_len_seq(sites_num_seqs),
sites_max_width(2 * nc),
columns(nc),
sites_num_seqs_with_sites(0),
sites_has_sites(sites_num_seqs)
{
	for(int i = 0; i < columns.size(); i++) {
		columns[i] = i;
	}
  for(int i = 0; i < sites_num_seqs; i++){
    sites_len_seq[i] = s.len_seq(i);
  }
	seq_cutoff = 0.00001;
	expr_cutoff = 0.70;
	iter = 0;
	dejavu = 0;
	spec = 0.0;
}

Sites::Sites(const Sites& s) :
sites_depth(s.sites_depth),
sites_num_seqs(s.sites_num_seqs),
sites_len_seq(s.sites_len_seq),
sites_max_width(s.sites_max_width),
columns(s.columns),
sites_num_seqs_with_sites(s.sites_num_seqs_with_sites),
sites_has_sites(s.sites_has_sites)
{
  mapsc = s.mapsc;
	spec = s.spec;
  iter = s.iter;
	dejavu = s.dejavu;
	seq_cutoff = s.seq_cutoff;
	expr_cutoff = s.expr_cutoff;
  *this=s;
}

Sites& Sites::operator= (const Sites& s) {
  if(this != &s){
    //assume that the same Seqset is referred to, so ignore some things
    sites_depth = s.sites_depth;
		sitelist.assign(s.sitelist.begin(), s.sitelist.end());
		sites_num_seqs_with_sites = s.sites_num_seqs_with_sites;
		sites_has_sites.assign(s.sites_has_sites.begin(), s.sites_has_sites.end());
    columns.assign(s.columns.begin(), s.columns.end());
		seq_cutoff = s.seq_cutoff;
		expr_cutoff = s.expr_cutoff;
		mapsc = s.mapsc;
		spec = s.spec;
		iter = s.iter;
		dejavu = s.dejavu;
  }
  return *this;
}

void Sites::clear_sites(){
  sitelist.clear();
	for(int i = 0; i < columns.size(); i++) {
		columns[i] = i;
	}
  sites_num_seqs_with_sites = 0;
	for(int i = 0; i < sites_num_seqs; i++) {
		sites_has_sites[i] = false;
	}
	mapsc = 0.0;
	spec = 0.0;
}

void Sites::remove_all_sites() {
	sitelist.clear();
	sites_num_seqs_with_sites = 0;
	for(int i = 0; i < sites_num_seqs; i++) {
		sites_has_sites[i] = false;
	}
}

void Sites::write(const Seqset& seqset, ostream& motout) const {
	map<char,char> nt;
  nt[0] = nt[5] = 'N';
  nt[1] = 'A';
	nt[2] = 'C';
	nt[3] = 'G';
	nt[4] = 'T';
  char** ss_seq = seqset.seq_ptr();
  for(int i = 0; i < number(); i++){
    int c = chrom(i);
    int p = posit(i);
    bool s = strand(i);
    for(int j = 0; j < width(); j++){
      if(s) {
				if(p + j >= 0 && p + j < seqset.len_seq(c))
					motout << nt[ss_seq[c][p + j]];
				else motout << ' ';
      }
      else {
				if(p + width() - 1 - j >= 0 && p + width()-1-j < seqset.len_seq(c))
					motout << nt[depth() - 1 - ss_seq[c][p + width() - 1 - j]];
				else motout << ' ';
      }
    }
    motout << '\t' << c << '\t' << p << '\t' << s << '\n';
  }
	int col = 0, prev_col = 0;
	vector<int>::const_iterator col_iter;
	for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter) {	
		col = *col_iter;
		motout << '*';
		for(int k = 0; k < (col - prev_col - 1); k++) motout << ' ';
		prev_col = col;
	}
  motout << endl;
	
	motout << "MAP Score: " << mapsc << endl;
	motout << "Specificity Score: " << spec << endl;
	motout << "Sequence cutoff: " << seq_cutoff << endl;
	motout << "Expression cutoff: " << expr_cutoff << endl;
	motout << "Iteration found: " << iter << endl;
	motout << "Dejavu: " << dejavu << endl << endl;
}

void Sites::read(istream& motin) {
	char* match;
	char line[200];
	vector<int> c;
	vector<int> p;
	vector<bool> s;
	
	// Read sites
	// (don't add yet, as they will get screwed up by the column changes)
	while(motin.getline(line, 200)) {
		if(line[0] == '*') break;
		match = strtok(line, "\t");
		c.push_back(atoi(strtok(NULL, "\t")));
		p.push_back(atoi(strtok(NULL, "\t")));
		s.push_back(atoi(strtok(NULL, "\0")));
	}
	
	int motwidth = strlen(line);
	columns.clear();
	for(int i = 0; i < motwidth; i++) {
		if(line[i] == '*') add_col(i);
	}
	
	// Add sites
	sitelist.clear();
	for(int i = 0; i < c.size(); i++) {
		assert(p[i] >= 0);
		add_site(c[i], p[i], s[i]);
	}
	
	char* heading;
	// Read MAP score
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_map(atof(strtok(NULL, "\0")));
	
	// Read specificity
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_spec(atof(strtok(NULL, "\0")));
	
	// Read sequence cutoff
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_seq_cutoff(atof(strtok(NULL, "\0")));
	
	// Read expression cutoff
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_expr_cutoff(atof(strtok(NULL, "\0")));
	
	// Read iteration found
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_iter(atoi(strtok(NULL, "\0")));
	
	// Read dejavu
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_dejavu(atoi(strtok(NULL, "\0")));
}

void Sites::destroy(){
  //for when the = operator can fail
  sites_max_width=0;
}

bool Sites::is_open_site(const int c, const int p){
	for(int i = 0; i < sitelist.size(); i++){
    if(chrom(i) == c) {
      int pp = posit(i);
      if(pp > p - width() && pp < p + width()) return false;
    }
  }
  return true;
}

void Sites::add_site(const int c, const int p, const bool s){
	Site st(c, p, s);
	sitelist.push_back(st);
	if(sites_has_sites[c] == 0) sites_num_seqs_with_sites++;
	sites_has_sites[c]++;
}

void Sites::calc_freq_matrix(const Seqset& b, int *fm){
	// fm will have allocation for depth() * ncols()
  char** ss_seq = b.seq_ptr();
  for(int i = 0; i < sites_depth * columns.size(); i++){
		fm[i] = 0;
  }
  for(int i = 0; i < number(); i++){ //i = site number
    int c = chrom(i);
    int p = posit(i);
    bool s = strand(i);
    int pos, matpos;
    if(s) {                              // forward strand
      matpos = 0;
			pos = 0;
			vector<int>::iterator col_iter;
			for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter) {
				pos = p + *col_iter;
				assert(p >= 0 && p < b.len_seq(c));
				int seq = ss_seq[c][pos];
				fm[matpos + seq]++;
				matpos += sites_depth;
      }
    } else {                             // reverse strand
      matpos = sites_depth - 1;
      pos = 0;
			vector<int>::iterator col_iter;
			for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter) {
				pos = p + width() - 1 - *col_iter;
				assert(p >= 0 && p < b.len_seq(c));
				int seq = ss_seq[c][p + width() - 1 - *col_iter];
				fm[matpos - seq]++;
				matpos += sites_depth;
      }
    }      
  }
}

bool Sites::column_freq(const int col, const Seqset& s, int *ret){
  char** ss_seq = s.seq_ptr();
  for(int i = 0; i < sites_depth; i++) ret[i] = 0;
  for(int i = 0; i < number(); i++){//i = site number
    int c = chrom(i);
    int p = posit(i);
    bool t = strand(i);
    if(t) {
      if( (p + col > s.len_seq(c) - 1) || (p + col < 0) ) return false;
      int seq = ss_seq[c][p + col];
      ret[seq]++;
    } else {
      if((p + width() - 1 - col > s.len_seq(c) - 1) || (p + width() - 1 - col < 0)) return false;
      int seq = ss_seq[c][p + width() - 1 - col];
      ret[sites_depth - seq - 1]++;
    }
  }
  return true;
}

int Sites::remove_col(const int c) {
  int ret = 0, ns = 0;         // return number of removed column in new numbering
	bool found = false;
	if(columns.size() == 0) {
		cerr << "remove_column called with no columns in motif" << endl;
		abort();
	}
	if(c == 0) {                   // if the column to be removed is the first column, we shift all columns so first is 0
    columns.erase(columns.begin());
		if(columns.size() > 0) {
			ns = columns[0];
			vector<int>::iterator col_iter;
			for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter)
				*col_iter -= ns;
		}
		ret = -ns;
		shift_sites(ns, 0);
		found = true;
  } else {
		vector<int>::iterator col_iter;
    for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter) {
      if(*col_iter == c) {
				columns.erase(col_iter);
				found = true;
				break;
			}
    }
	}
	if (! found) {
		cerr << "remove_column called for column " << c << " but it was not found!" << endl; 
		abort();
	}
  return ret;
}

void Sites::add_col(const int c) {
	int col, nxt, i;
	col = nxt = i = 0;
  if(c < 0){                   // column to the left of existing, shift existing to the right
		if(columns.size() > 0) {
			vector<int>::iterator col_iter;
			for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter){
				*col_iter -= c;
			}
		}
		columns.insert(columns.begin(), 0);
  } else {
		bool found = false;
		vector<int>::iterator col_iter;
		for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter){
			if(*col_iter > c) {
				columns.insert(col_iter, c);
				found = true;
				break;
			}
		}
		if(! found)
			columns.push_back(c);
  }
}

bool Sites::has_col(const int c) {
	return (find(columns.begin(), columns.end(), c) != columns.end());
}

void Sites::flip_sites(){
  int i;
  for(i = 0; i < number(); i++) {
    sitelist[i].flip();
  }
	int w = width();
	vector<int>::iterator col_iter;
	for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter) {
		*col_iter = w - *col_iter;
	}
	vector<int> temp;
	vector<int>::reverse_iterator rev_iter;
	for(rev_iter = columns.rbegin(); rev_iter != columns.rend(); ++rev_iter) {
		temp.push_back(*rev_iter);
	}
	columns.assign(temp.begin(), temp.end());
}

void Sites::shift_sites(const int l, const int r){
  // numbers for right movement of beg/end point of forward site
  // l: + for shorter/right
  // r: + for longer/right
  for(int i = 0; i < number(); i++) {
    if(strand(i))
			sitelist[i].shift(l);
    else 
			sitelist[i].shift(-r);
  }
}

int Sites::positions_available() const {
  int ret = 0;
  for(int i = 0; i < sites_num_seqs; i++){
    ret += sites_len_seq[i] - width() + 1;
  }
  return ret;
}

int Sites::positions_available(const bool* possible) const {
	int ret = 0;
	for(int i = 0; i < sites_num_seqs; i++) {
		if(possible[i])
			ret += sites_len_seq[i] - width() + 1;
	}
	return ret;
}

void Sites::columns_open(int &l, int &r){
  //input r/l are max values to be reduced
  //only works if site list is sorted
  int i;
  int c_prev=-1, p_prev;
  bool s_prev;
  for(i = 0; i < number(); i++) {
    int c = chrom(i);
    int p = posit(i);
    bool s = strand(i);
    if(c == c_prev){
      int d = p - p_prev - width();
      if(s == s_prev){
				if(l > d) l = d;
				if(r > d) r = d;
      } else {
				if(l > d/2) l = d/2;
				if(r > d/2) r = d/2;
      }
    } else {
      int f;
      if(c_prev != -1){
				f = sites_len_seq[c_prev] - width() - p_prev;
				if(s_prev){
					if(r > f) r = f;
				} else {
					if(l > f) l = f;
				}
      }
      f=p;
      if(s){
				if(l>f) l=f;
      }
      else{
				if(r>f) r=f;
      }
    }
  }
}

void Sites::freq_matrix_extended(const Seqset& b, double *fm) const {
  //assumes fm of dimension (sites_width+2*sites_num_cols)*depth, fills in edges with Ns
  char** ss_seq=b.seq_ptr();
  int i,col,j,seq;
  int fm_size = (width() + 2 * ncols()) * sites_depth;
  for(i = 0; i < fm_size; i++) fm[i] = 0.0;
  if(number() == 0) return;
  for(i = 0; i < number(); i++) {//i = site number
    int c = chrom(i);
    int p = posit(i);
    bool t = strand(i);
    for(j = 0, col = -ncols(); col < width() + ncols(); col++, j += sites_depth) {
      if(t) {
				if((p + col > b.len_seq(c) - 1) || (p + col < 0)) seq = 0;
				else seq = ss_seq[c][p + col];
				fm[j + seq] += 1.0;
      } else {
				if((p + width() - 1 - col > b.len_seq(c) - 1) || (p + width() - 1 - col < 0)) seq = 0;
				else seq = ss_seq[c][p + width() - 1 - col];
				fm[j + sites_depth - 1 - seq] += 1.0;
      }
    }
  }
  for(i = 0; i < fm_size; i++) fm[i] /= (double) number();
}

void Sites::print_columns(ostream& out) {
	out << "\t\t\t\t\t";
	vector<int>::iterator col_iter;
	for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter)
		out << " " << *col_iter;
	out << endl;
}
