#include "motif.h"

Motif::Motif(const Seqset& s, int nc, int dp, int np) :
seqset(s),
depth(dp),
npseudo(np),
num_seqs(s.num_seqs()),
len_seq(num_seqs),
max_width(2 * nc),
columns(nc),
num_seqs_with_sites(0),
has_sites(num_seqs)
{
	for(int i = 0; i < columns.size(); i++) {
		columns[i] = i;
	}
  for(int i = 0; i < num_seqs; i++){
    len_seq[i] = s.len_seq(i);
  }
	seq_cutoff = 0.00001;
	expr_cutoff = 0.70;
	iter = 0;
	dejavu = 0;
	spec = 0.0;
}

Motif::Motif(const Motif& m) :
seqset(m.seqset),
depth(m.depth),
npseudo(m.npseudo),
num_seqs(m.num_seqs),
len_seq(m.len_seq),
max_width(m.max_width),
columns(m.columns),
num_seqs_with_sites(m.num_seqs_with_sites),
has_sites(m.has_sites)
{
  mapsc = m.mapsc;
	spec = m.spec;
  iter = m.iter;
	dejavu = m.dejavu;
	seq_cutoff = m.seq_cutoff;
	expr_cutoff = m.expr_cutoff;
  *this = m;
}

Motif& Motif::operator= (const Motif& m) {
  if(this != &m){
    //assume that the same Seqset is referred to, so ignore some things
    depth = m.depth;
		npseudo = m.npseudo;
		sitelist.assign(m.sitelist.begin(), m.sitelist.end());
		num_seqs_with_sites = m.num_seqs_with_sites;
		has_sites.assign(m.has_sites.begin(), m.has_sites.end());
    columns.assign(m.columns.begin(), m.columns.end());
		seq_cutoff = m.seq_cutoff;
		expr_cutoff = m.expr_cutoff;
		mapsc = m.mapsc;
		spec = m.spec;
		iter = m.iter;
		dejavu = m.dejavu;
  }
  return *this;
}

void Motif::clear_sites(){
  sitelist.clear();
	for(int i = 0; i < columns.size(); i++) {
		columns[i] = i;
	}
  num_seqs_with_sites = 0;
	for(int i = 0; i < num_seqs; i++) {
		has_sites[i] = false;
	}
	mapsc = 0.0;
	spec = 0.0;
}

void Motif::remove_all_sites() {
	sitelist.clear();
	num_seqs_with_sites = 0;
	for(int i = 0; i < num_seqs; i++) {
		has_sites[i] = false;
	}
}

void Motif::write(ostream& motout) const {
	map<char,char> nt;
  nt[0] = nt[5] = 'N';
  nt[1] = 'A';
	nt[2] = 'C';
	nt[3] = 'G';
	nt[4] = 'T';
  char** ss_seq = seqset.seq_ptr();
	vector<Site>::const_iterator site_iter;
  for(site_iter = sitelist.begin(); site_iter != sitelist.end(); ++site_iter) {
    int c = site_iter->chrom();
    int p = site_iter->posit();
    bool s = site_iter->strand();
    for(int j = 0; j < width(); j++){
      if(s) {
				if(p + j >= 0 && p + j < seqset.len_seq(c))
					motout << nt[ss_seq[c][p + j]];
				else motout << ' ';
      }
      else {
				if(p + width() - 1 - j >= 0 && p + width()-1-j < seqset.len_seq(c))
					motout << nt[depth - 1 - ss_seq[c][p + width() - 1 - j]];
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

void Motif::read(istream& motin) {
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

void Motif::destroy(){
  //for when the = operator can fail
  max_width=0;
}

bool Motif::is_open_site(const int c, const int p){
	for(int i = 0; i < sitelist.size(); i++){
    if(chrom(i) == c) {
      int pp = posit(i);
      if(pp > p - width() && pp < p + width()) return false;
    }
  }
  return true;
}

void Motif::add_site(const int c, const int p, const bool s){
	assert(p >= 0 && p < seqset.len_seq(c));
	Site st(c, p, s);
	sitelist.push_back(st);
	if(has_sites[c] == 0) num_seqs_with_sites++;
	has_sites[c]++;
}

void Motif::calc_freq_matrix(int *fm){
	// fm will have allocation for depth() * ncols()
  char** ss_seq = seqset.seq_ptr();
  for(int i = 0; i < depth * ncols(); i++){
		fm[i] = 0;
  }
	int w = width();
	int c, p;
	bool s;
  vector<Site>::iterator site_iter;
  for(site_iter = sitelist.begin(); site_iter != sitelist.end(); ++site_iter) {
    c = site_iter->chrom();
    p = site_iter->posit();
    s = site_iter->strand();
    int pos, matpos;
    if(s) {                              // forward strand
      matpos = 0;
			pos = 0;
			vector<int>::iterator col_iter;
			for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter) {
				pos = p + *col_iter;
				assert(pos >= 0 && pos < seqset.len_seq(c));
				int seq = ss_seq[c][pos];
				fm[matpos + seq]++;
				matpos += depth;
      }
    } else {                             // reverse strand
      matpos = depth - 1;
      pos = 0;
			vector<int>::iterator col_iter;
			for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter) {
				pos = p + w - 1 - *col_iter;
				assert(pos >= 0 && pos < seqset.len_seq(c));
				int seq = ss_seq[c][pos];
				fm[matpos - seq]++;
				matpos += depth;
      }
    }      
  }
}

bool Motif::column_freq(const int col, int *ret){
  char** ss_seq = seqset.seq_ptr();
  for(int i = 0; i < depth; i++) ret[i] = 0;
	vector<Site>::iterator site_iter;
  for(site_iter = sitelist.begin(); site_iter != sitelist.end(); ++site_iter) {
    int c = site_iter->chrom();
    int p = site_iter->posit();
    bool t = site_iter->strand();
    if(t) {
      if( (p + col > seqset.len_seq(c) - 1) || (p + col < 0) ) return false;
      int seq = ss_seq[c][p + col];
      ret[seq]++;
    } else {
      if((p + width() - 1 - col > seqset.len_seq(c) - 1) || (p + width() - 1 - col < 0)) return false;
      int seq = ss_seq[c][p + width() - 1 - col];
      ret[depth - seq - 1]++;
    }
  }
  return true;
}

int Motif::remove_col(const int c) {
  int ret = 0, ns = 0;         // return number of removed column in new numbering
	bool found = false;
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

void Motif::add_col(const int c) {
  if(c < 0){                   // column to the left of existing, shift existing to the right
		if(columns.size() > 0) {
			vector<int>::iterator col_iter;
			for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter){
				*col_iter -= c;
			}
			shift_sites(-c, 0);
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

bool Motif::has_col(const int c) {
	return binary_search(columns.begin(), columns.end(), c);
}

void Motif::flip_sites() {
  int i;
  for(i = 0; i < number(); i++) {
    sitelist[i].flip();
  }
	int w = width();
	vector<int>::iterator col_iter;
	for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter) {
		*col_iter = w - *col_iter - 1;
	}
	vector<int> temp;
	vector<int>::reverse_iterator rev_iter;
	for(rev_iter = columns.rbegin(); rev_iter != columns.rend(); ++rev_iter) {
		temp.push_back(*rev_iter);
	}
	columns.assign(temp.begin(), temp.end());
}

void Motif::shift_sites(const int l, const int r){
  // numbers for right movement of beg/end point of forward site
  // l: + for shorter/right
  // r: + for longer/right
	vector<Site>::iterator site_iter;
  for(site_iter = sitelist.begin(); site_iter != sitelist.end(); ++site_iter) {
    if(site_iter->strand())
			site_iter->shift(l);
    else 
			site_iter->shift(-r);
		assert(site_iter->posit() >= 0);
		assert(site_iter->posit() < seqset.len_seq(site_iter->chrom()));
  }
}

int Motif::positions_available() const {
  int ret = 0;
  for(int i = 0; i < num_seqs; i++){
    ret += len_seq[i] - width() + 1;
  }
  return ret;
}

int Motif::positions_available(const bool* possible) const {
	int ret = 0;
	for(int i = 0; i < num_seqs; i++) {
		if(possible[i])
			ret += len_seq[i] - width() + 1;
	}
	return ret;
}

void Motif::columns_open(int &l, int &r){
	int w = width();
	l = r = (max_width - w)/2;
	int c, p, len;
	bool s;
	vector<Site>::iterator site_iter;
  for(site_iter = sitelist.begin(); site_iter != sitelist.end(); ++site_iter) {
    c = site_iter->chrom();
    p = site_iter->posit();
    s = site_iter->strand();
		len = seqset.len_seq(c);
		if(c == 1603)
			cerr << "\t\t\t\t\tbefore: (" << l << ", " << r << ")";
		if(s) {
			l = min(p, l);
			r = min(len - p - w, r);
		} else {
			l = min(len - p - w, l);
			r = min(p, r);
		}
		if(c == 1603)
			cerr << "   after: (" << l << ", " << r << ")" << endl;
  }
}

void Motif::freq_matrix_extended(double *fm) const {
  //assumes fm of dimension (sites_width+2*sites_num_cols)*depth, fills in edges with Ns
  char** ss_seq = seqset.seq_ptr();
  int i, col, j, seq;
  int fm_size = (width() + 2 * ncols()) * depth;
  for(i = 0; i < fm_size; i++) fm[i] = 0.0;
  if(number() == 0) return;
  for(i = 0; i < number(); i++) {//i = site number
    int c = chrom(i);
    int p = posit(i);
    bool t = strand(i);
    for(j = 0, col = -ncols(); col < width() + ncols(); col++, j += depth) {
      if(t) {
				if((p + col > seqset.len_seq(c) - 1) || (p + col < 0)) seq = 0;
				else seq = ss_seq[c][p + col];
				fm[j + seq] += 1.0;
      } else {
				if((p + width() - 1 - col > seqset.len_seq(c) - 1) || (p + width() - 1 - col < 0)) seq = 0;
				else seq = ss_seq[c][p + width() - 1 - col];
				fm[j + depth - 1 - seq] += 1.0;
      }
    }
  }
  for(i = 0; i < fm_size; i++) fm[i] /= (double) number();
}

void Motif::calc_score_matrix(double *sm, double* backfreq) {
	int* fm = new int[depth * ncols()];
  double tot = (double) number() + npseudo;
  calc_freq_matrix(fm);
  for(int i = 0; i < depth * ncols();i += depth){
    sm[i] = sm[i + 5] = 1.0;
    for(int j = 1; j <= 4; j++){
      int x = fm[i] + fm[i + 5];
      sm[i + j] = (fm[i + j] + (x + npseudo) * backfreq[j])/(tot * backfreq[j]);
    }
  }
	delete [] fm;
}

void Motif::orient() {
	int* freq_matrix = new int[depth * ncols()];
  double *freq = new double[6];
	for(int i = 0; i < 6; i++)
		freq[i] = 0.0;
	
  calc_freq_matrix(freq_matrix);
  for(int i = 0; i < depth * ncols(); i += depth) {
    for(int j = 1; j <= 4; j++)
      freq[j] += freq_matrix[i + j];
  }
	for(int i = 0; i < 6; i++)
		freq[i] /= number();
	
  double flip = 1.5 * freq[3] + 1.0 * freq[1] - 1.0 * freq[4] - 1.5 * freq[2];
  if(flip < 0.0) flip_sites();
	delete [] freq_matrix;
  delete [] freq;
}

void Motif::print_columns(ostream& out) {
	vector<int>::iterator col_iter;
	for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter)
		out << " " << *col_iter;
}

bool Motif::check_sites() {
	int c, p;
	bool s;
	int w = width();
	vector<Site>::iterator site_iter;
	for(site_iter = sitelist.begin(); site_iter != sitelist.end(); ++site_iter) {
		c = site_iter->chrom();
		p = site_iter->posit();
		if(p < 0 || p + w - 1 >= seqset.len_seq(c)) {
			cerr << "\t\t\t\t\tc: " << c << " p: " << p << " w: " << w << " len: " << seqset.len_seq(c) << endl;
			return false;
		}
	}
	return true;
}
