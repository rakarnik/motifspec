#include "motif.h"

Motif::Motif(const Seqset& s, const int nc, const vector<double>& p, const vector<double>& b) :
seqset(s),
init_nc(nc),
pseudo(p),
backfreq(b),
num_seqs(s.num_seqs()),
max_width(3 * init_nc),
columns(init_nc),
num_seqs_with_sites(0),
has_sites(num_seqs),
possible(num_seqs, false),
motif_score(0.0),
above_seqc(0),
ssp_size(0),
above_cutoffs(0),
seq_cutoff(0.00001),
ssp_cutoff(0.70),
dejavu(0)
{
	vector<int>::iterator cb = columns.begin();
	for(vector<int>::iterator ci = columns.begin(), ce = columns.end(); ci != ce; ++ci) {
		*ci = distance(cb, ci);
	}
	width = ((int) columns.back()) + 1;
	iter = "0";
}

Motif::Motif(const Motif& m) :
seqset(m.seqset),
init_nc(m.init_nc),
pseudo(m.pseudo),
backfreq(m.backfreq),
width(m.width),
num_seqs(m.num_seqs),
max_width(m.max_width),
columns(m.columns),
num_seqs_with_sites(m.num_seqs_with_sites),
has_sites(m.has_sites),
possible(m.possible),
motif_score(m.motif_score),
above_seqc(m.above_seqc),
ssp_size(m.ssp_size),
above_cutoffs(m.above_cutoffs),
seq_cutoff(m.seq_cutoff),
ssp_cutoff(m.ssp_cutoff),
iter(m.iter),
dejavu(m.dejavu)
{
	*this = m;
}

Motif& Motif::operator= (const Motif& m) {
	if(this != &m){
		//assume that the same Seqset is referred to, so ignore some things
		init_nc = m.init_nc;
		sitelist.assign(m.sitelist.begin(), m.sitelist.end());
		num_seqs_with_sites = m.num_seqs_with_sites;
		has_sites.assign(m.has_sites.begin(), m.has_sites.end());
		possible.assign(m.possible.begin(), m.possible.end());
		columns.assign(m.columns.begin(), m.columns.end());
		width = m.width;
		motif_score = m.motif_score;
		above_seqc = m.above_seqc;
		ssp_size = m.ssp_size;
		above_cutoffs = m.above_cutoffs;
		seq_cutoff = m.seq_cutoff;
		ssp_cutoff = m.ssp_cutoff;
		iter = m.iter;
		dejavu = m.dejavu;
	}
	return *this;
}

void Motif::clear_sites() {
	sitelist.clear();
	vector<int>(init_nc).swap(columns);
	vector<int>::iterator cb = columns.begin();
	for(vector<int>::iterator ci = columns.begin(), ce = columns.end(); ci != ce; ++ci)
		*ci = distance(cb, ci);
	width = ((int) columns.back()) + 1;
	num_seqs_with_sites = 0;
	has_sites.assign(num_seqs, false);
	motif_score = 0.0;
}

void Motif::remove_all_sites() {
	sitelist.clear();
	num_seqs_with_sites = 0;
	for(int i = 0; i < num_seqs; i++) {
		has_sites[i] = false;
	}
}

bool Motif::is_open_site(const int c, const int p){
	vector<Site>::iterator site_iter = sitelist.begin();
	for(vector<Site>::iterator si = sitelist.begin(), se = sitelist.end(); si != se; ++si) {
		if(si->chrom() == c) {
			int pp = si->posit();
			if(pp > p - width && pp < p + width) return false;
		}
	}
	return true;
}

void Motif::clear_search_space() {
	vector<bool> newp(num_seqs, false);
	swap(possible, newp);
	ssp_size = 0;
}

void Motif::add_to_search_space(const int g) {
	assert(ssp_size < num_seqs);
	assert(! possible[g]);
	possible[g] = true;
	ssp_size++;
	assert(ssp_size <= num_seqs);
}

void Motif::remove_from_search_space(const int g) {
	assert(ssp_size > 0);
	assert(possible[g]);
	possible[g] = false;
	ssp_size--;
	assert(ssp_size >= 0);
}

void Motif::add_site(const int c, const int p, const bool s){
	assert(p >= 0 && p < seqset.len_seq(c));
	Site st(c, p, s);
	sitelist.push_back(st);
	if(has_sites[c] == 0) num_seqs_with_sites++;
	has_sites[c]++;
}

void Motif::column_freq(const int col, int *ret){
	const vector<vector <char> >& seq = seqset.seq();
	for(int i = 0; i < 4; i++) ret[i] = 0;
	int c, p, len;
	bool s;
	vector<Site>::iterator si = sitelist.begin();
	vector<Site>::iterator se = sitelist.end();
	for(; si != se; ++si) {
		c = si->chrom();
		p = si->posit();
		s = si->strand();
		len = seqset.len_seq(c);
		if(s) {
			assert(p + col >= 0);
			assert(p + col < len);
			ret[seq[c][p + col]]++;
		} else {
			assert(p + width - 1 - col >= 0);
			assert(p + width - 1 - col < len);
			ret[3 - seq[c][p + width - 1 - col]]++;
		}
	}
}

double Motif::score_site(double* score_matrix, const int c, const int p, const bool s) const {
	const vector<vector <char> >& seq = seqset.seq();
	double L = 0.0;
	int matpos;
	vector<int>::const_iterator ci= columns.begin();
	vector<int>::const_iterator ce = columns.end();
	if(s) {
		matpos = 0;
		for(; ci != ce; ++ci) {
			assert(p + *ci >= 0);
			assert(p + *ci < seqset.len_seq(c));
			L += score_matrix[matpos + seq[c][p + *ci]];
			matpos += 4;
		}
	} else {
		matpos = 0;
		for(; ci != ce; ++ci) {
			assert(p + width - 1 - *ci >= 0);
			assert(p + width - 1 - *ci < seqset.len_seq(c));
			L += score_matrix[matpos + 3 - seq[c][p + width - 1 - *ci]];
			matpos += 4;
		}
	}
	return L;
}

void Motif::add_col(const int c) {
	if(c == 0) {
		assert(columns.size() == 0);
		columns.push_back(c);
	} else if(c < 0) {									 // column to the left of existing
		assert(columns.size() > 0);
		// Shift columns
		vector<int>::iterator ci = columns.begin();
		vector<int>::iterator ce = columns.end();
		for(; ci != ce; ++ci)
			*ci -= c;
		// Shift sites on Watson left
		vector<Site>::iterator si = sitelist.begin();
		vector<Site>::iterator se = sitelist.end();
		for(; si != se; ++si) 
			if(si->strand()) 
				si->shift(c);
		columns.insert(columns.begin(), 0);
	} else if(c > width - 1) {   // column to right of existing columns
		int shift = c - columns.back();
		// Shift sites on Crick strand left
		vector<Site>::iterator si = sitelist.begin();
		vector<Site>::iterator se = sitelist.end();
		for(; si != se; ++si)
			if(! si->strand())
				si->shift(-shift);
		columns.push_back(c);
	} else {
		bool found = false;
		vector<int>::iterator ci = columns.begin();
		vector<int>::iterator ce = columns.end();
		for(; ci != ce; ++ci){
			if(*ci > c) {
				columns.insert(ci, c);
				found = true;
				break;
			}
		}
		assert(found);
	}
	width = ((int) columns.back()) + 1;
	check_sites();
}

void Motif::remove_col(const int c) {
	if(c == 0) {									            // column to be removed is the first column
		columns.erase(columns.begin());
		if(columns.size() > 0) {
			int shift = columns.front();
			// Shift columns over to make first column 0
			vector<int>::iterator ci = columns.begin();
			vector<int>::iterator ce = columns.end();
			for(; ci != ce; ++ci)
				*ci -= shift;
			// Shift sites on Watson strand right
			vector<Site>::iterator si = sitelist.begin();
			vector<Site>::iterator se = sitelist.end();
			for(; si != se; ++si) 
				if(si->strand()) 
					si->shift(shift);
		}
	} else if(c == columns.back()) {          // column to be removed is the last column
		// Shift sites on Crick strand right
		int shift = columns.back() - columns[columns.size() - 2];
		vector<Site>::iterator si = sitelist.begin();
		vector<Site>::iterator se = sitelist.end();
		for(; si != se; ++si)
			if(! si->strand())
				si->shift(shift);
		columns.pop_back();
	} else {
		bool found = false;
		vector<int>::iterator ci = columns.begin();
		vector<int>::iterator ce = columns.end();
		for(; ci != ce; ++ci) {
			if(*ci == c) {
				columns.erase(ci);
				found = true;
				break;
			}
		}
		assert(found);
	}
	width = ((int) columns.back()) + 1;
	check_sites();
}

bool Motif::has_col(const int c) {
	return binary_search(columns.begin(), columns.end(), c);
}

bool Motif::column_sample(const bool add, const bool remove){
	bool changed = false;
	int freq[4];
	// Compute scores for current and surrounding columns
	int max_l, max_r;
	max_l = max_r = (max_width - width)/2;
	columns_open(max_l, max_r);
	int cs_span = max_l + max_r + width;
	
	// Compute information content for each column
	vector<struct idscore> wtx(cs_span);
	double best_wt = -DBL_MAX;
	for(int i = 0; i < cs_span; i++) {
		column_freq(i - max_l, freq);
		wtx[i].id = i;
		wtx[i].score = 0.0;
		for(int j = 0;j < 4; j++){
			wtx[i].score += gammaln(freq[j] + pseudo[j]);
			wtx[i].score -= (double)freq[j] * log(backfreq[j]);
		}
	}
	
	// Penalize outermost columns for length
	vector<struct idscore>::iterator wi = wtx.begin();
	int newwidth;
	for(; wi != wtx.end(); ++wi) {
		newwidth = width;
		if(wi->id < max_l)
			newwidth += max_l - wi->id;
		else if(wi->id > max_l + width - 1)
			newwidth += wi->id - max_l - width + 1;
		wi->score -= lnbico(newwidth - 2, ncols() - 2);
		if(wi->score > best_wt) best_wt = wi->score;
	}
	
	// Scale and sort columns by weight
	for(wi = wtx.begin(); wi != wtx.end(); ++wi)
		wi->score /= best_wt;
	sort(wtx.begin(), wtx.end(), isc);
	
	int nc = ncols();
	if(add) nc++;
	if(remove) nc--;
	int nadded = 0;
	vector<struct idscore>::iterator wi1 = wtx.begin();
	// Add top nc columns
	for(wi = wtx.begin(); wi != wtx.end() && nadded < nc; ++wi) {
		if(! has_col(wi->id - max_l)) {     // column is not already in motif
			add_col(wi->id - max_l);
			changed = true;
			if(wi->id - max_l < 0)                  // new column is to left of all existing columns, adjust column numbers
				for(wi1 = wi + 1; wi1 != wtx.end(); ++wi1)
					wi1->id -= wi->id - max_l;
		}
		nadded++;
	}
	
	// Remove any columns not in top nc but currently in motif
	for(; wi != wtx.end(); ++wi) {
		if(has_col(wi->id - max_l)) {
			if(wi->id - max_l == 0)
				for(wi1 = wi + 1; wi1 != wtx.end(); ++wi1)
					wi1->id -= column(1);
			remove_col(wi->id - max_l);
			changed = true;
		}
	}
	
	return changed;
}

void Motif::flip_sites() {
	int i;
	for(i = 0; i < number(); i++) {
		sitelist[i].flip();
	}
	vector<int>::iterator ci = columns.begin();
	vector<int>::iterator ce = columns.end();
	for(; ci != ce; ++ci) {
		*ci = width - *ci - 1;
	}
	vector<int> temp;
	vector<int>::reverse_iterator ri = columns.rbegin();
	vector<int>::reverse_iterator re = columns.rend();
	for(; ri != re; ++ri) {
		temp.push_back(*ri);
	}
	columns.assign(temp.begin(), temp.end());
}

void Motif::orient() {
	float* freq_matrix = new float[4 * ncols()];
	double *freq = new double[4];
	for(int i = 0; i < 4; i++)
		freq[i] = 0.0;
	
	calc_freq_matrix(freq_matrix);
	for(int i = 0; i < 4 * ncols(); i += 4) {
		for(int j = 0; j < 4; j++)
			freq[j] += freq_matrix[i + j];
	}
	for(int i = 0; i < 4; i++)
		freq[i] /= number();
	
	double flip = 1.5 * freq[2] + 1.0 * freq[0] - 1.0 * freq[3] - 1.5 * freq[1];
	if(flip < 0.0) {
		cerr << "\t\t\t\tFlipping motif...\n";
		flip_sites();
	}
	delete [] freq_matrix;
	delete [] freq;
}

int Motif::total_positions() const {
	int ret = 0;
	for(int i = 0; i < num_seqs; i++)
		ret += seqset.len_seq(i) - width + 1;
	return ret;
}

int Motif::positions_in_search_space() const {
	int ret = 0;
	for(int i = 0; i < num_seqs; i++)
		if(possible[i])
			ret += seqset.len_seq(i) - width + 1;
	return ret;
}

void Motif::columns_open(int &l, int &r){
	int c, p, len;
	bool s;
	vector<Site>::iterator si = sitelist.begin();
	vector<Site>::iterator se = sitelist.end();
	for(; si != se; ++si) {
		c = si->chrom();
		p = si->posit();
		s = si->strand();
		len = seqset.len_seq(c);
		if(s) {
			l = min(p, l);
			r = min(len - p - width, r);
		} else {
			l = min(len - p - width, l);
			r = min(p, r);
		}
		assert(l >= 0);
		assert(r >= 0);
	}
}

void Motif::calc_freq_matrix(float* fm) const {
	for(int i = 0; i < 4 * ncols(); i++){
		fm[i] = 0;
	}
	
	vector<float> w(num_seqs, 1.0);
	calc_freq_matrix(fm, w);
}

void Motif::calc_freq_matrix(float* fm, const vector<float>& w) const {
	const vector<vector <char> >& seq = seqset.seq();
	for(int i = 0; i < 4 * ncols(); i++) {
		fm[i] = 0.0;
	}
	
	int g, j, pos, matpos, len;
	bool s;
	vector<Site>::const_iterator si = sitelist.begin();
	vector<Site>::const_iterator se = sitelist.end();
	for(; si != se; ++si) {
		g = si->chrom();
		j = si->posit();
		s = si->strand();
		len = seqset.len_seq(g);
		if(s) {															 // forward strand
			matpos = 0;
			pos = 0;
			vector<int>::const_iterator ci = columns.begin();
			vector<int>::const_iterator ce = columns.end();
			for(; ci != ce; ++ci) {
				pos = j + *ci;
				assert(pos >= 0);
				assert(pos < len);
				fm[matpos + seq[g][pos]] += w[g];
				matpos += 4;
			}
		} else {														 // reverse strand
			matpos = 0;
			pos = 0;
			vector<int>::const_iterator ci = columns.begin();
			vector<int>::const_iterator ce = columns.end();
			for(; ci != ce; ++ci) {
				pos = j + width - 1 - *ci;
				assert(pos >= 0);
				assert(pos < len);
				fm[matpos + 3 - seq[g][pos]] += w[g];
				matpos += 4;
			}
		}
	}
}

void Motif::freq_matrix_extended(vector<float>& fm) const {
	const vector<vector <char> >& seq = seqset.seq();
	int i, col, j;
	int fm_size = fm.size();
	for(i = 0; i < fm_size; i++) fm[i] = 0.0;
	if(number() == 0) return;
	for(i = 0; i < number(); i++) {//i = site number
		int c = chrom(i);
		int p = posit(i);
		bool s = strand(i);
		for(j = 0, col = -ncols(); col < width + ncols(); col++, j += 4) {
			if(s) {
				if((p + col <= seqset.len_seq(c) - 1) && (p + col >= 0)) {
					assert(j + seq[c][p + col] < fm_size);
					fm[j + seq[c][p + col]] += 1.0;
				} else {
					for(int k = 0; k < 4; k++) {
						fm[j + k] += 0.25;
					}
				}
			} else {
				if((p + width - 1 - col <= seqset.len_seq(c) - 1) && (p + width - 1 - col >= 0)) {
					assert(j + 3 - seq[c][p + width - 1 - col] < fm_size);
					fm[j + 3 - seq[c][p + width - 1 - col]] += 1.0;
				} else {
					for(int k = 0; k < 4; k++) {
						fm[j + k] += 0.25;
					}
				}
			}
		}
	}
	for(i = 0; i < fm_size; i++) fm[i] /= (double) number();
}

void Motif::calc_score_matrix(double *sm) const {
	vector<float> w(num_seqs, 1.0);
	calc_score_matrix(sm, w);
}

void Motif::calc_score_matrix(double *sm, const vector<float>& w) const {
	float* fm = new float[4 * ncols()];
	double tot = (double) number();
	for(int j = 0; j < 4; j++) {
		tot += pseudo[j];
	}
	calc_freq_matrix(fm, w);
	for(int i = 0; i < 4 * ncols(); i += 4){
		for(int j = 0; j < 4; j++){
			sm[i + j] = log((fm[i + j] + pseudo[j])/tot);
		}
	}
	delete [] fm;
}

string Motif::consensus() const {
	int numsites = number();
	if(numsites < 1) return "";
	char nt[] = {'A', 'C', 'G', 'T'};
	const vector<vector <char> >& seq = seqset.seq();
	vector<string> hits;
	hits.reserve(numsites);
	vector<Site>::const_iterator si = sitelist.begin();
	vector<Site>::const_iterator se = sitelist.end();
	int c, p;
	bool s;
	for(; si != se; ++si) {
		c = si->chrom();
		p = si->posit();
		s = si->strand();
		if(s) {
			string hitseq = "";
			for(int k = p; k < p + width; k++)
				hitseq.append(1, nt[(int) seq[c][k]]);
			hits.push_back(hitseq);
		} else {
			string hitseq = "";
			for(int k = p + width - 1; k > p - 1; k--)
				hitseq.append(1, nt[3 - (int) seq[c][k]]);
			hits.push_back(hitseq);
		}
	}
	
	string cons;
	int wide1 = hits[0].size();
	int num1 = hits.size();
	int num1a, num1c, num1g, num1t;
	for(int i = 0; i < wide1; i++){
		num1a = num1c = num1g = num1t = 0;
		for(int j = 0; j < num1; j++){
			if(hits[j][i] == 'a' || hits[j][i] == 'A') num1a += 1;
			if(hits[j][i] == 'c' || hits[j][i] == 'C') num1c += 1;
			if(hits[j][i] == 'g' || hits[j][i] == 'G') num1g += 1;
			if(hits[j][i] == 't' || hits[j][i] == 'T') num1t += 1;
		}
		if(num1a > num1*0.7) cons += 'A';
		else if(num1c > num1*0.7) cons += 'C';
		else if(num1g > num1*0.7) cons += 'G';
		else if(num1t > num1*0.7) cons += 'T';
		else if((num1a + num1t) > num1 * 0.85) cons +='W';
		else if((num1a + num1c) > num1 * 0.85) cons +='M';
		else if((num1a + num1g) > num1 * 0.85) cons +='R';
		else if((num1t + num1c) > num1 * 0.85) cons +='Y';
		else if((num1t + num1g) > num1 * 0.85) cons +='K';
		else if((num1c + num1g) > num1 * 0.85) cons +='S';
		else cons += '-';
	}
	return cons;
}

void Motif::write(ostream& motout) const {
	char nt[] = {'A', 'C', 'G', 'T'};
	const vector<vector <char> >& seq = seqset.seq();
	vector<Site>::const_iterator si = sitelist.begin();
	vector<Site>::const_iterator se = sitelist.end();
	for(; si != se; ++si) {
		int c = si->chrom();
		int p = si->posit();
		bool s = si->strand();
		for(int j = 0; j < width; j++){
			if(s) {
				if(p + j >= 0 && p + j < seqset.len_seq(c))
					motout << nt[(int) (seq[c][p + j])];
				else motout << ' ';
			}
			else {
				if(p + width - 1 - j >= 0 && p + width - 1 - j < seqset.len_seq(c))
					motout << nt[3 - (int) seq[c][p + width - 1 - j]];
				else motout << ' ';
			}
		}
		motout << '\t' << c << '\t' << p << '\t' << s << '\n';
	}
	int col = 0, prev_col = 0;
	vector<int>::const_iterator ci = columns.begin();
	vector<int>::const_iterator ce = columns.end();
	for(; ci != ce; ++ci) {	
		col = *ci;
		for(int k = 0; k < (col - prev_col - 1); k++) motout << ' ';
		motout << '*';
		prev_col = col;
	}
	motout << endl;
	
	motout << "Score: " << motif_score << "\n";
	motout << "Sequences above sequence threshold: " << above_seqc << "\n";
	motout << "Size of search space: " << ssp_size << "\n";
	motout << "Sequence cutoff: " << seq_cutoff << "\n";
	motout << "Search space cutoff: " << ssp_cutoff << "\n";
	motout << "Iteration found: " << iter << "\n";
	motout << "Dejavu: " << dejavu << endl << "\n";
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
	int num_sites = c.size();
	for(int i = 0; i < num_sites; i++) {
		assert(p[i] >= 0);
		assert(p[i] + width - 1 < seqset.len_seq(c[i]));
		add_site(c[i], p[i], s[i]);
	}
	
	char* heading;
	
	// Read score
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_motif_score(atof(strtok(NULL, "\0")));
	
	// Read number of sequences above sequence threshold
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_above_seqc(atoi(strtok(NULL, "\0")));
	
	// Read size of search space
	motin.getline(line, 200);
	heading = strtok(line, ":");
	ssp_size = atoi(strtok(NULL, "\0"));
	
	// Read sequence cutoff
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_seq_cutoff(atof(strtok(NULL, "\0")));
	
	// Read expression cutoff
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_ssp_cutoff(atof(strtok(NULL, "\0")));
	
	// Read iteration found
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_iter(strtok(NULL, "\0"));
	
	// Read dejavu
	motin.getline(line, 200);
	heading = strtok(line, ":");
	set_dejavu(atoi(strtok(NULL, "\0")));
}

void Motif::print_columns(ostream& out) {
	vector<int>::iterator col_iter;
	for(col_iter = columns.begin(); col_iter != columns.end(); ++col_iter)
		out << " " << *col_iter;
}

void Motif::check_sites() {
	int c, p;
	bool s;
	vector<Site>::iterator si = sitelist.begin();
	vector<Site>::iterator se = sitelist.end();
	for(; si != se; ++si) {
		c = si->chrom();
		p = si->posit();
		s = si->strand();
		assert(p >= 0);
		assert(p + width - 1 < seqset.len_seq(c));
	}
}

void Motif::check_possible() {
	int poss_cnt = 0;
	vector<bool>::iterator pi = possible.begin();
	vector<bool>::iterator pe = possible.end();
	for(; pi != pe; ++pi)
		if(*pi)
			poss_cnt++;
	assert(poss_cnt == ssp_size);
}
