#include "semodel.h"

SEModel::SEModel(const vector<string>& seqs, float** exprtab, const int numexpr, const vector<string>& names, const int nc, const int bf, const double map_cut, const double sim_cut):
seqset(seqs),
sites(seqset, nc, 5 * nc),
select_sites(seqset, nc),
print_sites(seqset, nc, 5 * nc),
archive(sites, seqset, map_cut, sim_cut),
ngenes(names.size()),
seqscores(ngenes),
seqranks(ngenes),
expscores(ngenes),
expranks(ngenes)
{
	npossible = 0;
	possible = new bool[ngenes];
	clear_all_possible();
	nameset = names;
	verbose = false;
	
	set_default_params();
  max_motifs = bf;
  map_cutoff = map_cut;
  sim_cutoff = sim_cut;
  freq_matrix = new int[sites.depth() * sites.ncols()];
  score_matrix = new double[sites.depth() * sites.ncols()];
	
	expr = exprtab;
	npoints = numexpr;
	mean = new float[npoints];
	stdev = new float[npoints];
	pcorr = new float*[ngenes];
	for(int i = 0; i < ngenes; i++) {
		pcorr[i] = new float[ngenes];
		for(int j = 0; j < ngenes; j++) {
			pcorr[i][j] = -2;
		}
	}
}

SEModel::~SEModel(){
  delete [] possible;
	delete [] freq_matrix;
  delete [] score_matrix;
	delete [] mean;
	delete [] stdev;
  for(int i=0;i<seqset.num_seqs();i++){
		delete [] pcorr[i];
  }
	delete [] pcorr;
}

void SEModel::set_default_params(){
  separams.expect=10;
  separams.gcback=0.38;
  separams.minpass[0]=200;
  separams.seed=-1;
  separams.psfact=0.1;
  separams.weight=0.8; 
  separams.npass=1000000;
  separams.fragment=true;
  separams.flanking=0;
  separams.undersample=1;
  separams.oversample=1;
	separams.minsize=10;
	separams.mincorr=0.4;
}

void SEModel::set_final_params(){
  separams.npseudo = separams.expect * separams.psfact;
  separams.backfreq[0]=separams.backfreq[5] = 0.25;
  separams.backfreq[1]=separams.backfreq[4] = (1 - separams.gcback)/2.0;
  separams.backfreq[2]=separams.backfreq[3] = separams.gcback/2.0;
  for(int i = 0; i < 6; i++){
    separams.pseudo[i] = separams.npseudo * separams.backfreq[i];
  }
  separams.maxlen = 3 * sites.width();
  separams.minpass[1] = 2 * separams.minpass[0];
  separams.minpass[2] = 3 * separams.minpass[0];
	separams.nruns = sites.positions_available() / separams.expect / sites.ncols() / separams.undersample * separams.oversample;
	separams.select = 5.0;
}

void SEModel::ace_initialize(){
  ran_int.set_seed(separams.seed);
  separams.seed = ran_int.seed();
  ran_int.set_range(0, RAND_MAX);
  ran_dbl.set_seed(ran_int.rnum());
  ran_dbl.set_range(0.0, 1.0);
}

void SEModel::add_possible(const int gene) {
	if (! possible[gene]) {
		possible[gene] = true;
		npossible++;
	}
}

void SEModel::remove_possible(const int gene) {
	if (possible[gene]) {
		possible[gene] = false;
		npossible--;
	}
}

void SEModel::clear_all_possible() {
	for(int g = 0; g < ngenes; g++) {
		possible[g] = false;
	}
	npossible = 0;
}

bool SEModel::is_possible(const int gene) const {
	return possible[gene];
}

int SEModel::possible_size() const {
	return npossible;
}

int SEModel::possible_positions() const {
	return sites.positions_available(possible);
}

void SEModel::clear_sites() {
	sites.clear_sites();
	select_sites.clear_sites();
}

int SEModel::size() const {
	return sites.seqs_with_sites();
}

int SEModel::sites_size() const {
	return sites.number();
}

bool SEModel::is_member(const int gene) const {
	return sites.seq_has_site(gene);
}

void SEModel::genes(int* genes) const {
	int count = 0;
	for (int g = 0; g < ngenes; g++) {
		if (is_member(g)) {
			genes[count] = g;
			count++;
		}
	}
	assert(count == size());
}

void SEModel::seed_random_site() {
	if(npossible < 1) return;
	
	sites.clear_sites();
	int chosen_possible, chosen_seq, chosen_posit;
  bool watson;
	
	chosen_seq = chosen_posit = -1;
	
	/* First choose a sequence */
	ran_int.set_range(0, possible_size() - 1);
	chosen_possible = ran_int.rnum();
	int g;
	for(g = 0; g < ngenes; g++) {
		if(is_possible(g) && chosen_possible == 0) break; 
		if(is_possible(g)) chosen_possible--;
	}
	chosen_seq = g;
	
	/* Now choose a site */
  for(int j = 0; j < 50; j++) {
		ran_int.set_range(0, seqset.len_seq(chosen_seq) - sites.width() - 1);
		double db = ran_dbl.rnum();//random (0,1)
		watson = (db > 0.5);
		chosen_posit = ran_int.rnum();
		if(watson && (chosen_posit > seqset.len_seq(chosen_seq) - sites.width() - 1)) continue;
		if((! watson) && (chosen_posit < sites.width())) continue;
		if(sites.is_open_site(chosen_seq, chosen_posit)) {
      sites.add_site(chosen_seq, chosen_posit, watson);
      break;
    }
  }
}

void SEModel::calc_matrix() {
  int d = sites.depth();
  double tot = (double) sites.number() + separams.npseudo;
  sites.calc_freq_matrix(seqset, freq_matrix);
  for(int i = 0; i < d * sites.ncols();i += d){
    score_matrix[i] = score_matrix[i + 5] = 1.0;
    for(int j = 1;j <= 4; j++){
      int x = freq_matrix[i] + freq_matrix[i + 5];
      score_matrix[i + j] = (freq_matrix[i + j] + x * separams.backfreq[j] 
				+ separams.pseudo[j])/(tot * separams.backfreq[j]);
    }
  }
}

void SEModel::single_pass(const double minprob, bool greedy) {
	double ap = (separams.weight * separams.expect * 10
							+ (1 - separams.weight) * sites.number())
							/(2.0 * sites.positions_available(possible));
  calc_matrix();
	sites.remove_all_sites();
  select_sites.remove_all_sites();
  //will only update once per pass
	
  char **ss_seq;
  ss_seq = seqset.seq_ptr();
  double Lw, Lc, Pw, Pc, F;
	int matpos;
  int considered = 0;
  int gadd = -1,jadd = -1;
  for(int g = 0; g < seqset.num_seqs(); g++){
		if (! is_possible(g)) continue;
		for(int j = 0; j <= seqset.len_seq(g) - sites.width(); j++){
			Lw = 1.0;
			matpos = 0;
      for(int k = 0; k < sites.ncols(); k++){
				assert(j + sites.column(k) >= 0);
				assert(j + sites.column(k) <= seqset.len_seq(g));
				int seq = ss_seq[g][j + sites.column(k)];
				Lw *= score_matrix[matpos + seq];
				matpos += sites.depth();
      }
      Lc = 1.0;
			matpos = sites.depth() - 1;
      for(int k = 0; k < sites.ncols(); k++){
				assert(j + sites.width() - 1 - sites.column(k) >= 0);
				assert(j + sites.width() - 1 - sites.column(k) <= seqset.len_seq(g));
				int seq = ss_seq[g][j + sites.width() - 1 - sites.column(k)];
				Lc *= score_matrix[matpos - seq];
				matpos += sites.depth();
      }
      Pw = Lw * ap/(1.0 - ap + Lw * ap);
      Pc = Lc * ap/(1.0 - ap + Lc * ap);
      F = Pw + Pc - Pw * Pc;//probability of either
			if(F > (minprob * 0.025/separams.select)) {
				select_sites.add_site(g, j, true);
			}
			//strand irrelevant for select_sites
			if(g == gadd && j < jadd + sites.width()) continue;
			if(F < minprob) continue;
			considered++;
			Pw = F * Pw / (Pw + Pc);
			Pc = F - Pw;
			if(greedy) {
				if(Pw > Pc) {
					assert(j >= 0);
					assert(j <= seqset.len_seq(g) - sites.width());
					sites.add_site(g, j, true);
					gadd = g;
					jadd = j;
				} else {
					assert(j >= 0);
					assert(j <= seqset.len_seq(g) - sites.width());
					sites.add_site(g, j, false);
					gadd = g;
					jadd = j;
				}
			} else {
				double r = ran_dbl.rnum();
				if (r > F) continue;
				else if (r < Pw) {
					assert(j >= 0);
					assert(j <= seqset.len_seq(g) - sites.width());
					sites.add_site(g, j, true);
					gadd = g;
					jadd = j;
				} else {
					assert(j >= 0);
					assert(j <= seqset.len_seq(g) - sites.width());
					sites.add_site(g, j, false);
					gadd = g;
					jadd = j;
				}
			}
    }
  }
  if(verbose)cerr<<"sampling "<<considered<<"/"<<select_sites.number()<<"\n";
}

void SEModel::single_pass_select(const double minprob){
	int i,j,k,n;
  double ap= (separams.weight * separams.expect + (1-separams.weight) * sites.number())
							/ (2.0 * sites.positions_available(possible));
  calc_matrix();
	sites.remove_all_sites();
  //will only update once per pass
	
  char **ss_seq;
  ss_seq=seqset.seq_ptr();
  double Lw, Lc, Pw, Pc, F;
  int matpos;
  int iadd=-1,jadd=-1;
  for(n=0;n<select_sites.number();n++){
    i=select_sites.chrom(n);
		if(! is_possible(i)) continue;
		j=select_sites.posit(n);
    if(i==iadd&&j<jadd+sites.width()) continue;
    if(j<0||j>seqset.len_seq(i)-sites.width()) continue;
    //could be screwed up with column sampling
    Lw=1.0;matpos=0;
    for(k=0;k<sites.ncols();k++){
      int seq=ss_seq[i][j+sites.column(k)];
      Lw*=score_matrix[matpos+seq];
      matpos+=sites.depth();
    }
    Lc=1.0;matpos=sites.depth()-1;
    for(k=0;k<sites.ncols();k++){
      int seq=ss_seq[i][j+sites.width()-1-sites.column(k)];
      Lc*=score_matrix[matpos-seq];
      matpos+=sites.depth();
    }
		Pw = Lw * ap / (1.0 - ap + Lw * ap);
    Pc = Lc * ap / (1.0 - ap + Lc * ap);
    F = (Pw + Pc - Pw * Pc);//probability of either
		if(F < minprob) continue;
		Pw = F * Pw/(Pw + Pc);
		Pc = F - Pw;
		double r=ran_dbl.rnum();
		if (r>F)
			continue;
    else if (r<Pw) {
      sites.add_site(i,j,true);
      add_possible(i);
			iadd=i;jadd=j;
    } else {
      sites.add_site(i,j,false);
			add_possible(i);
      iadd=i;jadd=j;
    }
  }
}

void SEModel::compute_seq_scores() {
	double ap = (separams.weight * separams.expect * 10
							+ (1 - separams.weight) * sites.number())
							/(2.0 * sites.positions_available(possible));
  calc_matrix();
	
  char **ss_seq;
  ss_seq = seqset.seq_ptr();
  double Lw, Lc, Pw, Pc, F, bestF;
	int matpos;
	seqranks.clear();
	for(int g = 0; g < seqset.num_seqs(); g++) {
		bestF = 0.0;
		for(int j = 0; j <= seqset.len_seq(g) - sites.width(); j++) {
			Lw = 1.0;
			matpos = 0;
			for(int k = 0; k < sites.ncols(); k++) {
				assert(j + sites.column(k) >= 0);
				assert(j + sites.column(k) <= seqset.len_seq(g));
				int seq = ss_seq[g][j + sites.column(k)];
				Lw *= score_matrix[matpos + seq];
				matpos += sites.depth();
      }
      Lc = 1.0;
			matpos = sites.depth() - 1;
      for(int k = 0; k < sites.ncols(); k++){
				assert(j + sites.width() - 1 - sites.column(k) >= 0);
				assert(j + sites.width() - 1 - sites.column(k) <= seqset.len_seq(g));
				int seq = ss_seq[g][j + sites.width() - 1 - sites.column(k)];
				Lc *= score_matrix[matpos - seq];
				matpos += sites.depth();
      }
      Pw = Lw * ap/(1.0 - ap + Lw * ap);
      Pc = Lc * ap/(1.0 - ap + Lc * ap);
      F = Pw + Pc - Pw * Pc;//probability of either
			if(F > bestF) bestF = F;
		}
		seqscores[g] = bestF;
		struct idscore ids;
		ids.id = g;
		ids.score = bestF;
		seqranks.push_back(ids);
  }
	
	sort(seqranks.begin(), seqranks.end(), isc);
}

void SEModel::compute_expr_scores() {
	calc_mean();
	expranks.clear();
	for(int g = 0; g < ngenes; g++) {
		expscores[g] = get_corr_with_mean(expr[g]);
		struct idscore ids;
		ids.id = g;
		ids.score = expscores[g];
		expranks.push_back(ids);
	}
	sort(expranks.begin(), expranks.end(), isc);
}  

bool SEModel::column_sample(){
  int *freq = new int[sites.depth()];
	
	// Compute scores for current and surrounding columns
  int max_left, max_right;
  max_left = max_right = (sites.max_width() - sites.width())/2;
  sites.columns_open(max_left, max_right);
	int cs_span = max_left + max_right + sites.width();
  vector<struct idscore> wtx(cs_span);
	int x = max_left;
  //wtx[x + c] will refer to the weight of pos c in the usual numbering
	double wt, best_wt;
	for(int i = 0; i < cs_span; i++) {
		wtx[i].id = 1000;
		wtx[i].score = -DBL_MAX;
		if(sites.column_freq(i - x, seqset, freq)){
      wt = 0.0;
      for(int j = 0;j < sites.depth(); j++){
				wt += gammaln(freq[j] + separams.pseudo[j]);
				wt -= (double)freq[j] * log(separams.backfreq[j]);
      }
			wtx[i].id = i;
			wtx[i].score = wt;
			if(wt > best_wt) best_wt = wt;
    }
  }
	
	// Penalize outermost columns for length
  double scale = 0.0;
  if(best_wt > 100.0) scale = best_wt - 100.0; //keep exp from overflowing
	for(int i = 0; i < cs_span; i++){
		wtx[i].score -= scale;
		wtx[i].score = exp(wtx[i].score);
		int newwidth = sites.width();
		if(i < x)
			newwidth += (x - i);
		else if(i > (x + sites.width() - 1))
			newwidth += (i - x - sites.width() + 1);
		wtx[i].score /= bico(newwidth - 2, sites.ncols() - 2);
	}
	
	// Sort columns by their score in wtx
	sort(wtx.begin(), wtx.end(), isc);
	
	// Keep the top <ncol> columns
	int ncols = sites.ncols();
	int nseen = 0;
	int shift = 0;
	for(int i = 0; i < cs_span; ++i) {
		if(wtx[i].id - x > cs_span - 1)
			continue;
		if(nseen < ncols) {
			if(sites.has_col(wtx[i].id - x)) {        // column is ranked in top, and is already in motif
				nseen++;
			} else {
				sites.add_col(wtx[i].id - x);           // column is ranked in top, and is not in motif
				if(wtx[i].id - x < 0)                   // if column was to the left of the current columns, adjust the remaining columns
					for(int j = i + 1; j < cs_span; ++j)
						wtx[j].id -= wtx[i].id - x;
				nseen++;
			}
		} else {
			if(sites.has_col(wtx[i].id - x)) {        // column is not ranked in top, and is in motif
				if(wtx[i].id - x == 0)                  // we will remove the first column, adjust remaining columns
					for(int j = i + 1; j < cs_span; ++j)
						wtx[j].id -= sites.column(1);
				sites.remove_col(wtx[i].id - x);
			}
		}
	}
	assert(sites.ncols() == ncols);               // number of columns should not change
}

double SEModel::map_score() {
	int i,j,k;
  double ms = 0.0;
  double map_N = sites.positions_available(possible);  
  double w = separams.weight/(1.0-separams.weight);
  double map_alpha = (double) separams.expect * w;
  double map_beta = map_N * w - map_alpha;
  double map_success = (double)sites.number();
  ms += ( gammaln(map_success+map_alpha)+gammaln(map_N-map_success+map_beta) );
  ms -= ( gammaln(map_alpha)+gammaln(map_N + map_beta) );
	
  sites.calc_freq_matrix(seqset,freq_matrix);
  double sc[6]={0.0,0.0,0.0,0.0,0.0,0.0};
  int d=sites.depth();
  for(k=0;k!=d*sites.ncols();k+=d) {
    int x=freq_matrix[k]+freq_matrix[k+5];
    for(j=1;j<=4;j++){
      ms += gammaln((double)freq_matrix[k+j]+x*separams.backfreq[j]+separams.pseudo[j]);
      sc[j]+=freq_matrix[k+j]+x*separams.backfreq[j];
    }
  }
  ms-=sites.ncols()*gammaln((double)sites.number()+separams.npseudo);
  for (k=1;k<=4;k++) {
    ms -= sc[k]*log(separams.backfreq[k]);
  }
  /*This factor arises from a modification of the model of Liu, et al in which the background frequencies of DNA bases are taken to be constant for the organism under consideration*/
	double  vg;
  ms -= lnbico(sites.width()-2, sites.ncols()-2);
  for(vg=0.0,k=1;k<=4;k++) {
    vg += gammaln(separams.pseudo[k]);
  }
  vg-=gammaln((double)(separams.npseudo));
  ms-=((double)sites.ncols()*vg);
	return ms;
}

double SEModel::spec_score() {
	int x, s1, s2;
	x = s1 = s2 = 0;
	compute_seq_scores();
	compute_expr_scores();
	for(int g = 0; g < ngenes; g++) {
		if(seqscores[g] > sites.get_seq_cutoff()) s1++;
		if(expscores[g] > sites.get_expr_cutoff()) s2++;
		if(seqscores[g] > sites.get_seq_cutoff() && expscores[g] > sites.get_expr_cutoff()) x++;
	}
	
	double spec = prob_overlap(s1, s2, x, ngenes);
	if(spec > 0.99) return 0;
	else return -log10(spec);
}

void SEModel::optimize_sites(){
  int poss_sites = sites.positions_available(possible)+1;
  int *pos = new int[poss_sites];
  int *chr = new int[poss_sites];
  bool *str = new bool[poss_sites];
  Heap hp(poss_sites,3);
  double ap = (separams.weight * separams.expect
	            + (1 - separams.weight) * sites.number())
							/ (2.0 * sites.positions_available(possible));
	
  calc_matrix();
	calc_mean();
  char **ss_seq;
  ss_seq = seqset.seq_ptr();
  double Lw, Lc, Pw, Pc;
  int matpos;
  int h = 1;
  for(int i = 0; i < seqset.num_seqs(); i++) {
    if(! is_possible(i)) continue;
		for(int j = 0; j < seqset.len_seq(i) - sites.width() + 1; j++) {
      Lw = 1.0;
			matpos = 0;
      for(int k = 0; k < sites.ncols(); k++){
				int seq = ss_seq[i][j + sites.column(k)];
				Lw *= score_matrix[matpos + seq];
				matpos += sites.depth();
      }
      Lc = 1.0;
			matpos = sites.depth()-1;
      for(int k = 0; k < sites.ncols(); k++){
				int seq = ss_seq[i][j + sites.width() - 1 - sites.column(k)];
				Lc *= score_matrix[matpos - seq];
				matpos += sites.depth();
      }
      Pw = Lw * ap / (1.0 - ap + Lw * ap);
      Pc = Lc * ap / (1.0 - ap + Lc * ap);
      //F=(Pw+Pc-Pw*Pc);//probability of either
      if(Pw > sites.get_seq_cutoff() * 0.025 || Pc > sites.get_seq_cutoff() * 0.025) { 
				chr[h] = i;
				pos[h] = j;
				if(Pw >= Pc) {
					str[h] = true;
					hp.insert(h, -Pw);//lower is better
				} else {
					str[h] = false;
					hp.insert(h, -Pc);
				}
				h++;
      }
    }
  }
  
  sites.remove_all_sites();
  Sites best = sites;
  double ms, ms_best = -DBL_MAX;
	for(;;){
    h = hp.del_min();
		if(h == -1) break;
    if(sites.is_open_site(chr[h], pos[h])) {
      sites.add_site(chr[h], pos[h], str[h]);
      if(sites.number() <= 3) continue;
      if(hp.value(h) < -0.7) continue;
      ms = map_score();
      if(ms > ms_best){
				ms_best = ms;
				best = sites;
      }
    }
  }
	sites = best;
  delete [] pos;
  delete [] chr;
  delete [] str;
}

void SEModel::orient_motif(){
  double *info = new double[6];
  double *freq = new double[6];
  for(int i = 0; i < 6; i++) info[i] = 0.0;
  int d = sites.depth();
	
  double tot = (double)sites.number() + separams.npseudo;
  sites.calc_freq_matrix(seqset, freq_matrix);
  for(int i = 0; i < d * sites.ncols(); i += d){
    double ii = 0.0;
    for(int j = 1; j <= 4; j++) {
      int x = freq_matrix[i] + freq_matrix[i+5];
      freq[j] = (freq_matrix[i + j] + x * separams.backfreq[j] + separams.pseudo[j]) / tot;
      ii += freq[j]*log(freq[j]);
    }
    ii = 2 + ii;
    for(int j = 1; j <= 4; j++) info[j] += freq[j] * ii;
  }
	
  double flip = 1.5 * info[3] + 1.0 * info[1] - 1.0 * info[4] - 1.5 * info[2];
  //for(i=1;i<5;i++) cerr<<info[i]<<'\t';
  //cerr<<flip<<'\n';
  if(flip < 0.0) sites.flip_sites();
  delete [] info;
  delete [] freq;
}

void SEModel::orient_print_motif(){
  double *info = new double[6];
  double *freq = new double[6];
  for(int i = 0; i < 6; i++) info[i] = 0.0;
  int d = print_sites.depth();
	
  double tot = (double) print_sites.number() + separams.npseudo;
  print_sites.calc_freq_matrix(seqset, freq_matrix);
  for(int i = 0; i < d * print_sites.ncols(); i += d){
    double ii = 0.0;
    for(int j = 1; j <= 4; j++) {
      int x = freq_matrix[i] + freq_matrix[i+5];
      freq[j] = (freq_matrix[i + j] + x * separams.backfreq[j] + separams.pseudo[j]) / tot;
      ii += freq[j]*log(freq[j]);
    }
    ii = 2 + ii;
    for(int j = 1; j <= 4; j++) info[j] += freq[j] * ii;
  }
	
  double flip = 1.5 * info[3] + 1.0 * info[1] - 1.0 * info[4] - 1.5 * info[2];
  //for(i=1;i<5;i++) cerr<<info[i]<<'\t';
  //cerr<<flip<<'\n';
  if(flip < 0.0) print_sites.flip_sites();
  delete [] info;
  delete [] freq;
}

string SEModel::consensus() const {
	map<char,char> nt;
  nt[0]=nt[5]='N';
  nt[1]='A';nt[2]='C';nt[3]='G';nt[4]='T';
	char** ss_seq=seqset.seq_ptr();
	
	int numsites = sites.number();
	if(numsites < 1) return "";
	
	// cerr << "Computing consensus with " << numsites << " sites" << endl;
	vector<string> hits(numsites);
	for(int i = 0; i < numsites; i++){
		int c = sites.chrom(i);
    int p = sites.posit(i);
    bool s = sites.strand(i);
    for(int j = 0; j < sites.width(); j++){
      if(s) {
				if(p + j >= 0 && p+j < seqset.len_seq(c))
					 hits[i] += nt[ss_seq[c][p+j]];
				else hits[i] += ' ';
      } else {
				if(p + sites.width() - 1 - j >= 0 && p+sites.width() - 1 - j < seqset.len_seq(c))
					hits[i] += nt[sites.depth()-1-ss_seq[c][p+sites.width()-1-j]];
				else hits[i] += ' ';
      }
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

void SEModel::output_params(ostream &fout){
  fout<<" expect =      \t"<<separams.expect<<'\n';
  fout<<" gcback =      \t"<<separams.gcback<<'\n';
  fout<<" minpass =     \t"<<separams.minpass[0]<<'\n';
  fout<<" seed =        \t"<<separams.seed<<'\n';
  fout<<" numcols =     \t"<<sites.ncols()<<'\n';
  fout<<" undersample = \t"<<separams.undersample<<'\n';
  fout<<" oversample = \t"<<separams.oversample<<'\n';
}

void SEModel::modify_params(int argc, char *argv[]){
  GetArg2(argc,argv,"-expect",separams.expect);
  GetArg2(argc,argv,"-gcback",separams.gcback);
  GetArg2(argc,argv,"-minpass",separams.minpass[0]);
  GetArg2(argc,argv,"-seed",separams.seed);
  GetArg2(argc,argv,"-undersample",separams.undersample);
  GetArg2(argc,argv,"-oversample",separams.oversample);
}

void SEModel::calc_mean() {
	int i, j;
	for (i = 0; i < npoints; i++) {
		mean[i] = 0;
		for (j = 0; j < ngenes; j++)
			if(is_member(j))
				mean[i] += expr[j][i];
		mean[i] /= size();
	}
}

void SEModel::calc_stdev() {
	calc_mean();
	for(int i = 0; i < npoints; i++) {
		stdev[i] = 0;
		for(int j = 0; j < ngenes; j++)
			if(is_member(j))
				stdev[i] += (expr[j][i] - mean[i]) * (expr[j][i] - mean[i]);
		stdev[i] /= size() - 1;
		stdev[i] = sqrt(stdev[i]);
	} 
}

float SEModel::get_avg_pcorr() {
	float curr_avg_pcorr = 0.0;
	for(int g1 = 0; g1 < ngenes; g1++)
		for(int g2 = g1 + 1; g2 < ngenes; g2++)
			if(is_member(g1) && is_member(g2))
				curr_avg_pcorr += get_pcorr(g1, g2);
	curr_avg_pcorr = (2 * curr_avg_pcorr) / (size() * (size() - 1));
	return curr_avg_pcorr;
}

float SEModel::get_pcorr(const int g1, const int g2) {
	if(pcorr[g1][g2] == -2)
		pcorr[g1][g2] = pcorr[g2][g1] = corr(expr[g1], expr[g2], npoints);
	return pcorr[g1][g2];
}

float SEModel::get_corr_with_mean(const float* pattern) const {
	return corr(mean, pattern, npoints);
}

void SEModel::expand_search_around_mean(const double cutoff) {
	calc_mean();
	clear_all_possible();
	for(int g = 0; g < ngenes; g++)
		if(get_corr_with_mean(expr[g]) >= cutoff) add_possible(g);
}

void SEModel::expand_search_min_pcorr(const double cutoff) {
 	// First add all genes to list of candidates
 	list<int> candidates;
 	for(int g = 0; g < ngenes; g++)
 		candidates.push_back(g);
 	
 	// Now run through list of genes with sites, and only keep genes within minimum corr_cutoff
 	for(int g = 0; g < ngenes; g++) {
 		if(! is_member(g)) continue;
 		list<int> survivors;
 		for(list<int>::iterator iter = candidates.begin(); iter != candidates.end(); iter++)
 			if(get_pcorr(g, *iter) > cutoff) survivors.push_back(*iter);
 		candidates.assign(survivors.begin(), survivors.end());
 	}
 	
 	// Start with clean slate
 	clear_all_possible();
 	
 	// Add genes which already have sites to the search space
 	for(int g = 0; g < ngenes; g++)
 		if(is_member(g)) add_possible(g);
	for(list<int>::iterator iter = candidates.begin(); iter != candidates.end(); iter++)
		add_possible(*iter);
}

void SEModel::expand_search_avg_pcorr() {
	if(size() > 10) {
		sites.set_expr_cutoff(max(0.9 * get_avg_pcorr(), 0.4));
	}

	if(size() > 0) {
		float avg_pcorr;
		clear_all_possible();
		for(int g1 = 0; g1 < ngenes; g1++) {
			if(is_member(g1)) continue;
			avg_pcorr = 0;
			for(int g2 = 0; g2 < ngenes; g2++) {
				if(! is_member(g2)) continue;
				avg_pcorr += get_pcorr(g1, g2);
			}
			avg_pcorr /= size();
			if(avg_pcorr > sites.get_expr_cutoff()) add_possible(g1);
		}
		for(int g1 = 0; g1 < ngenes; g1++)
			if(is_member(g1)) add_possible(g1);
	}
}

void SEModel::search_for_motif(const int worker, const int iter) {
	sites.clear_sites();
	select_sites.clear_sites();
	sites.set_iter(iter);
	sites.set_expr_cutoff(0.85);
	sites.set_map(0);
	double sc_best_i, sp;
  int i_worse = 0;
	sc_best_i = map_cutoff;
	i_worse = 0;
	int phase = 0;
	int oldphase = 0;
	
	Sites best_sites = sites;
	for(int g = 0; g < ngenes; g++)
		add_possible(g);
	seed_random_site();
	if(size() < 1) {
		cerr << "\t\t\tSeeding failed -- restarting..." << endl;
		return;
	}
	
	expand_search_around_mean(sites.get_expr_cutoff());
	while(possible_size() < separams.minsize && sites.get_expr_cutoff() > 0.4) {
		print_status(cerr, 0, oldphase);
		sites.set_expr_cutoff(sites.get_expr_cutoff() - 0.05);
		expand_search_around_mean(sites.get_expr_cutoff());
	}
	if(possible_size() < 2) {
		cerr << "\t\t\tBad search start -- no genes within " << separams.mincorr << endl;
		return;
	}
	print_status(cerr, 0, oldphase);
	compute_seq_scores();
	
	double ap = (double) separams.expect/(2.0 * sites.positions_available(possible));
	sites.set_seq_cutoff(ap * 5.0);
	
	for(int i = 1; i <= separams.npass; i++){
		if(oldphase < phase) {
			print_status(cerr, i, oldphase);
			if(size() > 5) {
				compute_seq_scores();
				compute_expr_scores();
				set_seq_cutoffs(i);
				set_expr_cutoffs(i);
				expand_search_around_mean(sites.get_expr_cutoff());
			}
			oldphase = phase;
		}
		if(phase == 3) {
			if(size() < separams.minsize/2) {
				cerr << "\t\t\tReached phase " << phase << " with less than " << separams.minsize/2 << " sequences with sites. Restarting..." << endl;
				break;
			}
			double sc1 = map_score();
			sites.set_map(0.0);
			for(int z = 0; sites.get_map() < sc1 && z < 5; z++){
				column_sample();
				single_pass(sites.get_seq_cutoff(), true);
				sites.set_map(map_score());
				print_status(cerr, i, phase);
			}
			if(sites.get_map() < sc1) {
				sites = best_sites;
				sites.set_map(sc1);
			}
			print_status(cerr, i, phase);
			if(size() < separams.minsize) {
				cerr << "\t\t\tCompleted phase " << phase << " with less than " << separams.minsize << " sequences with sites. Restarting..." << endl;
				break;
			}
			sp = spec_score();
			sites.set_spec(sp);
			cerr << "\t\t\t\tSpecificity score was " << sp << endl;
			if(archive.check_motif(sites)) {
				char tmpfilename[30], motfilename[30];
				sprintf(tmpfilename, "%d.%d.mot.tmp", worker, iter);
				sprintf(motfilename, "%d.%d.mot", worker, iter);
				ofstream motout(tmpfilename);
				sites.write(seqset, motout);
				motout.close();
				rename(tmpfilename, motfilename);
				cerr << "\t\t\tCompleted phase " << phase << "! Restarting..." << endl;
			} else {
				cerr <<"\t\t\tToo similar! Restarting..." << endl;
			}
			break;
		}
		if(i_worse == 0) {
			single_pass(sites.get_seq_cutoff());
		} else {
			single_pass_select(sites.get_seq_cutoff());
		}
		if(sites_size() == 0) {
			if(sc_best_i == map_cutoff) {
				cerr << "\t\t\tNo sites and best score matched cutoff! Restarting..." << endl;
				break;
			}
			//if(best_sites.number()<4) break;
			sites = best_sites;
			select_sites = best_sites;
			phase++;
			i_worse = 0;
			continue;
		}
		if(i<=3) continue;
		if(phase < 3 && size() > separams.minsize)
			column_sample();
		sites.set_map(map_score());
		if(sites.get_map() - sc_best_i > 1e-3) {
			i_worse=0;
			if(size() > separams.minsize * 2) {
				if(! archive.check_motif(sites)) {
					print_status(cerr, i, phase);
					cerr <<"\t\t\tToo similar! Restarting..." << endl;
					break;
				}
			}
			sc_best_i = sites.get_map();
			best_sites = sites;
		}
		else i_worse++;
		if(i_worse > separams.minpass[phase]){
			if(sc_best_i == map_cutoff) {
				print_status(cerr, i, phase);
				cerr << "\t\t\ti_worse is greater than cutoff and best score at cutoff! Restarting..." << endl;
				break;
			}
			if(best_sites.number() < 2) {
				print_status(cerr, i, phase);
				cerr << "\t\t\ti_worse is greater than cutoff and only 1 site! Restarting..." << endl;
				break;
			}
			sites = best_sites;
			select_sites = best_sites;
			phase++;
			i_worse = 0;
		}
		
		if(i % 50 == 0 && size() > 5) {
			compute_seq_scores();
			compute_expr_scores();
			set_seq_cutoffs(i);
			set_expr_cutoffs(i);
			expand_search_around_mean(sites.get_expr_cutoff());
		}
	
		if(i % 50 == 0) {
			print_status(cerr, i, phase);
		}
	}
}

bool SEModel::consider_motif(const char* filename) {
	ifstream motin(filename);
	sites.clear_sites();
	sites.read(motin);
	bool ret = archive.consider_motif(sites);
	motin.close();
	return ret;
}

void SEModel::print_status(ostream& out, const int i, const int phase) {
	if(size() > 0) {
		out << "\t\t\t" << setw(5) << i;
		out << setw(3) << phase;
		int prec = cerr.precision(2);
		out << setw(5) << setprecision(2) << sites.get_expr_cutoff();
		out << setw(10) << setprecision(2) << sites.get_seq_cutoff();
		cerr.precision(prec);
		out << setw(5) << sites_size();
		out << setw(5) << size();
		out << setw(5) << possible_size();
		out << setw(40) << consensus();
		out << setw(10) << sites.get_map();
		out << endl;
	} else {
		out << "\t\t\tNo sites!" << endl;
	}
}

void SEModel::full_output(ostream &fout){
  Sites* s;
	for(int j = 0; j < archive.nmots(); j++){
    s = archive.return_best(j);
    if(s->get_map() > 0.0){
			fout << "Motif " << j + 1 << endl;
			s->write(seqset, fout);
		}
    else break;
  }
}

void SEModel::full_output(char *name){
  ofstream fout(name);
  full_output(fout);
}

void SEModel::set_seq_cutoffs(const int i) {
	int seqn, expn = 0, isect, next = 0;
	double best_po = 1;
	double best_sitecut = 0;
	double po, c;
	for(int i = 0; i < ngenes; i++)
		if(expscores[i] > sites.get_expr_cutoff()) expn++;
	for(c = seqranks[0].score; c > 0; c = seqranks[next].score) {
		seqn = isect = 0;
		for(int i = 0; i < ngenes; i++) {
			if(seqranks[i].score >= c) {
				seqn++;
				if(expscores[seqranks[i].id] >= sites.get_expr_cutoff())
					isect++;
			} else {
				next = i;
				break;
			}
		}
		po = prob_overlap(expn, seqn, isect, ngenes);
		// cerr << "\t\t\t\t\t\t" << c << ":\t" << isect << "/(" << seqn << "," << expn << ")/" << ngenes << "\t" << po << endl;
		if(po <= best_po) {
			best_po = po;
			best_sitecut = c;
		}
	}
	// cerr << "\t\t\t\t\tSetting sequence cutoff to " << best_sitecut << endl;
	sites.set_seq_cutoff(best_sitecut);
}

void SEModel::set_expr_cutoffs(const int i) {
	int seqn = 0, expn, isect;
	double best_po = 1;
	double best_exprcut = 0;
	double po, c;
	for(int i = 0; i < ngenes; i++)
		if(seqscores[i] >= sites.get_seq_cutoff()) seqn++;
	for(c = -1.0; c <= 1; c += 0.01) {
		expn = isect = 0;
		for(int i = 0; i < ngenes; i++) {
			if(expranks[i].score >= c) {
				expn++;
				if(seqscores[expranks[i].id] >= sites.get_seq_cutoff())
					isect++;
			} else {
				break;
			}
		}
		po = prob_overlap(seqn, expn, isect, ngenes);
		// cerr << "\t\t\t\t\t\t" << c << ":\t\t" << isect << "/(" << expn << "," << seqn << ")/" << ngenes << "\t\t" << po << endl; 
		if(po <= best_po) {
			best_po = po;
			best_exprcut = c;
		}
	}
	// cerr << "\t\t\t\t\tSetting expression cutoff to " << best_exprcut << endl;
	sites.set_expr_cutoff(best_exprcut);
}

void SEModel::print_possible(ostream& out) {
	for(int g = 0; g < ngenes; g++)
		if(is_possible(g)) out << nameset[g] << endl;
}

