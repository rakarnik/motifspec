#include "semodel.h"

SEModel::SEModel(const vector<string>& seqs, float** exprtab, const int numexpr, const vector<string>& names, const int nc, const int bf, const double map_cut, const double sim_cut):
nameset(names),
ngenes(names.size()),
seqset(seqs),
motif(seqset, nc),
archive(seqset, map_cut, sim_cut),
seqscores(ngenes),
bestpos(ngenes),
seqranks(ngenes),
expscores(ngenes),
expranks(ngenes)
{
	npossible = 0;
	possible = new bool[ngenes];
	clear_all_possible();
	verbose = false;
	
	set_default_params();
  max_motifs = bf;
  map_cutoff = map_cut;
  sim_cutoff = sim_cut;
  freq_matrix = new int[4 * motif.ncols()];
  score_matrix = new double[4 * motif.ncols()];
	
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
  separams.expect = 10;
  separams.minpass = 50;
  separams.seed = -1;
  separams.psfact = 0.1;
  separams.weight = 0.5; 
  separams.npass = 1000000;
  separams.fragment = true;
  separams.flanking = 0;
  separams.undersample = 1;
  separams.oversample = 1;
	separams.minsize = 5;
	separams.mincorr = 0.4;
}

void SEModel::set_final_params(){
  separams.npseudo = separams.expect * separams.psfact;
  separams.backfreq[0] = separams.backfreq[3] = (1 - seqset.gcgenome())/2.0;
  separams.backfreq[1] = separams.backfreq[2] = seqset.gcgenome()/2.0;
  for(int i = 0; i < 4; i++) {
		separams.pseudo[i] = separams.npseudo * separams.backfreq[i];
  }
  separams.maxlen = 3 * motif.get_width();
  separams.nruns = motif.positions_available() / separams.expect / motif.ncols() / separams.undersample * separams.oversample;
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

int SEModel::total_positions() const {
	return motif.positions_available();
}

int SEModel::possible_positions() const {
	return motif.positions_available(possible);
}

void SEModel::clear_sites() {
	motif.clear_sites();
}

int SEModel::size() const {
	return motif.seqs_with_sites();
}

int SEModel::motif_size() const {
	return motif.number();
}

bool SEModel::is_member(const int gene) const {
	return motif.seq_has_site(gene);
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
	
	motif.clear_sites();
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
	int width = motif.get_width();
  for(int j = 0; j < 50; j++) {
		ran_int.set_range(0, seqset.len_seq(chosen_seq) - width - 1);
		double db = ran_dbl.rnum();//random (0,1)
		watson = (db > 0.5);
		chosen_posit = ran_int.rnum();
		if(watson && (chosen_posit > seqset.len_seq(chosen_seq) - width - 1)) continue;
		if((! watson) && (chosen_posit < width)) continue;
		if(motif.is_open_site(chosen_seq, chosen_posit)) {
      cerr << "\t\t\tSeeding with (" << chosen_seq << "," << chosen_posit << "," << watson << ")\n";
			motif.add_site(chosen_seq, chosen_posit, watson);
      break;
    }
  }
}

void SEModel::calc_matrix() {
  motif.calc_score_matrix(score_matrix, separams.pseudo);
}

double SEModel::score_site(const int c, const int p, const bool s) {
	const vector<vector<int> >& ss_seq = seqset.seq();
	const vector<int>& bgpos = seqset.get_bgpos();
	const vector<float>& wbgscores = seqset.get_wbgscores();
	const vector<float>& cbgscores = seqset.get_cbgscores();
	double L = 0.0;
	int width = motif.get_width();
	int matpos, seq;
	vector<int>::iterator col_iter = motif.first_column();
	vector<int>::iterator last_col = motif.last_column();
	if(s) {
		matpos = 0;
		for(; col_iter != last_col; ++col_iter) {
			assert(p + *col_iter >= 0);
			assert(p + *col_iter < seqset.len_seq(c));
			L += score_matrix[matpos + ss_seq[c][p + *col_iter]];
			L -= wbgscores[bgpos[c] + p + *col_iter];
			matpos += 4;
		}
	} else {
		matpos = 0;
		for(; col_iter != last_col; ++col_iter) {
			assert(p + width - 1 - *col_iter >= 0);
			assert(p + width - 1 - *col_iter < seqset.len_seq(c));
			seq = ss_seq[c][p + width - 1 - *col_iter];
			L += score_matrix[matpos + 3 - ss_seq[c][p + width - 1 - *col_iter]];
			L -= cbgscores[bgpos[c] + p + width - 1 - *col_iter];
			matpos += 4;
		}
	}
	return exp(L);
}

void SEModel::single_pass(const double minprob, bool greedy) {
	double ap = (separams.weight * separams.expect * 10
							+ (1 - separams.weight) * motif.number())
							/(2.0 * motif.positions_available(possible));
  calc_matrix();
	motif.remove_all_sites();
  //will only update once per pass
	
  double Lw, Lc, Pw, Pc, F;
  int considered = 0;
  int gadd = -1, jadd = -1;
	int width = motif.get_width();
	for(int g = 0; g < seqset.num_seqs(); g++){
		if (! is_possible(g)) continue;
		for(int j = 0; j < seqset.len_seq(g) - width; j++){
			Lw = score_site(g, j, 1);
      Lc = score_site(g, j, 0);
      Pw = Lw * ap/(1.0 - ap + Lw * ap);
      Pc = Lc * ap/(1.0 - ap + Lc * ap);
      F = Pw + Pc - Pw * Pc;//probability of either
			if(g == gadd && j < jadd + width) continue;
			if(F < minprob) continue;
			// cerr << "\t\t\t\tGene:" << g << " Position:"<< j <<  " Score:" << F << "\n";
			considered++;
			Pw = F * Pw / (Pw + Pc);
			Pc = F - Pw;
			if(greedy) {                   // Always add if above minprob
				if(Pw > Pc) {
					assert(j >= 0);
					assert(j <= seqset.len_seq(g) - width);
					motif.add_site(g, j, true);
					gadd = g;
					jadd = j;
				} else {
					assert(j >= 0);
					assert(j <= seqset.len_seq(g) - width);
					motif.add_site(g, j, false);
					gadd = g;
					jadd = j;
				}
			} else {                       // Add with probability F
				double r = ran_dbl.rnum();
				if(r > F) continue;
				if (Pw > Pc) {
					assert(j >= 0);
					assert(j <= seqset.len_seq(g) - width);
					motif.add_site(g, j, true);
					gadd = g;
					jadd = j;
				} else {
					assert(j >= 0);
					assert(j <= seqset.len_seq(g) - width);
					motif.add_site(g, j, false);
					gadd = g;
					jadd = j;
				}
			}
    }
  }
}

void SEModel::compute_seq_scores() {
	double ap = (separams.weight * separams.expect * 10
							+ (1 - separams.weight) * motif.number())
							/(2.0 * motif.positions_available(possible));
  calc_matrix();
  double Lw, Lc, Pw, Pc, F, bestF;
	seqranks.clear();
	int width = motif.get_width();
	for(int g = 0; g < seqset.num_seqs(); g++) {
		bestF = 0.0;
		bestpos[g] = -1;
		for(int j = 0; j < seqset.len_seq(g) - width; j++) {
			Lw = score_site(g, j, 1);
			Lc = score_site(g, j, 0);
      Pw = Lw * ap/(1.0 - ap + Lw * ap);
      Pc = Lc * ap/(1.0 - ap + Lc * ap);
      F = Pw + Pc - Pw * Pc;
			if(F > bestF) {
				bestF = F;
				bestpos[g] = j;
			}
		}
		seqscores[g] = bestF;
		struct idscore ids;
		ids.id = g;
		ids.score = bestF;
		seqranks.push_back(ids);
  }
	sort(seqranks.begin(), seqranks.end(), isc);
}

void SEModel::compute_seq_scores_minimal() {
	double ap = (separams.weight * separams.expect * 10
							+ (1 - separams.weight) * motif.number())
							/(2.0 * motif.positions_available(possible));
  calc_matrix();
	int width = motif.get_width();
  double Lw, Lc, Pw, Pc, F;
	seqranks.clear();
	for(int g = 0; g < seqset.num_seqs(); g++) {
		// Some best positions might have been invalidated by column sampling
		// We mark these as invalid and don't score them
		if(bestpos[g] < 0 || bestpos[g] + width > seqset.len_seq(g)) {
			bestpos[g] = -1;
			F = 0;
		} else {
			Lw = score_site(g, bestpos[g], 1);
			Lc = score_site(g, bestpos[g], 0);
			Pw = Lw * ap/(1.0 - ap + Lw * ap);
			Pc = Lc * ap/(1.0 - ap + Lc * ap);
			F = Pw + Pc - Pw * Pc;
		}
		seqscores[g] = F;
		struct idscore ids;
		ids.id = g;
		ids.score = F;
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
	bool changed = false;
	int *freq = new int[4];
	int width = motif.get_width();
	// Compute scores for current and surrounding columns
  int max_left, max_right;
	max_left = max_right = (motif.get_max_width() - width)/2;
  motif.columns_open(max_left, max_right);
	int cs_span = max_left + max_right + width;
  vector<struct idscore> wtx(cs_span);
	int x = max_left;
  //wtx[x + c] will refer to the weight of pos c in the usual numbering
	double wt;
	double best_wt = -DBL_MAX;
	for(int i = 0; i < cs_span; i++) {
		wtx[i].id = 1000;
		wtx[i].score = -DBL_MAX;
		if(motif.column_freq(i - x, freq)){
      wt = 0.0;
      for(int j = 0;j < 4; j++){
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
		int newwidth = width;
		if(i < x)
			newwidth += (x - i);
		else if(i > (x + width - 1))
			newwidth += (i - x - width + 1);
		wtx[i].score /= bico(newwidth - 2, motif.ncols() - 2);
	}
	
	// Sort columns by their score in wtx
	sort(wtx.begin(), wtx.end(), isc);
	
	// Keep the top <ncol> columns
	int ncols = motif.ncols();
	int nseen = 0;
	for(int i = 0; i < cs_span; ++i) {
		if(wtx[i].id - x > cs_span - 1)
			continue;
		if(nseen < ncols) {
			if(motif.has_col(wtx[i].id - x)) {        // column is ranked in top, and is already in motif
				nseen++;
			} else {
				motif.add_col(wtx[i].id - x);           // column is ranked in top, and is not in motif
				if(wtx[i].id - x < 0)                   // if column was to the left of the current columns, adjust the remaining columns
					for(int j = i + 1; j < cs_span; ++j)
						wtx[j].id -= wtx[i].id - x;
				changed = true;
				nseen++;
			}
		} else {
			if(motif.has_col(wtx[i].id - x)) {        // column is not ranked in top, and is in motif
				if(wtx[i].id - x == 0)                  // we will remove the first column, adjust remaining columns
					for(int j = i + 1; j < cs_span; ++j)
						wtx[j].id -= motif.column(1);
				motif.remove_col(wtx[i].id - x);
				changed = true;
			}
		}
	}
	
	if(motif.ncols() != ncols) {                 // number of columns should not change
		cerr << "\t\t\t\t\tERROR: column sampling started with " << ncols << ", ended with " << motif.ncols() << endl;
		abort();
	}
	
	return changed;
}

double SEModel::matrix_score() {
	double ms = 0.0;
	motif.calc_freq_matrix(freq_matrix);
  int nc = motif.ncols();
	int w = motif.get_width();
	double sc[] = {0.0,0.0,0.0,0.0};
  for(int i = 0; i < 4 * nc; i += 4) {
    for(int j = 0; j < 4; j++) {
      ms += gammaln((double) freq_matrix[i + j] + separams.pseudo[j]);
      sc[j] += freq_matrix[i + j];
    }
  }
	ms -= nc * gammaln((double) motif.number() + separams.npseudo);
	for (int j = 0; j < 4; j++)
		ms -= sc[j] * log(separams.backfreq[j]);
	/* 
		This factor arises from a modification of the model of Liu, et al 
		in which the background frequencies of DNA bases are taken to be
		constant for the organism under consideration
	*/
	double vg = 0.0;
  ms -= lnbico(w - 2, nc - 2);
  for(int j = 0; j < 4; j++)
    vg += gammaln(separams.pseudo[j]);
  vg -= gammaln((double) (separams.npseudo));
  ms -= ((double) nc * vg);
	return ms;
}

double SEModel::entropy_score() {
	double es = 0.0;
	motif.calc_freq_matrix(freq_matrix);
	int nc = motif.ncols();
	double f;
	for(int i = 0; i < 4 * nc; i += 4) {
    for(int j = 1; j <= 4; j++) {
			f = (double) freq_matrix[i + j] + separams.pseudo[j];
			es += f * log(1/f);
    }
  }
	// Normalize for number of columns
	es /= nc;
	return es;
}

double SEModel::map_score() {
  double ms = 0.0;
  double map_N = motif.positions_available(possible);  
  double w = separams.weight/(1.0-separams.weight);
  double map_alpha = (double) separams.expect * w;
  double map_beta = map_N * w - map_alpha;
  double map_success = (double)motif.number();
  ms += ( gammaln(map_success+map_alpha)+gammaln(map_N-map_success+map_beta) );
  ms -= ( gammaln(map_alpha)+gammaln(map_N + map_beta) );
	ms += matrix_score();
	return ms;
}

double SEModel::spec_score() {
	int isect, seqn, expn;
	isect = seqn = expn = 0;
	compute_seq_scores_minimal();
	if(seqranks[0].score <= 0.85) compute_seq_scores();
	compute_expr_scores();
	for(int g = 0; g < ngenes; g++) {
		if(seqscores[g] >= motif.get_seq_cutoff()) seqn++;
		if(expscores[g] >= motif.get_expr_cutoff()) expn++;
		if(seqscores[g] >= motif.get_seq_cutoff() && expscores[g] >= motif.get_expr_cutoff()) isect++;
	}
	
	double spec = prob_overlap(expn, seqn, isect, ngenes);
	spec = (spec > 0.99)? 0 : -log10(spec);
	return spec;
}

void SEModel::output_params(ostream &fout){
  fout<<" expect =      \t"<<separams.expect<<'\n';
  fout<<" minpass =     \t"<<separams.minpass<<'\n';
  fout<<" seed =        \t"<<separams.seed<<'\n';
  fout<<" numcols =     \t"<<motif.ncols()<<'\n';
  fout<<" undersample = \t"<<separams.undersample<<'\n';
  fout<<" oversample = \t"<<separams.oversample<<'\n';
}

void SEModel::modify_params(int argc, char *argv[]){
  GetArg2(argc, argv, "-expect", separams.expect);
  GetArg2(argc, argv, "-minpass", separams.minpass);
  GetArg2(argc, argv, "-seed", separams.seed);
  GetArg2(argc, argv, "-undersample", separams.undersample);
  GetArg2(argc, argv, "-oversample", separams.oversample);
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
		motif.set_expr_cutoff(max(0.9 * get_avg_pcorr(), 0.4));
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
			if(avg_pcorr > motif.get_expr_cutoff()) add_possible(g1);
		}
		for(int g1 = 0; g1 < ngenes; g1++)
			if(is_member(g1)) add_possible(g1);
	}
}

int SEModel::search_for_motif(const int worker, const int iter) {
	motif.clear_sites();
	motif.set_iter(iter);
	motif.set_expr_cutoff(0.8);
	motif.set_map(0.0);
	motif.set_spec(0.0);
	int i_worse = 0;
	int phase = 0;
	
	for(int g = 0; g < ngenes; g++)
		add_possible(g);
	seed_random_site();
	if(size() < 1) {
		cerr << "\t\t\tSeeding failed -- restarting..." << endl;
		return BAD_SEED;
	}
	
	clear_all_possible();
	while(possible_size() < separams.minsize * 5 && motif.get_expr_cutoff() > 0.4) {
		motif.set_expr_cutoff(motif.get_expr_cutoff() - 0.05);
		expand_search_around_mean(motif.get_expr_cutoff());
	}
	if(possible_size() < 2) {
		cerr << "\t\t\tBad search start -- no genes within " << separams.mincorr << endl;
		return BAD_SEARCH_SPACE;
	}

	compute_seq_scores();
	compute_expr_scores();
	set_seq_cutoff();
	
	motif.set_map(map_score());
	motif.set_spec(spec_score());
	print_status(cerr, 0, phase);
	Motif best_motif = motif;

	int i, old_width;
	phase = 1;
	for(i = 1; i < 10000 && phase < 3; i++) {
		expand_search_around_mean(motif.get_expr_cutoff());
		single_pass(motif.get_seq_cutoff());
		motif.set_spec(spec_score());
		motif.set_map(map_score());
		print_status(cerr, i, phase);
		old_width = motif.get_width();
		column_sample();
		if(motif.get_width() > old_width)
			compute_seq_scores();
		if(size() > ngenes/3) {
			cerr << "\t\t\tToo many sites! Restarting..." << endl;
			return TOO_MANY_SITES;
		}
		if(size() < 2) {
			cerr << "\t\t\tZero or one sites, reloading best motif..." << endl;
			motif = best_motif;
			phase++;
			print_status(cerr, i, phase);
		}
		if(motif.get_spec() > best_motif.get_spec()) {
			if(! archive.check_motif(motif)) {
				cerr << "\t\t\tToo similar! Restarting..." << endl;
				return TOO_SIMILAR;
			}
			cerr << "\t\t\t\tNew best motif!" << endl;
			best_motif = motif;
			i_worse = 0;
		} else {
			i_worse++;
			if(i_worse > 100 * phase) {
				motif = best_motif;
				if(size() < 2) {
					cerr << "\t\t\tLess than 2 genes at bad move threshold! Restarting..." << endl;
					return TOO_FEW_SITES;
				}
				cerr << "\t\t\tReached bad move threshold, reloading best motif..." << endl;
				phase++;
				compute_seq_scores();
				compute_expr_scores();
				set_expr_cutoff();
				set_seq_cutoff();
			}
		}
	}
	
	motif = best_motif;
	cerr << "\t\t\tRunning final greedy pass...";
	single_pass(motif.get_seq_cutoff(), true);
	cerr << "done." << endl;
	motif.orient();
	compute_seq_scores();
	compute_expr_scores();
	motif.set_spec(spec_score());
	motif.set_map(map_score());
	print_status(cerr, i, phase);
	
	if(size() < separams.minsize) {
		cerr << "\t\t\tToo few sites! Restarting..." << endl;
		return TOO_FEW_SITES;
	}
	if(size() > ngenes/3) {
		cerr << "\t\t\tToo many sites! Restarting..." << endl;
		return TOO_MANY_SITES;
	}
	if(! archive.check_motif(motif)) {
		cerr << "\t\t\tToo similar! Restarting..." << endl;
		return TOO_SIMILAR;
	}
	
	char tmpfilename[30], motfilename[30];
	sprintf(tmpfilename, "%d.%d.mot.tmp", worker, iter);
	sprintf(motfilename, "%d.%d.mot", worker, iter);
	ofstream motout(tmpfilename);
	motif.write(motout);
	motout.close();
	rename(tmpfilename, motfilename);
	cerr << "\t\t\tWrote motif to " << motfilename << endl;
	return 0;
}

bool SEModel::consider_motif(const char* filename) {
	ifstream motin(filename);
	motif.clear_sites();
	motif.read(motin);
	bool ret = archive.consider_motif(motif);
	motin.close();
	return ret;
}

void SEModel::print_status(ostream& out, const int i, const int phase) {
	out << "\t\t\t"; 
	out << setw(5) << i;
	out << setw(3) << phase;
	int prec = cerr.precision(2);
	out << setw(5) << setprecision(2) << motif.get_expr_cutoff();
	out << setw(10) << setprecision(4) << motif.get_seq_cutoff();
	cerr.precision(prec);
	out << setw(5) << motif_size();
	out << setw(5) << size();
	out << setw(5) << possible_size();
	if(size() > 0) {
		out << setw(40) << motif.consensus();
		out << setw(15) << motif.get_spec();
		out << setw(15) << motif.get_map();
	} else {
		out << setw(40) << "-----------";
		out << setw(15) << "---";
		out << setw(15) << "---";
	}
	out << endl;
}

void SEModel::full_output(ostream &fout){
  Motif* s;
	for(int j = 0; j < archive.nmots(); j++){
    s = archive.return_best(j);
    if(s->get_spec() > 1){
			fout << "Motif " << j + 1 << endl;
			s->write(fout);
		}
    else break;
  }
}

void SEModel::full_output(char *name){
  ofstream fout(name);
  full_output(fout);
}

void SEModel::set_seq_cutoff() {
	int expn, seqn, isect;
	seqn = expn = isect = 0;
	double best_c, best_po;
	best_po = 1;
	best_c = 0;
	double po, c;
	vector<double>::iterator sc_iter = expscores.begin();
	for(; sc_iter != expscores.end(); ++sc_iter)
		if(*sc_iter >= motif.get_expr_cutoff()) expn++;
	vector<struct idscore>::iterator rank_iter = seqranks.begin();
	c = rank_iter->score;
	best_c = c;
	for(; rank_iter->score > 0.01 && rank_iter != seqranks.end(); ++rank_iter) {
		if(c > rank_iter->score) {
			po = prob_overlap(expn, seqn, isect, ngenes);
			// cerr << "\t\t\t\t" << c << ":\t" << isect << "/(" << seqn << "," << expn << ")/" << ngenes << "\t" << -log10(po) << endl;
			if(po <= best_po && isect >= 0) {
				best_c = c;
				best_po = po;
			}
		}
		seqn++;
		if(expscores[rank_iter->id] >= motif.get_expr_cutoff()) {
			isect++;
			c = rank_iter->score;
		}
  }
	
	cerr << "\t\t\tSetting sequence cutoff to " << best_c << endl;
	motif.set_seq_cutoff(best_c);
}

void SEModel::set_expr_cutoff() {
	set_expr_cutoff_spec();
}

void SEModel::set_expr_cutoff_slowexpand() {
	if(size() > separams.minsize) {
		double exprcut = 0.0;
		for(int i = 0; i < ngenes; i++)
			if(is_member(i))
				exprcut += seqscores[i];
		exprcut /= size();
		cerr << "\t\t\tAverage expression score was " << exprcut << endl;
		double expand = (1 - exprcut) * 1.05;
		exprcut = 1 - expand;
		cerr << "\t\t\tSetting expression cutoff to " << exprcut << endl;
		motif.set_expr_cutoff(exprcut);
	}
}

void SEModel::set_expr_cutoff_spec() {
	int seqn = 0, expn, isect;
	double best_po = 1;
	double best_exprcut = 0;
	double po, c;
	for(int i = 0; i < ngenes; i++)
		if(seqscores[i] >= motif.get_seq_cutoff()) seqn++;
	for(c = 0.0; c <= 1; c += 0.01) {
		expn = isect = 0;
		for(int i = 0; i < ngenes; i++) {
			if(expranks[i].score >= c) {
				expn++;
				if(seqscores[expranks[i].id] >= motif.get_seq_cutoff())
					isect++;
			} else {
				break;
			}
		}
		po = prob_overlap(seqn, expn, isect, ngenes);
		// cerr << "\t\t\t\t" << c << ":\t\t" << isect << "/(" << expn << "," << seqn << ")/" << ngenes << "\t\t" << po << endl; 
		if(po <= best_po && expn >= separams.minsize * 3) {
			best_po = po;
			best_exprcut = c;
		}
	}
	cerr << "\t\t\tSetting expression cutoff to " << best_exprcut << endl;
	motif.set_expr_cutoff(best_exprcut);
}

void SEModel::print_possible(ostream& out) {
	for(int g = 0; g < ngenes; g++)
		if(is_possible(g)) out << nameset[g] << endl;
}

