#include "semodel.h"

SEModel::SEModel(
const int s_type,
const vector<string>& names,
const vector<string>& seqs,
vector<vector <float> >& exprtab,
const vector<string>& sub,
vector<float>& sctab,
const int npts,
const int nc,
const int order,
const double sim_cut) :
search_type(s_type),
nameset(names),
ngenes(names.size()),
possible(ngenes),
seqset(seqs),
bgmodel(seqset, order),
motif(seqset, nc, separams.pseudo),
select_sites(seqset, nc, separams.pseudo),
archive(seqset, sim_cut, separams.pseudo),
expr(exprtab),
npoints(npts),
mean(npoints),
subset(sub),
scores(sctab),
mean_sc(0.0),
stdev_sc(0.0),
cumul_scores(ngenes),
seqscores(ngenes),
bestpos(ngenes),
seqranks(ngenes),
expscores(ngenes),
expranks(ngenes),
scranks(ngenes) {
	set_default_params();
	if(s_type == SCORE) {
		float ss = 0.0;
		for(int g = 0; g < ngenes; g++) {
			mean_sc += scores[g];
			ss += scores[g] * scores[g];
			cumul_scores[g] = g > 0? cumul_scores[g - 1] + scores[g] : scores[g];
		}
		mean_sc /= ngenes;
		stdev_sc = sqrt(ss/ngenes - mean_sc * mean_sc);
		for(int g = 0; g < ngenes; g++) {
			scores[g] = (scores[g] - mean_sc)/stdev_sc;
			scranks[g].id = g;
			scranks[g].score = scores[g];
		}
		sort(scranks.begin(), scranks.end(), isc);
	}
}

void SEModel::set_default_params(){
	separams.expect = 10;
	separams.minpass = 100;
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
	separams.backfreq[0] = separams.backfreq[3] = (1 - bgmodel.gcgenome())/2.0;
	separams.backfreq[1] = separams.backfreq[2] = bgmodel.gcgenome()/2.0;
	for(int i = 0; i < 4; i++) {
		separams.pseudo[i] = separams.npseudo * separams.backfreq[i];
	}
	separams.maxlen = 3 * motif.get_width();
	separams.nruns = motif.positions_available() / separams.expect / motif.ncols() / separams.undersample * separams.oversample;
	separams.select = 5.0;
	separams.minprob[0] = 0.000001;
	separams.minprob[1] = 0.001;
	separams.minprob[2] = 0.01;
	separams.minprob[3] = 0.1;
	separams.minscore[0] = 2.32;  // p ~ 0.01
	separams.minscore[1] = 2.32;  // p ~ 0.01
	separams.minscore[2] = 1.65;  // p ~ 0.05
	separams.minscore[3] = 1.28;  // p ~ 0.1
}

void SEModel::ace_initialize(){
	ran_int.set_seed(separams.seed);
	separams.seed = ran_int.seed();
	ran_int.set_range(0, RAND_MAX);
	ran_dbl.set_seed(ran_int.rnum());
	ran_dbl.set_range(0.0, 1.0);
}

void SEModel::reset_search_space() {
	clear_search_space();
	if(search_type == EXPRESSION) {
		for(int i = 0; i < ngenes; i++)
			add_to_search_space(i);
	} else if(search_type == SUBSET) {
		for(int i = 0; i < ngenes; i++)
			if(binary_search(subset.begin(), subset.end(), nameset[i]))
				add_to_search_space(i);
		assert((unsigned int) motif.get_possible() == subset.size());
	} else if(search_type == SCORE) {
		vector<struct idscore>::iterator scit = scranks.begin();
		for(; scit != scranks.end() && scit->score > 2.3; ++scit)
			add_to_search_space(scit->id);
	}
	assert(motif.get_possible() <= ngenes);
}

void SEModel::adjust_search_space() {
	int isect = 0;
	if(search_type == EXPRESSION) {
		calc_mean();
		clear_search_space();
		vector<struct idscore>::iterator exprit = expranks.begin();
		while(exprit != expranks.end() && exprit->score >= motif.get_expr_cutoff()) {
			add_to_search_space(exprit->id);
			if(seqscores[exprit->id] >= motif.get_seq_cutoff())
				isect++;
			++exprit;
		}
		motif.set_above_cutoffs(isect);
	} else if(search_type == SUBSET) {
		reset_search_space();
	} else if(search_type == SCORE) {
		clear_search_space();
		vector<struct idscore>::iterator scit = scranks.begin();
		while(scit != scranks.end() && scit->score >= motif.get_score_cutoff()) {
			add_to_search_space(scit->id);
			if(seqscores[scit->id] >= motif.get_seq_cutoff())
				isect++;
			++scit;
		}
		motif.set_above_cutoffs(isect);
	} else {
		cerr << "Invalid search type.\n";
	}
}

void SEModel::add_to_search_space(const int gene) {
	assert(motif.get_possible() < ngenes);
	assert(! possible[gene]);
	possible[gene] = true;
	motif.inc_possible();
	assert(motif.get_possible() <= ngenes);
}

void SEModel::remove_from_search_space(const int gene) {
	assert(motif.get_possible() > 0);
	assert(possible[gene]);
	possible[gene] = false;
	motif.dec_possible();
	assert(motif.get_possible() >= 0);
}

void SEModel::clear_search_space() {
	vector<bool> newp(ngenes, false);
	swap(possible, newp);
	motif.set_possible(0);
}

bool SEModel::is_in_search_space(const int gene) const {
	return possible[gene];
}

int SEModel::possible_size() const {
	return motif.get_possible();
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

void SEModel::update_seq_count() {
	int seqn = 0, isect = 0;
	vector<struct idscore>::iterator ids = seqranks.begin();
	for(; ids != seqranks.end() && ids->score >= motif.get_seq_cutoff(); ++ids) {
		seqn++;
		if((search_type == EXPRESSION && expscores[ids->id] >= motif.get_expr_cutoff())
				|| (search_type == SCORE && scores[ids->id] >= motif.get_score_cutoff()))
			isect++;
	}
	assert(isect <= seqn);
	assert(isect <= motif.get_possible());
	motif.set_above_seqc(seqn);
	motif.set_above_cutoffs(isect);
}

void SEModel::seed_random_site() {
	assert(motif.get_possible() > 0);

	int chosen_possible, chosen_seq, chosen_posit;
	bool watson;

	chosen_seq = chosen_posit = -1;

	/* First choose a sequence */
	ran_int.set_range(0, possible_size() - 1);
	chosen_possible = ran_int.rnum();
	int g;
	for(g = 0; g < ngenes; g++) {
		if(is_in_search_space(g) && chosen_possible == 0) break; 
		if(is_in_search_space(g)) chosen_possible--;
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

void SEModel::seed_high_scoring_site() {
	ran_dbl.set_range(0, cumul_scores[ngenes - 1]);
	double r = ran_dbl.rnum();
	int chosen_seq = 0;
	vector<float>::iterator csc_iter = cumul_scores.begin();
	for(; csc_iter != cumul_scores.end(); ++csc_iter) {
		if(*csc_iter >= r)
			break;
		chosen_seq++;
	}

	/* Now choose a site */
	int width = motif.get_width();
	int chosen_posit;
	bool watson;
	for(int j = 0; j < 50; j++) {
		ran_int.set_range(0, seqset.len_seq(chosen_seq) - width - 1);
		r = ran_dbl.rnum();//random (0,1)
		watson = (r > 0.5);
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

void SEModel::calc_matrix(double* score_matrix) {
	motif.calc_score_matrix(score_matrix, separams.pseudo);
}

double SEModel::score_site(double* score_matrix, const int c, const int p, const bool s) {
	double ms = motif.score_site(score_matrix, c, p, s);
	double bs = bgmodel.score_site(motif, c, p, s);
	return exp(ms - bs);
}

void SEModel::single_pass(bool greedy) {
	double ap = separams.weight * possible_size(); 
	ap += (1 - separams.weight) * motif.number();
	ap /= 2.0 * motif.positions_available(possible);
	double* score_matrix = new double[4 * motif.ncols()];
	calc_matrix(score_matrix);
	motif.remove_all_sites();
	select_sites.remove_all_sites();
	//will only update once per pass

	double Lw, Lc, Pw, Pc, F;
	int considered = 0;
	int gadd = -1, jadd = -1;
	int width = motif.get_width();
	for(int g = 0; g < seqset.num_seqs(); g++){
		if (! is_in_search_space(g)) continue;
		for(int j = 0; j < seqset.len_seq(g) - width; j++){
			Lw = score_site(score_matrix, g, j, 1);
			Lc = score_site(score_matrix, g, j, 0);
			Pw = Lw * ap/(1.0 - ap + Lw * ap);
			Pc = Lc * ap/(1.0 - ap + Lc * ap);
			F = Pw + Pc - Pw * Pc;//probability of either
			if(g == gadd && j < jadd + width) continue;
			if(F > motif.get_seq_cutoff()/5.0) select_sites.add_site(g, j, true);
			if(F < motif.get_seq_cutoff()) continue;
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
	delete [] score_matrix;
}

void SEModel::single_pass_select(bool greedy) {
	double ap = separams.weight * possible_size(); 
	ap += (1 - separams.weight) * motif.number();
	ap /= 2.0 * motif.positions_available(possible);
	double* score_matrix = new double[4 * motif.ncols()];
	calc_matrix(score_matrix);
	motif.remove_all_sites();

	double Lw, Lc, Pw, Pc, F;
	int g, j;
	int gadd = -1, jadd = -1;
	int width = motif.get_width();
	int num_sites = select_sites.number();
	for(int i = 0; i < num_sites; i++) {
		g = select_sites.chrom(i);
		j = select_sites.posit(i);
		if (! is_in_search_space(g)) continue;
		if(j < 0 || j + width > seqset.len_seq(g)) continue;
		Lw = score_site(score_matrix, g, j, 1);
		Lc = score_site(score_matrix, g, j, 0);
		Pw = Lw * ap/(1.0 - ap + Lw * ap);
		Pc = Lc * ap/(1.0 - ap + Lc * ap);
		F = Pw + Pc - Pw * Pc;//probability of either
		if(g == gadd && j < jadd + width) continue;
		if(F <= motif.get_seq_cutoff()) continue;
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
	delete [] score_matrix;
}

void SEModel::compute_seq_scores() {
	double ap = separams.weight * possible_size(); 
	ap += (1 - separams.weight) * motif.number();
	ap /= 2.0 * motif.positions_available(possible);
	double* score_matrix = new double[4 * motif.ncols()];
	calc_matrix(score_matrix);
	double Lw, Lc, Pw, Pc, F, bestF;
	seqranks.clear();
	int width = motif.get_width();
	int len;
	for(int g = 0; g < ngenes; g++) {
		bestF = 0.0;
		bestpos[g] = -1;
		len = seqset.len_seq(g);
		for(int j = 0; j < len - width; j++) {
			Lw = score_site(score_matrix, g, j, 1);
			Lc = score_site(score_matrix, g, j, 0);
			Pw = Lw * ap/(1.0 - ap + Lw * ap);
			Pc = Lc * ap/(1.0 - ap + Lc * ap);
			F = Pw + Pc - Pw * Pc;
			if(F > bestF) {
				bestF = F;
				bestpos[g] = j;
			}
		}
		seqscores[g] = bestF;
		struct idscore ids = {g, seqscores[g]};
		seqranks.push_back(ids);
	}
	delete [] score_matrix;
	sort(seqranks.begin(), seqranks.end(), isc);
	update_seq_count();
}

void SEModel::compute_seq_scores_minimal() {
	double ap = separams.weight * possible_size(); 
	ap += (1 - separams.weight) * motif.number();
	ap /= 2.0 * motif.positions_available(possible);
	double* score_matrix = new double[4 * motif.ncols()];
	calc_matrix(score_matrix);
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
			Lw = score_site(score_matrix, g, bestpos[g], 1);
			Lc = score_site(score_matrix, g, bestpos[g], 0);
			Pw = Lw * ap/(1.0 - ap + Lw * ap);
			Pc = Lc * ap/(1.0 - ap + Lc * ap);
			F = Pw + Pc - Pw * Pc;
		}
		seqscores[g] = F;
		struct idscore ids = {g, seqscores[g]};
		seqranks.push_back(ids);
	}
	delete [] score_matrix;
	sort(seqranks.begin(), seqranks.end(), isc);
	if(seqranks[0].score <= 0.85)
		compute_seq_scores();
	else
		update_seq_count();
}

void SEModel::compute_expr_scores() {
	calc_mean();
	expranks.clear();
	for(int g = 0; g < ngenes; g++) {
		expscores[g] = corr(mean, expr[g]);
		struct idscore ids = {g, expscores[g]};
		expranks.push_back(ids);
	}
	sort(expranks.begin(), expranks.end(), isc);
	adjust_search_space();
}

bool SEModel::column_sample(const bool add, const bool remove){
	bool changed = false;
	int freq[4];
	int width = motif.get_width();
	// Compute scores for current and surrounding columns
	int max_left, max_right;
	max_left = max_right = (motif.get_max_width() - width)/2;
	motif.columns_open(max_left, max_right);
	int cs_span = max_left + max_right + width;
	vector<double> wtx(cs_span);
	int x = max_left;
	//wtx[x + c] will refer to the weight of pos c in the usual numbering
	double wt;
	double best_wt = -DBL_MAX;
	for(int i = 0; i < cs_span; i++) {
		wtx[i] = -DBL_MAX;
		if(motif.column_freq(i - x, freq)){
			wt = 0.0;
			for(int j = 0;j < 4; j++){
				wt += gammaln(freq[j] + separams.pseudo[j]);
				wt -= (double)freq[j] * log(separams.backfreq[j]);
			}
			wtx[i] = wt;
			if(wt > best_wt) best_wt = wt;
		}
	}

	// Penalize outermost columns for length
	double scale = 0.0;
	if(best_wt > 100.0) scale = best_wt - 100.0; //keep exp from overflowing
	for(int i = 0; i < cs_span; i++){
		wtx[i]  -= scale;
		wtx[i] = exp(wtx[i]);
		int newwidth = width;
		if(i < x)
			newwidth += (x - i);
		else if(i > (x + width - 1))
			newwidth += (i - x - width + 1);
		wtx[i] /= bico(newwidth - 2, motif.ncols() - 2);
	}

	// Find worst column in the current set, best column outside current set
	best_wt = - DBL_MAX;
	double worst_wt = DBL_MAX;
	int best_col = 0, worst_col = 0;
	for(int i = 0; i < cs_span; i++) {
		if(motif.has_col(i - x)) { 
			if(wtx[i] < worst_wt) {
				worst_col = i - x;
				worst_wt = wtx[i];
			}
		} else if(wtx[i] > best_wt) {
			best_col = i - x;
			best_wt = wtx[i];
		}
	}

	// Swap if best column outside set is better than worst column inside set
	if(best_wt > worst_wt) {
		if(add) {
			motif.add_col(best_col);
			select_sites.add_col(best_col);
			if(best_col < 0)
				worst_col -= best_col;
			changed = true;
		}
		if(remove) {
			motif.remove_col(worst_col);
			select_sites.remove_col(worst_col);
			changed = true;
		}
	}

	return changed;
}

double SEModel::score() {
	double w = 1.0;
	double sc = spec_score() * w + matrix_score() * (1 - w);
	return sc;
}

double SEModel::matrix_score() {
	double ms = 0.0;
	int* freq_matrix = new int[4 * motif.ncols()];
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
	delete [] freq_matrix;
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

double SEModel::over_score() {
	double os = 0.0;
	double map_N = motif.positions_available(possible);  
	double w = separams.weight/(1.0-separams.weight);
	double map_alpha = (double) separams.expect * w;
	double map_beta = map_N * w - map_alpha;
	double map_success = (double)motif.number();
	os += ( gammaln(map_success+map_alpha)+gammaln(map_N-map_success+map_beta) );
	os -= ( gammaln(map_alpha)+gammaln(map_N + map_beta) );
	return os;
}

double SEModel::spec_score() {
	return -log_prob_overlap(motif.get_above_cutoffs(), motif.get_possible(), motif.get_above_seqc(), ngenes);
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

int SEModel::search_for_motif(const int worker, const int iter, const string outfile) {
	if(search_type == EXPRESSION) {
		return search_for_motif_expr(worker, iter, outfile);
	} else if(search_type == SUBSET) {
		return search_for_motif_subset(worker, iter, outfile);
	} else if(search_type == SCORE) {
		return search_for_motif_score(worker, iter, outfile);
	} else {
		cerr << "Invalid search type!\n";
		exit(1);
	}
}

int SEModel::search_for_motif_expr(const int worker, const int iter, const string outfile) {
	motif.clear_sites();
	select_sites.clear_sites();
	stringstream iterstr;
	iterstr << worker << '.' << iter;
	motif.set_iter(iterstr.str());
	int phase = 0;

	reset_search_space();
	seed_random_site();
	if(size() < 1) {
		cerr << "\t\t\tSeeding failed -- restarting...\n";
		return BAD_SEED;
	}

	motif.set_expr_cutoff(0.8);
	adjust_search_space();
	while(possible_size() < separams.minsize * 5 && motif.get_expr_cutoff() > 0.4) {
		motif.set_expr_cutoff(motif.get_expr_cutoff() - 0.05);
		adjust_search_space();
	}
	if(possible_size() < 2) {
		cerr << "\t\t\tBad search start -- no genes within " << separams.mincorr << '\n';
		return BAD_SEARCH_SPACE;
	}

	compute_seq_scores();
	compute_expr_scores();
	set_seq_cutoff(phase);
	motif.set_motif_score(score());
	print_status(cerr, 0, phase);
	Motif best_motif = motif;

	int i, i_worse = 0;
	double r = 0.0;
	phase = 1;
	for(i = 1; i < 10000 && phase < 3; i++) {
		adjust_search_space();
		if(i_worse == 0)
			single_pass(false);
		else
			single_pass_select(false);
		for(int j = 0; j < 3; j++)
			if(! column_sample(true, true)) break;
		compute_seq_scores_minimal();
		compute_expr_scores();
		motif.set_motif_score(score());
		print_status(cerr, i, phase);
		r = ran_dbl.rnum();
		if(r > 0.5) {
			column_sample(true, false);
		} else {
			column_sample(false, true);
		}
		if(size() > ngenes/3) {
			cerr << "\t\t\tToo many sites! Restarting...\n";
			return TOO_MANY_SITES;
		}
		if(size() < 2) {
			cerr << "\t\t\tZero or one sites, reloading best motif...\n";
			phase++;
			best_motif.set_seq_cutoff(separams.minprob[phase]);
			motif = best_motif;
			select_sites = best_motif;
			compute_seq_scores_minimal();
			compute_expr_scores();
			set_seq_cutoff(phase);
			set_expr_cutoff();
			print_status(cerr, i, phase);
			i_worse = 0;
			continue;
		}
		if(motif.get_motif_score() > best_motif.get_motif_score()) {
			if(! archive.check_motif(motif)) {
				cerr << "\t\t\tToo similar! Restarting...\n";
				return TOO_SIMILAR;
			}
			cerr << "\t\t\t\tNew best motif!\n";
			best_motif = motif;
			i_worse = 0;
		} else {
			i_worse++;
			if(i_worse > separams.minpass * phase) {
				if(size() < 2) {
					cerr << "\t\t\tLess than 2 genes at bad move threshold! Restarting...\n";
					return TOO_FEW_SITES;
				}
				cerr << "\t\t\tReached bad move threshold, reloading best motif...\n";
				phase++;
				best_motif.set_seq_cutoff(separams.minprob[phase]);
				motif = best_motif;
				select_sites = best_motif;
				compute_seq_scores_minimal();
				compute_expr_scores();
				set_seq_cutoff(phase);
				set_expr_cutoff();
				print_status(cerr, i, phase);
			}
		}
	}

	cerr << "\t\t\tRunning final greedy pass...";
	single_pass(true);
	cerr << "done.\n";
	motif.orient();
	compute_seq_scores();
	compute_expr_scores();
	motif.set_motif_score(score());
	print_status(cerr, i, phase);

	if(size() < separams.minsize) {
		cerr << "\t\t\tToo few sites! Restarting...\n";
		return TOO_FEW_SITES;
	}
	if(size() > ngenes/3) {
		cerr << "\t\t\tToo many sites! Restarting...\n";
		return TOO_MANY_SITES;
	}
	if(! archive.check_motif(motif)) {
		cerr << "\t\t\tToo similar! Restarting...\n";
		return TOO_SIMILAR;
	}

	archive.consider_motif(motif);
	char tmpfilename[30], motfilename[30];
	sprintf(tmpfilename, "%s.%d.%d.mot.tmp", outfile.c_str(), worker, iter);
	sprintf(motfilename, "%s.%d.%d.mot", outfile.c_str(), worker, iter);
	ofstream motout(tmpfilename);
	motif.write(motout);
	motout.close();
	rename(tmpfilename, motfilename);
	cerr << "\t\t\tWrote motif to " << motfilename << '\n';
	return 0;
}

int SEModel::search_for_motif_subset(const int worker, const int iter, const string outfile) {
	motif.clear_sites();
	select_sites.clear_sites();
	stringstream iterstr;
	iterstr << worker << '.' << iter;
	motif.set_iter(iterstr.str());
	int phase = 0;

	reset_search_space();
	seed_random_site();
	if(size() < 1) {
		cerr << "\t\t\tSeeding failed -- restarting...\n";
		return BAD_SEED;
	}

	compute_seq_scores();
	set_seq_cutoff(phase);
	motif.set_motif_score(score());
	print_status(cerr, 0, phase);
	Motif best_motif = motif;

	int i, i_worse = 0;
	double r = 0.0;
	phase = 1;
	for(i = 1; i < 10000 && phase < 3; i++) {
		if(i_worse == 0)
			single_pass(false);
		else
			single_pass_select(false);
		for(int j = 0; j < 3; j++)
			if(! column_sample(true, true)) break;
		compute_seq_scores_minimal();
		motif.set_motif_score(score());
		print_status(cerr, i, phase);
		r = ran_dbl.rnum();
		if(r > 0.5) {
			column_sample(true, false);
		} else {
			column_sample(false, true);
		}
		if(size() > ngenes/3) {
			cerr << "\t\t\tToo many sites! Restarting...\n";
			return TOO_MANY_SITES;
		}
		if(size() < 2) {
			cerr << "\t\t\tZero or one sites, reloading best motif...\n";
			phase++;
			best_motif.set_seq_cutoff(separams.minprob[phase]);
			motif = best_motif;
			select_sites = best_motif;
			compute_seq_scores_minimal();
			set_seq_cutoff(phase);
			print_status(cerr, i, phase);
			i_worse = 0;
			continue;
		}
		if(motif.get_motif_score() > best_motif.get_motif_score()) {
			if(! archive.check_motif(motif)) {
				cerr << "\t\t\tToo similar! Restarting...\n";
				return TOO_SIMILAR;
			}
			cerr << "\t\t\t\tNew best motif!\n";
			best_motif = motif;
			i_worse = 0;
		} else {
			i_worse++;
			if(i_worse > separams.minpass * phase) {
				if(size() < 2) {
					cerr << "\t\t\tLess than 2 genes at bad move threshold! Restarting...\n";
					return TOO_FEW_SITES;
				}
				cerr << "\t\t\tReached bad move threshold, reloading best motif...\n";
				phase++;
				best_motif.set_seq_cutoff(separams.minprob[phase]);
				motif = best_motif;
				select_sites = best_motif;
				compute_seq_scores_minimal();
				set_seq_cutoff(phase);
				print_status(cerr, i, phase);
			}
		}
	}

	cerr << "\t\t\tRunning final greedy pass...";
	single_pass(true);
	cerr << "done.\n";
	motif.orient();
	compute_seq_scores();
	motif.set_motif_score(score());
	print_status(cerr, i, phase);

	if(size() < separams.minsize) {
		cerr << "\t\t\tToo few sites! Restarting...\n";
		return TOO_FEW_SITES;
	}
	if(size() > ngenes/3) {
		cerr << "\t\t\tToo many sites! Restarting...\n";
		return TOO_MANY_SITES;
	}
	if(! archive.check_motif(motif)) {
		cerr << "\t\t\tToo similar! Restarting...\n";
		return TOO_SIMILAR;
	}

	archive.consider_motif(motif);
	char tmpfilename[30], motfilename[30];
	sprintf(tmpfilename, "%s.%d.%d.mot.tmp", outfile.c_str(), worker, iter);
	sprintf(motfilename, "%s.%d.%d.mot", outfile.c_str(), worker, iter);
	ofstream motout(tmpfilename);
	motif.write(motout);
	motout.close();
	rename(tmpfilename, motfilename);
	cerr << "\t\t\tWrote motif to " << motfilename << '\n';
	return 0;
}

int SEModel::search_for_motif_score(const int worker, const int iter, const string outfile) {
	motif.clear_sites();
	select_sites.clear_sites();
	stringstream iterstr;
	iterstr << worker << '.' << iter;
	motif.set_iter(iterstr.str());
	int phase = 0;

	reset_search_space();
	/*
		 seed_random_site();
		 if(size() < 1) {
		 cerr << "\t\t\tSeeding failed -- restarting...\n";
		 return BAD_SEED;
		 }
	 */
	motif.add_site(63, 19, 0);

	compute_seq_scores();
	set_seq_cutoff(phase);
	set_score_cutoff(phase);
	motif.set_motif_score(score());
	print_status(cerr, 0, phase);
	Motif best_motif = motif;

	int i, i_worse = 0;
	double r = 0.0;
	phase = 1;
	for(i = 1; i < 10000 && phase < 3; i++) {
		adjust_search_space();
		if(i_worse == 0)
			single_pass(false);
		else
			single_pass_select(false);
		for(int j = 0; j < 3; j++)
			if(! column_sample(true, true)) break;
		compute_seq_scores_minimal();
		motif.set_motif_score(score());
		print_status(cerr, i, phase);
		r = ran_dbl.rnum();
		if(r > 0.5) {
			column_sample(true, false);
		} else {
			column_sample(false, true);
		}
		if(size() > ngenes/3) {
			cerr << "\t\t\tToo many sites! Restarting...\n";
			return TOO_MANY_SITES;
		}
		if(size() < 2) {
			cerr << "\t\t\tZero or one sites, reloading best motif...\n";
			phase++;
			best_motif.set_seq_cutoff(separams.minprob[phase]);
			motif = best_motif;
			select_sites = best_motif;
			compute_seq_scores();
			set_seq_cutoff(phase);
			set_score_cutoff(phase);
			motif.set_motif_score(score());
			print_status(cerr, i, phase);
			i_worse = 0;
			continue;
		}
		if(motif.get_motif_score() > best_motif.get_motif_score()) {
			if(! archive.check_motif(motif)) {
				cerr << "\t\t\tToo similar! Restarting...\n";
				return TOO_SIMILAR;
			}
			cerr << "\t\t\t\tNew best motif!\n";
			best_motif = motif;
			i_worse = 0;
		} else {
			i_worse++;
			if(i_worse > separams.minpass * phase) {
				if(size() < 2) {
					cerr << "\t\t\tLess than 2 genes at bad move threshold! Restarting...\n";
					return TOO_FEW_SITES;
				}
				cerr << "\t\t\tReached bad move threshold, reloading best motif...\n";
				phase++;
				best_motif.set_seq_cutoff(separams.minprob[phase]);
				motif = best_motif;
				select_sites = best_motif;
			}
		}
	}

	cerr << "\t\t\tRunning final greedy pass...";
	single_pass(true);
	cerr << "done.\n";
	motif.orient();
	print_status(cerr, i, phase);

	if(size() < separams.minsize) {
		cerr << "\t\t\tToo few sites! Restarting...\n";
		return TOO_FEW_SITES;
	}
	if(size() > ngenes/3) {
		cerr << "\t\t\tToo many sites! Restarting...\n";
		return TOO_MANY_SITES;
	}
	if(! archive.check_motif(motif)) {
		cerr << "\t\t\tToo similar! Restarting...\n";
		return TOO_SIMILAR;
	}

	archive.consider_motif(motif);
	char tmpfilename[30], motfilename[30];
	sprintf(tmpfilename, "%s.%d.%d.mot.tmp", outfile.c_str(), worker, iter);
	sprintf(motfilename, "%s.%d.%d.mot", outfile.c_str(), worker, iter);
	ofstream motout(tmpfilename);
	motif.write(motout);
	motout.close();
	rename(tmpfilename, motfilename);
	cerr << "\t\t\tWrote motif with score " << motif.get_motif_score() << " to " << motfilename << '\n';
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
	out << setw(5) << i;
	out << setw(3) << phase;
	int prec = cerr.precision(2);
	out << setw(10) << setprecision(3) << motif.get_seq_cutoff();
	if(search_type == EXPRESSION)
		out << setw(10) << setprecision(2) << motif.get_expr_cutoff();
	else if(search_type == SUBSET)
		out << setw(10) << "N/A";
	else if(search_type == SCORE)
		out << setw(10) << motif.get_score_cutoff();
	out << setw(7) << motif.seqs_with_sites();
	out << setw(7) << motif.get_above_cutoffs();
	out << setw(7) << motif.get_above_seqc();
	out << setw(7) << motif.get_possible();
	if(size() > 0) {
		out << setw(35) << motif.consensus();
		out << setw(15) << setprecision(10) << motif.get_motif_score();
	} else {
		out << setw(50) << "---------------------------------------";
	}
	out << '\n';
	cerr.precision(prec);
}

void SEModel::full_output(ostream &fout){
	Motif* s;
	for(int j = 0; j < archive.nmots(); j++){
		s = archive.return_best(j);
		if(s->get_motif_score() > 1){
			fout << "Motif " << j + 1 << '\n';
			s->write(fout);
		}
		else break;
	}
}

void SEModel::full_output(char *name){
	ofstream fout(name);
	full_output(fout);
}

void SEModel::set_seq_cutoff(const int phase) {
	int seqn = 0, isect = 0;
	float seqcut, best_seqcut = separams.minprob[phase];
	double po, best_po = DBL_MAX;
	vector<struct idscore>::const_iterator sr_iter = seqranks.begin();
	for(seqcut = sr_iter->score; sr_iter->score >= separams.minprob[phase] && sr_iter != seqranks.end(); ++sr_iter) {
		if(seqcut >= sr_iter->score) {
			assert(isect <= motif.get_possible());
			assert(isect <= seqn);
			po = log_prob_overlap(isect, motif.get_possible(), seqn, ngenes);
			// cerr << "\t\t\t" << motif.get_score_cutoff() << '\t' << seqcut << '\t' << motif.get_possible() << '\t' << seqn << '\t' << isect << '\t' << po << '\n';
			if(po < best_po) {
				best_seqcut = seqcut;
				best_po = po;
			}
		}
		seqn++;
		if(is_in_search_space(sr_iter->id))
			isect++;
		seqcut = sr_iter->score;
	}
	// cerr << "\t\t\tSetting sequence cutoff to " << best_seqcut << "(minimum " << separams.minprob[phase] << ")\n";
	motif.set_seq_cutoff(best_seqcut);
	update_seq_count();
}

void SEModel::set_expr_cutoff() {
	int expn = 0, isect = 0;
	float exprcut, best_exprcut = 0.0;
	double po, best_po = DBL_MAX;
	vector<struct idscore>::const_iterator er_iter = expranks.begin();
	for(exprcut = 0.95; exprcut >= 0; exprcut -= 0.01) {
		while(er_iter->score >= exprcut) {
			expn++;
			if(seqscores[er_iter->id] > motif.get_seq_cutoff())
				isect++;
			++er_iter;
		}
		assert(isect <= expn);
		assert(isect <= motif.get_above_seqc());
		po = log_prob_overlap(isect, expn, motif.get_above_seqc(), ngenes);
		// cerr << "\t\t\t" << exprcut << '\t' << motif.get_seq_cutoff() << '\t' << expn << '\t' << motif.get_above_seqc() << '\t' << isect << '\t' << po << '\n';
		if(po <= best_po) {
			best_exprcut = exprcut;
			best_po = po;
		}
	}
	// cerr << "\t\t\tSetting expression cutoff to " << best_exprcut << '\n';
	motif.set_expr_cutoff(best_exprcut);
	adjust_search_space();
}

void SEModel::set_score_cutoff(const int phase) {
	int scn = 0, isect = 0;
	float sccut, best_sccut = 0.0;
	double po, best_po = DBL_MAX;
	vector<struct idscore>::const_iterator sc_iter = scranks.begin();
	for(sccut = sc_iter->score; sc_iter != scranks.end() && sc_iter->score > separams.minscore[phase]; ++sc_iter) {
		if(sccut >= sc_iter->score) {
			assert(isect <= motif.get_above_seqc());
			assert(isect <= scn);
			po = log_prob_overlap(isect, scn, motif.get_above_seqc(), ngenes);
			// cerr << "\t\t\t" << sccut << '\t' << motif.get_seq_cutoff() << '\t' << scn << '\t' << motif.get_above_seqc() << '\t' << isect << '\t' << po << '\n';
			if(po < best_po) {
				best_sccut = sccut;
				best_po = po;
			}
		}
		scn++;
		if(seqscores[sc_iter->id] >= motif.get_seq_cutoff())
			isect++;
		sccut = sc_iter->score;
	}
	// cerr << "\t\t\tSetting score cutoff to " << best_sccut << " (minimum " << separams.minscore[phase] << ")\n";
	motif.set_score_cutoff(best_sccut);
	adjust_search_space();
}

void SEModel::print_possible(ostream& out) {
	for(int g = 0; g < ngenes; g++)
		if(is_in_search_space(g)) out << nameset[g] << endl;
}


