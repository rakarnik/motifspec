#include "motifsearch.h"

MotifSearch::MotifSearch(const vector<string>& names,
		const vector<string>& seqs,
		const int nc,
		const int order,
		const double sim_cut,
		const int maxm) :
nameset(names),
ngenes(names.size()),
seqset(seqs),
bgmodel(seqset, order),
motif(seqset, nc, params.pseudo, params.backfreq),
select_sites(seqset, nc, params.pseudo, params.backfreq),
archive(seqset, sim_cut, maxm, params.pseudo, params.backfreq),
seqscores(ngenes),
seqranks(ngenes),
bestpos(ngenes),
beststrand(ngenes){
	set_default_params();
}

void MotifSearch::set_default_params(){
	params.expect = 10;
	params.minpass = 100;
	params.seed = -1;
	params.psfact = 0.1;
	params.weight = 0.5; 
	params.npass = 1000000;
	params.fragment = true;
	params.flanking = 0;
	params.undersample = 1;
	params.oversample = 1;
	params.minsize = 5;
	params.mincorr = 0.4;
}

void MotifSearch::set_final_params(){
	params.npseudo = params.expect * params.psfact;
	params.backfreq.push_back((1 - bgmodel.gcgenome())/2.0);
	params.backfreq.push_back(bgmodel.gcgenome()/2.0);
	params.backfreq.push_back(params.backfreq[0]);
	params.backfreq.push_back(params.backfreq[1]);
	for(int i = 0; i < 4; i++) {
		params.pseudo.push_back(params.npseudo * params.backfreq[i]);
	}
	params.maxlen = 3 * motif.get_width();
	params.nruns = motif.total_positions() / params.expect / motif.ncols() / params.undersample * params.oversample;
	params.select = 5.0;
	params.minprob[0] = 0.000001;
	params.minprob[1] = 0.0002;
	params.minprob[2] = 0.01;
	params.minprob[3] = 0.2;
}

void MotifSearch::ace_initialize(){
	ran_int.set_seed(params.seed);
	params.seed = ran_int.seed();
	ran_int.set_range(0, RAND_MAX);
	ran_dbl.set_seed(ran_int.rnum());
	ran_dbl.set_range(0.0, 1.0);
}

int MotifSearch::total_positions() const {
	return motif.total_positions();
}

int MotifSearch::positions_in_search_space() const {
	return motif.positions_in_search_space();
}

void MotifSearch::clear_sites() {
	motif.clear_sites();
}

int MotifSearch::size() const {
	return motif.seqs_with_sites();
}

int MotifSearch::motif_size() const {
	return motif.number();
}

bool MotifSearch::is_member(const int gene) const {
	return motif.seq_has_site(gene);
}

void MotifSearch::genes(int* genes) const {
	int count = 0;
	for (int g = 0; g < ngenes; g++) {
		if (is_member(g)) {
			genes[count] = g;
			count++;
		}
	}
	assert(count == size());
}

void MotifSearch::update_seq_count() {
	int seqn = 0, isect = 0;
	vector<struct idscore>::iterator ids = seqranks.begin();
	for(; ids != seqranks.end() && ids->score >= motif.get_seq_cutoff(); ++ids) {
		seqn++;
		if(motif.in_search_space(ids->id))
			isect++;
	}
	assert(isect <= seqn);
	assert(isect <= motif.get_search_space_size());
	motif.set_above_seqc(seqn);
	motif.set_above_cutoffs(isect);
}

void MotifSearch::seed_random_site() {
	assert(motif.get_search_space_size() > 0);

	int chosen_possible, chosen_seq, chosen_posit;
	bool watson;

	chosen_seq = chosen_posit = -1;

	/* First choose a sequence */
	ran_int.set_range(0, motif.get_search_space_size() - 1);
	chosen_possible = ran_int.rnum();
	int g;
	for(g = 0; g < ngenes; g++) {
		if(motif.in_search_space(g) && chosen_possible == 0) break; 
		if(motif.in_search_space(g)) chosen_possible--;
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

void MotifSearch::calc_matrix(double* score_matrix) {
	motif.calc_score_matrix(score_matrix);
}

double MotifSearch::score_site(double* score_matrix, const int c, const int p, const bool s) {
	double ms = motif.score_site(score_matrix, c, p, s);
	double bs = bgmodel.score_site(motif.first_column(), motif.last_column(), motif.get_width(), c, p, s);
	return fastexp(ms - bs);
}

void MotifSearch::single_pass(bool greedy) {
	double ap = params.weight * motif.get_search_space_size(); 
	ap += (1 - params.weight) * motif.number();
	ap /= 2.0 * motif.positions_in_search_space();
	double* score_matrix = new double[4 * motif.ncols()];
	calc_matrix(score_matrix);
	motif.remove_all_sites();
	select_sites.remove_all_sites();

	double Lw, Lc, Pw, Pc, F;
	int considered = 0;
	int gadd = -1, jadd = -1;
	int width = motif.get_width();
	for(int g = 0; g < seqset.num_seqs(); g++){
		if (! motif.in_search_space(g)) continue;
		for(int j = 0; j < seqset.len_seq(g) - width; j++){
			Lw = score_site(score_matrix, g, j, 1);
			Lc = score_site(score_matrix, g, j, 0);
			Pw = Lw * ap/(1.0 - ap + Lw * ap);
			Pc = Lc * ap/(1.0 - ap + Lc * ap);
			F = Pw + Pc - Pw * Pc;
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

void MotifSearch::single_pass_select(bool greedy) {
	double ap = params.weight * motif.get_search_space_size(); 
	ap += (1 - params.weight) * motif.number();
	ap /= 2.0 * motif.positions_in_search_space();
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
		if (! motif.in_search_space(g)) continue;
		if(j < 0 || j + width > seqset.len_seq(g)) continue;
		Lw = score_site(score_matrix, g, j, 1);
		Lc = score_site(score_matrix, g, j, 0);
		Pw = Lw * ap/(1.0 - ap + Lw * ap);
		Pc = Lc * ap/(1.0 - ap + Lc * ap);
		F = Pw + Pc - Pw * Pc;
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

void MotifSearch::compute_seq_scores() {
	double ap = params.weight * motif.get_search_space_size(); 
	ap += (1 - params.weight) * motif.number();
	ap /= 2.0 * motif.positions_in_search_space();
	double* score_matrix = new double[4 * motif.ncols()];
	calc_matrix(score_matrix);
	double Lw, Lc, Pw, Pc, F, bestF;
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
				beststrand[g] = Pw > Pc? 1 : 0;
			}
		}
		seqscores[g] = bestF;
		seqranks[g].id = g;
		seqranks[g].score = seqscores[g];
	}
	delete [] score_matrix;
	sort(seqranks.begin(), seqranks.end(), isc);
}

void MotifSearch::compute_seq_scores_minimal() {
	double ap = params.weight * motif.get_search_space_size(); 
	ap += (1 - params.weight) * motif.number();
	ap /= 2.0 * motif.positions_in_search_space();
	double* score_matrix = new double[4 * motif.ncols()];
	calc_matrix(score_matrix);
	int width = motif.get_width();
	double Lw, Lc, Pw, Pc, F;
	for(int g = 0; g < seqset.num_seqs(); g++) {
		// Some best positions might have been invalidated by column sampling
		// We mark these as invalid and don't score them
		if(bestpos[g] >= 0 && bestpos[g] + width <= seqset.len_seq(g)) {
			Lw = score_site(score_matrix, g, bestpos[g], 1);
			Lc = score_site(score_matrix, g, bestpos[g], 0);
			Pw = Lw * ap/(1.0 - ap + Lw * ap);
			Pc = Lc * ap/(1.0 - ap + Lc * ap);
			F = Pw + Pc - Pw * Pc;
		} else {
			bestpos[g] = -1;
			F = 0.0;
		}
		seqscores[g] = F;
		seqranks[g].id = g;
		seqranks[g].score = seqscores[g];
	}
	delete [] score_matrix;
	sort(seqranks.begin(), seqranks.end(), isc);
	if(seqranks[0].score < 0.85)
		compute_seq_scores();
}

double MotifSearch::score() {
	double w = 1.0;
	double sc = spec_score() * w + matrix_score() * (1 - w);
	return sc;
}

double MotifSearch::matrix_score() {
	double ms = 0.0;
	float* freq_matrix = new float[4 * motif.ncols()];
	motif.calc_freq_matrix(freq_matrix);
	int nc = motif.ncols();
	int w = motif.get_width();
	double sc[] = {0.0,0.0,0.0,0.0};
	for(int i = 0; i < 4 * nc; i += 4) {
		for(int j = 0; j < 4; j++) {
			ms += gammaln((double) freq_matrix[i + j] + params.pseudo[j]);
			sc[j] += freq_matrix[i + j];
		}
	}
	delete [] freq_matrix;
	ms -= nc * gammaln((double) motif.number() + params.npseudo);
	for (int j = 0; j < 4; j++)
		ms -= sc[j] * log(params.backfreq[j]);
	/* 
		 This factor arises from a modification of the model of Liu, et al 
		 in which the background frequencies of DNA bases are taken to be
		 constant for the organism under consideration
	 */
	double vg = 0.0;
	ms -= lnbico(w - 2, nc - 2);
	for(int j = 0; j < 4; j++)
		vg += gammaln(params.pseudo[j]);
	vg -= gammaln((double) (params.npseudo));
	ms -= ((double) nc * vg);
	return ms;
}

double MotifSearch::over_score() {
	double os = 0.0;
	double map_N = motif.positions_in_search_space();  
	double w = params.weight/(1.0-params.weight);
	double map_alpha = (double) params.expect * w;
	double map_beta = map_N * w - map_alpha;
	double map_success = (double)motif.number();
	os += ( gammaln(map_success+map_alpha)+gammaln(map_N-map_success+map_beta) );
	os -= ( gammaln(map_alpha)+gammaln(map_N + map_beta) );
	return os;
}

double MotifSearch::spec_score() {
	return -log_prob_overlap(motif.get_above_cutoffs(), motif.get_search_space_size(), motif.get_above_seqc(), ngenes);
}

void MotifSearch::output_params(ostream &fout){
	fout<<" expect =      \t"<<params.expect<<'\n';
	fout<<" minpass =     \t"<<params.minpass<<'\n';
	fout<<" seed =        \t"<<params.seed<<'\n';
	fout<<" numcols =     \t"<<motif.ncols()<<'\n';
	fout<<" undersample = \t"<<params.undersample<<'\n';
	fout<<" oversample = \t"<<params.oversample<<'\n';
}

void MotifSearch::modify_params(int argc, char *argv[]){
	GetArg2(argc, argv, "-expect", params.expect);
	GetArg2(argc, argv, "-minpass", params.minpass);
	GetArg2(argc, argv, "-seed", params.seed);
	GetArg2(argc, argv, "-undersample", params.undersample);
	GetArg2(argc, argv, "-oversample", params.oversample);
}

bool MotifSearch::consider_motif(const char* filename) {
	ifstream motin(filename);
	motif.clear_sites();
	motif.read(motin);
	bool ret = archive.consider_motif(motif);
	motin.close();
	return ret;
}

void MotifSearch::full_output(ostream &fout){
	fout << "Parameter values:\n";
	output_params(fout);
	fout << "\nInput sequences:\n";
	for(unsigned int x = 0; x < nameset.size(); x++)
		fout << "#" << x << '\t' << nameset[x] << endl;
	fout << '\n';
	archive.write(fout);
}

void MotifSearch::full_output(char *name){
	ofstream fout(name);
	full_output(fout);
}

void MotifSearch::set_seq_cutoff(const int phase) {
	int seqn = 0, isect = 0;
	double seqcut, best_seqcut = params.minprob[phase];
	double po, best_po = DBL_MAX;
	vector<struct idscore>::const_iterator sr_iter = seqranks.begin();
	for(seqcut = 0.999; seqcut >= params.minprob[phase]; seqcut -= 0.001) {
		while(sr_iter->score >= seqcut) {
			seqn++;
			if(motif.in_search_space(sr_iter->id))
				isect++;
			++sr_iter;
		}
		assert(isect <= seqn);
		assert(isect <= motif.get_search_space_size());
		po = log_prob_overlap(isect, seqn, motif.get_search_space_size(), ngenes);
		// cerr << "\t\t\t" << seqcut << '\t' << motif.get_seq_cutoff() << '\t' << seqn << '\t' << motif.get_search_space_size() << '\t' << isect << '\t' << po << '\n';
		if(po <= best_po) {
			best_seqcut = seqcut;
			best_po = po;
		}
	}
	// cerr << "\t\t\tSetting sequence cutoff to " << best_seqcut << "(minimum " << params.minprob[phase] << ")\n";
	motif.set_seq_cutoff(best_seqcut);
	update_seq_count();
}

void MotifSearch::print_status(ostream& out, const int i, const int phase) {
	out << setw(5) << i;
	out << setw(3) << phase;
	int prec = cerr.precision(2);
	out << setw(10) << setprecision(3) << motif.get_seq_cutoff();
	out << setw(10) << setprecision(2) << motif.get_ssp_cutoff();
	out << setw(7) << motif.seqs_with_sites();
	out << setw(7) << motif.get_above_cutoffs();
	out << setw(7) << motif.get_above_seqc();
	out << setw(7) << motif.get_search_space_size();
	if(size() > 0) {
		out << setw(50) << motif.consensus();
		out << setw(15) << setprecision(10) << motif.get_motif_score();
	} else {
		out << setw(50) << "---------------------------------------";
	}
	out << '\n';
	cerr.precision(prec);
}

