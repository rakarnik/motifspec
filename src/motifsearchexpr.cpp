#include "motifsearchexpr.h"

MotifSearchExpr::MotifSearchExpr(const vector<string>& names,
																 const vector<string>& seqs,
																 const int nc,
																 const int order,
																 const double sim_cut,
																 const int maxm,
																 vector<vector <float> >& exprtab,
																 const int npts) :
MotifSearch(names, seqs, nc, order, sim_cut, maxm),
expr(exprtab),
npoints(npts),
mean(npoints),
expscores(ngenes),
expranks(ngenes) {
	reset_search_space();
}

void MotifSearchExpr::set_final_params() {
	MotifSearch::set_final_params();
	params.minprob[0] = 0.00001;
	params.minprob[1] = 0.002;
	params.minprob[2] = 0.01;
	params.minprob[3] = 0.2;
	params.minpass = 200;
}

void MotifSearchExpr::calc_mean() {
	int i, j;
	for (i = 0; i < npoints; i++) {
		mean[i] = 0;
		for (j = 0; j < ngenes; j++)
			if(is_member(j))
				mean[i] += expr[j][i];
		mean[i] /= size();
	}
}

void MotifSearchExpr::compute_expr_scores() {
	calc_mean();
	for(int g = 0; g < ngenes; g++) {
		expscores[g] = corr(mean, expr[g]);
		expranks[g].id = g;
		expranks[g].score = expscores[g];
	}
	sort(expranks.begin(), expranks.end(), isc);
	adjust_search_space();
}

void MotifSearchExpr::reset_search_space() {
	motif.clear_search_space();
	for(int i = 0; i < ngenes; i++)
			motif.add_to_search_space(i);
	assert(motif.get_search_space_size() > 0);
	assert(motif.get_search_space_size() <= ngenes);
}

void MotifSearchExpr::adjust_search_space() {
	int isect = 0;
	calc_mean();
	motif.clear_search_space();
	vector<struct idscore>::iterator exprit = expranks.begin();
	while(exprit != expranks.end() && exprit->score >= motif.get_ssp_cutoff()) {
		motif.add_to_search_space(exprit->id);
		if(seqscores[exprit->id] >= motif.get_seq_cutoff())
			isect++;
		++exprit;
	}
	motif.set_above_cutoffs(isect);
}

void MotifSearchExpr::set_search_space_cutoff(const int phase) {
	(void) phase;
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
	motif.set_ssp_cutoff(best_exprcut);
	adjust_search_space();
}

int MotifSearchExpr::search_for_motif(const int worker, const int iter, const string outfile) {
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
	
	motif.set_ssp_cutoff(0.8);
	while(motif.get_search_space_size() < params.minsize * 5 && motif.get_ssp_cutoff() > 0.4) {
		motif.set_ssp_cutoff(motif.get_ssp_cutoff() - 0.05);
		adjust_search_space();
	}
	if(motif.get_search_space_size() < 2) {
		cerr << "\t\t\tBad search start -- no genes within " << params.mincorr << '\n';
		return BAD_SEARCH_SPACE;
	}
	
	compute_expr_scores();
	compute_seq_scores();
	set_seq_cutoff(phase);
	update_seq_count();
	motif.set_motif_score(score());
	print_status(cerr, 0, phase);
	Motif best_motif = motif;
	
	int i, i_worse = 0;
	phase = 1;
	for(i = 1; i < 10000 && phase < 3; i++) {
		adjust_search_space();
		if(i_worse == 0) {
			single_pass(false);
		} else {
			single_pass_select(false);
		}
		for(int j = 0; j < 3; j++)
			if(! motif.column_sample()) break;
		if(i % 2 == 0 && motif.ncols() < 2 * motif.init_ncols()) {
			motif.column_sample(true, false);
		} else if(motif.ncols() > 6) {
			motif.column_sample(false, true);
		}
		compute_seq_scores_minimal();
		compute_expr_scores();
		motif.set_motif_score(score());
		print_status(cerr, i, phase);
		if(size() > ngenes/3) {
			cerr << "\t\t\tToo many sites! Restarting...\n";
			return TOO_MANY_SITES;
		}
		if(size() < 2) {
			cerr << "\t\t\tZero or one sites, reloading best motif...\n";
			phase++;
			best_motif.set_seq_cutoff(params.minprob[phase]);
			motif = best_motif;
			select_sites = best_motif;
			compute_seq_scores_minimal();
			compute_expr_scores();
			update_seq_count();
			set_seq_cutoff(phase);
			set_search_space_cutoff(phase);
			print_status(cerr, i, phase);
			i_worse = 0;
			continue;
		}
		if(motif.get_motif_score() > best_motif.get_motif_score() ||
				(motif.get_motif_score() > best_motif.get_motif_score() * 0.99 && motif.number() > best_motif.number())) {
			if(! archive.check_motif(motif)) {
				cerr << "\t\t\tToo similar! Restarting...\n";
				return TOO_SIMILAR;
			}
			cerr << "\t\t\t\tNew best motif!\n";
			best_motif = motif;
			i_worse = 0;
		} else {
			i_worse++;
			if(i_worse > params.minpass * phase) {
				if(size() < 2) {
					cerr << "\t\t\tLess than 2 genes at bad move threshold! Restarting...\n";
					return TOO_FEW_SITES;
				}
				cerr << "\t\t\tReached bad move threshold, reloading best motif...\n";
				phase++;
				best_motif.set_seq_cutoff(params.minprob[phase]);
				motif = best_motif;
				select_sites = best_motif;
				compute_seq_scores_minimal();
				compute_expr_scores();
				update_seq_count();
				set_seq_cutoff(phase);
				set_search_space_cutoff(phase);
				print_status(cerr, i, phase);
				i_worse = 0;
			}
		}
	}
	
	cerr << "\t\t\tRunning final greedy pass...";
	single_pass(true);
	cerr << "done.\n";
	motif.orient();
	compute_seq_scores();
	compute_expr_scores();
	update_seq_count();
	motif.set_motif_score(score());
	print_status(cerr, i, phase);
	
	if(size() < params.minsize) {
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


