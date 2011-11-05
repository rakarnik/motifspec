#include "motifsearchscore.h"

MotifSearchScore::MotifSearchScore(const vector<string>& names,
																 const vector<string>& seqs,
																 const int nc,
																 const int order,
																 const double sim_cut,
																 vector<float>& sctab) :
MotifSearch(names, seqs, nc, order, sim_cut),
scores(sctab),
scranks(ngenes) {
	for(int i = 0; i < ngenes; i++) {
		scranks[i].id = i;
		scranks[i].score = scores[i];
	}
	sort(scranks.begin(), scranks.end(), isc);
	reset_search_space();
}

void MotifSearchScore::reset_search_space() {
	motif.clear_search_space();
	vector<struct idscore>::iterator scit = scranks.begin();
	for(; scit != scranks.end() && motif.get_search_space_size() <= 100; ++scit)
		motif.add_to_search_space(scit->id);
	motif.set_score_cutoff(scit->score);
	assert(motif.get_search_space_size() > 0);
	assert(motif.get_search_space_size() <= ngenes);
}

void MotifSearchScore::adjust_search_space() {
	int isect = 0;
	motif.clear_search_space();
	vector<struct idscore>::iterator scit = scranks.begin();
	while(scit != scranks.end() && scit->score >= motif.get_score_cutoff()) {
		motif.add_to_search_space(scit->id);
		if(seqscores[scit->id] >= motif.get_seq_cutoff())
			isect++;
		++scit;
	}
	motif.set_above_cutoffs(isect);
}

void MotifSearchScore::set_search_space_cutoff(const int phase) {
	int scn = 0, isect = 0;
	float sccut, best_sccut = 0.0;
	double po, best_po = DBL_MAX;
	vector<struct idscore>::const_iterator sc_iter = scranks.begin();
	for(sccut = sc_iter->score; sc_iter != scranks.end(); ++sc_iter) {
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
	// cerr << "\t\t\tSetting score cutoff to " << best_sccut << " (minimum " << params.minscore[phase] << ")\n";
	motif.set_score_cutoff(best_sccut);
	adjust_search_space();
}

void MotifSearchScore::calc_matrix(double* score_matrix) {
	motif.calc_score_matrix(score_matrix, scores);
}

int MotifSearchScore::search_for_motif(const int worker, const int iter, const string outfile) {
	motif.clear_sites();
	select_sites.clear_sites();
	stringstream iterstr;
	iterstr << worker << '.' << iter;
	motif.set_iter(iterstr.str());
	int phase = 0;
	
	motif.set_score_cutoff(2.3);
	reset_search_space();
	seed_random_site();
	if(size() < 1) {
		cerr << "\t\t\tSeeding failed -- restarting...\n";
		return BAD_SEED;
	}
	
	compute_seq_scores();
	set_seq_cutoff(phase);
	set_search_space_cutoff(phase);
	motif.set_motif_score(score());
	print_status(cerr, 0, phase);
	Motif best_motif = motif;
	
	int i, i_worse = 0;
	phase = 1;
	for(i = 1; i < 10000 && phase < 3; i++) {
		adjust_search_space();
		if(i_worse == 0)
			single_pass(false);
		else
			single_pass_select(false);
		for(int j = 0; j < 3; j++)
			if(! motif.column_sample()) break;
		compute_seq_scores_minimal();
		motif.set_motif_score(score());
		print_status(cerr, i, phase);
		if(size() > ngenes/3) {
			cerr << "\t\t\tToo many sites! Restarting...\n";
			return TOO_MANY_SITES;
		}
		if(size() < 2) {
			cerr << "\t\t\tZero or one sites, reloading best motif...\n";
			phase++;
			motif = best_motif;
			select_sites = best_motif;
			compute_seq_scores();
			set_seq_cutoff(phase);
			set_search_space_cutoff(phase);
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
			if(i_worse > params.minpass * phase) {
				if(size() < 2) {
					cerr << "\t\t\tLess than 2 genes at bad move threshold! Restarting...\n";
					return TOO_FEW_SITES;
				}
				cerr << "\t\t\tReached bad move threshold, reloading best motif...\n";
				phase++;
				motif = best_motif;
				select_sites = best_motif;
				if(phase < 3) {
					compute_seq_scores_minimal();
					set_seq_cutoff(phase);
					set_search_space_cutoff(phase);
				}
			}
		}
	}
	
	cerr << "\t\t\tRunning final greedy pass...";
	single_pass(true);
	cerr << "done.\n";
	motif.orient();
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
	cerr << "\t\t\tWrote motif with score " << motif.get_motif_score() << " to " << motfilename << '\n';
	return 0;
}

void MotifSearchScore::print_status(ostream& out, const int i, const int phase) {
	out << setw(5) << i;
	out << setw(3) << phase;
	int prec = cerr.precision(2);
	out << setw(10) << setprecision(3) << motif.get_seq_cutoff();
	out << setw(10) << motif.get_score_cutoff();
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
