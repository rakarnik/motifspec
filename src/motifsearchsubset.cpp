#include "motifsearchsubset.h"

MotifSearchSubset::MotifSearchSubset(const vector<string>& names,
																 const vector<string>& seqs,
																 const int nc,
																 const int order,
																 const double sim_cut,
																 const int maxm,
																 const vector<string>& sub) :
MotifSearch(names, seqs, nc, order, sim_cut, maxm),
subset(sub) {
	reset_search_space();
}

void MotifSearchSubset::set_final_params() {
	MotifSearch::set_final_params();
	params.minprob[0] = 0.00001;
	params.minprob[1] = 0.002;
	params.minprob[2] = 0.1;
	params.minprob[3] = 0.2;
}

void MotifSearchSubset::reset_search_space() {
	motif.clear_search_space();
	for(int i = 0; i < ngenes; i++)
		if(binary_search(subset.begin(), subset.end(), nameset[i]))
			motif.add_to_search_space(i);
	assert(motif.get_search_space_size() > 0);
	assert(motif.get_search_space_size() <= ngenes);
}

void MotifSearchSubset::adjust_search_space() {
	reset_search_space();
}

void MotifSearchSubset::set_search_space_cutoff(const int phase) {
	(void) phase;
	return;
}

double MotifSearchSubset::score() {
	return -log_prob_overlap(motif.get_above_cutoffs(), motif.get_search_space_size(), motif.get_above_seqc(), ngenes);
}

int MotifSearchSubset::search_for_motif(const int worker, const int iter, const string outfile) {
	motif.clear_sites();
	select_sites.clear_sites();
	stringstream iterstr;
	iterstr << worker << '.' << iter;
	motif.set_iter(iterstr.str());
	int phase = 0;
	
	seed_random_site();
	if(size() < 1) {
		cerr << "\t\t\tSeeding failed -- restarting...\n";
		return BAD_SEED;
	}
	
	compute_seq_scores();
	set_seq_cutoff(phase);
	update_seq_count();
	motif.set_motif_score(score());
	print_status(cerr, 0, phase);
	Motif best_motif = motif;
	
	int i, i_worse = 0;
	phase = 1;
	for(i = 1; i < 10000 && phase < 3; i++) {
		if(i_worse == 0)
			single_pass();
		else
			single_pass_select();
		for(int j = 0; j < 3; j++)
			if(! motif.column_sample()) break;
		compute_seq_scores_minimal();
		motif.set_motif_score(score());
		print_status(cerr, i, phase);
		if(size() < 2) {
			cerr << "\t\t\tZero or one sites, reloading best motif...\n";
			phase++;
			best_motif.set_seq_cutoff(params.minprob[phase]);
			motif = best_motif;
			select_sites = best_motif;
			compute_seq_scores_minimal();
			update_seq_count();
			set_seq_cutoff(phase);
			update_seq_count();
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
				update_seq_count();
				set_seq_cutoff(phase);
				update_seq_count();
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
	update_seq_count();
	motif.set_motif_score(score());
	print_status(cerr, i, phase);
	
	if(size() < params.minsize) {
		cerr << "\t\t\tToo few sites! Restarting...\n";
		return TOO_FEW_SITES;
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

