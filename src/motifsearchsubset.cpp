#include "motifsearchsubset.h"

MotifSearchSubset::MotifSearchSubset(const vector<string>& names,
																 const vector<string>& seqs,
																 const int nc,
																 const int order,
																 const double sim_cut,
																 const vector<string>& sub) :
MotifSearch(names, seqs, nc, order, sim_cut),
subset(sub) {
	reset_search_space();
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
	motif.set_motif_score(score());
	print_status(cerr, 0, phase);
	Motif best_motif = motif;
	
	int i, i_worse = 0;
	phase = 1;
	for(i = 1; i < 10000 && phase < 3; i++) {
		if(i_worse == 0)
			single_pass(false);
		else
			single_pass_select(false);
		for(int j = 0; j < 3; j++)
			if(! column_sample()) break;
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
			best_motif.set_seq_cutoff(params.minprob[phase]);
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

void MotifSearchSubset::print_status(ostream& out, const int i, const int phase) {
	out << setw(5) << i;
	out << setw(3) << phase;
	int prec = cerr.precision(2);
	out << setw(10) << setprecision(3) << motif.get_seq_cutoff();
	out << setw(10) << "N/A";
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
