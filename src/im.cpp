#include "im.h"

int main(int argc, char *argv[]) {
	set_new_handler(out_of_memory);

	string seqfile;                       // file with sequences
	string exprfile;                      // file with expression data
	string clusfile;
	string outfile;
	if(argc < 6) {
    print_usage(cout);
    exit(0);
  }
	
	if ((! GetArg2(argc, argv, "-s", seqfile)) 
				|| (! GetArg2(argc, argv, "-e", exprfile)) 
				|| (! GetArg2(argc, argv, "-c", clusfile))
				|| (! GetArg2(argc, argv, "-o", outfile))) {
		print_usage(cout);
    exit(0);
	}
	
	// Read parameters
	if(! GetArg2(argc, argv, "-k", k)) k = 1;
	int nc;                               // number of columns in motif
	if(! GetArg2(argc,argv,"-numcols", nc)) nc = 10;
	
  vector<string> seqs, nameset1;
  cerr << "Reading sequence data from '" << seqfile << "'... ";
	get_fasta_fast(seqfile.c_str(), seqs, nameset1);
	cerr << "done.\n";
	
	vector<string> nameset2;
	cerr << "Reading expression data from '" << exprfile << "'... ";
  expr = get_expr(exprfile.c_str(), &npoints, nameset2);
	cerr << "done.\n";
	
	if(nameset1.size() == nameset2.size()) {
		ngenes = nameset1.size();
	} else {
		cerr << "Inconsistent sizes for sequence and expression datasets\n";
		exit(0);
	}
	
	bool match = true;
	for (int i = 0; i < ngenes; i++) {
		match = match && (nameset1[i] == nameset2[i]);
		if (nameset1[i] != nameset2[i]) {
			cerr << "Row " << i << ": Sequence: '" << nameset1[i] << "'  " << "Expression: '" << nameset2[i] << "'" << endl;
		}
	}
	if(! match) {
		cerr << "Inconsistent gene names for sequence and expression datasets\n";
		exit(0);
	}
	
	cerr << "Successfully read input files -- dataset size is " << ngenes << " genes X " << npoints << " timepoints" << endl;

	
	jcorr = new float*[ngenes];
	for(int i = 0; i < ngenes; i++) {
		jcorr[i] = new float[ngenes];
		for(int j = 0; j < ngenes; j++) {
			jcorr[i][j] = - 2;
		}
	}
	
	float jc = 0;
	jcorr = new float*[ngenes];
	for(int i = 0; i < ngenes; i++) {
		jcorr[i] = new float[ngenes];
	}
	/*
	ofstream jcorrout("jcorr.out");
	for(int i = 0; i < ngenes; i++) {
		cerr << "Computing correlation for gene " << i << endl;
		for(int j = i + 1; j < ngenes; j++) {
			jc = jack_corr(expr[i], expr[j], npoints);
			jcorr[i][j] = jcorr[j][i] = jc;
			jcorrout << i << '\t' << j << '\t' << jcorr[i][j] << endl;
			jcorrout << j << '\t' << i << '\t' << jcorr[i][j] << endl;
		}
	}
	*/
	
	ifstream corrin("jcorr.out");
	int corrcount = 0;
	string line;
	vector<string> values;
	while(getline(corrin, line)) {
		values = split(line, '\t');
		int i = atoi(values[0].c_str());
		int j = atoi(values[1].c_str());
		float jc = atof(values[2].c_str());
		jcorr[i][j] = jcorr[j][i] = jc;
		corrcount++;
	}
	cerr << "Read " << corrcount << " correlation values" << endl;
	
	cerr << "Setting up AlignACE... ";
	AlignACE a;
	a.init(seqs, nc);
	a.modify_params(argc, argv);
	a.set_final_params();
	a.ace_initialize();
	cerr << "done." << endl;
	
	string tmpstr(outfile);
	tmpstr.append(".tmp.ace");
	doit(tmpstr.c_str(), a, nameset1);
	
	string outstr(outfile);
	outstr.append(".adj.ace");
	ofstream out(outstr.c_str(), ios::trunc);
	print_ace(out, a, nameset1);
	out.close();
}

void doit(const char* outfile, AlignACE& a, vector<string>& nameset) {
	double corr_cutoff[] = {0.65, 0.60, 0.50, 0.50};
  double sc, cmp, sc_best_i;
  int i_worse;
  Sites best_sites = a.ace_sites;
	
	int nruns = a.ace_sites.positions_available()
	            / a.ace_params.ap_expect
							/ a.ace_sites.ncols()
							/ a.ace_params.ap_undersample
							* a.ace_params.ap_oversample;
	
	for(int j = 1; j <= nruns; j++) {
		cerr << "\t\tSearch restart #" << j << "/" << nruns << endl;
		
		sc_best_i = a.ace_map_cutoff;
    i_worse = 0;
    int phase = 0;
    int old_phase = 0;
		
		a.ace_sites.clear_sites();
		a.ace_select_sites.clear_sites();
		for(int g = 0; g < ngenes; g++)
			a.add_possible(g);
		a.seed_random_site();
		cerr << "Just seeded: ace_members is " << a.ace_members << endl;
		expand_ace_search(a, corr_cutoff[phase]);
		if(a.ace_members < 2) {
			cerr << "\t\t\tSeed was too far from other genes!" << endl;
			continue;
		}
		
		for(int i = 1; i <= a.ace_params.ap_npass; i++){
			sc = a.map_score();
			print_ace_status(cerr, a, i, old_phase, sc);
			if(old_phase < phase) {
				expand_ace_search(a, corr_cutoff[phase]);
				old_phase = phase;
			}
			if(phase == 3) {
				if(a.ace_sites.sites_num < 5) {
					cerr << "\t\t\tReached phase 3 with less than 5 sites. Restarting..." << endl;
					break;
				}
				double sc1 = a.map_score();
				sc = 0.0;
				for(int z = 0; sc < sc1 && z < 5; z++){
					a.optimize_columns();
					a.optimize_sites();
					sc = a.map_score();
					print_ace_status(cerr, a, i, phase, sc);
				}
				if(sc < sc1) {
					a.ace_sites = best_sites;
					sc = sc1;
				}
				sc = a.map_score();
				print_ace_status(cerr, a, i, phase, sc);
				a.ace_archive.consider_motif(a.ace_sites, sc);
				cerr << "\t\t\tCompleted phase 3! Restarting..." << endl;
				break;
      }
      if(i_worse == 0)
				a.single_pass(a.ace_params.ap_sitecut[phase]);
      else 
				a.single_pass_select(a.ace_params.ap_sitecut[phase]);
      if(a.ace_sites.number() == 0) {
				if(sc_best_i == a.ace_map_cutoff) {
					cerr << "\t\t\tNo sites and best score matched cutoff! Restarting..." << endl;
					break;
				}
				//if(best_sites.number()<4) break;
				a.ace_sites=best_sites;
				a.ace_select_sites=best_sites;
				phase++;
				i_worse = 0;
				continue;
      }
      if(i<=3) continue;
			if(phase < 2) {
				if(a.column_sample(0)) {}
				if(a.column_sample(a.ace_sites.width()-1)) {}
				for(int m = 0; m < 3; m++) {
					if(!(a.column_sample())) break;
				}
			}
      sc = a.map_score();
      if(sc - sc_best_i > 1e-3){
				i_worse=0;
				cmp = a.ace_archive.check_motif(a.ace_sites, sc);
				if(cmp > a.ace_sim_cutoff) {
					print_ace_status(cerr, a, i, phase, sc);
					cerr <<"\t\t\tToo similar! Restarting..." << endl;
					break;
				}
				sc_best_i = sc;
				best_sites = a.ace_sites;
      }
      else i_worse++;
      if(i_worse > a.ace_params.ap_minpass[phase]){
				if(sc_best_i == a.ace_map_cutoff) {
					print_ace_status(cerr, a, i, phase, sc);
					cerr << "\t\t\ti_worse is greater than cutoff and best score at cutoff! Restarting..." << endl;
					break;
				}
				if(best_sites.number() < 2) {
					print_ace_status(cerr, a, i, phase, sc);
					cerr << "\t\t\ti_worse is greater than cutoff and only 1 site! Restarting..." << endl;
					break;
				}
				a.ace_sites = best_sites;
				a.ace_select_sites = best_sites;
				phase++;
				i_worse = 0;
      }
			
		
			print_ace_status(cerr, a, i, phase, sc);
			expand_ace_search(a, corr_cutoff[phase]);
		}
		
		if(j % 100 == 0) {
			ofstream out(outfile, ios::trunc);
			print_ace(out, a, nameset);
		}
	}
}

void print_full_ace(ostream& out, AlignACE& a, const vector <string>& nameset) {
	out << "Parameter values:\n";
  a.output_params(out);
  out << "\nInput sequences:\n";
  for(int x = 0; x < nameset.size(); x++) out << "#" << x << '\t' << nameset[x] << endl;
  out << '\n';
  a.full_output(out);
}

void print_ace(ostream& out, AlignACE& a, const vector <string>& nameset) {
	out << "Parameter values:\n";
  a.output_params(out);
  out << "\nInput sequences:\n";
  for(int x = 0; x < nameset.size(); x++) out << "#" << x << '\t' << nameset[x] << endl;
  out << endl;
	a.full_output(out);
}

void print_ace_status(ostream& out, AlignACE& a, const int i, const int phase, const double sc) {
	if(a.ace_sites.number() > 0) {
		out << "\t\t\t" << setw(5) << i;
		out << setw(3) << phase; 
		out << setw(5) << a.ace_sites.number();
		out << setw(5) << a.ace_members;
		out << setw(40) << a.consensus();
		out << setw(10) << sc;
		out << endl;
	} else {
		out << "\t\t\tNo sites!" << endl;
	}
}

void print_usage(ostream& fout) {
	fout<<"Usage: im -s seqfile -e exprfile (options)\n";
  fout<<" Seqfile must be in FASTA format.\n";
	fout<<" Exprfile must be in tab-delimited format.\n";
  fout<<"Options:\n";
	fout<<" -k          \tnumber of clusters (10)\n";
  fout<<" -numcols    \tnumber of columns to align (10)\n";
  fout<<" -expect     \tnumber of sites expected in model (10)\n";
  fout<<" -gcback     \tbackground fractional GC content of input sequence (0.38)\n";
  fout<<" -minpass    \tminimum number of non-improved passes in phase 1 (200)\n";
  fout<<" -seed       \tset seed for random number generator (time)\n";
  fout<<" -undersample\tpossible sites / (expect * numcols * seedings) (1)\n"; 
  fout<<" -oversample\t1/undersample (1)\n";
}

void expand_ace_search(AlignACE& a, double mincorr) {
	if(a.ace_members == 1) {
		for(int g1 = 0; g1 < ngenes; g1++) {
			if(! a.is_possible(g1)) continue;
			cerr << "\t\t\tJust seeded(" << a.ace_members << "): adding all genes within 0.7 of gene " << g1 << endl;
			for(int g2 = 0; g2 < ngenes; g2++)
				if(jcorr_lookup(g1, g2) > 0.7)
					a.add_possible(g2);
			break;
		}
	} else {
		cerr << "\t\t\tMore than one gene in search space: adding closest gene to those already in space..." << endl;
		for(int g1 = 0; g1 < ngenes; g1++) {
			if(! a.ace_sites.sites_has_sites[g1]) continue;
			int closest = -1;
			float maxjc = -1;
			float jc;
			cerr << "\t\t\t\tLooking for gene closest to gene " << g1 << endl;
			for(int g2 = 0; g2 < ngenes; g2++) {
				if(a.is_possible(g2)) continue;
				jc = jcorr_lookup(g1, g2);
				if(jc > maxjc) {
					closest = g2;
					maxjc = jc;
				}
			}
			if(closest != -1 && maxjc > 0.7) {
				cerr << "\t\t\t\tClosest was " << closest << endl;
				assert(! a.is_possible(closest));
				a.add_possible(closest);
			}
		}
	}
}

float jcorr_lookup(const int g1, const int g2) {
	float jc;
	if(jcorr[g1][g2] == -2) {
		jc = jack_corr(expr[g1], expr[g2], npoints);
		jcorr[g1][g2] = jcorr[g2][g1] = jc;
	}
	return jcorr[g1][g2];
}
