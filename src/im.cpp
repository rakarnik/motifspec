#include "im.h"

int main(int argc, char *argv[]) {
	set_new_handler(out_of_memory);

	string seqfile;                       // file with sequences
	string exprfile;                      // file with expression data
	string clusfile;
	if(argc < 6) {
    print_usage(cout);
    exit(0);
  }
	
	if ((! GetArg2(argc, argv, "-s", seqfile)) || (! GetArg2(argc, argv, "-e", exprfile)) || (! GetArg2(argc, argv, "-c", clusfile))) {
		print_usage(cout);
    exit(0);
	}
	
	// Read parameters
	if(! GetArg2(argc, argv, "-k", k)) k = 1;
	int nc;                               // number of columns in motif
	if(! GetArg2(argc,argv,"-numcols", nc)) nc = 10;
	
	string outfile = "";
	GetArg2(argc, argv, "-o", outfile);
	
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
	/*
	vector<string> clusnames;
	cerr << "Reading cluster data from '" << clusfile << "'... ";
  get_cluster(clusfile.c_str(), k, clusnames);
	cerr << "done.\n";
	
	int nclusgenes = clusnames.size();
	cerr << nclusgenes << " genes read for cluster " << k << endl;
	
	cerr << "Checking cluster genes..." << endl;
	vector<string>::iterator result;
	int notfound = 0;
	for (int g = 0; g < nclusgenes; g++) {
		if(find(nameset1.begin(), nameset1.end(), clusnames[g]) == nameset1.end()) {
			cerr << "\t\tGene '" << clusnames[g] << "' in cluster not found in dataset" << endl;
			notfound++;
		}
	}
	cerr << "done." << endl;
	if(notfound > 0) cerr << notfound << " genes were not found in the dataset" << endl;
	
	nclusgenes -= notfound;
	int clusgenes[nclusgenes];
	int index = 0;
	for (int g = 0; g < ngenes; g++) {
		if(find(clusnames.begin(), clusnames.end(), nameset1[g]) != clusnames.end()) {
			clusgenes[index] = g;
			index++;
		}
	}
	*/
	
	cerr << "Setting up adjustment cluster... ";
	Cluster c;
	c.init(expr, npoints, nameset2);
	/*
	c.add_genes(clusgenes, nclusgenes); 
	c.calc_mean();
	cerr << "done." << endl;
	*/
	
	cerr << "Setting up AlignACE... ";
	AlignACE a;
	a.init(seqs, nc);
	a.modify_params(argc, argv);
	a.set_final_params();
	a.ace_initialize();
	cerr << "done." << endl;
	
	cerr << "Adjusting clusters using sequence information... " << endl;
	doit(c, a, nameset1);
	cerr << "done." << endl;
	
	string outstr(outfile);
	if(outfile != "") {
		outstr.append(".adj.ace");
		ofstream out(outstr.c_str(), ios::trunc);
		print_ace(out, a, nameset1);
		out.close();
	} else {
		print_ace(cerr, a, nameset1);
	}
}

void doit(Cluster& c, AlignACE& a, vector<string>& nameset) {
	double corr_cutoff[] = {0.65, 0.60, 0.50, 0.50};
  double sc, cmp, sc_best_i;
  int i_worse;
  Sites best_sites = a.ace_sites;
	
	// Use the final search neighborhood to decide the number of runs
	// sync_ace_neighborhood(c, a, corr_cutoff[3]);
	int nruns = a.ace_sites.positions_available()
	            / a.ace_params.ap_expect
							/ a.ace_sites.ncols()
							/ a.ace_params.ap_undersample
							* a.ace_params.ap_oversample;
	
	// Reset the search area strictly to cluster members
	// sync_ace_members(c, a);
	
	for(int j = 1; j <= nruns; j++) {
		cerr << "\t\tSearch restart #" << j << "/" << nruns << endl;
		// Create a copy of the cluster which we can modify
		//Cluster* c1 = new Cluster(c);
		Cluster* c1 = new Cluster();
		c1->init(expr, npoints, nameset);
		
		sc_best_i = a.ace_map_cutoff;
    i_worse = 0;
    int phase = 0;
    int old_phase = 0;
		
		a.ace_sites.clear_sites();
		a.ace_select_sites.clear_sites();
		for(int g = 0; g < ngenes; g++)
			a.add_possible(g);
		a.seed_random_sites(1);
		sync_cluster(*c1, a);
		sync_ace_neighborhood(*c1, a, corr_cutoff[phase]);
		
		for(int i = 1; i <= a.ace_params.ap_npass; i++){
			if(old_phase < phase) {
				sync_cluster(*c1, a);
				sync_ace_neighborhood(*c1, a, corr_cutoff[phase]);
				print_ace_status(cerr, a, i, phase, sc);
				old_phase = phase;
			}
			if(phase == 4) {
				if (a.ace_sites.number() < 5) {
					cerr << "\t\t\tReached phase 5 and not enough sites! Restarting..." << endl;
					break;
				}
				double sc1 = a.map_score();
				sc = 0.0;
				for(int z = 0; sc < sc1 && z < 5; z++){
					a.optimize_columns();
					a.optimize_sites();
					sc = a.map_score();
				}
				if(sc < sc1) {
					a.ace_sites = best_sites;
					sc = sc1;
				}
				a.ace_archive.consider_motif(a.ace_sites, sc);
				print_ace_status(cerr, a, i, phase, sc);
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
			
			if(i == 1 || i % 50 == 0) print_ace_status(cerr, a, i, phase, sc);
		}
		
		delete c1;
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

void sync_ace_members(const Cluster& c, AlignACE& a) {
	// cerr << "\t\t\t\tSyncing AlignACE with cluster members... ";
	int size = c.size();
	int* genes = new int[size];
	c.genes(genes);
	a.clear_possible();
	for(int g = 0; g < size; g++)
		a.add_possible(genes[g]);
	// cerr << "done." << endl;
}

void sync_ace_neighborhood(const Cluster& c, AlignACE& a, double mincorr) {
	// cerr << "\t\t\t\tSyncing AlignACE with cluster neighborhood... ";
	int count = 0;
	a.clear_possible();
	for(int g = 0; g < ngenes; g++) {
		if(c.corr(expr[g]) > mincorr) {
			a.add_possible(g);
			count++;
		}
	}
	// cerr << "\t\t\t\tSearch neighborhood with corr > " << mincorr << " contains " << count << " genes" << endl; 
	// cerr << "done." << endl;
}

void sync_cluster(Cluster&c, const AlignACE& a) {
	c.remove_all_genes();
	for(int s = 0; s < a.ace_sites.number(); s++) {
		c.add_gene(a.ace_sites.sites_chrom[s]);
	}
	c.calc_mean();
}
