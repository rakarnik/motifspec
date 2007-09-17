#include "im.h"

int main(int argc, char *argv[]) {
	string seqfile;                       // file with sequences
	string exprfile;                      // file with expression data
	if(argc < 4) {
    print_usage(cout);
    exit(0);
  }
	
	if ((! GetArg2(argc, argv, "-s", seqfile)) || (! GetArg2(argc, argv, "-e", exprfile))) {
		print_usage(cout);
    exit(0);
	}
	
	// Read parameters
	if(! GetArg2(argc, argv, "-k", k)) k = 10;
	int rounds;                           // number of rounds
	if(! GetArg2(argc, argv, "-rounds", rounds)) rounds = 50;
	int nc;                               // number of columns in motif
	if(! GetArg2(argc,argv,"-numcols", nc)) nc = 10;
	float mthresh;                        // how close clusters are allowed to get
	if(! GetArg2(argc, argv, "-mthresh", mthresh)) mthresh = 0.9999;
	float cthresh;                        // how far genes can be from the cluster mean
	if(! GetArg2(argc, argv, "-cthresh", cthresh)) cthresh = 0.65;
	int minsize;                          // minimum size for a cluster
	if(! GetArg2(argc, argv, "-minsize", minsize)) minsize = 10;
	string bootstrap;
	if(! GetArg2(argc, argv, "-bootstrap", bootstrap)) bootstrap = "no";
	int maxmots;
	if(! GetArg2(argc, argv, "-maxmots", maxmots)) maxmots = 10;
	
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
	
	gene_cluster = new int[ngenes];
	for (int g = 0; g < ngenes; g++) {
		gene_cluster[g] = -1;                             // -1 means not assigned to a cluster
	}
	
	cerr << "Successfully read input files -- dataset size is " << ngenes << " genes X " << npoints << " timepoints" << endl;
	
	
	cerr << "Setting up kmeans clusters... ";
	clusters = new Cluster[k];
	for (int c = 0; c < k; c++) {
		clusters[c].init(expr, npoints, nameset2);
	}
	cerr << "done.\n";
	
	cerr << "Assigning cluster centers to random genes... ";
	srand(100);
	for (int c = 0; c < k; c++) {
		recenter_cluster_random(c);
	}
	cerr << "done" << endl;
	
	// First run a complete round of k-means clustering
	int nchanges = 1;
	int toosmall;
	int assigned;
	for (int i = 0; (i < rounds && nchanges > 0) || toosmall > 0; i++) {
		cerr << "Running iteration " << i + 1 << endl;
		nchanges = assign_genes_to_nearest_cluster(cthresh);
		cerr << "\t" << nchanges << " changes in this iteration" << endl;
		
		sort_clusters();
		
		assigned = 0;
		cerr << "\t";
		for (int j = 0; j < k; j++) {
			cerr << clusters[j].size() << " ";
			assigned += clusters[j].size();
		}
		cerr << "(" << ngenes - assigned << " unassigned)" << endl;
		
		check_cluster_distances(mthresh);
		
		toosmall = 0;
		cerr << "\tChecking for small clusters... " << endl; 
		toosmall = check_cluster_sizes(minsize);
		cerr << "\tdone!" << endl;
		cerr << "\tThere were " << toosmall << " small clusters" << endl;
		
		if ((i+1) % 5 == 0) {
			int reset = 0;
			cerr << "\tRecentering small clusters... " << endl;
			reset = reset_small_clusters(minsize);
			cerr << "\tdone!" << endl;
			cerr << "\t" << reset << " small clusters were recentered." << endl;
			assigned = 0;
			cerr << "\t";
			for (int j = 0; j < k; j++) {
				cerr << clusters[j].size() << " ";
				assigned += clusters[j].size();
			}
			cerr << "(" << ngenes - assigned << " unassigned)" << endl;
		}
		
		if(i >= rounds && toosmall > 0) {
			cerr << "\tRemoving cluster " << k << endl;
			int size = clusters[k - 1].size();
			int* genes = new int[size];
			clusters[k - 1].genes(genes);
			for(int g = 0; g < size; g++) {
				gene_cluster[genes[g]] = -1;
			}
			k--;
			delete[] genes;
		}
		
		if (outfile != "") {
			ofstream out(outfile.c_str());
			print_clusters(out, nameset2);
			out.close();
		}
	}
	
	cerr << "Setting up adjustment models... " << endl;
	aces = new AlignACE[k];
	string line;
	for (int c = 0; c < k; c++) {
		aces[c].init(seqs, nc);
		aces[c].modify_params(argc, argv);
		aces[c].set_final_params();
		aces[c].ace_initialize();
	}
	cerr << "done." << endl;

	cerr << "Adjusting clusters using sequence information... " << endl;
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	pthread_t tids[k];
	int m[k];
	int rc, status;
	for(int c = 0; c < k; c++) {
		m[c] = c;
		rc = pthread_create(&tids[c], NULL, doit, (void*)(&m[c]));
		if(rc) {
			cerr << endl << "Error occurred during pthread_create() for cluster " << c + 1 << ". Return code was " << rc << endl;
			exit(1);
		}
	}
	pthread_attr_destroy(&attr);
	for(int c = 0; c < k; c++) {
		rc = pthread_join(tids[c], (void **)&status);
		if(rc) {
			cerr << "Error occurred during pthread_join() for cluster " << c + 1 << ". Return code was " << rc << endl;
			exit(1);
		}
	}
	cerr << "done." << endl;
	
	if (outfile == "") {
		print_clusters(cout, nameset2);
	} else {
		string outfile2 = outfile + ".adj";
		ofstream outf(outfile2.c_str());
		print_clusters(outf, nameset2);
	}
	
	for(int c = 0; c < k; c++) {
		stringstream outstream;
		outstream << c + 1 << ".adj.ace";
		string outstr;
		outstream >> outstr;
		ofstream out(outstr.c_str());
		print_ace(out, aces[c], nameset1);
	}
	
	delete [] clusters;
	delete [] aces;

	pthread_exit(NULL);
}

void* doit(void* m) {
	int c = *((int*) m); 
	cerr << "Optimizing cluster " << c + 1 << endl;
	
	double corr_cutoff[] = {0.65, 0.60, 0.50};
  double sc, cmp, sc_best_i;
  int i_worse;
  Sites best_sites = aces[c].ace_sites;
	
	for(int j = 1; j <= aces[c].ace_params.ap_nruns; j++) {
		// cerr << "\t\tSearch restart #" << j << "/" << a.ace_params.ap_nruns << endl;
		// Create a copy of the cluster which we can modify
		Cluster c1 = clusters[c];
		c1.calc_mean();
		
		sc_best_i = aces[c].ace_map_cutoff;
    i_worse = 0;
    int phase = 0;
    
		sync_ace_members(c1, aces[c]);
		aces[c].seed_random_sites_restricted(1);
    aces[c].ace_select_sites.clear_sites();
		
		for(int i = 1; i <= aces[c].ace_params.ap_npass; i++){
			if(phase == 3) {
				double sc1 = aces[c].map_score();
				sc = 0.0;
				for(int z = 0; sc < sc1 && z < 5; z++){
					aces[c].optimize_columns();
					aces[c].optimize_sites();
					sc = aces[c].map_score();
				}
				if(sc < sc1) {
					aces[c].ace_sites = best_sites;
					sc = sc1;
				}
				aces[c].ace_archive.consider_motif(aces[c].ace_sites, sc);
				// cerr << "\t\t\tReached phase 3! Restarting..." << endl;
				break;
      }
      if(i_worse == 0)
				aces[c].single_pass_restricted(aces[c].ace_params.ap_sitecut[phase]);
      else 
				aces[c].single_pass_select(aces[c].ace_params.ap_sitecut[phase]);
      if(aces[c].ace_sites.number() == 0) {
				if(sc_best_i == aces[c].ace_map_cutoff) {
					// cerr << "\t\t\tNo sites and best score matched cutoff! Restarting..." << endl;
					break;
				}
				//if(best_sites.number()<4) break;
				aces[c].ace_sites=best_sites;
				aces[c].ace_select_sites=best_sites;
				phase++;
				i_worse = 0;
				continue;
      }
      if(i<=3) continue;
      if(aces[c].column_sample(0)) {}
      if(aces[c].column_sample(aces[c].ace_sites.width()-1)) {}
      for(int m = 0; m < 3; m++) {
				if(!(aces[c].column_sample())) break;
      }
      sc = aces[c].map_score();
      if(sc - sc_best_i > 1e-3){
				i_worse=0;
				cmp = aces[c].ace_archive.check_motif(aces[c].ace_sites, sc);
				if(cmp > aces[c].ace_sim_cutoff) {
					// cerr <<"\t\t\tToo similar! Restarting..." << endl;
					break;
				}
				sc_best_i = sc;
				best_sites = aces[c].ace_sites;
      }
      else i_worse++;
      if(i_worse > aces[c].ace_params.ap_minpass[phase]){
				if(sc_best_i == aces[c].ace_map_cutoff) {
					// cerr << "\t\t\ti_worse is greater than cutoff and best score matched cutoff! Restarting..." << endl;
					break;
				}
				if(best_sites.number() < 2) {
					// cerr << "\t\t\ti_worse is greater than cutoff and best score matched cutoff! Restarting..." << endl;
					break;
				}
				aces[c].ace_sites = best_sites;
				aces[c].ace_select_sites = best_sites;
				phase++;
				i_worse = 0;
      }
			/*
			if((i == 1 || i % 50 == 0) && aces[c].ace_sites.number() > 0) {
				cerr << "\t\t\t" << setw(5) << i;
				cerr << setw(3) << phase; 
				cerr << setw(5) << aces[c].ace_sites.number();
				cerr << setw(5) << aces[c].poss_count;
				cerr << setw(40) << aces[c].consensus();
				cerr << setw(10) << sc;
				cerr << endl;
			}
			*/
			if(phase > 0 and i % 50 == 0) {
				sync_cluster(c1, aces[c]);
				sync_ace_neighborhood(c1, aces[c], corr_cutoff[phase]);
			}
		}
	}
	cerr << "Done with cluster " << c + 1 << endl;
	
	return m;
}

int assign_genes_to_nearest_cluster (float threshold) {
	float corr, max_corr;
	int max_clus, nchanges = 0;
	for (int g = 0; g < ngenes; g++) {
		max_clus = -1;
		max_corr = 0;
		for(int c = 0; c < k; c++) {
			corr = clusters[c].corr(expr[g]);
			// cerr << "\tCorrelation of gene " << (g + 1) << " with cluster " << (c + 1) << " is " << corr << endl;
			if (corr > max_corr && corr > threshold) {
				max_clus = c;
				max_corr = corr;
			}
		}
		if (max_clus != gene_cluster[g]) {
			nchanges++;
			reassign_gene(g, max_clus);
		}
	}
	update_cluster_means();
	return nchanges;
}

int check_cluster_distances(float threshold) {
	float c;
	int big, small, count;
	for (int i = 0; i < k; i++) {
		for (int j = i+1; j < k; j++) {
			c = clusters[i].corr(clusters[j].get_mean());
			//cerr << "Cluster " << (i + 1) << " and cluster " << (j + 1) << " have correlation " << c << endl;
			if (c > threshold) {
				if (clusters[i].size() > clusters[j].size()) {
					big = i;
					small = j;
				} else {
					big = j;
					small = i;
				}
				cerr << "\tMerging cluster " << (small + 1) << " into cluster " << (big + 1) << "(correlation was " << c << ")" << endl;
				merge_clusters(big, small);
				recenter_cluster(small, find_outlier());
				count++;
			}
		}
	}
	update_cluster_means();
	return count;
}

int check_cluster_sizes(int minsize) {
	int smallcount = 0;
	for (int i = 0; i < k; i++) {
		if (clusters[i].size() < minsize) {
			smallcount++;
		}
	}
	return smallcount;
}

int reset_small_clusters(int minsize) {
	int delcount = 0;
	for (int i = 0; i < k; i++) {
		if (clusters[i].size() < minsize) {
			cerr << "\tResetting cluster " << i + 1 << "(" << clusters[i].size() << " < " << minsize <<")" << endl;
			recenter_cluster_random(i);
			delcount++;
		}
	}
	return delcount;
}

void sort_clusters() {
	for (int i = 0; i < k; i++) {
		for (int j = i + 1; j < k; j++) {
			if (clusters[i].size() < clusters[j].size()) {
				swap(i, j);
			}
		}
	}
}

void swap(int c1, int c2) {
	int size1 = clusters[c1].size();
	int* genes1 = new int[size1];
	clusters[c1].genes(genes1);
	for (int i = 0; i < size1; i++) {
		gene_cluster[genes1[i]] = c2;
	}
	int size2 = clusters[c2].size();
	int* genes2 = new int[size2];
	clusters[c2].genes(genes2);
	for (int i = 0; i < size2; i++) {
		gene_cluster[genes2[i]] = c1;
	}
	delete [] genes1;
	delete [] genes2;
	swap(clusters[c1], clusters[c2]);
}

void merge_clusters(int c1, int c2) {
	int size2 = clusters[c2].size();
	int* genes2 = new int[size2];
	clusters[c2].genes(genes2);
	for (int i = 0; i < size2; i++) {
		reassign_gene(genes2[i], c1);
	}
	clusters[c1].calc_mean();
	clusters[c2].calc_mean();
}

void reassign_gene (int g, int c) {
	if(gene_cluster[g] == c) return;
	/*
	if(gene_cluster[g] == -1) {
		cerr << "\tMoving unassigned gene " << (g + 1) << " to cluster " << (c + 1) << endl;
	} else if(c == -1) {
		cerr << "\tMoving gene " << (g + 1) << " from cluster " << (gene_cluster[g] + 1) << " to unassigned" << endl;
	} else { 
		cerr << "\tMoving gene " << (g + 1) << " from cluster " << (gene_cluster[g] + 1) << " to cluster " << (c + 1) << endl;
	}
	*/
	if (c != -1) {
		clusters[c].add_gene(g);
	}
	if (gene_cluster[g] != -1) {
		clusters[gene_cluster[g]].remove_gene(g);
	}
	gene_cluster[g] = c;
}

void update_cluster_means() {
	for (int c = 0; c < k; c++)
		clusters[c].calc_mean();
}

void recenter_cluster(int c, int g) {
	int size = clusters[c].size();
	int* genes = new int[size];
	clusters[c].genes(genes);
	for (int i = 0; i < size; i++) {
		clusters[c].remove_gene(genes[i]);
		gene_cluster[genes[i]] = -1;
	}
	reassign_gene(g, c);
	clusters[c].calc_mean();
	assert(clusters[c].corr(expr[g]) > 0.95);
}

void recenter_cluster_random(int c) {
	int new_center = (int) ((float) rand()/RAND_MAX * ngenes);
	recenter_cluster(c, new_center);
}

int find_outlier() {
	float corr, min_corr = 1.0;
	int outlier;
	for (int g = 0; g < ngenes; g++) {
		for(int c = 0; c < k; c++) {
			corr = clusters[c].corr(expr[g]);
			if (corr < min_corr) {
				outlier = g;
				min_corr = corr;
			}
		}
	}
	return outlier;
}

void print_clusters(ostream& out, const vector<string>& nameset) {
	for(int c = 0; c < k; c++) {
		int size = clusters[c].size();
		int* genes = new int[size];
		clusters[c].genes(genes);
		out << "Cluster " << c + 1 << ", " << size << " orfs" << endl;
		for(int j = 0; j < size; j++) {
			out << nameset[genes[j]] << endl;
		}
		delete [] genes;
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

void debug_check_membership() {
	cerr << "\tChecking cluster membership arrays" << endl;
	for(int c = 0; c < k; c++) {
		cerr << "\t\tChecking cluster " << c + 1 << " ...";
		for (int g = 0; g < ngenes; g++) {
			if(gene_cluster[g] == c) {
				assert(clusters[c].is_member(g));
			} else {
				assert(! clusters[c].is_member(g));
			}
		}
		cerr << "done." << endl;
	}
}

void sync_ace_members(Cluster& c, AlignACE& a) {
	int size = c.size();
	int* genes = new int[size];
	c.genes(genes);
	a.set_possible(genes, size);
}

void sync_ace_neighborhood(Cluster& c, AlignACE& a, double mincorr) {
	vector<int> possibles;
	for(int g = 0; g < ngenes; g++) {
		if(c.corr(expr[g]) > mincorr) {
			possibles.push_back(g);
		}
	}
	if(possibles.size() > 0)
		a.set_possible(&possibles[0], possibles.size());
}

void sync_cluster(Cluster&c, AlignACE& a) {
	c.remove_all_genes();
	for(int s = 0; s < a.ace_sites.number(); s++) {
		c.add_gene(a.ace_sites.sites_chrom[s]);
	}
	c.calc_mean();
}
