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
	bsclusters = new Cluster[k];
	for (int c = 0; c < k; c++) {
		bsclusters[c].init(expr, npoints, nameset2);
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
			cerr << bsclusters[j].size() << " ";
			assigned += bsclusters[j].size();
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
				cerr << bsclusters[j].size() << " ";
				assigned += bsclusters[j].size();
			}
			cerr << "(" << ngenes - assigned << " unassigned)" << endl;
		}
		
		if(i >= rounds && toosmall > 0) {
			cerr << "\tRemoving cluster " << k << endl;
			int size = bsclusters[k - 1].size();
			int* genes = new int[size];
			bsclusters[k - 1].genes(genes);
			for(int g = 0; g < size; g++) {
				gene_cluster[genes[g]] = -1;
			}
			k--;
			delete[] genes;
		}
		
		if (outfile != "") {
			ofstream out(outfile.c_str());
			print_bsclusters(out, nameset2);
			out.close();
		}
	}
	
	// Setting up and running bootstrap AlignACE models
	cerr << "Running bootstrap AlignACE models... " << endl;
	bsaces = new AlignACE[k];
	vector<string> aceoutstr(k);
	for (int c = 0; c < k; c++) {
		stringstream aceoutstrm;
		aceoutstrm << c + 1;
		aceoutstrm << ".ace";
		aceoutstrm >> aceoutstr[c];
		cerr << "\tRunning on cluster " << c + 1 << "... " << endl;
		int size = bsclusters[c].size();
		int* genes = new int[size];
		bsclusters[c].genes(genes);
		vector<string> clus_names;
		vector<string> clus_seqs;
		for(int g = 0; g < size; g++) {
			clus_names.push_back(nameset1[genes[g]]);
			clus_seqs.push_back(seqs[genes[g]]);
		}
		bsaces[c].init(clus_seqs, nc);
		bsaces[c].modify_params(argc, argv);
		if(bootstrap == "yes") {
			bsaces[c].doit();
			ofstream aceout(aceoutstr[c].c_str());
			print_full_ace(aceout, bsaces[c], clus_names);
		}
		cerr << "\tdone." << endl;
	}
	cerr << "done." << endl;
	
	// Set up adjustment AlignACE models, seeding with sites found above
	cerr << "Setting up adjustment models... " << endl;
	nmots = new int[k];
	aces = new AlignACE*[k];
	clusters = new Cluster*[k];
	int** nsites;
	int** additions;
	int** subtractions;
	double** score;
	nsites = new int*[k];
	additions = new int*[k];
	subtractions = new int*[k];
	score = new double*[k];
	string line;
	for (int c = 0; c < k; c++) {
		int size = bsclusters[c].size();
		int* genes = new int[size];
		bsclusters[c].genes(genes);
		
		// First figure out the number of motifs in this cluster 
		ifstream acein1(aceoutstr[c].c_str());
		nmots[c] = 0;
		while(! acein1.eof()) {
			while((! acein1.eof()) && line.find("Motif") == line.npos) getline(acein1, line);
			if(! acein1.eof()) (nmots[c])++;
			line.clear();
		}
		cerr << "\t\t" << nmots[c] << " motifs read for cluster " << c + 1 << " from file " << aceoutstr[c] << endl;
		if(nmots[c] > maxmots) nmots[c] = maxmots;
		aces[c] = new AlignACE[nmots[c]];
		clusters[c] = new Cluster[nmots[c]];
		nsites[c] = new int[nmots[c]];
		additions[c] = new int[nmots[c]];
		subtractions[c] = new int[nmots[c]];
		score[c] = new double[nmots[c]];
		
		ifstream acein2(aceoutstr[c].c_str());
		for(int d = 0; d < nmots[c] && ! acein2.eof(); d++) {
			// cerr << "\t\tSetting up AlignACE model for motif " << c + 1 << "." << d + 1 << endl;
			aces[c][d].init(seqs, nc);
			aces[c][d].modify_params(argc, argv);
			aces[c][d].set_final_params();
			aces[c][d].ace_initialize();
			nsites[c][d] = 0; 
			while((! acein2.eof()) && line.find("Motif") == line.npos) getline(acein2, line); // Discard everything before "Motif..."
			
			while(! acein2.eof()) {                                                          // Read everything until "*"
				getline(acein2, line);
				if(line.find("*") != line.npos) break;
				vector<string> fields = split(line, '\t');
				aces[c][d].ace_sites.add_site(genes[str_to_int(fields[1])], str_to_int(fields[2]), str_to_int(fields[3]));
				(nsites[c][d])++;
			}
			int width = line.length();
			int prevpos = line.find("*");
			int currpos = prevpos;
			while(prevpos < width - 1) {
				currpos = line.find("*", prevpos + 1);
				aces[c][d].ace_sites.sites_active_fwd[prevpos] = currpos;
				prevpos = currpos;
			}
			aces[c][d].ace_sites.sites_active_fwd[width - 1] = width - 1;
			aces[c][d].ace_sites.sites_width = width;
			
			// cerr << "\t\tSetting up adjustment cluster for motif " << c + 1 << "." << d + 1 << endl;
			
			clusters[c][d].init(expr, npoints, nameset2);
			for(int g = 0; g < size; g++) {
				clusters[c][d].add_gene(genes[g]);
			}
			clusters[c][d].calc_mean();
			
			// cerr << "\t\tSet up adjustment model for motif " << c + 1 << "." << d + 1 << endl;
		}
	}
	cerr << "done." << endl;
	
	cerr << "Initial state of clusters:" << endl;
	for (int c = 0; c < k; c++) {
		for(int d = 0; d < nmots[c]; d++) {
			cerr << setw(3) << right << c + 1;
			cerr << ".";
			cerr << setw(4) << left << d + 1;
			cerr << setw(30) << left << aces[c][d].consensus();
			cerr << setw(6) << nsites[c][d];
			cerr << setw(6) << clusters[c][d].size();
			cerr << setw(12) << aces[c][d].map_score();
			cerr << endl;
		}
	}
	cerr << endl;
	
	delete [] bsaces;
	delete [] bsclusters;
	
	cerr << "Adjusting clusters using sequence information... " << endl;
	nchanges = 1;
	for (int r = 0; r < 50 && nchanges > 0; r++) {
		cerr << "\tRound " << r + 1 << "/" << rounds << ":" << endl;
		nchanges = 0;
		for (int c = 0; c < k; c++) {
			for(int d = 0; d < nmots[c]; d++) {
				additions[c][d] = 0;
				subtractions[c][d] = 0;
				vector<int> possibles;
				aces[c][d].calc_matrix();
				aces[c][d].ace_sites.remove_all_sites();
				// If gene is in or near this cluster, check if we have a high-scoring site
				for(int g = 0; g < ngenes; g++) {
					float corr = clusters[c][d].corr(expr[g]);
					if(corr > 0.50) {
						possibles.push_back(g);
						if(aces[c][d].consider_site(g, corr, 0.2)) {
							// cerr << "\t\t\tAdding gene " << g + 1 << " to cluster " << c + 1 << endl;
							if(! clusters[c][d].is_member(g)) {
								clusters[c][d].add_gene(g);
								nchanges++;
								additions[c][d]++;
							}
						} else {
							// cerr << "\t\t\tRemoving gene " << g + 1 << " to cluster " << c + 1 << endl;
							if(clusters[c][d].is_member(g)) {
								clusters[c][d].remove_gene(g);
								nchanges++;
								subtractions[c][d]++;
							}
						}
					}
				}
				
				score[c][d] = aces[c][d].map_score_restricted(&possibles[0], possibles.size());
			}
		}
		
		update_cluster_means();
		
		for (int c = 0; c < k; c++) {
			for(int d = 0; d < nmots[c]; d++) {
				cerr << setw(3) << right << c + 1;
				cerr << ".";
				cerr << setw(4) << left << d + 1;
				cerr << setw(30) << left << aces[c][d].consensus();
				cerr << setw(6) << right << clusters[c][d].size();
				cerr << setw(12) <<  right << score[c][d];
				cerr << setw(3) << right << "+";
				cerr << setw(3) << left << additions[c][d];
				cerr << setw(3) << right << "-";
				cerr << setw(3) << left << subtractions[c][d];
				cerr << endl;
			}
		}
		
		cerr << "\t" << nchanges << " changes in total" << endl << endl;
	}
	
	if (outfile == "") {
		print_clusters(cout, nameset2);
	} else {
		string outfile2 = outfile + ".adj";
		ofstream outf(outfile2.c_str());
		print_clusters(outf, nameset2);
	}
	
	for(int c = 0; c < k; c++) {
		for(int d = 0; d < nmots[c]; d++) {
			stringstream outstream;
			outstream << c + 1 << "." << setfill('0') << setw(2) << d + 1 << ".ace";
			string outstr;
			outstream >> outstr;
			ofstream out(outstr.c_str());
			print_ace(out, aces[c][d], score[c][d], nameset1);
		}
	}
	
	for(int c = 0; c < k; c++) {
		delete [] clusters[c];
		delete [] aces[c];
		delete [] nsites[c];
		delete [] additions[c];
		delete [] subtractions[c];
		delete [] score[c];
	}
	delete [] nmots;
	delete [] clusters;
	delete [] aces;
	delete [] nsites;
	delete [] additions;
	delete [] subtractions;
}

int assign_genes_to_nearest_cluster (float threshold) {
	float corr, max_corr;
	int max_clus, nchanges = 0;
	for (int g = 0; g < ngenes; g++) {
		max_clus = -1;
		max_corr = 0;
		for(int c = 0; c < k; c++) {
			corr = bsclusters[c].corr(expr[g]);
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
	update_bscluster_means();
	return nchanges;
}

int check_cluster_distances(float threshold) {
	float c;
	int big, small, count;
	for (int i = 0; i < k; i++) {
		for (int j = i+1; j < k; j++) {
			c = bsclusters[i].corr(bsclusters[j].get_mean());
			//cerr << "Cluster " << (i + 1) << " and cluster " << (j + 1) << " have correlation " << c << endl;
			if (c > threshold) {
				if (bsclusters[i].size() > bsclusters[j].size()) {
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
	update_bscluster_means();
	return count;
}

int check_cluster_sizes(int minsize) {
	int smallcount = 0;
	for (int i = 0; i < k; i++) {
		if (bsclusters[i].size() < minsize) {
			smallcount++;
		}
	}
	return smallcount;
}

int reset_small_clusters(int minsize) {
	int delcount = 0;
	for (int i = 0; i < k; i++) {
		if (bsclusters[i].size() < minsize) {
			cerr << "\tResetting cluster " << i + 1 << "(" << bsclusters[i].size() << " < " << minsize <<")" << endl;
			recenter_cluster_random(i);
			delcount++;
		}
	}
	return delcount;
}

void sort_clusters() {
	for (int i = 0; i < k; i++) {
		for (int j = i + 1; j < k; j++) {
			if (bsclusters[i].size() < bsclusters[j].size()) {
				swap(i, j);
			}
		}
	}
}

void swap(int c1, int c2) {
	int size1 = bsclusters[c1].size();
	int* genes1 = new int[size1];
	bsclusters[c1].genes(genes1);
	for (int i = 0; i < size1; i++) {
		gene_cluster[genes1[i]] = c2;
	}
	int size2 = bsclusters[c2].size();
	int* genes2 = new int[size2];
	bsclusters[c2].genes(genes2);
	for (int i = 0; i < size2; i++) {
		gene_cluster[genes2[i]] = c1;
	}
	delete [] genes1;
	delete [] genes2;
	swap(bsclusters[c1], bsclusters[c2]);
}

void merge_clusters(int c1, int c2) {
	int size2 = bsclusters[c2].size();
	int* genes2 = new int[size2];
	bsclusters[c2].genes(genes2);
	for (int i = 0; i < size2; i++) {
		reassign_gene(genes2[i], c1);
	}
	bsclusters[c1].calc_mean();
	bsclusters[c2].calc_mean();
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
		bsclusters[c].add_gene(g);
	}
	if (gene_cluster[g] != -1) {
		bsclusters[gene_cluster[g]].remove_gene(g);
	}
	gene_cluster[g] = c;
}

void update_bscluster_means() {
	for (int i = 0; i < k; i++) {
		bsclusters[i].calc_mean();
	}
}

void update_cluster_means() {
	for (int c = 0; c < k; c++) {
		for(int d = 0; d < nmots[c]; d++) {
			clusters[c][d].calc_mean();
		}
	}
}

void recenter_cluster(int c, int g) {
	int size = bsclusters[c].size();
	int* genes = new int[size];
	bsclusters[c].genes(genes);
	for (int i = 0; i < size; i++) {
		bsclusters[c].remove_gene(genes[i]);
		gene_cluster[genes[i]] = -1;
	}
	reassign_gene(g, c);
	bsclusters[c].calc_mean();
	assert(bsclusters[c].corr(expr[g]) > 0.95);
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
			corr = bsclusters[c].corr(expr[g]);
			if (corr < min_corr) {
				outlier = g;
				min_corr = corr;
			}
		}
	}
	return outlier;
}

void print_bsclusters(ostream& out, const vector<string>& nameset) {
	for(int c = 0; c < k; c++) {
		int size = bsclusters[c].size();
		int* genes = new int[size];
		bsclusters[c].genes(genes);
		out << "Cluster " << c + 1 << ", " << size << " orfs" << endl;
		for(int j = 0; j < size; j++) {
			out << nameset[genes[j]] << endl;
		}
		delete [] genes;
	}
}

void print_clusters(ostream& out, const vector<string>& nameset) {
	for(int c = 0; c < k; c++) {
		for(int d = 0; d < nmots[c]; d++) {
			int size = clusters[c][d].size();
			int* genes = new int[size];
			clusters[c][d].genes(genes);
			out << "Cluster " << c + 1 << "." << d + 1 << ", " << size << " orfs" << endl;
			for(int j = 0; j < size; j++) {
				out << nameset[genes[j]] << endl;
			}
			delete [] genes;
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

void print_ace(ostream& out, AlignACE& a, const double score, const vector <string>& nameset) {
	out << "Parameter values:\n";
  a.output_params(out);
  out << "\nInput sequences:\n";
  for(int x = 0; x < nameset.size(); x++) out << "#" << x << '\t' << nameset[x] << endl;
  out << endl;
	out << "Motif 1" << endl;
  a.output(out);
	out << "MAP Score: " << score << endl << endl;
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
				assert(bsclusters[c].is_member(g));
			} else {
				assert(! bsclusters[c].is_member(g));
			}
		}
		cerr << "done." << endl;
	}
}

