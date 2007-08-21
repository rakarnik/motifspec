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
	
	int nchanges = 1;
	int toosmall;
	int assigned;
	for (int i = 0; i < rounds || nchanges > 0 || toosmall > 0; i++) {
		cerr << "Running iteration " << i + 1 << endl;
		nchanges = assign_genes_to_nearest_cluster(cthresh);
		cerr << "\t" << nchanges << " changes in this iteration" << endl;
		
		cerr << "\tSorting clusters by size... ";
		sort_clusters();
		cerr << "done." << endl;
		
		assigned = 0;
		cerr << "\t";
		for (int j = 0; j < k; j++) {
			cerr << clusters[j].size() << " ";
			assigned += clusters[j].size();
		}
		cerr << "(" << ngenes - assigned << " unassigned)" << endl;
		
		
		check_cluster_distances(mthresh);
		
		toosmall = 0;
		if ((i+1) % 5 == 0) {
			cerr << "\tResetting small clusters... " << endl; 
			toosmall = check_cluster_sizes(minsize);
			cerr << "\tdone!" << endl;
			cerr << "\tThere were " << toosmall << " small clusters" << endl;
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
			vector<int> genes;
			clusters[k - 1].genes(genes);
			for(int g = 0; g < genes.size(); g++) {
				gene_cluster[genes[g]] = -1;
			}
			k--;
		}
		
		if (outfile != "") {
			ofstream out(outfile.c_str());
			print_clusters(out, nameset2);
			out.close();
		}
	}
	
	cerr << "Searching for overrepresented sequences within clusters... " << endl;
	pthread_t threads[k];
	int rc;
	struct aa_init inits[k];
	for (int c = 0; c < k; c++) {
		vector<string> clus_seqs;
		vector<int> genes;
		clusters[c].genes(genes);
		for (int g = 0; g < genes.size(); g++) {
			clus_seqs.push_back(seqs[genes[g]]);
		}
		inits[c].c = c;
		inits[c].seqset = clus_seqs;
		inits[c].nc = nc;
		inits[c].argc = argc;
		inits[c].argv = argv;
		rc = pthread_create(&threads[c], NULL, doit, &inits[c]); 
		if(rc != 0) {
			cerr << "Unable to create thread for cluster " << c + 1 << ", error code is " << rc << endl;
			exit(-1);
		}
	}
	
	int status;
	for (int c = 0; c < k; c++) {
		rc = pthread_join(threads[c], (void **) &status);
		if (rc != 0) {
			cerr << "Error while joining thread, return code from pthread_join() is " << rc << endl;
			exit(-1);
		}
		cerr << "Completed join with thread " << c << ", status is "<< status << endl;
	}
	
	cerr << "done." << endl;
	
	if (outfile == "") {
		print_clusters(cout, nameset2);
	}
	
	delete [] clusters;
	pthread_exit(NULL);
}

void* doit (void* a) {
	struct aa_init* init = (aa_init*) a;
	cerr << "Starting search for cluster " << (init->c + 1) << endl;
	AlignACE ace(init->seqset, init->nc);
	ace.modify_params(init->argc, init->argv);
	ace.doit();
	string acefn;
	stringstream acefnstr;
	acefnstr << (init->c + 1);
	acefnstr << ".ace";
	acefnstr >> acefn;
	ofstream motifout(acefn.c_str());
	ace.full_output(motifout);
	cerr << "Done with cluster " << (init->c + 1) << endl;
	return a;
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
	vector<int> genes1;
	clusters[c1].genes(genes1);
	for (int i = 0; i < genes1.size(); i++) {
		gene_cluster[genes1[i]] = c2;
	}
	vector<int> genes2;
	clusters[c2].genes(genes2);
	for (int i = 0; i < genes2.size(); i++) {
		gene_cluster[genes2[i]] = c1;
	}
	swap(clusters[c1], clusters[c2]);
}

void merge_clusters(int c1, int c2) {
	vector<int> genes2;
	clusters[c2].genes(genes2);
	for (int i = 0; i < genes2.size(); i++) {
		reassign_gene(genes2[i], c1);
	}
	clusters[c1].calc_mean();
	clusters[c2].calc_mean();
}

void reassign_gene (int g, int c) {
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
	for (int i = 0; i < k; i++) {
		clusters[i].calc_mean();
	}
}

void recenter_cluster(int c, int g) {
	vector<int> genes;
	clusters[c].genes(genes);
	int numgenes = genes.size();
	for (int i = 0; i < numgenes; i++) {
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
	for(int i = 0; i < k; i++) {
		vector<int> genes;
		clusters[i].genes(genes);
		out << "Cluster " << i + 1 << ", " << genes.size() << " orfs" << endl;
		for(int j = 0; j < genes.size(); j++) {
			out << nameset[genes[j]] << endl;
		}
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