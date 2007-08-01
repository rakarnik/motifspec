#include "im.h"

int main(int argc, char *argv[]) {
	string seqfile;                       // file with sequences
	string exprfile;                      // file with expression data
	if(argc < 4) {
    print_usage(cout);
    exit(0);
  }
	
  GetArg2(argc, argv, "-s", seqfile);
	GetArg2(argc, argv, "-e", exprfile);
	if (seqfile == "" || exprfile == "") {
		print_usage(cout);
    exit(0);
	}
	
	// Read parameters
	int k;                                // number of clusters
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
	GetArg2(argc, argv, "-minsize", minsize);
	
	string outfile = "";
	GetArg2(argc, argv, "-o", outfile);
	
  vector<string> seqs, nameset1, nameset2;
  cerr << "Reading sequence data from '" << seqfile << "'... ";
	get_fasta_fast(seqfile.c_str(), seqs, nameset1);
	cerr << "done.\n";
	
	vector<vector<float> > expr;          // the expression data
	cerr << "Reading expression data from '" << exprfile << "'... ";
  get_expr(exprfile.c_str(), expr, nameset2);
	cerr << "done.\n";
	
	int ngenes = nameset2.size();                           // number of genes
	/*
	if(nameset1.size() == nameset2.size()) {
		ngenes = nameset1.size();
	} else {
		cerr << "Inconsistent sizes for sequence and expression datasets\n";
		exit(0);
	}
	*/
	
	cerr << "Setting up AlignACE and kmeans models... ";
	vector<AlignACE> gibbs(k, AlignACE(seqs, nc));
	for (int i = 0; i < k; i++) {
		gibbs[i].modify_params (argc, argv);
	}
	vector<Cluster> clusters(k, Cluster(expr, nameset2));
	cerr << "done.\n";
	
	vector<int> gene_cluster(ngenes, (int) -1);             // cluster assignments for genes
	
	cerr << "Assigning cluster centers to random genes... ";
	srand(1);
	for (int i = 0; i < k; i++) {
		recenter_cluster_random(clusters, gene_cluster, i);
	}
	cerr << "done.\n";
	
	int nchanges = 1;
	int toosmall;
	for (int i = 0; i < rounds || nchanges > 0 || toosmall > 0; i++) {
		cerr << "Running iteration " << i + 1 << endl;
		nchanges = assign_genes_to_nearest_cluster(clusters, gene_cluster, expr, cthresh);
		cerr << "\t" << nchanges << " changes in this iteration" << endl;
		cerr << "\t";
		for (int j = 0; j < clusters.size(); j++) {
			cerr << clusters[j].size() << " ";
		}
		cerr << endl;
		check_cluster_distances(clusters, gene_cluster, expr, mthresh);
		
		toosmall = 0;
		if ((i+1) % 10 == 0 || i >= rounds) {
			cerr << "\tResetting small clusters... " << endl; 
			toosmall = check_cluster_sizes(clusters, gene_cluster, expr, minsize);
			cerr << "done!" << endl;
			cerr << "\tThere were " << toosmall << " small clusters" << endl;
		}
		
		sort_clusters(clusters, gene_cluster);
				
		if(i >= rounds && toosmall > 0) {
			cerr << "\tRemoving cluster " << clusters.size() << endl;
			vector<int> genes;
			clusters[clusters.size() - 1].genes(genes);
			for(int g = 0; g < genes.size(); g++) {
				gene_cluster[genes[g]] = -1;
			}
			clusters.erase(clusters.end());
		}
		
		if (outfile != "") {
			ofstream out(outfile.c_str());
			print_clusters(out, clusters, nameset2);
			out.close();
		}
	}
	
	if (outfile == "") {
		print_clusters(cout, clusters, nameset2);
	}
}

int assign_genes_to_nearest_cluster(vector<Cluster>& clusters, vector<int>& gene_cluster, const vector<vector<float> >& expr, float threshold) {
	float corr, max_corr;
	int max_clus, nchanges = 0;
	for (int g = 0; g < gene_cluster.size(); g++) {
		max_clus = -1;
		max_corr = 0;
		for(int i = 0; i < clusters.size(); i++) {
			corr = clusters[i].corr(expr[g]);
			if (corr > max_corr && corr > threshold) {
				max_clus = i;
				max_corr = corr;
			}
		}
		if (max_clus != gene_cluster[g]) {
			nchanges++;
			reassign_gene(clusters, gene_cluster, g, max_clus);
		}
	}
	update_cluster_means(clusters);
	return nchanges;
}

int check_cluster_distances(vector<Cluster>& clusters, vector<int>& gene_cluster, const vector<vector<float> >& expr, float threshold) {
	float c;
	int big, small, count;
	for (int i = 0; i < clusters.size(); i++) {
		for (int j = i+1; j < clusters.size(); j++) {
			c = clusters[i].corr(clusters[j].get_mean());
			if (c > threshold) {
				if (clusters[i].size() > clusters[j].size()) {
					big = i;
					small = j;
				} else {
					big = j;
					small = i;
				}
				cerr << "\tMerging cluster " << (small + 1) << " into cluster " << (big + 1) << "(correlation was " << c << ")" << endl;
				merge_clusters(clusters, gene_cluster, big, small);
				recenter_cluster(clusters, gene_cluster, small, find_outlier(clusters, gene_cluster, expr));
				count++;
			}
		}
	}
	update_cluster_means(clusters);
	return count;
}

int check_cluster_sizes(vector<Cluster>& clusters, vector<int>& gene_cluster, const vector<vector<float> >& expr, int minsize) {
	int delcount = 0;
	for (int i = 0; i < clusters.size(); i++) {
		if (clusters[i].size() < minsize) {
			cerr << "\tResetting cluster " << i + 1 << "(" << clusters[i].size() << ")" << endl;
			recenter_cluster(clusters, gene_cluster, i, find_outlier(clusters, gene_cluster, expr));
			delcount++;
		}
	}
	return delcount;
}

void sort_clusters(vector<Cluster>& clusters, vector<int>& gene_cluster) {
	for (int i = 0; i < clusters.size(); i++) {
		for (int j = i + 1; j < clusters.size(); j++) {
			if (clusters[i].size() < clusters[j].size()) {
				swap(clusters, gene_cluster, i, j);
			}
		}
	}
}

void swap(vector<Cluster>& clusters, vector<int>& gene_cluster, int c1, int c2) {
	swap(clusters[c1], clusters[c2]);
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
}

void merge_clusters(vector<Cluster>& clusters, vector<int>& gene_cluster, int c1, int c2) {
	vector<int> genes2;
	clusters[c2].genes(genes2);
	for (int i = 0; i < genes2.size(); i++) {
		reassign_gene(clusters, gene_cluster, genes2[i], c1);
	}
	clusters[c1].calc_mean();
	clusters[c2].calc_mean();
}

void reassign_gene (vector<Cluster>& clusters, vector<int>& gene_cluster, int g, int c) {
		if (c != -1) {
			clusters[c].add_gene(g);
		}
		if (gene_cluster[g] != -1) {
			clusters[gene_cluster[g]].remove_gene(g);
		}
		gene_cluster[g] = c;
}

void update_cluster_means(vector<Cluster>& clusters) {
	for (int i = 0; i < clusters.size(); i++) {
		clusters[i].calc_mean();
	}
}

void recenter_cluster(vector<Cluster>& clusters, vector<int>& gene_cluster, int c, int g) {
	vector<int> genes;
	clusters[c].genes(genes);
	for (int i = 0; i < genes.size(); i++) {
		clusters[c].remove_gene(genes[i]);
		gene_cluster[genes[i]] = -1;
	}
	reassign_gene(clusters, gene_cluster, g, c);
	clusters[c].calc_mean();
}

void recenter_cluster_random(vector<Cluster>& clusters, vector<int>& gene_cluster, int c) {
	int new_center = (int) ((float) rand()/RAND_MAX * gene_cluster.size());
	recenter_cluster(clusters, gene_cluster, c, new_center);
}

int find_outlier(const vector<Cluster>& clusters, const vector<int>& gene_cluster, const vector<vector<float> >& expr) {
	float corr, min_corr = 1.0;
	int outlier;
	for (int g = 0; g < gene_cluster.size(); g++) {
		for(int i = 0; i < clusters.size(); i++) {
			corr = clusters[i].corr(expr[g]);
			if (corr < min_corr) {
				outlier = g;
				min_corr = corr;
			}
		}
	}
	return outlier;
}

void print_clusters(ostream& out, const vector<Cluster>& clusters, const vector<string>& nameset) {
	for(int i = 0; i < clusters.size(); i++) {
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