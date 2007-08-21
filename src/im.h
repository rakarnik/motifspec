/*
 *  im.cpp
 */
#include <iostream>
#include <pthread.h>
#include "alignace.h"
#include "cluster.h"
#include "standard.h"

int k;                     // the number of clusters
int ngenes;                // the number of genes
int npoints;               // the number of data points
Cluster* clusters;         // the k-means clusters
int* gene_cluster;         // mapping of genes to clusters
float** expr;              // the expression data

struct aa_init {
	int c;
	vector<string> seqset;
	int nc;
	int argc;
	char** argv;
};

void* doit(void* a);
int assign_genes_to_nearest_cluster(float threshold);
int check_cluster_distances(float threshold);
int check_cluster_sizes(int minsize);
void sort_clusters();
void swap(int c1, int c2);
void merge_clusters(int c1, int c2);
void reassign_gene(int g, int c);
void update_cluster_means();
void recenter_cluster(int c, int g);
void recenter_cluster_random(int c);
int find_outlier();
void print_clusters(ostream& out, const vector<string>& nameset);
void print_usage(ostream& fout);
void debug_check_membership();
