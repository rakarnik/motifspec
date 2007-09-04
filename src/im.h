/*
 *  im.cpp
 */
#include <iostream>
#include "alignace.h"
#include "cluster.h"
#include "standard.h"

int k;                     // the number of clusters
int ngenes;                // the number of genes
int npoints;               // the number of data points
int* nmots;								 // the number of motifs found in each cluster
Cluster* bsclusters;       // the bootstrap k-means clusters
Cluster** clusters;        // the adjustment k-means clusters
AlignACE* bsaces;          // the bootstrap AlignACE models
AlignACE** aces;           // the adjustment AlignACE models
int* gene_cluster;         // mapping of genes to clusters
float** expr;              // the expression data

int assign_genes_to_nearest_cluster(float threshold);
int check_cluster_distances(float threshold);
int check_cluster_sizes(int minsize);
int reset_small_clusters(int minsize);
void sort_clusters();
void swap(int c1, int c2);
void merge_clusters(int c1, int c2);
void reassign_gene(int g, int c);
void update_bscluster_means();
void update_cluster_means();
void recenter_cluster(int c, int g);
void recenter_cluster_random(int c);
int find_outlier();
void print_bsclusters(ostream& out, const vector<string>& nameset);
void print_clusters(ostream& out, const vector<string>& nameset);
void print_full_ace(ostream& out, AlignACE& a, const vector <string>& nameset);
void print_ace(ostream& out, AlignACE& a, const double score, const vector <string>& nameset);
void print_motifs(ostream& out, const vector<string>& nameset);
void print_usage(ostream& fout);
void debug_check_membership();
