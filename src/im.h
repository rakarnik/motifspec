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
Cluster* clusters;         // the adjustment k-means clusters
AlignACE* aces;            // the adjustment AlignACE models
int* gene_cluster;         // mapping of genes to clusters
float** expr;              // the expression data

void doit(Cluster& c, AlignACE& a);
int assign_genes_to_nearest_cluster(float threshold);
int check_cluster_distances(float threshold);
int check_cluster_sizes(int minsize);
int reset_small_clusters(int minsize);
void sort_clusters();
void swap(int c1, int c2);
void merge_clusters(int c1, int c2);
void reassign_gene(int g, int c);
void update_cluster_means();
void recenter_cluster(int c, int g);
void recenter_cluster_random(int c);
int find_outlier();
void print_clusters(ostream& out, const vector<string>& nameset);
void print_full_ace(ostream& out, AlignACE& a, const vector <string>& nameset);
void print_ace(ostream& out, AlignACE& a, const vector <string>& nameset);
void print_motifs(ostream& out, const vector<string>& nameset);
void print_usage(ostream& fout);
void debug_check_membership();
void sync_ace_members(Cluster& c, AlignACE& a);
void sync_ace_neighborhood(Cluster& c, AlignACE& a, double mincorr);
void sync_cluster(Cluster& c, AlignACE& a);
