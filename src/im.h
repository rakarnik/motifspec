/*
 *  im.cpp
 */
#include <iostream>
#include "alignace.h"
#include "cluster.h"
#include "standard.h"



int assign_genes_to_nearest_cluster(vector<Cluster>& clusters, vector<int>& gene_cluster, const vector<vector<float> >& expr, float threshold);
int check_cluster_distances(vector<Cluster>& clusters, vector<int>& gene_cluster, const vector<vector<float> >& expr, float threshold);
int check_cluster_sizes(vector<Cluster>& clusters, vector<int>& gene_cluster, const vector<vector<float> >& expr, int minsize);
void sort_clusters(vector<Cluster>& clusters, vector<int>& gene_cluster);
void swap(vector<Cluster>& clusters, vector<int>& gene_cluster, int c1, int c2);
void merge_clusters(vector<Cluster>& clusters, vector<int>& gene_cluster, int c1, int c2);
void reassign_gene(vector<Cluster>& clusters, vector<int>& gene_cluster, int g, int c);
void update_cluster_means(vector<Cluster>& clusters);
void recenter_cluster(vector<Cluster>& clusters, vector<int>& gene_cluster, int c, int g);
void recenter_cluster_random(vector<Cluster>& clusters, vector<int>& gene_cluster, int c);
int find_outlier(const vector<Cluster>& clusters, const vector<int>& gene_cluster, const vector<vector<float> >& expr);
void print_clusters(ostream& out, const vector<Cluster>& clusters, const vector<string>& nameset);
void print_usage(ostream& fout);
