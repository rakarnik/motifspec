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
float** expr;              // the expression data
float** jcorr;             // the pairwise jackknife correlation value
int minsize;               // minimum number of sequences required for motif to be considered
double mincorr;               // minimum correlation required to be considered a hit

void doit(const char* outfile, AlignACE& a, vector<string>& nameset);
void expand_ace_search_around_mean(AlignACE& a, Cluster& c, const double corr_cutoff);
void expand_ace_search_all_pairs(AlignACE& a, const double corr_cutoff);
void expand_ace_search_pairs_avg(AlignACE& a, const double corr_cutoff);
float jcorr_lookup(const int g1, const int g2);
float avg_jcorr(AlignACE& a);
void print_clusters(ostream& out, const vector<string>& nameset);
void print_full_ace(ostream& out, AlignACE& a, const vector <string>& nameset);
void print_ace(ostream& out, AlignACE& a, const vector <string>& nameset);
void print_ace_status(ostream& out, AlignACE& a, const int i, const int phase, const double cutoff, const double sc);
void print_motifs(ostream& out, const vector<string>& nameset);
void print_usage(ostream& fout);
