/*
 *  im.cpp
 */
#include <iostream>
#include "semodel.h"
#include "standard.h"

int k;                     // the number of clusters
int ngenes;                // the number of genes
int npoints;               // the number of data points
int ncol;
float** expr;              // the expression data
float** jcorr;             // the pairwise jackknife correlation value
int minsize;               // minimum number of sequences required for motif to be considered
double mincorr;               // minimum correlation required to be considered a hit

void doit(const char* outfile, SEModel& se);
void print_clusters(ostream& out, const vector<string>& nameset);
void print_full_ace(ostream& out, SEModel& se);
void print_ace(ostream& out, SEModel& se);
void print_ace_status(ostream& out, SEModel& a, const int i, const int phase, const double cutoff, const double sc);
void print_motifs(ostream& out, const vector<string>& nameset);
void print_usage(ostream& fout);
