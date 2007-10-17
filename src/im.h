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

struct model {
	int num;
	Cluster& c;
	AlignACE& a;
};

void doit(const char* outfile, AlignACE& a, vector<string>& nameset);
void expand_ace_search(AlignACE& a, double mincorr);
float jcorr_lookup(const int g1, const int g2);
float avg_jcorr(AlignACE& a);
void print_clusters(ostream& out, const vector<string>& nameset);
void print_full_ace(ostream& out, AlignACE& a, const vector <string>& nameset);
void print_ace(ostream& out, AlignACE& a, const vector <string>& nameset);
void print_ace_status(ostream& out, AlignACE& a, const int i, const int phase, const double sc);
void print_motifs(ostream& out, const vector<string>& nameset);
void print_usage(ostream& fout);
