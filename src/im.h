#include <iostream>
#include <dirent.h>
#include <signal.h>
#include <errno.h>
#include "semodel.h"
#include "standard.h"

// Search types
#define UNDEFINED 0
#define EXPRESSION 1
#define SUBSET 2

int worker;                                // worker ID if worker, -1 if archive
bool archive;                              // archive mode

int ngenes;                                // number of genes
int npoints;                               // number of expression data points
int nsubset;                               // number of sequences in search subset
int ncol;                                  // number of columns
int order;                                 // order of background model
double simcut;                             // similarity cutoff for motifs
vector<vector <float> > expr;              // the expression data
string outfile;                            // name of output file

int read_motifs(SEModel& se);
void output(SEModel& se);
void print_clusters(ostream& out, const vector<string>& nameset);
void print_full_ace(ostream& out, SEModel& se);
void print_ace(ostream& out, SEModel& se);
void print_ace_status(ostream& out, SEModel& a, const int i, const int phase, const double cutoff, const double sc);
void print_motifs(ostream& out, const vector<string>& nameset);
void print_usage(ostream& fout);
