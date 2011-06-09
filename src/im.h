#include <iostream>
#include <dirent.h>
#include <signal.h>
#include <errno.h>
#include "standard.h"
#include "motifsearch.h"
#include "motifsearchexpr.h"
#include "motifsearchscore.h"
#include "motifsearchsubset.h"

// Search types
#define UNDEFINED 0
#define EXPRESSION 1
#define SUBSET 2
#define SCORE 3

int worker;                                // worker ID if worker, -1 if archive
bool archive;                              // archive mode

int ngenes;                                // number of genes
int npoints;                               // number of expression data points
int nsubset;                               // number of sequences in search subset
int ncol;                                  // number of columns
int order;                                 // order of background model
double simcut;                             // similarity cutoff for motifs
string outfile;                            // name of output file

int read_motifs(MotifSearch* se);
void output(MotifSearch* se);
void print_clusters(ostream& out, const vector<string>& nameset);
void print_full_ace(ostream& out, MotifSearch* se);
void print_ace(ostream& out, MotifSearch* se);
void print_ace_status(ostream& out, MotifSearch* a, const int i, const int phase, const double cutoff, const double sc);
void print_motifs(ostream& out, const vector<string>& nameset);
void print_usage(ostream& fout);
