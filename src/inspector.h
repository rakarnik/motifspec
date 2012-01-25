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

vector<string> seq_nameset;                // names associated with sequences
vector<string> data_nameset;               // names associated with expression values/binding scores
vector<vector <float> > expr;              // expression data
vector<float> scores;                      // binding scores
vector<string> subset;                     // names of sequences to be searched
int ngenes;                                // number of sequences
int npoints;                               // number of expression data points
int nsubset;                               // number of sequences in search subset
int ncol;                                  // number of columns
int order;                                 // order of background model
double simcut;                             // similarity cutoff for motifs
int maxm;                                  // maximum number of motifs
string outfile;                            // name of output file

void order_data_expr(vector<vector <float> >& newexpr);
void order_data_scores(vector <float>& newscores);
int read_motifs(MotifSearch* se);
void output(MotifSearch* se);
void print_clusters(ostream& out, const vector<string>& nameset);
void print_full_ace(ostream& out, MotifSearch* se);
void print_ace(ostream& out, MotifSearch* se);
void print_ace_status(ostream& out, MotifSearch* a, const int i, const int phase, const double cutoff, const double sc);
void print_motifs(ostream& out, const vector<string>& nameset);
void print_usage(ostream& fout);

