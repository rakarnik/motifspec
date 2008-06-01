/*
 *  im.cpp
 */
#include <iostream>
#include <dirent.h>
#include <signal.h>
#include "semodel.h"
#include "standard.h"

int worker;                // worker ID if worker, -1 if archive
bool archive;              // archive mode

int ngenes;                // number of genes
int npoints;               // number of data points
int ncol;                  // number of columns
float** expr;              // the expression data
SEModel se;                // the sequence-expression model

string outfile;

int read_motifs();
void output();
void final_output(int param);
void print_clusters(ostream& out, const vector<string>& nameset);
void print_full_ace(ostream& out, SEModel& se);
void print_ace(ostream& out, SEModel& se);
void print_ace_status(ostream& out, SEModel& a, const int i, const int phase, const double cutoff, const double sc);
void print_motifs(ostream& out, const vector<string>& nameset);
void print_usage(ostream& fout);
