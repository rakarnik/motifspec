//Copyright 1998 President and Fellows of Harvard University
//standard.h

#ifndef _standard
#define _standard
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <iterator>
#include <cmath>
#include <sstream>
#include <fcntl.h>
#include "limits.h"
#include "time.h"
#include "stdlib.h"
#include "string.h"
#include "float.h"
using namespace std;

void get_fasta_fast(const char* filename, vector<string>& seqset, 
			   vector<string>& nameset);
void get_fasta_fast(const char* filename, vector<string>& seqset); 
void get_fasta_fast(istream &test,vector<string>& seqset, vector<string>& nameset);
float** get_expr(const char* filename, int* npoints, vector<string>& nameset);
float** get_expr(istream &is, int* npoints, vector<string>& nameset);
void get_cluster(const char* filename, const int num, vector<string>& nameset);
void get_cluster(istream &is, const int num, vector<string>& nameset);
bool GetArg2(int argc, char *argv[], const char *c, int &cval);
bool GetArg2(int argc, char *argv[], const char *c, float &cval);
bool GetArg2(int argc, char *argv[], const char *c, double &cval);
bool GetArg2(int argc, char *argv[], const char *c, string &cval);
bool GetArg2(int argc, char *argv[], const char *c);
string reverse_comp(const string &forward);
vector<string> split(string s, char c, bool skipall=false);
int convert_roman(string s);
string capitalize(string s);
string lower_case(string s);
int str_to_int(const string &s);
string int_to_str(int x);
string clip_white(const string &s);
double str_to_dbl(const string &s);
string random_dna(int len);
double gammaln(double x);
double stirlingln(int n);
double  lnfact(int n);
double  bico(int N, int k);
double  lnbico(int N, int k);
float corr(const float* expr1, const float* expr2, const int num, const int jindex = -1);
float jack_corr(const float* expr1, const float* expr2, const int num);
double find_cutoff(double sum, double sumsq, int num, int num_sdevs_below);
int number_motifs(const char* file);
int number_lines_beg(const char* file, string k);
bool is_number(string s);
double prob_overlap(int x, int y, int i, int t);
string ace_consensus(const char* file, int mot_num);
double ace_mapscore(const char* file, int mot_num);
void alloc_error();
#endif






