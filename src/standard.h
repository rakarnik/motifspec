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

#define MAX_LN_FACT 50000

void get_fasta_fast(const char* filename, vector<string>& seqset, 
			   vector<string>& nameset);
void get_fasta_fast(const char* filename, vector<string>& seqset); 
void get_fasta_fast(istream &test,vector<string>& seqset, vector<string>& nameset);
void get_expr(const char* filename, vector<vector <float> >& expr, vector<string>& nameset);
void get_expr(istream &is, vector<vector <float> >& expr, vector<string>& nameset);
void get_list(const char* filename, vector<string>& listset);
void get_scores(const char* filename, vector<float>& sc, vector<string>& nameset);
void get_scores(istream& scfile, vector<float>& sc, vector<string>& nameset);
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
double log_prob_overlap(int x, int s1, int s2, int n);
double prob_overlap(int x, int s1, int s2, int n);
double bico(int n, int k);
double lnbico(int n, int k);
double lnfact(int n);
double gammaln(double x);
double stirlingln(int n);
double logsum(double x, double y);
float corr(const vector<float>& expr1, const vector<float>& expr2);
float corr(const vector<float>& expr1, const vector<float>& expr2, const unsigned int start1, const unsigned int start2, const unsigned int length);
double find_cutoff(double sum, double sumsq, int num, int num_sdevs_below);
int number_motifs(const char* file);
int number_lines_beg(const char* file, string k);
bool is_number(string s);
string ace_consensus(const char* file, int mot_num);
double ace_mapscore(const char* file, int mot_num);
void alloc_error();
#endif






