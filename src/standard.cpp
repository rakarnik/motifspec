//Copyright 1999 President and Fellows of Harvard University
//standard.cpp

#include "standard.h"

void get_fasta_fast(istream &test, vector<string>& seqset, vector<string>& nameset) {
/*This is a routine that reads in a fasta formatted file fname
and returns a vector of strings seqset corresponding to the sequences
in the file.  Fasta format allows for a set of sequences separated by
titles beginning with '>' and ending with '\n'.*/
	string line, name, seq;
	while(getline(test, line)) {
		if(line[0] == '>') {
			if(! name.empty()) {
				string namen = name;
				nameset.push_back(namen);
				string seqn = seq;
				seqset.push_back(seqn);
				seq.clear();
			}
			name = line.substr(1);
		} else {
			seq.append(line);
		}
	}
	if(! name.empty()) {
		if(seq.size() < 1) {
			cerr << "Invalid sequence for '" << name << "'!\n";
			exit(1);
		}
		nameset.push_back(name);
		string seqn = seq;
		seqset.push_back(seqn);
	}
	assert(nameset.size() == seqset.size());
}

void get_fasta_fast(const char* fname, vector<string>& seqset, vector<string>& nameset) {
	ifstream test(fname);
	if(! test) {
	  cerr << "No such file '" << fname << "'\n";
	  exit(0);
	}
	get_fasta_fast(test, seqset, nameset);
}

void get_fasta_fast(const char* filename, vector<string>& seqset){
	//when you don't care about titles
	vector<string> nameset;
	get_fasta_fast(filename, seqset, nameset);
}

void get_expr(const char* filename, vector<vector <float> >& expr, vector<string>& nameset){
	ifstream exprfile(filename);
	if(! exprfile){
	  cerr << "No such file '" << filename << "'\n";
	  exit(0);
	}
	return get_expr(exprfile, expr, nameset);
}

void get_expr(istream& exprfile, vector<vector <float> >& expr, vector<string>& nameset) {
	string line;
	int npoints = 0;
	
	vector<string> values;
	vector<string>::iterator val_iter;
	while(getline(exprfile, line)) {
		vector<float> row;
		row.reserve(npoints);
		values = split(line, '\t');
		nameset.push_back(values[0]);
		val_iter = values.begin() + 1;
		for (; val_iter != values.end(); ++val_iter) {
			row.push_back(atof(val_iter->c_str()));
		}
		expr.push_back(row);
	}
}

void get_list(const char* filename, vector<string>& listset) {
	ifstream listfile(filename);
	if(! listfile){
	  cerr<<"No such file "<< filename<<'\n';
	  exit(0);
	}
	string line;
	while(getline(listfile, line)) {
		if(line != "") {
			listset.push_back(line);
		}
	}
}

void get_scores(const char* filename, vector<float>& sc, vector<string>& nameset){
	ifstream scfile(filename);
	if(! scfile){
	  cerr << "No such file '" << filename << "'\n";
	  exit(0);
	}
	return get_scores(scfile, sc, nameset);
}

void get_scores(istream& scfile, vector<float>& sc, vector<string>& nameset) {
	string line;
	vector<string> fields;
	vector<string>::iterator val_iter;
	while(getline(scfile, line)) {
		fields = split(line, '\t');
		assert(fields.size() >= 2);
		nameset.push_back(fields[0]);
		sc.push_back(atof(fields[1].c_str()));
	}
}

void get_cluster(const char* filename, const int num, vector<string>& nameset) {
	ifstream clusfile(filename);
	if(! clusfile){
	  cerr<<"No such file "<< filename<<'\n';
	  exit(0);
	}
	return get_cluster(clusfile, num, nameset);
}

void get_cluster(istream &clusfile, const int num, vector<string>& nameset) {
	bool in_cluster = false;;
	string line;
	string clus1name, clus2name;
	stringstream clus1strm, clus2strm;
	clus1strm << "Cluster " << num << ",";
	getline(clus1strm, clus1name);
	clus2strm << "Cluster " << (num + 1) << ",";
	getline(clus2strm, clus2name);
	while(getline(clusfile, line)) {
		if(line.find(clus1name.c_str()) != line.npos) {
			in_cluster = true;
			continue;
		} else if(line.find (clus2name.c_str()) != line.npos) {
			in_cluster = false;
			break;
		}
		if(in_cluster && line != "") {
			nameset.push_back(line);
		}
	}
}


//The following GetArg2 commands change the command line option parsing
//to expect spaces, i.e., -n 2 instead of -n2.  This allows for more
//flexible parameter input and allows for tab completion from unix
//command line for filename inputs.


bool GetArg2(int argc, char *argv[], const char *c, int &cval){
  if(argc < 2) return false;
  for (int arg=1;arg<argc;arg++) {
    if (strcmp(argv[arg],c)==0) {
      string s=argv[++arg];
      cval=str_to_int(s);
      return true;
    }
  }
  return false;
}

bool GetArg2(int argc, char *argv[], const char *c, float &cval){
  if(argc < 2) return false;
  for (int arg=1;arg<argc;arg++) {
    if (strcmp(argv[arg],c)==0) {
      string s=argv[++arg];
      cval=(float)str_to_dbl(s);
      return true;
    }
  }
  return false;
}

bool GetArg2(int argc, char *argv[], const char *c, double &cval){
  if(argc < 2) return false;
  for (int arg=1;arg<argc;arg++) {
    if (strcmp(argv[arg],c)==0) {
      string s=argv[++arg];
      cval=str_to_dbl(s);
      return true;
    }
  }
  return false;
}

bool GetArg2(int argc, char *argv[], const char *c, string &cval){
  if(argc < 2) return false;
  for (int arg=1;arg<argc;arg++) {
    if (strcmp(argv[arg],c)==0) {
      cval=argv[++arg];
      return true;
    }
  }
  return false;
}

bool GetArg2(int argc, char *argv[], const char *c){
  if(argc < 2) return false;
  for (int arg=1;arg<argc;arg++) {
    if (strcmp(argv[arg],c)==0) {
      return true;
    }
  }
  return false;
}

string reverse_comp(const string &forward){
	static map<char,char> comp;
	static bool init=false;
	if(!init){
		comp['a']='t';
		comp['A']='T';
		comp['c']='g';
		comp['C']='G';
		comp['g']='c';
		comp['G']='C';
		comp['t']='a';
		comp['T']='A';
		comp[' ']=' ';
		init=true;
	}
	string rev;
	string::const_reverse_iterator p;
	for(p=forward.rbegin();p!=forward.rend();p++){
		rev+=comp[*p];
	}
	return rev;
}

float mean(const vector<float>& v) {
	float m = 0.0;
	vector<float>::const_iterator viter = v.begin();
	for(; viter != v.end(); ++viter)
		m += *viter;
	m /= v.size();
	return m;
}

float stdev(const vector<float>& v) {
	float m = mean(v);
	float std = 0.0;
	vector<float>::const_iterator viter = v.begin();
	for(; viter != v.end(); ++viter)
		std += (*viter - m) * (*viter - m);
	std /= v.size();
	std = sqrt(std);
	return std;
}

float corr(const vector<float>& expr1, const vector<float>& expr2) {
	return corr(expr1, expr2, 0, 0, min(expr1.size(), expr2.size()));
}

float corr(const vector<float>& expr1, const vector<float>& expr2, const unsigned int start1, const unsigned int start2, const unsigned int length) {
	if(length == 0) return 0.0;
	float u1, u2, u11, u22, u12, s1, s2, c;
	u1 = u2 = u11 = u22 = u12 = s1 = s2 = c = 0.0;
	assert(start1 + length - 1 <= expr1.size());
	assert(start2 + length - 2 <= expr2.size());
	unsigned int i, j;
	for (i = start1, j = start2; i < start1 + length && j < start2 + length; i++, j++) {
		u1 += expr1[i];
		u2 += expr2[j];
		u11 += expr1[i] * expr1[i];
		u22 += expr2[j] * expr2[j];
		u12 += expr1[i] * expr2[j];
	}
	u1 /= length;
	u2 /= length;
	u11 /= length;
	u22 /= length;
	u12 /= length;
	if(u11 > u1 * u1) {
		s1 = sqrt(u11 - u1 * u1);
	}
	if(u22 > u2 * u2) {
		s2 = sqrt(u22 - u2 * u2);
	}
	if (s1 > 0 && s2 > 0) {
		c = (u12 - u1 * u2)/(s1*s2);
	}
	return c;
}

double log_prob_overlap(int x, int s1, int s2, int n) {
	assert(x <= s1);
	assert(x <= s2);
	assert(s1 <= n);
	assert(s2 <= n);

	if(x == 0)
		return 0;

	int m1 = min(s1, s2);
	int m2 = max(s1, s2);
	double ret = 0.0, lt;
	for(int i = x; i <= m1; i++) {
		lt = lnbico(m1, i) + lnbico(n - m1, m2 - i) - lnbico(n, m2);
		if(! ret)
			ret = lt;
		else
			ret = logsum(ret, lt);
		if(lt < ret - 5) break;
	}
	return ret;
}

double prob_overlap(int x, int s1, int s2, int n) {
  assert(x <= s1);
	assert(x <= s2);
	assert(s1 <= n);
	assert(s2 <= n);

	if(x == 0)
		return 1;

	int m1 = min(s1, s2);
	int m2 = max(s1, s2);
	double ret, lt;
	if(x > m1 - x) {
		ret = 0.0;
		for(int i = x; i <= m1; i++) {
			lt = lnbico(m1, i) + lnbico(n - m1, m2 - i) - lnbico(n, m2);
			ret += lt;
		}
		return ret;
	} else {
		ret = 1.0;
		for(int i = 0; i < x; i++) {
			lt = lnbico(m1, i) + lnbico(n - m1, m2 - i) - lnbico(n, m2);
			ret -= lt;
		}
	}
	return ret;
}

double  bico(int n, int k) {
	return floor(0.5 + exp(lnfact(n) - lnfact(k) - lnfact(n-k)));
}

double  lnbico(int n, int k) {
	return lnfact(n) - lnfact(k) - lnfact(n-k);
}

double	lnfact(int n)
/* static variables are guaranteed to be initialized to zero */
{
	static double lnft[MAX_LN_FACT];
	assert(n < MAX_LN_FACT);

	if (n <= 1) return 0.0;
	if (n <= 50) return lnft[n] ? lnft[n] : (lnft[n] = gammaln(n + 1.0));
	else return lnft[n]? lnft[n] : (lnft[n] = stirlingln(n));
}

double gammaln(double xx) {
	double x, y, tmp, ser;
	static double cof[6]={76.18009172947146, -86.50532032941677,
		24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
		-0.5395239384953e-5};
	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.000000000190015;
	for (int j = 0; j <= 5; j++)
		ser += cof[j]/++y;
	return -tmp + log(2.5066282746310005*ser/x);
}

double stirlingln(int n) {
	return (0.91984 + (n + 0.5) * log(n) - n);
}

double logsum(double x, double y) {
	return x + log(1 + exp(y - x));
}

vector<string> split(string s, char c, bool skipall){
	//skipall default to false, i.e. for delimited files
	istringstream st(s.c_str());
	vector<string> vs;
	string w;
	if(s[s.length()-1]!=c) s+=c;//get last field if no final delimiter
	while(getline(st,w,c)){
		if(skipall && w.size()>0) vs.push_back(w);
		if(!skipall) vs.push_back(w);
	}
	return vs;
}

string capitalize(string s){
	//make all caps
	string::iterator p=s.begin();
	for(;p!=s.end();p++){
		*p=toupper(*p);
	}
	return s;
}

string lower_case(string s){
	//make all lower-case
	string::iterator p=s.begin();
	for(;p!=s.end();p++){
		*p=tolower(*p);
	}
	return s;
}

int convert_roman(string s){
	static bool init=false;
	static map<char,int> roman;
	if(!init){
		roman['I']=1;
		roman['V']=5;
		roman['X']=10;
		roman['L']=50;
		roman['C']=100;
		roman['D']=500;
		roman['M']=1000;
		init=true;
	}
	s=capitalize(s);
	int res=0;
	int last=0;
	int val;
	string::reverse_iterator i=s.rbegin();
	for(;i!=s.rend();i++){
		val=roman[*i];
		if(val>=last) res+=val;
		else res-=val;
		last=val;
	}
	return res;
}

int str_to_int(const string &s){
	//assume no overflow
	//assume that the input is only digits, everything else will be ignored
	double ret=str_to_dbl(s);
	if(ret>INT_MAX) return -1;
	return (int)ret;
	//default behavior is truncation
}

double str_to_dbl(const string &s){
	//assume no overflow
	//assume that the input is only digits, everything else will be ignored
	static bool init=false;
	static map<char,int> digit;
	if(!init){
		digit['0']=0;
		digit['1']=1;
		digit['2']=2;
		digit['3']=3;
		digit['4']=4;
		digit['5']=5;
		digit['6']=6;
		digit['7']=7;
		digit['8']=8;
		digit['9']=9;
		init=true;
	}
	double ret=0.0;double dec=0.0;int pnx=1;int xp=0;double pn=1.0;
	unsigned int i=0;
	if(s[i]=='-') {pn=-1.0;i++;}
	if(s[i]=='+') {pn=1.0;i++;}
	for(;i<s.size();i++){
		if(!isdigit(s[i])) break;
		ret*=10.0;
		ret+=(double)digit[s[i]];
	}
	if(s[i]=='.'){
		i++;
		double dp=0.1;
		for(;i<s.size();i++){
			if(!isdigit(s[i])) break;
			dec+=digit[s[i]]*dp;
			dp/=10.0;
		}
	}
	ret+=dec;
	ret*=pn;
	if(s[i]=='e'||s[i]=='E'){
		i++;
		if(s[i]=='-') {pnx=-1;i++;}
		if(s[i]=='+') {pnx=1;i++;}
		for(;i<s.size();i++){
			xp*=10;
			xp+=digit[s[i]];
		}
	}
	xp*=pnx;
	ret*=pow(10.0,(double)xp);
	return ret;
}

string random_dna(int len){
	static bool seeded=false;
	if(!seeded){
		seeded=true;
		srand( (unsigned)time( NULL ) );
		rand();rand();rand();rand();rand();
	}
	string s;
	for(int i=0;i<len;i++){
		double r=(double)rand()/RAND_MAX;
		if(r<0.25) s+='A';
		else if(r<0.50) s+='C';
		else if(r<0.75) s+='G';
		else s+='T';
	}
	return s;
}

string int_to_str(int x){
	static bool init=false;
	static map<int,char> conv;
	if(!init){
		conv[0]='0';
		conv[1]='1';
		conv[2]='2';
		conv[3]='3';
		conv[4]='4';
		conv[5]='5';
		conv[6]='6';
		conv[7]='7';
		conv[8]='8';
		conv[9]='9';
		init=true;
	}
	string s,t;
	if(x<0){
		s+='-';
		x*=-1;
	}
	for(;x>0;){
		int y=x%10;
		x=x/10;
		t+=conv[y];//reverse order
	}
	for(int j=t.size()-1;j>=0;j--){
		s+=t[j];
	}
	return s;
}

string clip_white(const string &s){
	string t;
	unsigned int first=0;
	unsigned int i;
	for(i = 0; i < s.size(); i++){
		if(!isspace(s[i])) {
			first = i;
			break;
		}
	}
	unsigned int last = s.size() - 1;
	for(i = s.size() - 1; i > first; i--){
		if(!isspace(s[i])) {
			last = i;
			break;
		}
	}
	for(i = first; i <= last; i++){
		t += s[i];
	}
	return t;
}

double find_cutoff(double sum, double sumsq, int num, int num_sdevs_below){
	double sd=sqrt( sumsq/(num-1) - (sum/num) * (sum/(num-1)) );
	double avg=sum/num;
	return avg-num_sdevs_below*sd;
}

int number_motifs(const char* file){
  ifstream fin(file);
  string s;
  int x=0;
  while (getline(fin,s,'\n')){
    if(s.size()<5) continue;
    if(s.substr(0,5)=="Motif"){
      x++;
    }
  }
  fin.close();
  return x;
}

bool is_number(string s){
  for(unsigned int i=0;i<s.size();i++){
    if(!isdigit(s[i])&&s[i]!='.') return false;
  }
  return true;
}

int number_lines_beg(const char* file, string k){
  ifstream fin(file);
  string s;
  int x=0;
  int n=k.size();
  while (getline(fin,s,'\n')){
    if(s.substr(0,n)==k){
      x++;
    }
  }
  fin.close();
  return x;
}

string ace_consensus(const char* file, int mot_num){
  ifstream fin(file);
  string s;
  int x=0;
  vector<string> a,s1;
  while (getline(fin,s,'\n')){
    if(s.size()<5) continue;
    if(s.substr(0,5)=="Motif"){
      a=split(s,' ');
      x=str_to_int(a[a.size()-1]);
      continue;
    }
    if(s.find('*')!=string::npos) x=0;
    if(x<mot_num) continue;
    else if(x>mot_num) break;
    else if(x==mot_num){
      string t=s.substr(0,s.find_first_of('\t'));
      s1.push_back(t);
      //      cerr<<t<<'\n';
    }
  }
  string ret;
  int wide1=s1[0].size();
  int num1=s1.size();
  int i,j;
  int num1a,num1c,num1g,num1t;
  for(i=0;i<wide1;i++){
    num1a=num1c=num1g=num1t=0;
    for(j=0;j<num1;j++){
      if(s1[j][i]=='a'||s1[j][i]=='A') num1a+=1;
      if(s1[j][i]=='c'||s1[j][i]=='C') num1c+=1;
      if(s1[j][i]=='g'||s1[j][i]=='G') num1g+=1;
      if(s1[j][i]=='t'||s1[j][i]=='T') num1t+=1;
    }
    if(num1a>num1*0.7) ret+='A';
    else if(num1c>num1*0.7) ret+='C';
    else if(num1g>num1*0.7) ret+='G';
    else if(num1t>num1*0.7) ret+='T';
    else if((num1a+num1t)>num1*0.85) ret+='W';
    else if((num1a+num1c)>num1*0.85) ret+='M';
    else if((num1a+num1g)>num1*0.85) ret+='R';
    else if((num1t+num1c)>num1*0.85) ret+='Y';
    else if((num1t+num1g)>num1*0.85) ret+='K';
    else if((num1c+num1g)>num1*0.85) ret+='S';
    else ret+='-';
  }
  return ret;
}

double ace_mapscore(const char* file, int mot_num){
  ifstream fin(file);
  string s;
  int x=0;
  double ret=-1.0;
  vector<string> a;
  while (getline(fin,s,'\n')){
    if(s.size()<5) continue;
    if(s.substr(0,5)=="Motif"){
      a=split(s,' ');
      x=str_to_int(a[a.size()-1]);
      continue;
    }
    if(s.substr(0,3)=="MAP"&&mot_num==x){
      a=split(s,' ');
      ret=str_to_dbl(a[a.size()-1]);
      break;
    }
  }
  return ret;
}

void alloc_error() {
	cerr << "Memory allocation failed!\n";
	abort();
}

