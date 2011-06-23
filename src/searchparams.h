#ifndef _searchparams
#define _searchparams

struct SearchParams {
	int expect;						           // number of expected sites
	double weight;				           // fractional weight on priors
	double psfact;				           // psfact * numsites=npseudo
	double npseudo;				           // number of pseudo counts
	vector<double> backfreq;         // array for gc content
	vector<double> pseudo;           // pseudocounts for any frequency calculations
	int maxlen;						           // maximum length of sites
	int npass;
	int minpass;
	int nruns;
	float minprob[4];                // minimum P(m|s) cutoff in each phase
	float minscore[4];               // minimum score cutoff in each phase
	bool fragment;
	int seed;
	double select;
	int flanking;
	int undersample;
	int oversample;
	int minsize;
	float mincorr;
};

#endif

