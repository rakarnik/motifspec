//Copyright 1999 President and Fellows of Harvard University
//sites.h

#ifndef _sites
#define _sites
#include "standard.h"
#include "seqset.h"

class Site {
	int chr;
	int pos;
	bool str;
public:
	Site() {};
	Site(const int c, const int p, const bool s);
	int chrom() const { return chr; };
	int posit() const { return pos; };
	bool strand() const { return str; };
	void shift(const int shift) { pos += shift; };
	void flip() { str =! str; };
};

#endif
