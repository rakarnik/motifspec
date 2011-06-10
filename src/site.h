#ifndef _sites
#define _sites

class Site {
	int chr;
	int pos;
	bool str;
public:
	Site() {};
	Site(const int c, const int p, const bool s);
	Site(const Site& s);
	Site& operator= (const Site& s);
	int chrom() const { return chr; };
	void chrom(const int c) { chr = c; };
	int posit() const { return pos; };
	void posit(const int p) { pos = p; };
	bool strand() const { return str; };
	void strand(const bool strand) { str = strand; };
	void shift(const int shift) { pos += shift; };
	void flip() { str =! str; };
};

#endif
