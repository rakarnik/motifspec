//Copyright 1999 President and Fellows of Harvard University
//sites.cpp

#include "site.h"

Site::Site(const int c, const int p, const bool s) :
chr(c),
pos(p),
str(s) {
}

Site::Site(const Site& s) :
chr(s.chr),
pos(s.pos),
str(s.str) {
	*this = s;
}

Site& Site::operator= (const Site& s) {
	if(this != &s) {
		chr = s.chr;
		pos = s.pos;
		str = s.str;
	}
	return *this;
}

