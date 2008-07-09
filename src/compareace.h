//Copyright 1999 President and Fellows of Harvard University
//compareace.h

#ifndef _compareace
#define _compareace
#include "standard.h"
#include "motif.h"

class CompareACESites{
 	Motif comp_motif;
  bool comp_ready;
  double *comp_fm1;
  double *comp_fm2;
  int comp_fmsize;
  int comp_l1,comp_r1,comp_l2,comp_r2;

public:
	CompareACESites(const Motif &m);
	CompareACESites(const CompareACESites& s);
  ~CompareACESites();
  CompareACESites& operator=(const CompareACESites& s);
  double compare(const CompareACESites &c);
  Motif* motif() { return &comp_motif; };
  void kill_motif();//useful for some copy operations when comp_sites not needed
  bool ready() const {return comp_ready;};
};

CompareACESites read_motif(const string fname, const int mot);
void cCompareACE(int argc, char *argv[]);

#endif
