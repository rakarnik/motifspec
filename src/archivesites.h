//Copyright 1999 President and Fellows of Harvard University
//archivesites.h

#ifndef _archivesites
#define _archivesites
#include "standard.h"
#include "compareace.h"

class ArchiveSites{
  double arch_map_cutoff;
  double arch_sim_cutoff;
  int arch_min_visits;
  const Seqset& arch_seqset;
	vector<CompareACESites> arch_sites;

 public:
  ArchiveSites(Seqset& seq, double map_cut, double sim_cut);
	int nmots() const { return arch_sites.size(); }
  bool check_motif(const Motif& m);               // Returns true if no better motif, false otherwise
  bool consider_motif(const Motif& m);            // Returns true if motif was added, false otherwise
  Motif* return_best(const int i=0);
	void clear();
	void read(istream& archin);
	void write(ostream& archout);
};


#endif
