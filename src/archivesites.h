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
	vector<CompareACESites> arch_sites;
  const Seqset& arch_seqset;

 public:
  ArchiveSites(Sites& s, Seqset& seq, double map_cut, double sim_cut);
	int motifcount() const { return arch_sites.size(); }
  bool check_motif(const Sites& s);               // Returns true if better motif not seen, false otherwise
  bool consider_motif(const Sites& s, bool fnl=true);   // Returns true if motif was added, false otherwise
  Sites* return_best(const int i=0);
	void clear();
	void read(istream& archin);
	void write(ostream& archout);
};


#endif
