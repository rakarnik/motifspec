//Copyright 1999 President and Fellows of Harvard University
//archivesites.h

#ifndef _archivesites
#define _archivesites
#include "standard.h"
#include "compareace.h"

class ArchiveSites{
	int arch_num;
  int arch_max_num;
  double arch_map_cutoff;
  double arch_sim_cutoff;
  int arch_min_visits;
  Sites* arch_base_site;
	CompareACESites *arch_sites;
  Seqset* arch_seqset;

 public:
 	ArchiveSites();
  ArchiveSites(Sites& s, Seqset& seq, int max, double map_cut, double sim_cut);
	~ArchiveSites();
	void init(Sites& s, Seqset& seq, int max, double map_cut, double sim_cut);
	void clear();                                         // Remove all motifs from the archive
  bool check_motif(const Sites& s);                     // Returns true if better motif not seen, false otherwise
  bool consider_motif(const Sites& s, bool fnl=true);   // Returns true if motif was added, false otherwise
  double return_best(Sites& s, int i=0);
	int motifcount() const { return arch_num; };
	void read(istream& archin);
	void write(ostream& archout);
};


#endif
