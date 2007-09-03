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
  CompareACESites *arch_sites;
  double *arch_score;
  int *arch_dejavu;
  Seqset* arch_seqset;

 public:
 	ArchiveSites();
  ArchiveSites(const Sites& s, Seqset& seq, int max, double map_cut, double sim_cut);
	~ArchiveSites();
	void init(const Sites& s, Seqset& seq, int max, double map_cut, double sim_cut);
  double check_motif(const Sites& s, double sc);
  double consider_motif(const Sites& s, double sc, bool fnl=true);
  double return_best(Sites& s, int i=0);
  
};


#endif
