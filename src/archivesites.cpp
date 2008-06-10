//Copyright 1999 President and Fellows of Harvard University
//archivesites.cpp

#include "archivesites.h"

ArchiveSites::ArchiveSites(Sites& s, Seqset& seq, double map_cut, double sim_cut) : 
arch_seqset(seq) {
  arch_map_cutoff = 0.0;
  arch_sim_cutoff = sim_cut;
  arch_min_visits = 3;
}

bool ArchiveSites::check_motif(const Sites& s) {
  double cmp;
	CompareACESites c;
  c.init(s, arch_seqset);
  for(int i = 0; i < arch_sites.size(); i++){
    if(s.get_map() <= arch_sites[i].sites()->get_map()){
      cmp = c.compare(arch_sites[i]);
			if(cmp > arch_sim_cutoff && arch_sites[i].sites()->get_dejavu() >= arch_min_visits) {
				return false;
      }
    }
  }
  return true;
}

bool ArchiveSites::consider_motif(const Sites& s, bool fnl) {
  bool ret = false;
	//increment dejavu only for final motif
  //fnl=false assumes that no better motif is similar, so always add
  if(s.get_map() <= arch_map_cutoff) return false;
  CompareACESites c;
  c.init(s, arch_seqset);
  double cmp;

	// Check if similar to better motif.
	// If so, increment dejavu for better motif and return false
	for(int i = 0; i < arch_sites.size(); i++) {
		if(s.get_map() <= arch_sites[i].sites()->get_map()) {
			cmp = c.compare(arch_sites[i]);
			if(cmp > arch_sim_cutoff) {
				arch_sites[i].sites()->inc_dejavu();
				return false;
			}
		}
	}
	
	// There are no better motifs similar to this one, so we add
	// Step 1: Delete similar motifs with lower scores and shift up
	vector<CompareACESites> new_arch_sites;
	for(int i = 0; i < arch_sites.size(); i++) {
		if(s.get_map() <= arch_sites[i].sites()->get_map())
			new_arch_sites.push_back(arch_sites[i]);
		else {
			cmp = c.compare(arch_sites[i]);
			if(cmp <= arch_sim_cutoff) {
				new_arch_sites.push_back(arch_sites[i]);
			}
		}
	}
	arch_sites.swap(new_arch_sites);
	// Step 2: Add the new motif at the correct position by score
	vector<CompareACESites>::iterator iter;
	for(iter = arch_sites.begin(); iter != arch_sites.end(); iter++) {
		if(iter->sites()->get_map() < s.get_map()) break;
	}
	arch_sites.insert(iter, c);
	return true;
}

Sites* ArchiveSites::return_best(const int i) {
  return arch_sites[i].sites();
}

void ArchiveSites::clear() {
	arch_sites.clear();
}

void ArchiveSites::read(istream& archin) {
	int arch_num = 0;
	char line[200];
	while(archin.getline(line, 200)) {
		if(strstr(line, "Motif")) {
			Sites s(arch_seqset);
			s.read(archin);
			CompareACESites c;
			c.init(s, arch_seqset);
			arch_sites.push_back(c);
		}
	}
}

void ArchiveSites::write(ostream& archout) {
	for(int i = 0; i < arch_sites.size(); i++) {
		if(arch_sites[i].sites()->get_map() > arch_map_cutoff) {
			archout << "Motif " << i + 1 << endl;
			arch_sites[i].sites()->write(arch_seqset, archout);
		}
	}
}

