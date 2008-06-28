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
	CompareACESites c(s, arch_seqset);
  for(int i = 0; i < arch_sites.size(); i++){
    if(s.get_spec() <= arch_sites[i].sites()->get_spec()){
      cmp = c.compare(arch_sites[i]);
			if(cmp > arch_sim_cutoff && arch_sites[i].sites()->get_dejavu() >= arch_min_visits) {
				return false;
      }
    }
  }
  return true;
}

bool ArchiveSites::consider_motif(const Sites& s) {
  if(s.get_spec() < 1) return false; 
  
	CompareACESites c(s, arch_seqset);
  double cmp;

	// Check if similar to better motif.
	// If so, increment dejavu for better motif and return false
	vector<CompareACESites>::iterator iter;
	for(iter = arch_sites.begin(); iter != arch_sites.end(); ++iter) {
		if(s.get_spec() <= iter->sites()->get_spec()) {
			cmp = c.compare(*iter);
			if(cmp > arch_sim_cutoff) {
				iter->sites()->inc_dejavu();
				return false;
			}
		} else {
			break;
		}
	}
	
	// There are no better motifs similar to this one, so we add
	// Step 1: Delete similar motifs with lower scores
	while(iter != arch_sites.end()) {
		assert(s.get_spec() > iter->sites()->get_spec());
		cmp = c.compare(*iter);
		if(cmp > arch_sim_cutoff)
			iter = arch_sites.erase(iter);
		else
			++iter;
	}
	
	// Step 2: Add the new motif at the correct position by score
	for(iter = arch_sites.begin(); iter != arch_sites.end(); ++iter)
		if(iter->sites()->get_spec() < s.get_spec()) break;
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
	char line[200];
	while(archin.getline(line, 200)) {
		if(strstr(line, "Motif")) {
			Sites s(arch_seqset, 12, 5 * 12);
			s.read(archin);
			CompareACESites c(s, arch_seqset);
			arch_sites.push_back(c);
		}
	}
}

void ArchiveSites::write(ostream& archout) {
	for(int i = 0; i < arch_sites.size(); i++) {
		if(arch_sites[i].sites()->get_spec() > 1) {
			archout << "Motif " << i + 1 << endl;
			arch_sites[i].sites()->write(archout);
		}
	}
}

