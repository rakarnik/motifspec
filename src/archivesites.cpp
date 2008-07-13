//Copyright 1999 President and Fellows of Harvard University
//archivesites.cpp

#include "archivesites.h"

ArchiveSites::ArchiveSites(Seqset& seq, double map_cut, double sim_cut) : 
arch_seqset(seq),
arch_sites(0, Motif(seq)) {
  arch_map_cutoff = map_cut;
  arch_sim_cutoff = sim_cut;
  arch_min_visits = 3;
}

bool ArchiveSites::check_motif(const Motif& m) {
  double cmp;
	CompareACESites c(m);
  for(int i = 0; i < arch_sites.size(); i++){
    if(m.get_spec() <= arch_sites[i].motif()->get_spec()){
      cmp = c.compare(arch_sites[i]);
			if(cmp > arch_sim_cutoff && arch_sites[i].motif()->get_dejavu() >= arch_min_visits) {
				return false;
      }
    } else {
			break;
		}
  }
  return true;
}

bool ArchiveSites::consider_motif(const Motif& m) {
  if(m.get_spec() < 1) return false; 
  
	CompareACESites c(m);
  double cmp;

	// Check if similar to better motif.
	// If so, increment dejavu for better motif and return false
	vector<CompareACESites>::iterator iter;
	for(iter = arch_sites.begin(); iter != arch_sites.end(); ++iter) {
		if(m.get_spec() <= iter->motif()->get_spec()) {
			cmp = c.compare(*iter);
			if(cmp > arch_sim_cutoff) {
				iter->motif()->inc_dejavu();
				return false;
			}
		} else {
			break;
		}
	}
	
	// There are no better motifs similar to this one, so we add
	// Step 1: Delete similar motifs with lower scores
	while(iter != arch_sites.end()) {
		assert(m.get_spec() > iter->motif()->get_spec());
		cmp = c.compare(*iter);
		if(cmp > arch_sim_cutoff)
			iter = arch_sites.erase(iter);
		else
			++iter;
	}
	
	// Step 2: Add the new motif at the correct position by score
	for(iter = arch_sites.begin(); iter != arch_sites.end(); ++iter)
		if(iter->motif()->get_spec() < m.get_spec()) break;
	arch_sites.insert(iter, c);
	return true;
}

Motif* ArchiveSites::return_best(const int i) {
  return arch_sites[i].motif();
}

void ArchiveSites::clear() {
	arch_sites.clear();
}

void ArchiveSites::read(istream& archin) {
	char line[200];
	while(archin.getline(line, 200)) {
		if(strstr(line, "Motif")) {
			Motif m(arch_seqset, 12);
			m.read(archin);
			CompareACESites c(m);
			arch_sites.push_back(c);
		}
	}
}

void ArchiveSites::write(ostream& archout) {
	for(int i = 0; i < arch_sites.size(); i++) {
		if(arch_sites[i].motif()->get_spec() > 1) {
			archout << "Motif " << i + 1 << endl;
			arch_sites[i].motif()->write(archout);
		}
	}
}

