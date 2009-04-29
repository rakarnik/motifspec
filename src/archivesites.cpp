//Copyright 1999 President and Fellows of Harvard University
//archivesites.cpp

#include "archivesites.h"

ArchiveSites::ArchiveSites(Seqset& seq, double map_cut, double sim_cut) : 
arch_seqset(seq),
arch_comp(0) {
  arch_map_cutoff = map_cut;
  arch_sim_cutoff = sim_cut;
  arch_min_visits = 3;
}

bool ArchiveSites::check_motif(const Motif& m) {
  double cmp;
	CompareACESites c(m);
	vector<CompareACESites>::iterator iter = arch_comp.begin();
  for(; iter != arch_comp.end() && m.get_spec() <= iter->motif()->get_spec(); ++iter){
    cmp = c.compare(*iter);
		if(cmp > arch_sim_cutoff && iter->motif()->get_dejavu() >= arch_min_visits)
			return false;
  }
  return true;
}

bool ArchiveSites::consider_motif(const Motif& m) {
  if(m.get_spec() < 1) return false; 
  
	CompareACESites c(m);
  double cmp;

	// Check if similar to better motif.
	// If so, increment dejavu for better motif and return false
	vector<CompareACESites>::iterator iter = arch_comp.begin();
	for(; iter != arch_comp.end() && m.get_spec() <= iter->motif()->get_spec(); ++iter) {
		cmp = c.compare(*iter);
		if(cmp > arch_sim_cutoff) {
			iter->motif()->inc_dejavu();
			return false;
		}
	}
	
	// There are no better motifs similar to this one, so we add
	// Step 1: Delete similar motifs with lower scores
	while(iter != arch_comp.end()) {
		assert(m.get_spec() > iter->motif()->get_spec());
		cmp = c.compare(*iter);
		if(cmp > arch_sim_cutoff)
			iter = arch_comp.erase(iter);
		else
			++iter;
	}
	
	// Step 2: Add the new motif at the correct position by score
	for(iter = arch_comp.begin(); iter != arch_comp.end(); ++iter)
		if(iter->motif()->get_spec() < m.get_spec()) break;
	arch_comp.insert(iter, c);
	return true;
}

Motif* ArchiveSites::return_best(const int i) {
  return arch_comp[i].motif();
}

void ArchiveSites::clear() {
	vector<CompareACESites>().swap(arch_comp);
}

void ArchiveSites::read(istream& archin) {
	char line[200];
	while(archin.getline(line, 200)) {
		if(strstr(line, "Motif")) {
			Motif m(arch_seqset, 12);
			m.read(archin);
			CompareACESites c(m);
			arch_comp.push_back(c);
		}
	}
}

void ArchiveSites::write(ostream& archout) {
	int i = 0;
	vector<CompareACESites>::iterator iter = arch_comp.begin();
	for(; iter != arch_comp.end(); ++iter) {
		if(iter->motif()->get_spec() > 1) {
			archout << "Motif " << i + 1 << "\n";
			iter->motif()->write(archout);
			i++;
		}
	}
}

