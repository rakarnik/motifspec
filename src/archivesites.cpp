//Copyright 1999 President and Fellows of Harvard University
//archivesites.cpp

#include "archivesites.h"

ArchiveSites::ArchiveSites(Seqset& seq, double sim_cut, double* pseudo) : 
arch_seqset(seq),
archive(0),
arch_sim_cutoff(sim_cut),
pseudo(pseudo),
arch_min_visits(3) {
}

bool ArchiveSites::check_motif(const Motif& m) {
	float cmp;
	vector<Motif>::iterator iter = archive.begin();
  for(; iter != archive.end() && m.get_motif_score() <= iter->get_motif_score(); ++iter){
    cmp = iter->compare(m);
		if(cmp > arch_sim_cutoff && iter->get_dejavu() >= arch_min_visits)
			return false;
  }
  return true;
}

bool ArchiveSites::consider_motif(const Motif& m) {
  if(m.get_motif_score() < 1) return false;
	float cmp;

	// Check if similar to better motif.
	// If so, increment dejavu for better motif and return false
	int motnum = 0;
	vector<Motif>::iterator iter = archive.begin();
	for(; iter != archive.end() && m.get_motif_score() <= iter->get_motif_score(); ++iter) {
		cmp = iter->compare(m);
		cerr << "Comparing with motif " << motnum << ", score was " << cmp << "\n";
		if(cmp > arch_sim_cutoff) {
			iter->inc_dejavu();
			return false;
		}
		motnum++;
	}
	
	Motif m1(m);

	// There are no better motifs similar to this one, so we add
	// Step 1: Delete similar motifs with lower scores
	int delcount = 0;
	while(iter != archive.end()) {
		assert(m.get_motif_score() > iter->get_motif_score());
		cmp = iter->compare(m1);
		cerr << "Comparing with motif " << motnum << ", score was " << cmp << "\n";
		if(cmp > arch_sim_cutoff) {
	 		iter = archive.erase(iter);
			m1.inc_dejavu();
			delcount++;
		} else {
			++iter;
		}
		motnum++;
	}
	cerr << "Deleted " << delcount << " similar motifs\n";

	
	// Step 2: Add the new motif at the correct position by score
	for(iter = archive.begin(); iter != archive.end(); ++iter)
		if(iter->get_motif_score() < m1.get_motif_score()) break;
	archive.insert(iter, m1);
	return true;
}

Motif* ArchiveSites::return_best(const int i) {
  return &archive[i];
}

void ArchiveSites::clear() {
	vector<Motif>().swap(archive);
}

void ArchiveSites::read(istream& archin) {
	char line[200];
	while(archin.getline(line, 200)) {
		if(strstr(line, "Motif")) {
			Motif m(arch_seqset, 12, pseudo);
			m.read(archin);
			archive.push_back(m);
		}
	}
}

void ArchiveSites::write(ostream& archout) {
	int i = 0;
	vector<Motif>::iterator iter = archive.begin();
	for(; iter != archive.end(); ++iter) {
		if(iter->get_motif_score() > 1) {
			archout << "Motif " << i + 1 << "\n";
			iter->write(archout);
			i++;
		}
	}
}

