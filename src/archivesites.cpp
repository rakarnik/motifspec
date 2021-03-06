//Copyright 1999 President and Fellows of Harvard University
//archivesites.cpp

#include "archivesites.h"

ArchiveSites::ArchiveSites(const Seqset& seq, const double sim_cut, const int maxm,
		const vector<double>& p, const vector<double>& b) : 
seqset(seq),
mc(),
archive(0, Motif(seq, 12, p, b)),
sim_cutoff(sim_cut),
max_motifs(maxm),
pseudo(p),
backfreq(b),
min_visits(3) {
}

bool ArchiveSites::check_motif(const Motif& m) {
	vector<Motif>::iterator iter = archive.begin();
	float cmp1, cmp2;
	Motif rm(m);
	rm.flip_sites();
	for(; iter != archive.end() && m.get_motif_score() <= 0.9 * iter->get_motif_score(); ++iter){
		cmp1 = mc.compare(*iter, m);
		cmp2 = mc.compare(*iter, rm);
		if((cmp1 >= 0.95 || cmp2 >= 0.95) && iter->get_dejavu() >= min_visits)
			return false;
	}
	return true;
}

bool ArchiveSites::consider_motif(const Motif& m) {
	if(m.get_motif_score() < 1) return false;
	float cmp1, cmp2;
	Motif rm(m);
	rm.flip_sites();
	
	// Check if similar to better motif.
	// If so, increment dejavu for better motif and return false
	int motnum = 0;
	vector<Motif>::iterator iter = archive.begin();
	for(; iter != archive.end() && m.get_motif_score() <= iter->get_motif_score(); ++iter) {
		cmp1 = mc.compare(*iter, m);
		cmp2 = mc.compare(*iter, rm);
		// cerr << "Comparing with motif " << motnum << ", scores were " << cmp1 << " and " << cmp2 << '\n';
		if(cmp1 >= sim_cutoff || cmp2 >= sim_cutoff) {
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
		cmp1 = mc.compare(*iter, m);
		cmp2 = mc.compare(*iter, rm);
		// cerr << "Comparing with motif " << motnum << ", scores were " << cmp1 << " and " << cmp2 << '\n';
		if(cmp1 >= sim_cutoff || cmp2 >= sim_cutoff) {
	 		iter = archive.erase(iter);
			m1.inc_dejavu();
			delcount++;
		} else {
			++iter;
		}
		motnum++;
	}
	// cerr << "Deleted " << delcount << " motifs\n";
	
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
	archive.clear();
}

void ArchiveSites::read(istream& archin) {
	char line[200];
	while(archin.getline(line, 200)) {
		if(strstr(line, "Motif")) {
			Motif m(seqset, 12, pseudo, backfreq);
			m.read(archin);
			archive.push_back(m);
		}
	}
}

void ArchiveSites::write(ostream& archout) {
	int i = 1;
	vector<Motif>::iterator iter = archive.begin();
	for(; i <= max_motifs && iter != archive.end(); ++iter) {
		if(iter->get_motif_score() > 1) {
			archout << "Motif " << i << "\n";
			iter->write(archout);
		}
		i++;
	}
}
