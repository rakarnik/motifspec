//Copyright 1999 President and Fellows of Harvard University
//archivesites.cpp

#include "archivesites.h"

ArchiveSites::ArchiveSites() {
}

ArchiveSites::ArchiveSites(Sites& s, Seqset& seq, int max, double map_cut, double sim_cut) : arch_seqset(&seq) {
	arch_base_site = &s;
	arch_max_num = max;
  arch_num = 0;
  arch_map_cutoff = 0.0;
  arch_sim_cutoff = sim_cut;
  arch_min_visits = 3;
	arch_sites = new CompareACESites[max];
  for(int i = 0; i < max; i++){
    arch_sites[i].init(*arch_base_site, *arch_seqset);
  }
}

ArchiveSites::~ArchiveSites(){
	delete [] arch_sites;
}

void ArchiveSites::init(Sites& s, Seqset& seq, int max, double map_cut, double sim_cut) {
	arch_base_site = &s;
	arch_seqset = &seq;
	arch_max_num = max;
  arch_num = 0;
  arch_map_cutoff = 0.0;
  arch_sim_cutoff = sim_cut;
  arch_min_visits = 3;
  arch_sites = new CompareACESites[max];
  for(int i = 0; i < max; i++){
    arch_sites[i].init(*arch_base_site, *arch_seqset);
  }
}

void ArchiveSites::clear() {
	delete [] arch_sites;
	arch_sites = new CompareACESites[arch_max_num];
  for(int i = 0; i < arch_max_num; i++)
    arch_sites[i].init(*arch_base_site, *arch_seqset);
	arch_num = 0;
}

bool ArchiveSites::check_motif(const Sites& s){
  int i;
  CompareACESites c;
  c.init(s, *arch_seqset);
  for(i = 0; i < arch_num; i++){
    if(s.get_map() <= arch_sites[i].sites()->get_map()){
      if(c.compare(arch_sites[i]) > arch_sim_cutoff && arch_sites[i].sites()->get_dejavu() >= arch_min_visits) {
				return false;
      }
    }
  }
  return true;
}


bool ArchiveSites::consider_motif(const Sites& s, bool fnl){
  bool ret = false;
	//increment dejavu only for final motif
  //fnl=false assumes that no better motif is similar, so always add
  if(s.get_map() <= arch_map_cutoff) return false;
  CompareACESites c;
  c.init(s, *arch_seqset);
	int j;
  double cmp;
  for(int i = 0; i < arch_num; i++){
    if(s.get_map() <= arch_sites[i].sites()->get_map()){
      if(! fnl) continue;
      cmp = c.compare(arch_sites[i]);
      //if very similar to better motif, ++dejavu and return false
      if(cmp > arch_sim_cutoff) {
				if(fnl) arch_sites[i].sites()->inc_dejavu();
				break;
      }
    } else {
      //first delete similar motifs with lower scores and shift up
      for(j = i; j < arch_max_num; j++) {
				if(arch_sites[j].sites()->get_map() == arch_map_cutoff) break;
				cmp = c.compare(arch_sites[j]);
				if(cmp > arch_sim_cutoff) {
					arch_sites[j].sites()->set_map(arch_map_cutoff);
				}
			}
      for(int k = i, m = i + 1; k < j && m < j; )	{
				if(arch_sites[k].sites()->get_map() == arch_map_cutoff){
					if(arch_sites[m].sites()->get_map() > arch_map_cutoff){
						arch_sites[k] = arch_sites[m];
						arch_sites[m].sites()->set_map(arch_map_cutoff);
					} else m++;
				} else { 
					k++;
					m=k+1;
				}
			}
      //then insert the new motif and shift the list down
      for(j = arch_num - 1; j >= i + 1; j--) {
				arch_sites[j] = arch_sites[j - 1];
      }
      arch_sites[i] = c;
      if(fnl) arch_sites[i].sites()->set_dejavu(1);
      else arch_sites[i].sites()->set_dejavu(0);
      ret = true;
			arch_num++;
			break;
    }
  }
	return ret;
}

double ArchiveSites::return_best(Sites& s, int i){
  s = *(arch_sites[i].sites());
  //  cerr<<i<<'\t'<<arch_dejavu[i]<<'\n';
  return arch_sites[i].sites()->get_map();
}

void ArchiveSites::read(istream& archin) {
	arch_num = 0;
	char line[200];
	while(archin.getline(line, 200)) {
		if(strstr(line, "Motif")) {
			// Read list of sites in motif
			Sites s;
			s.init(*arch_seqset);
			s.read(archin);
			CompareACESites c;
			c.init(s, *arch_seqset);
			arch_sites[arch_num] = c;
			arch_num++;
		}
	}
}

void ArchiveSites::write(ostream& archout) {
	for(int i = 0; i < arch_num; i++) {
		archout << "Motif " << i + 1 << endl;
		arch_sites[i].sites()->write(*arch_seqset, archout);
	}
}

