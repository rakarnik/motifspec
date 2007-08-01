//Copyright 1999 President and Fellows of Harvard University
//archivesites.cpp

#include "archivesites.h"

ArchiveSites::ArchiveSites(const Sites& s, const Seqset& seq, int max, double map_cut, double sim_cut):arch_seqset(seq){
  arch_max_num=max;
  arch_num=0;
  arch_map_cutoff=0.0;
  arch_sim_cutoff=sim_cut;
  arch_min_visits=3;
  arch_sites=new CompareACESites[max];
  arch_score=new double[max];
  arch_dejavu=new int[max];
  int i;
  for(i=0;i<max;i++){
    arch_sites[i].init(s,seq);
    arch_score[i]=map_cut;
    arch_dejavu[i]=0;
  }
}

ArchiveSites::~ArchiveSites(){
  delete [] arch_sites;
  delete [] arch_score;
  delete [] arch_dejavu;
}


double ArchiveSites::check_motif(const Sites& s, double sc){
  int i;
  CompareACESites c;
  c.init(s,arch_seqset);
  double cmp, ret=-1.0;
  bool csd=true;
  for(i=0;i<arch_max_num;i++){
    if(sc<=arch_score[i]){
      //if very similar to better motif with min visits, return that score
      cmp=c.compare(arch_sites[i]);
      if(cmp>arch_sim_cutoff) csd=false;
      if(arch_dejavu[i]>=arch_min_visits){
	if(cmp>ret) ret=cmp;
      }
    }
    else break;
  }
  if(csd) cmp=consider_motif(s,sc,false);
  //returns best similarity to min_visited motif with greater map score
  return ret;
}


double ArchiveSites::consider_motif(const Sites& s, double sc, bool fnl){
  //increment dejavu only for final motif
  //fnl=false assumes that no better motif is similar, so always add
  int i,j,k,m;
  if(sc<=arch_map_cutoff) return -1.0;
  CompareACESites c;
  c.init(s,arch_seqset);
  double cmp;
  for(i=0;i<arch_max_num;i++){
    if(sc<=arch_score[i]){
      if(!fnl) continue;
      cmp=c.compare(arch_sites[i]);
      //if very similar to better motif, ++dejavu and return that score
      if(cmp>arch_sim_cutoff) {
	if(fnl) arch_dejavu[i]++;
	return cmp;
      }
    }
    else {
      //first delete similar motifs with lower scores and shift up
      for(j=i;j<arch_max_num;j++){
	if(arch_score[j]==arch_map_cutoff) break;
	cmp=c.compare(arch_sites[j]);
	if(cmp>arch_sim_cutoff) {
	  arch_score[j]=arch_map_cutoff;
	}
      }
      for(k=i,m=i+1;k<j&&m<j;){
	if(arch_score[k]==arch_map_cutoff){
	  if(arch_score[m]>arch_map_cutoff){
	    arch_sites[k]=arch_sites[m];
	    arch_score[k]=arch_score[m];
	    arch_dejavu[k]=arch_dejavu[m];
	    arch_score[m]=arch_map_cutoff;
	  }
	  else m++;
	}
	else{ 
	  k++;
	  m=k+1;
	}
      }
      //then insert the new motif and shift the list down
      for(j=arch_max_num-1;j>=i+1;j--){
	arch_sites[j]=arch_sites[j-1];
	arch_score[j]=arch_score[j-1];
	arch_dejavu[j]=arch_dejavu[j-1];
      }
      arch_sites[i]=c;
      arch_score[i]=sc;
      if(fnl) arch_dejavu[i]=1;
      else arch_dejavu[i]=0;
      break;
    }
  }
  return 0.0;
}

double ArchiveSites::return_best(Sites& s, int i){
  s = *(arch_sites[i].sites());
  //  cerr<<i<<'\t'<<arch_dejavu[i]<<'\n';
  return arch_score[i];
}


