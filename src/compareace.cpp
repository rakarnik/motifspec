//Copyright 1999 President and Fellows of Harvard University
//compareace.cpp

#include "compareace.h"

CompareACESites::CompareACESites() {
	comp_ready = false;
}

CompareACESites::CompareACESites(const Motif& m): 
comp_motif(m) 
{
  comp_ready=true;
	int i,j,k;
  int min_num=6;
  int *b1=new int[min_num];
  double *b1v=new double[min_num];
  for(i=0;i<6;i++){
    b1[i]=-1;
    b1v[i]=0.0;
  }
  comp_fmsize = (m.width() + 2 * m.ncols()) * m.get_depth();
  comp_fm1 = new double[comp_fmsize];
  comp_fm2 = new double[comp_fmsize];//2=reverse complement of 1
  m.freq_matrix_extended(comp_fm1);
  //assuming nacgtn here
  for(i = 0; i < comp_fmsize; i += 6){
    comp_fm1[i+1] += (comp_fm1[i]+comp_fm1[i+5])/4.0;
    comp_fm1[i+2] += (comp_fm1[i]+comp_fm1[i+5])/4.0;
    comp_fm1[i+3] += (comp_fm1[i]+comp_fm1[i+5])/4.0;
    comp_fm1[i+4] += (comp_fm1[i]+comp_fm1[i+5])/4.0;
  } 
  for(i = 0; i < comp_fmsize; i++){
    comp_fm2[i] = comp_fm1[comp_fmsize - 1 - i];
  }
  
  double i1;
  double lg2=log(2.0);
  for(i = 0; i < comp_fmsize; i += 6) {
    i1 = 2.0;
    for(j = 1;j <= 4; j++) {
      if(comp_fm1[i + j] > 0.0) i1 += comp_fm1[i + j] *log(comp_fm1[i + j])/lg2;
    }
    for(j = 0; j < min_num; j++){
      if (i1 > b1v[j]){
				for(k = min_num - 1; k > j; k--) {
					b1v[k] = b1v[k - 1];
					b1[k] = b1[k - 1];
				}
				b1v[j] = i1;
				b1[j] = i;
				break;
      }
    }
  }
	
  comp_l1 = comp_fmsize;
  comp_r1 = -1;
  for(i = 0; i < min_num; i++){
    if(b1[i] > comp_r1) comp_r1 = b1[i];
    if(b1[i] < comp_l1) comp_l1 = b1[i];
  }
  
  //also consider ties to the worst considered position
  for(i = comp_l1 - 6; i >= 0;i -= 6) {
    i1 = 2.0;
    for(j = 1; j <= 4; j++){
      if(comp_fm1[i + j] > 0.0) i1 += comp_fm1[i + j] * log(comp_fm1[i + j])/lg2;
    }
    if(fabs(i1 - b1v[min_num - 1]) < .0001) comp_l1 = i;
  }
  for(i = comp_r1 + 6; i < comp_fmsize;i += 6){
    i1 = 2.0;
    for(j = 1;j <= 4; j++){
      if(comp_fm1[i + j] > 0.0) i1 += comp_fm1[i + j] * log(comp_fm1[i + j])/lg2;
    }
    if(fabs(i1 - b1v[min_num - 1]) < .0001) comp_r1 = i;
  }
	
  comp_r2 = comp_fmsize - 6 - comp_l1;
  comp_l2 = comp_fmsize - 6 - comp_r1;
	
  delete [] b1;
  delete [] b1v;
}

CompareACESites::CompareACESites(const CompareACESites& s):comp_motif(s.comp_motif) {
  comp_ready = false;
  if(!s.comp_ready) return;
  comp_ready = true;
  comp_fmsize = s.comp_fmsize;
  comp_l1 = s.comp_l1;
  comp_r1 = s.comp_r1;
  comp_l2 = s.comp_l2;
  comp_r2 = s.comp_r2;
  comp_fm1 = new double[comp_fmsize];
  comp_fm2 = new double[comp_fmsize];
  for(int i = 0; i < comp_fmsize; i++){
    comp_fm1[i] = s.comp_fm1[i];
    comp_fm2[i] = s.comp_fm2[i];
  }
}

CompareACESites::~CompareACESites(){
  if(comp_ready){
    delete [] comp_fm1;
    delete [] comp_fm2;
  }
}

CompareACESites& CompareACESites::operator= (const CompareACESites& s){
  if(this!=&s&&s.comp_ready){
    if(comp_ready){
      delete [] comp_fm1;
      delete [] comp_fm2;
    }
    comp_ready = true;
    comp_motif = s.comp_motif;
    comp_fmsize = s.comp_fmsize;
    comp_l1 = s.comp_l1;
    comp_r1 = s.comp_r1;
    comp_l2 = s.comp_l2;
    comp_r2 = s.comp_r2;
    comp_fm1 = new double[comp_fmsize];
    comp_fm2 = new double[comp_fmsize];
    for(int i = 0; i < comp_fmsize; i++){
      comp_fm1[i]=s.comp_fm1[i];
      comp_fm2[i]=s.comp_fm2[i];
    }
  }
  return *this;
}

double CompareACESites::compare(const CompareACESites &c){
	
  int i,j;
  int ref,ref1,ref2;
  int ref1rc,ref2rc;
  ref1=max(0,comp_l1-c.comp_r1+c.comp_l1);
  ref1=max(ref1,comp_r1-c.comp_fmsize+6+c.comp_l1);
  ref2=min(comp_fmsize-6-c.comp_r1+c.comp_l1,comp_r1);
  ref2=min(ref2,c.comp_l1+comp_l1);
  ref1rc=max(0,comp_l2-c.comp_r1+c.comp_l1);
  ref1rc=max(ref1rc,comp_r2-c.comp_fmsize+6+c.comp_l1);
  ref2rc=min(comp_fmsize-6-c.comp_r1+c.comp_l1,comp_r2);
  ref2rc=min(ref2rc,c.comp_l1+comp_l2);
  //these limits put ref frame within array bounds and require at least one overlap pos
  //ref1,ref2, etc. refer to the position of the second seq aligned with c.comp_l1
	
  double best_coeff=-1.1, coeff;
  double avg_1,avg_2,avg_12,avg_1sq,avg_2sq,dp;
  int x1,y1,x2,y2;
	
  for(ref=ref1;ref<=ref2;ref+=6){
    x1=c.comp_l1,y1=c.comp_r1,x2=ref,y2=ref+c.comp_r1-c.comp_l1;
    if(x2>comp_l1){
      x1-=(x2-comp_l1);
      x2-=(x2-comp_l1);
    }
    if(y2<comp_r1){
      y1+=(comp_r1-y2);
      y2+=(comp_r1-y2);
    }
    //x1..y1 in seq1, x2..y2 in seq2
    dp=4.0*((y1-x1)/6+1);
    avg_1=avg_2=avg_12=avg_1sq=avg_2sq=0.0;
    for(i=0;i<=(y1-x1);i+=6){
      for(j=1;j<=4;j++){
				avg_1+=c.comp_fm1[x1+j+i];
				avg_2+=comp_fm1[x2+j+i];
				avg_12+=c.comp_fm1[x1+j+i] * comp_fm1[x2+j+i];
				avg_1sq+=c.comp_fm1[x1+j+i] * c.comp_fm1[x1+j+i];
				avg_2sq+=comp_fm1[x2+j+i] * comp_fm1[x2+j+i];
      }
    }
    avg_1/=dp;
    avg_2/=dp;
    avg_12/=dp;
    avg_1sq/=dp;
    avg_2sq/=dp;
    coeff=(avg_12-avg_1*avg_2)/sqrt((avg_1sq-avg_1*avg_1)*(avg_2sq-avg_2*avg_2));
    if(coeff>best_coeff) {
      best_coeff=coeff;
    }
  }
	
  for(ref=ref1rc;ref<=ref2rc;ref+=6){
    x1=c.comp_l1,y1=c.comp_r1,x2=ref,y2=ref+c.comp_r1-c.comp_l1;
    if(x2>comp_l2){
      x1-=(x2-comp_l2);
      x2-=(x2-comp_l2);
    }
    if(y2<comp_r2){
      y1+=(comp_r2-y2);
      y2+=(comp_r2-y2);
    }
    //x1..y1 in seq1, x2..y2 in seq2
    dp=4.0*((y1-x1)/6+1);
    avg_1=avg_2=avg_12=avg_1sq=avg_2sq=0.0;
    for(i=0;i<=(y1-x1);i+=6){
      for(j=1;j<=4;j++){
				avg_1+=c.comp_fm1[x1+j+i];
				avg_2+=comp_fm2[x2+j+i];
				avg_12+=c.comp_fm1[x1+j+i] * comp_fm2[x2+j+i];
				avg_1sq+=c.comp_fm1[x1+j+i] * c.comp_fm1[x1+j+i];
				avg_2sq+=comp_fm2[x2+j+i] * comp_fm2[x2+j+i];
      }
    }
    avg_1/=dp;
    avg_2/=dp;
    avg_12/=dp;
    avg_1sq/=dp;
    avg_2sq/=dp;
    coeff=(avg_12-avg_1*avg_2)/sqrt((avg_1sq-avg_1*avg_1)*(avg_2sq-avg_2*avg_2));
    if(coeff>best_coeff) {
      best_coeff=coeff;
    }
  }
  return best_coeff;
}

void CompareACESites::kill_motif(){
  comp_motif.destroy();
}
