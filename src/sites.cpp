//Copyright 1999 President and Fellows of Harvard University
//sites.cpp

#include "sites.h"

Sites::Sites(const vector<string>& v, int nc, int mx, int dp){
  sites_num=0;
  sites_width=nc;
  sites_num_cols=nc;
  sites_depth=dp;
  sites_num_seqs=v.size();
  sites_len_seq=new int[sites_num_seqs];
  for(int i=0;i<v.size();i++){
    sites_len_seq[i]=v[i].length();
  }
  sites_max_num_sites=0;
  for(int i=0;i<sites_num_seqs;i++){
    sites_max_num_sites+=sites_len_seq[i]/mx;
  }
  sites_max_width=3*sites_width;
	allocate_mem();
  clear_sites();
}

Sites::Sites(const Sites& s){
  if(!s.sites_alloc) return;
  sites_num_seqs=s.sites_num_seqs;
  sites_len_seq=new int[sites_num_seqs];
  for(int i=0;i<sites_num_seqs;i++){
    sites_len_seq[i]=s.sites_len_seq[i];
  }
  sites_max_num_sites=s.sites_max_num_sites;
  sites_max_width=s.sites_max_width;
  allocate_mem();
  *this=s;
}

void Sites::sites_init(const Sites& s){
  if(!s.sites_alloc) return;
  sites_num_seqs=s.sites_num_seqs;
  sites_len_seq=new int[sites_num_seqs];
  for(int i=0;i<sites_num_seqs;i++){
    sites_len_seq[i]=s.sites_len_seq[i];
  }
  sites_max_num_sites=s.sites_max_num_sites;
  sites_max_width=s.sites_max_width;
  allocate_mem();
  *this=s;
}

void Sites::init(const vector<string>& v, int nc, int mx, int dp){
  sites_num = 0;
  sites_width = nc;
  sites_num_cols = nc;
  sites_depth = dp;
  sites_num_seqs = v.size();
  sites_len_seq = new int[sites_num_seqs];
  for(int i = 0; i < v.size(); i++){
		sites_len_seq[i] = v[i].length();
	}
  sites_max_num_sites = 0;
  for(int i = 0;i < sites_num_seqs; i++){
    sites_max_num_sites += sites_len_seq[i]/mx;
  }
  sites_max_width = 3*sites_width;
  allocate_mem();
	sites_num_seqs_with_sites = 0;
	for(int i = 0; i < v.size(); i++){
		sites_has_sites[i] = 0;
	}
  clear_sites();
	corr_cutoff = 0.70;
}

Sites& Sites::operator= (const Sites& s){
  if(this != &s && s.sites_alloc){
    //assume that the same Seqset is referred to, so ignore some things
    sites_num = s.sites_num;
    sites_width = s.sites_width;
    sites_num_cols = s.sites_num_cols;
    sites_depth = s.sites_depth;
    for(int i = 0; i < sites_num; i++){
      sites_chrom[i] = s.sites_chrom[i];
      sites_posit[i] = s.sites_posit[i];
      sites_strand[i] = s.sites_strand[i];
    }
		sites_num_seqs_with_sites = s.sites_num_seqs_with_sites;
		for(int i = 0; i < sites_num_seqs; i++) {
			sites_has_sites[i] = s.sites_has_sites[i];
		}
    for(int i = 0; i < sites_max_width; i++){
      sites_active_fwd[i] = s.sites_active_fwd[i];
    }
		corr_cutoff = s.corr_cutoff;
  }
  return *this;
}

void Sites::allocate_mem(){
  sites_alloc = true;
  sites_chrom = new int[sites_max_num_sites];
  sites_posit = new int[sites_max_num_sites];
  sites_strand = new bool[sites_max_num_sites];
	sites_has_sites = new int[sites_num_seqs];
  sites_active_fwd = new int[sites_max_width];
}

Sites::~Sites(){
	if(sites_alloc){
    delete [] sites_len_seq;
    delete [] sites_chrom;
    delete [] sites_posit;
    delete [] sites_strand;
    delete [] sites_active_fwd;
  }
}

void Sites::clear_sites(){
  sites_num = 0;
  sites_width=sites_num_cols;
  for(int i = 0; i < sites_num_cols - 1; i++){
    sites_active_fwd[i]=i+1;
  }
  sites_active_fwd[sites_num_cols-1] = sites_num_cols - 1;
	sites_num_seqs_with_sites = 0;
	for(int i = 0; i < sites_num_seqs; i++) {
		sites_has_sites[i] = false;
	}
}

void Sites::remove_all_sites(){
  sites_num = 0;
	sites_num_seqs_with_sites = 0;
	for(int i = 0; i < sites_num_seqs; i++) {
		sites_has_sites[i] = false;
	}
}

void Sites::destroy(){
  //for when the = operator can fail
  sites_num=0;
  sites_max_width=0;
}

bool Sites::is_open_site(const int c, const int p){
  for(int i=0;i<sites_num;i++){
    if(sites_chrom[i]==c){
      int pp=sites_posit[i];
      if(pp>p-sites_width&&pp<p+sites_width) return false;
    }
  }
  return true;
}

void Sites::add_site(const int c, const int p, const bool s){
	sites_chrom[sites_num]=c;
  sites_posit[sites_num]=p;
  sites_strand[sites_num]=s;
	if(sites_has_sites[c] == 0) sites_num_seqs_with_sites++;
	sites_has_sites[c]++;
	sites_num++;
}

void Sites::remove_site(const int c, const int p){
  int i;
	for(i = 0; i < sites_num; i++){
    if(sites_chrom[i] != c || sites_posit[i] != p) continue;
    break;
  }
  // i now equal to index of site to be removed
  sites_num--;
  sites_chrom[i] = sites_chrom[sites_num];
  sites_posit[i] = sites_posit[sites_num];
  sites_strand[i] = sites_strand[sites_num];
  
	if(i != sites_num) {               // site was found
		sites_has_sites[c]--;
		if(sites_has_sites[c] == 0) sites_num_seqs_with_sites--;
	}
}

void Sites::calc_freq_matrix(const Seqset& b, int *fm){
	//fm will have allocation for depth()*ncols()
  char** ss_seq = b.seq_ptr();
  for(int i = 0; i < sites_depth * sites_num_cols; i++){
		fm[i] = 0;
  }
  for(int i = 0; i < sites_num; i++){ //i = site number
    int c = sites_chrom[i];
    int p = sites_posit[i];
    bool s = sites_strand[i];
    int col, pos, matpos;
    if(s) {                              // forward strand
      matpos = 0;
			col = 0;
      pos = 0;
			for(int j = 0; j < sites_num_cols; j++){ //j = position number
				pos = p + col;
				assert(p >= 0 && p < b.len_seq(c));
				int seq = ss_seq[c][pos];
				fm[matpos + seq]++;
				col = sites_active_fwd[col];
				matpos += sites_depth;
      }
    } else {                             // reverse strand
      matpos = sites_depth - 1;
      col = 0;
      pos = 0;
			for(int j = 0; j < sites_num_cols; j++){
				pos = p + sites_width - 1 - col;
				assert(p >= 0 && p < b.len_seq(c));
				int seq = ss_seq[c][p + sites_width - 1 - col];
				fm[matpos - seq]++;
				col = sites_active_fwd[col];
				matpos += sites_depth;
      }
    }      
  }
}

bool Sites::column_freq(const int col, const Seqset& s, int *ret){
  char** ss_seq=s.seq_ptr();
  int i,j;
  for(i=0;i<sites_depth;i++) ret[i]=0;
  for(i=0;i<sites_num;i++){//i = site number
    int c=sites_chrom[i];
    int p=sites_posit[i];
    bool t=sites_strand[i];
    if(t){
      if( (p+col > s.len_seq(c)-1) || (p+col <0) ){
	return false;
      }
      int seq=ss_seq[c][p+col];
      ret[seq]++;
    }
    else{
      if((p+sites_width-1-col>s.len_seq(c)-1)||(p+sites_width-1-col<0)){
	return false;
      }
      int seq=ss_seq[c][p+sites_width-1-col];
      ret[sites_depth-seq-1]++;
    }
  }
  return true;
}

int Sites::remove_col(const int c) {
  int col = 0, nxt = 0, ret = 0, ns = 0;         //return number of removed column in new numbering
	bool found = false;
	if(c == 0) {                   // if the column to be removed is the first column
    ns = sites_active_fwd[0];
    ret = -ns;
    for(col = ns; col < sites_width - 1;){
      nxt = sites_active_fwd[col];
      sites_active_fwd[col - ns] = nxt - ns;
      col = nxt;
    }
    shift_sites(ns, 0);
    sites_width -= ns;
    sites_active_fwd[sites_width - 1] = sites_width - 1;
		found = true;
  }
  else if(c == (sites_width - 1)) {   // if the column to be removed is that last column
    ret = c;
    for(col = 0;;) {
      nxt = sites_active_fwd[col];
			if(nxt == (sites_width - 1)) {
				sites_width = col + 1;
				sites_active_fwd[col] = col;
				shift_sites(0, col - nxt);
				found = true;
				break;
      } else {
				col = nxt;
			}
		}
  } else {                            // somewhere in the middle, so find right insertion point
    ret = c;
		col = 0;
    for(col = 0;;) {
      nxt = sites_active_fwd[col];
			if(nxt == col) break;
      if(nxt == c) {
				sites_active_fwd[col] = sites_active_fwd[nxt];
				found = true;
				break;
      }
      col = nxt;
    }
  }
  sites_num_cols--;
	if (! found) {
		cerr << "remove_column called for column " << c << " but it was not found!" << endl; 
		abort();
	}
  return ret;
}

void Sites::add_col(const int c){
	int col, nxt, i;
	col = nxt = i = 0;
  if(c < 0){
    for(i=sites_width-1;i>=0;i--){
      sites_active_fwd[i-c]=sites_active_fwd[i]-c;
    }
    sites_active_fwd[0]=-c;
    sites_width+=-c;
    shift_sites(c,0);
  }
  else if(c<sites_width){
    for(col=0;;){
      nxt=sites_active_fwd[col];
      if(nxt>c){
				sites_active_fwd[col]=c;
				sites_active_fwd[c]=nxt;
				break;
      }
      else col=nxt;
    }
  }
  else{
		shift_sites(0,c-sites_width+1);
    sites_active_fwd[sites_width-1]=c;
    sites_active_fwd[c]=c;
    sites_width=c+1;
  }
  sites_num_cols++;
}

void Sites::flip_sites(){
  int i;
  for(i=0;i<sites_num;i++){
    sites_strand[i]=!(sites_strand[i]);
  }
  int *temp=new int[sites_width];
  int l=0,n=0;
  while(n!=(sites_width-1)){
    n=sites_active_fwd[l];
    temp[sites_width-1-n]=sites_width-1-l;
    l=n;
  }
  for(i=0;i<sites_width;i++){
    sites_active_fwd[i]=temp[i];
  }
  delete [] temp;
}

void Sites::shift_sites(const int l, const int r){
  // numbers for right movement of beg/end point of forward site
  // l: + for shorter/right
  // r: + for longer/right
  for(int i=0;i<sites_num;i++){
    int c=sites_chrom[i];
    int p=sites_posit[i];
    bool s=sites_strand[i];
    int newp;
    if(s) newp=p+l;
    else newp=p-r;
    sites_posit[i]=newp;
    //    sites_stat[c][newp]=sites_stat[c][p];
    //sites_stat[c][p]=0;
  }
}

int Sites::positions_available() const {
  int ret = 0;
  for(int i = 0; i < sites_num_seqs; i++){
    ret += sites_len_seq[i] - sites_width + 1;
  }
  return ret;
}

int Sites::positions_available(const bool* possible) const {
	int ret = 0;
	for(int i = 0; i < sites_num_seqs; i++) {
		if(possible[i])
			ret += sites_len_seq[i] - sites_width + 1;
	}
	return ret;
}

void Sites::columns_open(int &l, int &r){
  //input r/l are max values to be reduced
  //only works if site list is sorted
  int i;
  int c_prev=-1, p_prev;
  bool s_prev;
  for(i=0;i<sites_num;i++){
    int c=sites_chrom[i];
    int p=sites_posit[i];
    bool s=sites_strand[i];
    if(c==c_prev){
      int d=p-p_prev-sites_width;
      if(s==s_prev){
				if(l>d) l=d;
				if(r>d) r=d;
      }
      else{
				if(l>d/2) l=d/2;
				if(r>d/2) r=d/2;
      }
    }
    else{
      int f;
      if(c_prev!=-1){
				f=sites_len_seq[c_prev]-sites_width-p_prev;
				if(s_prev){
					if(r>f) r=f;
				}
				else{
					if(l>f) l=f;
				}
      }
      f=p;
      if(s){
				if(l>f) l=f;
      }
      else{
				if(r>f) r=f;
      }
    }
  }
}

void Sites::freq_matrix_extended(const Seqset& b, double *fm) const {
  //assumes fm of dimension (sites_width+2*sites_num_cols)*depth, fills in edges with Ns
  char** ss_seq=b.seq_ptr();
  int i,col,j,seq;
  int fm_size=(sites_width+2*sites_num_cols)*sites_depth;
  for(i=0;i<fm_size;i++) fm[i]=0.0;
  if(sites_num==0) return;
  for(i=0;i<sites_num;i++){//i = site number
    int c=sites_chrom[i];
    int p=sites_posit[i];
    bool t=sites_strand[i];
    for(j=0,col=-sites_num_cols;col<sites_width+sites_num_cols;col++,j+=sites_depth){
      if(t){
	if( (p+col > b.len_seq(c)-1) || (p+col <0) ) seq=0;
	else seq=ss_seq[c][p+col];
	fm[j+seq]+=1.0;
      }
      else{
	if((p+sites_width-1-col>b.len_seq(c)-1)||(p+sites_width-1-col<0)) seq=0;
	else seq=ss_seq[c][p+sites_width-1-col];
	fm[j+sites_depth-1-seq]+=1.0;
      }
    }
  }
  for(i=0;i<fm_size;i++) fm[i]/=(double)sites_num;
}

double CompareACE(const Sites &s1, const Sites &s2, const Seqset &t1, const Seqset &t2){

  //  if(s1.number()==0||s2.number()==0) return 0.0;
  int i,j,k;
  int min_num=6;
  int *b1=new int[min_num];
  int *b2=new int[min_num];
  int *b3=new int[min_num];
  double *b1v=new double[min_num];
  double *b2v=new double[min_num];
  double *b3v=new double[min_num];//3=reverse complement of 2
  for(i=0;i<6;i++){
    b1[i]=b2[i]=b3[i]=-1;
    b1v[i]=b2v[i]=b3v[i]=0.0;
  }
  int fm_size1=(s1.width()+2*s1.ncols())*s1.depth();
  int fm_size2=(s2.width()+2*s2.ncols())*s2.depth();
  double *f1=new double[fm_size1];
  double *f2=new double[fm_size2];
  double *f3=new double[fm_size2];
  s1.freq_matrix_extended(t1,f1);
  s2.freq_matrix_extended(t2,f2);
  //assuming nacgtn here
  for(i=0;i<fm_size1;i+=6){
    f1[i+1]+=(f1[i]+f1[i+5])/4.0;
    f1[i+2]+=(f1[i]+f1[i+5])/4.0;
    f1[i+3]+=(f1[i]+f1[i+5])/4.0;
    f1[i+4]+=(f1[i]+f1[i+5])/4.0;
  }
  for(i=0;i<fm_size2;i+=6){
    f2[i+1]+=(f2[i]+f2[i+5])/4.0;
    f2[i+2]+=(f2[i]+f2[i+5])/4.0;
    f2[i+3]+=(f2[i]+f2[i+5])/4.0;
    f2[i+4]+=(f2[i]+f2[i+5])/4.0;
  }
 
 for(i=0;i<fm_size2;i++){
    f3[i]=f2[fm_size2-1-i];
  }

 double i1,i2,i3;
  for(i=0;i<fm_size1;i+=6){
    i1=2.0;
    for(j=1;j<=4;j++){
      if(f1[i+j]>0.0) i1+=f1[i+j]*log(f1[i+j])/log(2.0);
    }
    for(j=0;j<min_num;j++){
      if (i1>b1v[j]){
	for(k=min_num-1;k>j;k--){
	  b1v[k]=b1v[k-1];
	  b1[k]=b1[k-1];
	}
	b1v[j]=i1;
	b1[j]=i;
	break;
      }
    }
  }

  for(i=0;i<fm_size2;i+=6){
    i2=i3=2.0;
    for(j=1;j<=4;j++){
      if(f2[i+j]>0.0) i2+=f2[i+j]*log(f2[i+j])/log(2.0);
      if(f3[i+j]>0.0) i3+=f3[i+j]*log(f3[i+j])/log(2.0);
    }
    //    cerr<<"i i2 i3 "<<i<<'\t'<<i2<<'\t'<<i3<<'\n';
    for(j=0;j<min_num;j++){
      if (i2>b2v[j]){
	for(k=min_num-1;k>j;k--){
	  b2v[k]=b2v[k-1];
	  b2[k]=b2[k-1];
	}
	b2v[j]=i2;
	b2[j]=i;
	break;
      }
    }
    for(j=0;j<min_num;j++){
      if (i3>b3v[j]){
	for(k=min_num-1;k>j;k--){
	  b3v[k]=b3v[k-1];
	  b3[k]=b3[k-1];
	}
	b3v[j]=i3;
	b3[j]=i;
	break;
      }
    }
  }
  int l1,l2,l3,r1,r2,r3;
  l1=fm_size1;l2=l3=fm_size2;
  r1=r2=r3=-1;
  for(i=0;i<min_num;i++){
    if(b1[i]>r1) r1=b1[i];
    if(b1[i]<l1) l1=b1[i];
    if(b2[i]>r2) r2=b2[i];
    if(b2[i]<l2) l2=b2[i];
    if(b3[i]>r3) r3=b3[i];
    if(b3[i]<l3) l3=b3[i];
  }

  //also consider ties to the worst considered position
  for(i=l1-6;i>=0;i-=6){
    i1=2.0;
    for(j=1;j<=4;j++){
      if(f1[i+j]>0.0) i1+=f1[i+j]*log(f1[i+j])/log(2.0);
    }
    if(fabs(i1-b1v[min_num-1])<.0001) l1=i;
  }
  for(i=r1+6;i<fm_size1;i+=6){
    i1=2.0;
    for(j=1;j<=4;j++){
      if(f1[i+j]>0.0) i1+=f1[i+j]*log(f1[i+j])/log(2.0);
    }
    if(fabs(i1-b1v[min_num-1])<.0001) r1=i;
  }

  for(i=l2-6;i>=0;i-=6){
    i2=2.0;
    for(j=1;j<=4;j++){
      if(f2[i+j]>0.0) i2+=f2[i+j]*log(f2[i+j])/log(2.0);
    }
    if(fabs(i2-b2v[min_num-1])<.0001) l2=i;
  }
  for(i=r2+6;i<fm_size2;i+=6){
    i2=2.0;
    for(j=1;j<=4;j++){
      if(f2[i+j]>0.0) i2+=f2[i+j]*log(f2[i+j])/log(2.0);
    }
    if(fabs(i2-b2v[min_num-1])<.0001) r2=i;
  }

  for(i=l3-6;i>=0;i-=6){
    i3=2.0;
    for(j=1;j<=4;j++){
      if(f3[i+j]>0.0) i3+=f3[i+j]*log(f3[i+j])/log(2.0);
    }
    if(fabs(i3-b3v[min_num-1])<.0001) l3=i;
  }
  for(i=r3+6;i<fm_size2;i+=6){
    i3=2.0;
    for(j=1;j<=4;j++){
      if(f3[i+j]>0.0) i3+=f3[i+j]*log(f3[i+j])/log(2.0);
    }
    if(fabs(i3-b3v[min_num-1])<.0001) r3=i;
  }

  //Up to this point is all initialization

  int ref,ref1,ref2;
  int refrc,ref1rc,ref2rc;
  ref1=max(0,l2-r1+l1);
  ref1=max(ref1,r2-fm_size1+6+l1);
  ref2=min(fm_size2-6-r1+l1,r2);
  ref2=min(ref2,l1+l2);
  ref1rc=max(0,l3-r1+l1);
  ref1rc=max(ref1rc,r3-fm_size1+6+l1);
  ref2rc=min(fm_size2-6-r1+l1,r3);
  ref2rc=min(ref2rc,l1+l3);
  //these limits put ref frame within array bounds and require at least one overlap pos
  //ref1,ref2, etc. refer to the position of the second seq aligned with l1

  double best_coeff=-1.1, coeff;
  double avg_1,avg_2,avg_12,avg_1sq,avg_2sq;
  int x1,y1,x2,y2,dp;

  for(ref=ref1;ref<=ref2;ref+=6){
    x1=l1,y1=r1,x2=ref,y2=ref+r1-l1;
    if(x2>l2){
      x1-=(x2-l2);
      x2-=(x2-l2);
    }
    if(y2<r2){
      y1+=(r2-y2);
      y2+=(r2-y2);
    }
    //x1..y1 in seq1, x2..y2 in seq2
    dp=(y1-x1)/6+1;
    avg_1=avg_2=avg_12=avg_1sq=avg_2sq=0.0;
    for(i=0;i<=(y1-x1);i+=6){
      for(j=1;j<=4;j++){
	avg_1+=f1[x1+j+i];
	avg_2+=f2[x2+j+i];
	avg_12+=f1[x1+j+i] * f2[x2+j+i];
	avg_1sq+=f1[x1+j+i] * f1[x1+j+i];
	avg_2sq+=f2[x2+j+i] * f2[x2+j+i];
      }
    }
    avg_1/=(4.0*dp);
    avg_2/=(4.0*dp);
    avg_12/=(4.0*dp);
    avg_1sq/=(4.0*dp);
    avg_2sq/=(4.0*dp);
    coeff=(avg_12-avg_1*avg_2)/sqrt((avg_1sq-avg_1*avg_1)*(avg_2sq-avg_2*avg_2));
    if(coeff>best_coeff) {
      best_coeff=coeff;
      cerr<<"\rw"<<ref<<" ";
    }
  }

  for(ref=ref1rc;ref<=ref2rc;ref+=6){
    x1=l1,y1=r1,x2=ref,y2=ref+r1-l1;
    if(x2>l3){
      x1-=(x2-l3);
      x2-=(x2-l3);
    }
    if(y2<r3){
      y1+=(r3-y2);
      y2+=(r3-y2);
    }
    //x1..y1 in seq1, x2..y2 in seq2
    dp=(y1-x1)/6+1;
    avg_1=avg_2=avg_12=avg_1sq=avg_2sq=0.0;
    for(i=0;i<=(y1-x1);i+=6){
      for(j=1;j<=4;j++){
	avg_1+=f1[x1+j+i];
	avg_2+=f3[x2+j+i];
	avg_12+=f1[x1+j+i] * f3[x2+j+i];
	avg_1sq+=f1[x1+j+i] * f1[x1+j+i];
	avg_2sq+=f3[x2+j+i] * f3[x2+j+i];
      }
    }
    avg_1/=(4.0*dp);
    avg_2/=(4.0*dp);
    avg_12/=(4.0*dp);
    avg_1sq/=(4.0*dp);
    avg_2sq/=(4.0*dp);
    coeff=(avg_12-avg_1*avg_2)/sqrt((avg_1sq-avg_1*avg_1)*(avg_2sq-avg_2*avg_2));
    if(coeff>best_coeff) {
      best_coeff=coeff;
      cerr<<"\rc"<<ref<<" ";
    }
  }
  return best_coeff;
  delete [] f1;
  delete [] f2;
  delete [] f3;
  delete [] b1;
  delete [] b2;
  delete [] b3;
  delete [] b1v;
  delete [] b2v;
  delete [] b3v;
}
