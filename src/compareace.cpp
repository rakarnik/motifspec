//Copyright 1999 President and Fellows of Harvard University
//compareace.cpp

#include "compareace.h"

CompareACESites::~CompareACESites(){
  if(comp_ready){
    delete [] comp_fm1;
    delete [] comp_fm2;
  }
}

CompareACESites::CompareACESites(const CompareACESites& s):comp_sites(s.comp_sites){
  comp_ready=false;
  if(!s.comp_ready) return;
  comp_ready=true;
  comp_fmsize=s.comp_fmsize;
  comp_l1=s.comp_l1;
  comp_r1=s.comp_r1;
  comp_l2=s.comp_l2;
  comp_r2=s.comp_r2;
  comp_fm1=new double[comp_fmsize];
  comp_fm2=new double[comp_fmsize];
  for(int i=0;i<comp_fmsize;i++){
    comp_fm1[i]=s.comp_fm1[i];
    comp_fm2[i]=s.comp_fm2[i];
  }
}

CompareACESites& CompareACESites::operator= (const CompareACESites& s){
  if(this!=&s&&s.comp_ready){
    if(comp_ready){
      delete [] comp_fm1;
      delete [] comp_fm2;
    }
    comp_ready=true;
    comp_sites=s.comp_sites;
    comp_fmsize=s.comp_fmsize;
    comp_l1=s.comp_l1;
    comp_r1=s.comp_r1;
    comp_l2=s.comp_l2;
    comp_r2=s.comp_r2;
    comp_fm1=new double[comp_fmsize];
    comp_fm2=new double[comp_fmsize];
    for(int i=0;i<comp_fmsize;i++){
      comp_fm1[i]=s.comp_fm1[i];
      comp_fm2[i]=s.comp_fm2[i];
    }
  }
  return *this;
}


void CompareACESites::init(const Sites &s, const Seqset &t){
  comp_ready=true;

  comp_sites.sites_init(s);
  comp_sites=s;
  int i,j,k;
  int min_num=6;
  int *b1=new int[min_num];
  double *b1v=new double[min_num];
  for(i=0;i<6;i++){
    b1[i]=-1;
    b1v[i]=0.0;
  }
  comp_fmsize=(s.width()+2*s.ncols())*s.depth();
  comp_fm1=new double[comp_fmsize];
  comp_fm2=new double[comp_fmsize];//2=reverse complement of 1
  s.freq_matrix_extended(t,comp_fm1);
  //assuming nacgtn here
  for(i=0;i<comp_fmsize;i+=6){
    comp_fm1[i+1]+=(comp_fm1[i]+comp_fm1[i+5])/4.0;
    comp_fm1[i+2]+=(comp_fm1[i]+comp_fm1[i+5])/4.0;
    comp_fm1[i+3]+=(comp_fm1[i]+comp_fm1[i+5])/4.0;
    comp_fm1[i+4]+=(comp_fm1[i]+comp_fm1[i+5])/4.0;
  } 
  for(i=0;i<comp_fmsize;i++){
    comp_fm2[i]=comp_fm1[comp_fmsize-1-i];
  }
  
  double i1;
  double lg2=log(2.0);
  for(i=0;i<comp_fmsize;i+=6){
    i1=2.0;
    for(j=1;j<=4;j++){
      if(comp_fm1[i+j]>0.0) i1+=comp_fm1[i+j]*log(comp_fm1[i+j])/lg2;
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

  comp_l1=comp_fmsize;
  comp_r1=-1;
  for(i=0;i<min_num;i++){
    if(b1[i]>comp_r1) comp_r1=b1[i];
    if(b1[i]<comp_l1) comp_l1=b1[i];
  }
  
  //also consider ties to the worst considered position
  for(i=comp_l1-6;i>=0;i-=6){
    i1=2.0;
    for(j=1;j<=4;j++){
      if(comp_fm1[i+j]>0.0) i1+=comp_fm1[i+j]*log(comp_fm1[i+j])/lg2;
    }
    if(fabs(i1-b1v[min_num-1])<.0001) comp_l1=i;
  }
  for(i=comp_r1+6;i<comp_fmsize;i+=6){
    i1=2.0;
    for(j=1;j<=4;j++){
      if(comp_fm1[i+j]>0.0) i1+=comp_fm1[i+j]*log(comp_fm1[i+j])/lg2;
    }
    if(fabs(i1-b1v[min_num-1])<.0001) comp_r1=i;
  }

  comp_r2=comp_fmsize-6-comp_l1;
  comp_l2=comp_fmsize-6-comp_r1;

  delete [] b1;
  delete [] b1v;
}


double CompareACESites::compare(const CompareACESites &c){

  int i,j,k;
  int ref,ref1,ref2;
  int refrc,ref1rc,ref2rc;
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

void CompareACESites::kill_sites(){
  comp_sites.destroy();
}

CompareACESites read_motif(const string fname, const int mot){
  string s;
  int x;
  ifstream fin(fname.c_str());
  if(!fin){
    cerr<<"No such file "<<fname<<'\n';
    exit(0);
  }
  if(mot>0){
    while (getline(fin,s,'\n')){
      if(s.size()<5) continue;
      if(s.substr(0,5)=="Motif"){
	if(s.size()>=12&&s.substr(0,12)=="Motif number"){
	  x=str_to_int(s.substr(13));
	}
	else{
	  x=str_to_int(s.substr(6));
	}
	if(x==mot){
	  break;
	}
      }
    }
  }
  CompareACESites cs;
  if(mot>0&&x!=mot) return cs;
  vector<string> seqset;
  int wd;
  while(getline(fin,s,'\n')){
    if(s.size()<5) continue;
    if(s[0]=='#') continue;
    if(s.find('*')!=string::npos) break;
    s = s.substr(0,s.find('\t'));
    //this sequence is pre- and post- padded with spaces for equal length sequences
    seqset.push_back(s);
    wd=s.length();
  }
  //  cout<<wd<<"\n";
  
  Seqset sq;
	sq.init(seqset);
  Sites st(sq,wd);
	for(x=0;x<seqset.size();x++){
    st.add_site(x,0,true);
  }
  cs.init(st,sq);
  cs.kill_sites();
  return cs;
}


void cCompareACE(int argc, char *argv[]){
  if(argc<2){
    cerr<<"This code compares motifs found by AlignACE.\n";
    cerr<<"Usage 1a: CompareACE ace_file1 mot1 ace_file2 mot2\n";
    cerr<<" This compares motif mot1 of ace_file1 with motif number mot2 of ace_file2.\n";
    cerr<<"Usage 1b: CompareACE ace_file1 ace_file2\n";
    cerr<<" This performs all pairwise comparisons between the motifs in the specified files (assumes all filenames have some non-digit characters).\n";
    cerr<<"Usage 2: CompareACE -all file ace_col mot_col (-c [cutoff])\n";
    cerr<<" This usage performs all pairwise comparisons between all motifs in the specified file and returns those scoring better than the given cutoff (default all).\n";
    cerr<<" The file should be tab-delimited with the name of the AlignACE file in column ace_col and the motif number in column mot_col, specifying one motif per line.\n";
    cerr<<"Usage 3: CompareACE -form file1 file2\n";
    cerr<<" These files should contain a set of aligned sites, optionally with asterisks underneath for selection of active columns, and any comment lines (starting with #).\n";
    cerr<<"The score returned by CompareACE is the maximum value of the Pearson correlation coefficient between motif base frequencies, considering all possible alignments spanning the most informative six positions of each motif.\n";
    //cerr<<"The -n option may be used to change the number of required best columns to be considered.\n";
    //cerr<<"The -s option may be used to only require that the strong region of the first motif be included in the scored overlap.\n";
    //cerr<<"Any -n or -s options must be added at the end of the command line.\n";
    //cerr<<"For example: CompareACE file1.ace 5 file2.ace 3 -n 8 -s\n";
    exit(0);
  }

  int i,j;
  double ans;
  double cutoff=-999.0;
  GetArg2(argc,argv,"-c",cutoff);
  if(!GetArg2(argc,argv,"-all")&&!GetArg2(argc,argv,"-form")){
    string ace1,ace2;
    int mot1_a,mot1_z,mot2_a,mot2_z;
    ace1=argv[1];
    if(!is_number(argv[2])){
      ace2=argv[2];
      mot1_a=1;mot2_a=1;
      mot1_z=number_motifs(ace1.c_str());
      mot2_z=number_motifs(ace2.c_str());
      if(mot1_z==0) mot1_a=mot1_z=0;
      if(mot2_z==0) mot2_a=mot2_z=0;
    }
    else{
      ace2=argv[3];
      mot1_a=mot1_z=str_to_int(argv[2]);
      mot2_a=mot2_z=str_to_int(argv[4]);
    }
    //cout<<mot1_a<<'\t'<<mot1_z<<'\t'<<mot2_a<<'\t'<<mot2_z<<'\n';
    for(i=mot1_a;i<=mot1_z;i++){
      CompareACESites csa1=read_motif(ace1,i);
      if(!csa1.ready()) exit(0);
      for(j=mot2_a;j<=mot2_z;j++){
	if(ace1==ace2&&i>=j) continue;
	CompareACESites csa2=read_motif(ace2,j);
	if(!csa2.ready()) exit(0);
	ans=csa1.compare(csa2);
	if(ans<cutoff) continue;
	if(mot1_a==mot1_z&&mot2_a==mot2_z) cout<<ans<<'\n';
	else cout<<ace1<<'\t'<<i<<'\t'<<ace2<<'\t'<<j<<'\t'<<ans<<'\n';
      }
    }
  }
  else if(GetArg2(argc,argv,"-all")){
    string fn=argv[2];
    int ace_col=str_to_int(argv[3]);
    int mot_col=str_to_int(argv[4]);
    
    ifstream fin(fn.c_str());
    string s;
    int num=0;
    while(getline(fin,s,'\n')){
      num++;
    }
    fin.close();
    
    CompareACESites* ca=new CompareACESites[num];
    for(i=0;i<num;i++) ca[i].kill_sites();
    string* fa=new string[num];
    int* ma=new int[num];
    ifstream fin2(fn.c_str());
    i=0;
    while(getline(fin2,s,'\n')){
      vector<string> v=split(s,'\t');
      string ace=v[ace_col-1];
      int mot=str_to_int(v[mot_col-1]);
      CompareACESites cs=read_motif(ace,mot);
      if(cs.ready()){
	ca[i]=cs;
	fa[i]=ace;
	ma[i]=mot;
	i++;
      }
    }
    for(i=0;i<num;i++){
      for(j=i+1;j<num;j++){
	double ans=(ca[i]).compare(ca[j]);
	if(ans>=cutoff){
	  cout<<fa[i]<<'\t'<<ma[i]<<'\t'<<fa[j]<<'\t'<<ma[j]<<'\t'<<ans<<'\n';
	}
      }
    }
  }
  else if(GetArg2(argc,argv,"-form")){
    string file1=argv[2];
    string file2=argv[3];
    CompareACESites csa1=read_motif(file1,0);
    if(!csa1.ready()) exit(0);
    CompareACESites csa2=read_motif(file2,0);
    if(!csa2.ready()) exit(0);
    ans=csa1.compare(csa2);
    cout<<ans<<'\n';    
  }  
}

