//Copyright 1998 President and Fellows of Harvard University
//alignace.cpp

#include "alignace.h"

AlignACE::AlignACE() {
}

AlignACE::~AlignACE(){
  delete [] ace_freq_matrix;
  delete [] ace_score_matrix;
  delete [] ace_membership;
	for(int i=0;i<ace_seqset.num_seqs();i++){
    delete [] ace_site_bias[i];
  }
}

void AlignACE::init(const vector<string>& v, const int nc, const int bf, const double map_cut, const double sim_cut) {
  ace_seqset.init(v);
	ace_sites.init(v, nc, nc);
	ace_select_sites.init(v, nc);
	ace_print_sites.init(v, nc, nc);
	ace_archive.init(ace_sites,ace_seqset,bf,map_cut,sim_cut);
	set_default_params();
  ace_max_motifs=bf;
  ace_map_cutoff=map_cut;
  ace_sim_cutoff=sim_cut;
  ace_freq_matrix=new int[ace_sites.depth()*ace_sites.ncols()];
  ace_score_matrix=new double[ace_sites.depth()*ace_sites.ncols()];
  ace_membership = new bool[ace_seqset.num_seqs()];
	clear_possible();
	ace_site_bias=new double*[ace_seqset.num_seqs()];
	for(int i=0;i<ace_seqset.num_seqs();i++){
		ace_site_bias[i]=new double[ace_seqset.len_seq(i)];
  }
  ace_verbose=false;
}

void AlignACE::set_default_params(){
  ace_params.ap_expect=10;
  ace_params.ap_gcback=0.38;
  ace_params.ap_minpass[0]=200;
  ace_params.ap_seed=-1;
  ace_params.ap_psfact=0.1;
  ace_params.ap_weight=0.8; 
  ace_params.ap_npass=1000000;
  ace_params.ap_fragment=true;
  ace_params.ap_flanking=0;
  ace_params.ap_undersample=1;
  ace_params.ap_oversample=1;
}

void AlignACE::set_final_params(){
  ace_params.ap_npseudo=ace_params.ap_expect * ace_params.ap_psfact;
  ace_params.ap_backfreq[0]=ace_params.ap_backfreq[5]=0.25;
  ace_params.ap_backfreq[1]=ace_params.ap_backfreq[4]=(1-ace_params.ap_gcback)/2.0;
  ace_params.ap_backfreq[2]=ace_params.ap_backfreq[3]=ace_params.ap_gcback/2.0;
  for(int i=0;i<6;i++){
    ace_params.ap_pseudo[i]=ace_params.ap_npseudo * ace_params.ap_backfreq[i];
  }
  ace_params.ap_maxlen=3*ace_sites.width();
  double ap=(double)ace_params.ap_expect/(2.0*ace_sites.positions_available());
  ace_params.ap_minpass[1]=2*ace_params.ap_minpass[0];
  ace_params.ap_minpass[2]=3*ace_params.ap_minpass[0];
	ace_params.ap_minpass[3]=4*ace_params.ap_minpass[0];
	ace_params.ap_minpass[4]=5*ace_params.ap_minpass[0];
  ace_params.ap_sitecut[0]=ap*5.0;//10.0*ap;
	ace_params.ap_sitecut[4]=0.2;
	ace_params.ap_sitecut[2]=sqrt(ace_params.ap_sitecut[0]*ace_params.ap_sitecut[4]);
	ace_params.ap_sitecut[3]=sqrt(ace_params.ap_sitecut[2]*ace_params.ap_sitecut[4]);
	ace_params.ap_sitecut[1]=sqrt(ace_params.ap_sitecut[0]*ace_params.ap_sitecut[2]);
	ace_params.ap_nruns=ace_sites.positions_available()/ace_params.ap_expect/ace_sites.ncols()/ace_params.ap_undersample*ace_params.ap_oversample;
	ace_params.ap_select=5.0;
}

void AlignACE::ace_initialize(){
  ace_ran_int.set_seed(ace_params.ap_seed);
  ace_params.ap_seed=ace_ran_int.seed();
  ace_ran_int.set_range(0,RAND_MAX);
  ace_ran_dbl.set_seed(ace_ran_int.rnum());
  ace_ran_dbl.set_range(0.0,1.0);
}

void AlignACE::add_possible(int poss) {
	if(! ace_membership[poss]) {
		ace_membership[poss] = true;
		ace_members++;
	}
}

void AlignACE::remove_possible(int poss) {
	if(ace_membership[poss]) {
		ace_membership[poss] = false;
		ace_members--;
	}
}

void AlignACE::clear_possible() {
	for(int i = 0; i < ace_seqset.num_seqs(); i++) {
		ace_membership[i] = false;
	}
	ace_members = 0;
}

bool AlignACE::is_possible(int poss) {
	return ace_membership[poss];
}

void AlignACE::seed_random_site() {
	if(ace_members < 1) return;
	
	ace_sites.clear_sites();
  int chosen_seq, chosen_posit;
  bool watson;
	
	chosen_seq = chosen_posit = -1;
	
  for(int j = 0; j < 50; j++) {
		ace_ran_int.set_range(0, ace_seqset.num_seqs() - 1);
		chosen_seq = ace_ran_int.rnum();
		if(ace_membership[chosen_seq] == 0)  continue;
		if(ace_seqset.len_seq(chosen_seq) < ace_sites.width()) continue;
		ace_ran_int.set_range(0, ace_seqset.len_seq(chosen_seq) - ace_sites.width() - 1);
		double db = ace_ran_dbl.rnum();//random (0,1)
		watson = (db > 0.5);
		chosen_posit = ace_ran_int.rnum();
		if(watson && (chosen_posit > ace_seqset.len_seq(chosen_seq) - ace_sites.width() - 1)) continue;
		if((! watson) && (chosen_posit < ace_sites.width())) continue;
		if(ace_sites.is_open_site(chosen_seq, chosen_posit)) {
      ace_sites.add_site(chosen_seq, chosen_posit, watson);
      clear_possible();
			add_possible(chosen_seq);
			break;
    }
  }
}


void AlignACE::seed_biased_site(){
  int i,j,k,m;
  double ap=(double)ace_params.ap_expect/(2.0*ace_sites.positions_available());
  char **ss_seq;
  ss_seq=ace_seqset.seq_ptr();
  double Lw,Lc,Pw,Pc,T,sc;
  int matpos,col;
  for(i=0;i<ace_seqset.num_seqs();i++){
    for(j=0;j<ace_seqset.len_seq(i)-ace_sites.ncols()+1;j++){
      ace_site_bias[i][j]=1.0;
    }
  }
  for(m=0;m<ace_max_motifs;m++){
    sc=get_best_motif(m);
    if(sc<=0.0) break;
    calc_matrix();
    for(i=0;i<ace_seqset.num_seqs();i++){
      for(j=0;j<ace_seqset.len_seq(i)-ace_sites.width()+1;j++){
				Lw=1.0;matpos=0;col=0;
				for(k=0;k<ace_sites.ncols();k++){
					int seq=ss_seq[i][j+col];
					Lw*=ace_score_matrix[matpos+seq];
					col=ace_sites.next_column(col);
					matpos+=ace_sites.depth();
				}
				Lc=1.0;matpos=ace_sites.depth()-1;col=0;
				for(k=0;k<ace_sites.ncols();k++){
					int seq=ss_seq[i][j+ace_sites.width()-1-col];
					Lc*=ace_score_matrix[matpos-seq];
					col=ace_sites.next_column(col);
					matpos+=ace_sites.depth();
				}
				Pw=Lw*ap/(1.0-ap+Lw*ap);
				Pc=Lc*ap/(1.0-ap+Lc*ap);
				//if(Pw<0.8&&Pc<0.8) continue;
				for(k=j;k<=(j+ace_sites.width()-ace_sites.ncols());k++){
					ace_site_bias[i][k]*=(1-Pw)*(1-Pc);
					//ace_site_bias[i][k]=0.0;
				}
      }
    }
  }
  if(m==0){
    seed_random_site();
    return;
  }
	
  for(T=0.0,i=0;i<ace_seqset.num_seqs();i++){
    for(j=0;j<ace_seqset.len_seq(i)-ace_sites.ncols()+1;j++){
      T+=ace_site_bias[i][j];
      //      if(m>3) cerr<<ace_site_bias[i][j]<<'\t'<<T<<'\n';
    }
  }
	
  cerr<<"seeding "<<T<<" biased\n";
	
  ace_sites.clear_sites();
	
  T*=ace_ran_dbl.rnum();
  for(i=0;i<ace_seqset.num_seqs();i++){
    for(j=0;j<ace_seqset.len_seq(i)-ace_sites.ncols()+1;j++){
      T-=ace_site_bias[i][j];
      if(T<=0.0) break;
    }
    if(T<=0.0) break;
  }
	
  double db=ace_ran_dbl.rnum();//random (0,1)
		bool watson=true;
		if(db>0.5) watson=false;
		ace_sites.add_site(i,j,watson);
}


void AlignACE::calc_matrix() {
  int d = ace_sites.depth();
  double tot = (double) ace_sites.number() + ace_params.ap_npseudo;
  ace_sites.calc_freq_matrix(ace_seqset, ace_freq_matrix);
  for(int i = 0; i < d * ace_sites.ncols();i += d){
    ace_score_matrix[i] = ace_score_matrix[i + 5] = 1.0;
    for(int j = 1;j <= 4; j++){
      int x = ace_freq_matrix[i] + ace_freq_matrix[i + 5];
      ace_score_matrix[i + j] = (ace_freq_matrix[i + j] + x * ace_params.ap_backfreq[j] 
				+ ace_params.ap_pseudo[j])/(tot * ace_params.ap_backfreq[j]);
    }
  }
}

void AlignACE::single_pass(const double minprob) {
	double ap = (ace_params.ap_weight * ace_params.ap_expect
							+ (1 - ace_params.ap_weight) * ace_sites.number())
							/(2.0 * ace_sites.positions_available(ace_membership));
  calc_matrix();
	ace_sites.remove_all_sites();
  ace_select_sites.remove_all_sites();
  //will only update once per pass
	
  char **ss_seq;
  ss_seq = ace_seqset.seq_ptr();
  double Lw, Lc, Pw, Pc, F;
  int matpos, col;
  int considered = 0;
  int gadd = -1,jadd = -1;
  for(int g = 0; g < ace_seqset.num_seqs(); g++){
		if (! ace_membership[g]) continue;
		for(int j = 0; j <= ace_seqset.len_seq(g) - ace_sites.width(); j++){
			Lw = 1.0;
			matpos = 0;
			col = 0;
      for(int k = 0; k < ace_sites.ncols(); k++){
				assert(j + col >= 0);
				assert(j + col <= ace_seqset.len_seq(g));
				int seq = ss_seq[g][j + col];
				Lw *= ace_score_matrix[matpos + seq];
				col = ace_sites.next_column(col);
				matpos += ace_sites.depth();
      }
      Lc = 1.0;
			matpos = ace_sites.depth() - 1;
			col = 0;
      for(int k = 0; k < ace_sites.ncols(); k++){
				assert(j + ace_sites.width() - 1 - col >= 0);
				assert(j + ace_sites.width() - 1 - col <= ace_seqset.len_seq(g));
				int seq = ss_seq[g][j + ace_sites.width() - 1 - col];
				Lc *= ace_score_matrix[matpos - seq];
				col = ace_sites.next_column(col);
				matpos += ace_sites.depth();
      }
      Pw = Lw * ap/(1.0 - ap + Lw * ap);
      Pc = Lc * ap/(1.0 - ap + Lc * ap);
      F = Pw + Pc - Pw * Pc;//probability of either
			if(F > (minprob/ace_params.ap_select)) {
				ace_select_sites.add_site(g, j, true);
			}
			//strand irrelevant for ace_select_sites
			if(g == gadd && j < jadd + ace_sites.width()) continue;
			if(F < minprob) continue;
			considered++;
			Pw = F * Pw / (Pw + Pc);
			Pc = F - Pw;
			double r = ace_ran_dbl.rnum();
			if (r > F)
				continue;
			else if (r < Pw) {
				assert(j >= 0);
				assert(j <= ace_seqset.len_seq(g) - ace_sites.width());
				ace_sites.add_site(g, j, true);
				gadd = g;
				jadd = j;
			} else {
				assert(j >= 0);
				assert(j <= ace_seqset.len_seq(g) - ace_sites.width());
				ace_sites.add_site(g, j, false);
				gadd = g;
				jadd = j;
			}
    }
  }
  if(ace_verbose)cerr<<"sampling "<<considered<<"/"<<ace_select_sites.number()<<"\n";
}

void AlignACE::single_pass_select(const double minprob){
	int i,j,k,n;
  double ap=(ace_params.ap_weight*ace_params.ap_expect+(1-ace_params.ap_weight)*ace_sites.number())/(2.0*ace_sites.positions_available());
  calc_matrix();
	ace_sites.remove_all_sites();
  //will only update once per pass
	
  char **ss_seq;
  ss_seq=ace_seqset.seq_ptr();
  double Lw,Lc,Pw,Pc,F;
  int matpos,col;
  int iadd=-1,jadd=-1;
  for(n=0;n<ace_select_sites.number();n++){
    i=ace_select_sites.chrom(n);
    j=ace_select_sites.posit(n);
    if(i==iadd&&j<jadd+ace_sites.width()) continue;
    if(j<0||j>ace_seqset.len_seq(i)-ace_sites.width()) continue;
    //could be screwed up with column sampling
    Lw=1.0;matpos=0;col=0;
    for(k=0;k<ace_sites.ncols();k++){
      int seq=ss_seq[i][j+col];
      Lw*=ace_score_matrix[matpos+seq];
      col=ace_sites.next_column(col);
      matpos+=ace_sites.depth();
    }
    Lc=1.0;matpos=ace_sites.depth()-1;col=0;
    for(k=0;k<ace_sites.ncols();k++){
      int seq=ss_seq[i][j+ace_sites.width()-1-col];
      Lc*=ace_score_matrix[matpos-seq];
      col=ace_sites.next_column(col);
      matpos+=ace_sites.depth();
    }
    Pw=Lw*ap/(1.0-ap+Lw*ap);
    Pc=Lc*ap/(1.0-ap+Lc*ap);
    F=(Pw+Pc-Pw*Pc);//probability of either
		if(F<minprob) continue;
		Pw=F*Pw/(Pw+Pc);
		Pc=F-Pw;
		double r=ace_ran_dbl.rnum();
		if (r>F)
			continue;
    else if (r<Pw) {
      ace_sites.add_site(i,j,true);
      iadd=i;jadd=j;
    } else {
      ace_sites.add_site(i,j,false);
      iadd=i;jadd=j;
    }
  }
}

bool AlignACE::column_sample(const int c, const bool sample){
	//cerr << "\t\t\tRunning column sample with c = " << c << "... ";
	//sample default to true, if false then replace with best column
  //just consider throwing out the worst column, sample for replacement, no need to sample both ways, unless column is specified
  int i, j;
  int col = 0, col_worst = 0, col_removed = 0;
  int *freq = new int[ace_sites.depth()];
  double wt, wt_worst = DBL_MAX;
	
  if(c!=-1)
		col_worst = c;   // user chosen, hopefully a real column
  else {
    for(i = 0; i < ace_sites.ncols(); i++){
      ace_sites.column_freq(col, ace_seqset, freq);
      wt = 0.0;
      for(j = 0; j < ace_sites.depth(); j++) {
				wt += gammaln(freq[j] + ace_params.ap_pseudo[j]);
				wt -= (double)freq[j] * log(ace_params.ap_backfreq[j]);
      }
      if(wt < wt_worst) {
				col_worst = col;
				wt_worst = wt;
      }
      col = ace_sites.next_column(col);
    }
  }
	
  col_removed = ace_sites.remove_col(col_worst);
  i = ace_select_sites.remove_col(col_worst);
  if(i != col_removed) {
		cerr << i << "  " << col_removed << " wrong assumption in column_sample\n";
		abort();
	}
	
  int max_left, max_right;
  max_left = max_right = (ace_sites.max_width()-ace_sites.width())/2;
  ace_sites.columns_open(max_left,max_right);
  int cs_span=max_left+max_right+ace_sites.width();
  double *wtx=new double[cs_span];
  int x=max_left;
  //wtx[x+c] will refer to the weight of pos c in the usual numbering
  col=0;
  double best_wt=0.0;
  for(i=0;i<cs_span;i++){
    wtx[i]=0.0;
    if((i-x)==col){
      col=ace_sites.next_column(col);
      continue;
    }
    if(ace_sites.column_freq(i-x,ace_seqset,freq)){
      wt=0.0;
      for(j=0;j<ace_sites.depth();j++){
				wt+=gammaln(freq[j]+ace_params.ap_pseudo[j]);
				wt-=(double)freq[j]*log(ace_params.ap_backfreq[j]);
      }
      wtx[i]=wt;
      if(wt>best_wt) best_wt=wt;
    }
  }
	
  double scale=0.0;
  if(best_wt>100.0) scale=best_wt-100.0;//keep exp from overflowing
	double tot2=0.0;
	for(i=0;i<cs_span;i++){
		if(wtx[i]==0.0) continue;
		wtx[i]-=scale;
		wtx[i]=exp(wtx[i]);
		int newwidth=ace_sites.width();
		if(i<x) newwidth+=(x-i);
		else if(i>(x+ace_sites.width()-1)) newwidth+=(i-x-ace_sites.width()+1);
		wtx[i]/=bico(newwidth-2,ace_sites.ncols()-2);
		tot2+=wtx[i];
	}
	
	double pick;
	int col_pick;
	double cutoff=.01;
	if(sample){
		if(1-(wtx[x+col_removed]/tot2)<cutoff){
			ace_sites.add_col(col_removed);
			ace_select_sites.add_col(col_removed);
			delete [] wtx;
			delete [] freq;
			return false;
		} 
		pick=ace_ran_dbl.rnum()*tot2;
		col_pick=373;
		for(i=0;i<cs_span;i++){
			if(wtx[i]==0.0) continue;
			pick-=wtx[i];
			if(pick<=0.0) {
				col_pick=i-x;
				break;
			}
		}
	}
	else{//select best
		pick=-DBL_MAX;
		for(i=0;i<cs_span;i++){
			if(wtx[i]>pick){
				col_pick=i-x;
				pick=wtx[i];
			}
		}
	}
	if(col_pick==373){
		cout<<tot2<<'\t'<<"373 reached.\n";
		abort();
	}
	ace_sites.add_col(col_pick);
	ace_select_sites.add_col(col_pick);
	
	delete [] wtx;
	delete [] freq;
	
	// cerr << "done." <<  endl;
	
	if(col_removed==col_pick) return false;
  else return true;
}

double AlignACE::map_score(){
  int i,j,k;
  double ms=0.0;
  double map_N=ace_sites.positions_available(ace_membership);  
  double w=ace_params.ap_weight/(1.0-ace_params.ap_weight);
  double map_alpha=(double) ace_params.ap_expect*w;
  double map_beta=map_N*w - map_alpha;
  double map_success=(double)ace_sites.number();
  ms +=( gammaln(map_success+map_alpha)+gammaln(map_N-map_success+map_beta) );
  ms -=( gammaln(map_alpha)+gammaln(map_N + map_beta) );
	
  ace_sites.calc_freq_matrix(ace_seqset,ace_freq_matrix);
  double sc[6]={0.0,0.0,0.0,0.0,0.0,0.0};
  int d=ace_sites.depth();
  for(k=0;k!=d*ace_sites.ncols();k+=d) {
    int x=ace_freq_matrix[k]+ace_freq_matrix[k+5];
    for(j=1;j<=4;j++){
      ms += gammaln((double)ace_freq_matrix[k+j]+x*ace_params.ap_backfreq[j]+ace_params.ap_pseudo[j]);
      sc[j]+=ace_freq_matrix[k+j]+x*ace_params.ap_backfreq[j];
    }
  }
  ms-=ace_sites.ncols()*gammaln((double)ace_sites.number()+ace_params.ap_npseudo);
  for (k=1;k<=4;k++) {
    ms -= sc[k]*log(ace_params.ap_backfreq[k]);
  }
  /*This factor arises from a modification of the model of Liu, et al in which the background frequencies of DNA bases are taken to be constant for the organism under consideration*/
  double  vg;
  ms -= lnbico(ace_sites.width()-2, ace_sites.ncols()-2);
  for(vg=0.0,k=1;k<=4;k++) {
    vg += gammaln(ace_params.ap_pseudo[k]);
  }
  vg-=gammaln((double)(ace_params.ap_npseudo));
  ms-=((double)ace_sites.ncols()*vg);
  return ms;
}

double AlignACE::get_best_motif(int i){
  return ace_archive.return_best(ace_sites,i);
}

void AlignACE::optimize_columns(){
  while(column_sample(false)){}
  //replaces worst column with best column until these are the same
}

void AlignACE::optimize_sites(){
  int possible = ace_sites.positions_available()+1;
  int *pos = new int[possible];
  int *chr = new int[possible];
  bool *str = new bool[possible];
  Heap hp(possible,3);
  double ap = (ace_params.ap_weight * ace_params.ap_expect
	            + (1 - ace_params.ap_weight) * ace_sites.number())
							/ (2.0 * ace_sites.positions_available(ace_membership));
	
  calc_matrix();
  char **ss_seq;
  ss_seq = ace_seqset.seq_ptr();
  double Lw, Lc, Pw, Pc;
  int matpos, col;
  int h = 1;
  double cutoff = 0.2;
  for(int i = 0; i < ace_seqset.num_seqs(); i++) {
    if(ace_membership[i] == 0) continue;
		for(int j = 0; j < ace_seqset.len_seq(i) - ace_sites.width() + 1; j++) {
      Lw = 1.0;
			matpos = 0;
			col = 0;
      for(int k = 0;k < ace_sites.ncols(); k++){
				int seq = ss_seq[i][j + col];
				Lw *= ace_score_matrix[matpos + seq];
				col = ace_sites.next_column(col);
				matpos += ace_sites.depth();
      }
      Lc = 1.0;
			matpos = ace_sites.depth()-1;
			col = 0;
      for(int k = 0; k < ace_sites.ncols(); k++){
				int seq = ss_seq[i][j + ace_sites.width()-1-col];
				Lc *= ace_score_matrix[matpos - seq];
				col = ace_sites.next_column(col);
				matpos += ace_sites.depth();
      }
      Pw = Lw * ap / (1.0 - ap + Lw * ap);
      Pc = Lc * ap / (1.0 - ap + Lc * ap);
      //F=(Pw+Pc-Pw*Pc);//probability of either
      if(Pw > cutoff || Pc > cutoff) {
				chr[h] = i;
				pos[h] = j;
				if(Pw >= Pc) {
					str[h] = true;
					hp.insert(h, -Pw);//lower is better
				} else {
					str[h] = false;
					hp.insert(h, -Pc);
				}
				h++;
      }
    }
  }
  
  ace_sites.remove_all_sites();
  Sites best = ace_sites;
  double ms, ms_best = -DBL_MAX;
  for(;;){
    h = hp.del_min();
    if(h == -1) break;
    if(ace_sites.is_open_site(chr[h], pos[h])) {
      ace_sites.add_site(chr[h], pos[h], str[h]);
      if(ace_sites.number() <= 3) continue;
      if(hp.value(h) < -0.7) continue;
      ms = map_score();
      if(ms > ms_best){
				ms_best = ms;
				best = ace_sites;
      }
    }
  }
  ace_sites = best;
  delete [] pos;
  delete [] chr;
  delete [] str;
}

void AlignACE::orient_motif(){
  double *info = new double[6];
  double *freq = new double[6];
  for(int i = 0; i < 6; i++) info[i] = 0.0;
  int d = ace_sites.depth();
	
  double tot = (double)ace_sites.number() + ace_params.ap_npseudo;
  ace_sites.calc_freq_matrix(ace_seqset, ace_freq_matrix);
  for(int i = 0; i < d * ace_sites.ncols(); i += d){
    double ii = 0.0;
    for(int j = 1; j <= 4; j++) {
      int x = ace_freq_matrix[i] + ace_freq_matrix[i+5];
      freq[j] = (ace_freq_matrix[i + j] + x * ace_params.ap_backfreq[j] + ace_params.ap_pseudo[j]) / tot;
      ii += freq[j]*log(freq[j]);
    }
    ii = 2 + ii;
    for(int j = 1; j <= 4; j++) info[j] += freq[j] * ii;
  }
	
  double flip = 1.5 * info[3] + 1.0 * info[1] - 1.0 * info[4] - 1.5 * info[2];
  //for(i=1;i<5;i++) cerr<<info[i]<<'\t';
  //cerr<<flip<<'\n';
  if(flip < 0.0) ace_sites.flip_sites();
  delete [] info;
  delete [] freq;
}

void AlignACE::orient_print_motif(){
  double *info = new double[6];
  double *freq = new double[6];
  for(int i = 0; i < 6; i++) info[i] = 0.0;
  int d = ace_print_sites.depth();
	
  double tot = (double) ace_print_sites.number() + ace_params.ap_npseudo;
  ace_print_sites.calc_freq_matrix(ace_seqset, ace_freq_matrix);
  for(int i = 0; i < d * ace_print_sites.ncols(); i += d){
    double ii = 0.0;
    for(int j = 1; j <= 4; j++) {
      int x = ace_freq_matrix[i] + ace_freq_matrix[i+5];
      freq[j] = (ace_freq_matrix[i + j] + x * ace_params.ap_backfreq[j] + ace_params.ap_pseudo[j]) / tot;
      ii += freq[j]*log(freq[j]);
    }
    ii = 2 + ii;
    for(int j = 1; j <= 4; j++) info[j] += freq[j] * ii;
  }
	
  double flip = 1.5 * info[3] + 1.0 * info[1] - 1.0 * info[4] - 1.5 * info[2];
  //for(i=1;i<5;i++) cerr<<info[i]<<'\t';
  //cerr<<flip<<'\n';
  if(flip < 0.0) ace_print_sites.flip_sites();
  delete [] info;
  delete [] freq;
}

string AlignACE::consensus() {
	map<char,char> nt;
  nt[0]=nt[5]='N';
  nt[1]='A';nt[2]='C';nt[3]='G';nt[4]='T';
	char** ss_seq=ace_seqset.seq_ptr();
	
	int numsites = ace_sites.number();
	if(numsites < 1) return "";
	
	// cerr << "Computing consensus with " << numsites << " sites" << endl;
	vector<string> hits(numsites);
	for(int i = 0; i < numsites; i++){
		int c = ace_sites.chrom(i);
    int p = ace_sites.posit(i);
    bool s = ace_sites.strand(i);
    for(int j = 0; j < ace_sites.width(); j++){
      if(s) {
				if(p + j >= 0 && p+j < ace_seqset.len_seq(c))
					 hits[i] += nt[ss_seq[c][p+j]];
				else hits[i] += ' ';
      } else {
				if(p + ace_sites.width() - 1 - j >= 0 && p+ace_sites.width() - 1 - j < ace_seqset.len_seq(c))
					hits[i] += nt[ace_sites.depth()-1-ss_seq[c][p+ace_sites.width()-1-j]];
				else hits[i] += ' ';
      }
    }
  }
	
	string cons;
  int wide1 = hits[0].size();
  int num1 = hits.size();
  int num1a, num1c, num1g, num1t;
  for(int i = 0; i < wide1; i++){
    num1a = num1c = num1g = num1t = 0;
    for(int j = 0; j < num1; j++){
      if(hits[j][i] == 'a' || hits[j][i] == 'A') num1a += 1;
      if(hits[j][i] == 'c' || hits[j][i] == 'C') num1c += 1;
      if(hits[j][i] == 'g' || hits[j][i] == 'G') num1g += 1;
      if(hits[j][i] == 't' || hits[j][i] == 'T') num1t += 1;
    }
    if(num1a > num1*0.7) cons += 'A';
    else if(num1c > num1*0.7) cons += 'C';
    else if(num1g > num1*0.7) cons += 'G';
    else if(num1t > num1*0.7) cons += 'T';
    else if((num1a + num1t) > num1 * 0.85) cons +='W';
    else if((num1a + num1c) > num1 * 0.85) cons +='M';
    else if((num1a + num1g) > num1 * 0.85) cons +='R';
    else if((num1t + num1c) > num1 * 0.85) cons +='Y';
    else if((num1t + num1g) > num1 * 0.85) cons +='K';
    else if((num1c + num1g) > num1 * 0.85) cons +='S';
    else cons += '-';
  }
  return cons;
}

void AlignACE::output(ostream &fout){
  map<char,char> nt;
  nt[0] = nt[5] = 'N';
  nt[1] = 'A';
	nt[2] = 'C';
	nt[3] = 'G';
	nt[4] = 'T';
  char** ss_seq = ace_seqset.seq_ptr();
  int x = ace_params.ap_flanking;
  for(int i = 0; i < ace_print_sites.number(); i++){
    int c = ace_print_sites.chrom(i);
    int p = ace_print_sites.posit(i);
    bool s = ace_print_sites.strand(i);
    for(int j = -x; j < ace_print_sites.width() + x; j++){
      if(s) {
				if(p + j >= 0 && p + j < ace_seqset.len_seq(c))
					fout << nt[ss_seq[c][p + j]];
				else fout << ' ';
      }
      else {
				if(p + ace_print_sites.width() - 1 - j >= 0 && p + ace_print_sites.width()-1-j < ace_seqset.len_seq(c))
					fout << nt[ace_print_sites.depth() - 1 - ss_seq[c][p + ace_print_sites.width() - 1 - j]];
				else fout << ' ';
      }
    }
    fout << '\t' << c << '\t' << p << '\t' << s << '\n';
  }
  for(int i = 0; i < x; i++) fout << ' ';
  int j = 0;
	for(int i = 0;;){
    j = ace_print_sites.next_column(i);
    fout << '*';
    if(i == ace_print_sites.width() - 1) break;
    for(int k = 0; k < (j - i - 1); k++) fout << ' ';
    i = j;
  }
  fout << "\n";
}

void AlignACE::full_output(ostream &fout){
  for(int j = 0; j < ace_max_motifs; j++){
    double sc = ace_archive.return_best(ace_print_sites, j);
    orient_print_motif();
    if(sc > 0.0){
      fout << "Motif " << j + 1 << '\n';
      output(fout);
      fout << "MAP Score: " << sc << "\n\n";
    }
    else break;
  }
}

void AlignACE::full_output(char *name){
  ofstream fout(name);
  full_output(fout);
}

void AlignACE::output_params(ostream &fout){
  fout<<" expect =      \t"<<ace_params.ap_expect<<'\n';
  fout<<" gcback =      \t"<<ace_params.ap_gcback<<'\n';
  fout<<" minpass =     \t"<<ace_params.ap_minpass[0]<<'\n';
  fout<<" seed =        \t"<<ace_params.ap_seed<<'\n';
  fout<<" numcols =     \t"<<ace_sites.ncols()<<'\n';
  fout<<" undersample = \t"<<ace_params.ap_undersample<<'\n';
  fout<<" oversample = \t"<<ace_params.ap_oversample<<'\n';
}

void AlignACE::modify_params(int argc, char *argv[]){
  GetArg2(argc,argv,"-expect",ace_params.ap_expect);
  GetArg2(argc,argv,"-gcback",ace_params.ap_gcback);
  GetArg2(argc,argv,"-minpass",ace_params.ap_minpass[0]);
  GetArg2(argc,argv,"-seed",ace_params.ap_seed);
  GetArg2(argc,argv,"-undersample",ace_params.ap_undersample);
  GetArg2(argc,argv,"-oversample",ace_params.ap_oversample);
}

void AlignACE::print_usage(ostream &fout){
  print_version(fout);
  fout<<"Usage: AlignACE -i seqfile (options)\n";
  fout<<" Seqfile must be in FASTA format.\n";
  fout<<"Options:\n";
  fout<<" -numcols    \tnumber of columns to align (10)\n";
  fout<<" -expect     \tnumber of sites expected in model (10)\n";
  fout<<" -gcback     \tbackground fractional GC content of input sequence (0.38)\n";
  fout<<" -minpass    \tminimum number of non-improved passes in phase 1 (200)\n";
  fout<<" -seed       \tset seed for random number generator (time)\n";
  fout<<" -undersample\tpossible sites / (expect * numcols * seedings) (1)\n"; 
  fout<<" -oversample\t1/undersample (1)\n"; 
  cout<<"Output format:\n";
  cout<<" column 1:\tsite sequence\n";
  cout<<" column 2:\tsequence number\n";
  cout<<" column 3:\tposition of site within sequence\n";
  cout<<" column 4:\tstrand of site (1=forward, 0=reverse)\n";
}

void AlignACE::print_version(ostream &fout){
  fout<<"AlignACE 4.0 05/13/04\n";
}

void AlignACE::debug_check_columns() {
	assert(ace_sites.width() == ace_select_sites.width());
	int col1 = 0, nxt1 = 0;
	int col2 = 0, nxt2 = 0;
	while(true) {
		nxt1 = ace_sites.next_column(col1);
		nxt2 = ace_select_sites.next_column(col2);
		if(nxt1 == col1 && nxt2 == col2) break;
		assert(nxt1 == nxt2);
		col1 = nxt1;
		col2 = nxt2;
	}
}

