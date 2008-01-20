#include "semodel.h"

SEModel::~SEModel(){
  delete [] possible;
	delete [] freq_matrix;
  delete [] score_matrix;
  for(int i=0;i<seqset.num_seqs();i++){
    delete [] site_bias[i];
  }
	delete [] mean;
}

void SEModel::init(const vector<string>& seqs, float** exprtab, const int numexpr, const vector<string>& names, const int nc, const int bf, const double map_cut, const double sim_cut) {
	npossible = 0;
	ngenes = names.size();
	possible = new bool[ngenes];
	clear_all_possible();
	nameset = names;
	verbose = false;
	
	seqset.init(seqs);
	sites.init(seqs, nc, nc);
	select_sites.init(seqs, nc);
	print_sites.init(seqs, nc, nc);
	archive.init(sites, seqset,bf, map_cut, sim_cut);
	set_default_params();
  max_motifs = bf;
  map_cutoff = map_cut;
  sim_cutoff = sim_cut;
  freq_matrix = new int[sites.depth() * sites.ncols()];
  score_matrix = new double[sites.depth() * sites.ncols()];
  site_bias = new double*[seqset.num_seqs()];
	for(int i = 0;i < seqset.num_seqs();i++){
		site_bias[i]=new double[seqset.len_seq(i)];
  }
  
	expr = exprtab;
	npoints = numexpr;
	mean = new float[npoints];
	stdev = new float[npoints];
	pcorr = new float*[ngenes];
	for(int i = 0; i < ngenes; i++) {
		pcorr[i] = new float[ngenes];
		for(int j = 0; j < ngenes; j++) {
			pcorr[i][j] = -2;
		}
	}	
}

void SEModel::set_default_params(){
  separams.expect=10;
  separams.gcback=0.38;
  separams.minpass[0]=200;
  separams.seed=-1;
  separams.psfact=0.1;
  separams.weight=0.8; 
  separams.npass=1000000;
  separams.fragment=true;
  separams.flanking=0;
  separams.undersample=1;
  separams.oversample=1;
}

void SEModel::set_final_params(){
  separams.npseudo = separams.expect * separams.psfact;
  separams.backfreq[0]=separams.backfreq[5] = 0.25;
  separams.backfreq[1]=separams.backfreq[4] = (1 - separams.gcback)/2.0;
  separams.backfreq[2]=separams.backfreq[3] = separams.gcback/2.0;
  for(int i = 0; i < 6; i++){
    separams.pseudo[i] = separams.npseudo * separams.backfreq[i];
  }
  separams.maxlen = 3 * sites.width();
  separams.minpass[1] = 2 * separams.minpass[0];
  separams.minpass[2] = 3 * separams.minpass[0];
	set_cutoffs();
	separams.nruns = sites.positions_available() / separams.expect / sites.ncols() / separams.undersample * separams.oversample;
	separams.select = 5.0;
}

void SEModel::ace_initialize(){
  ran_int.set_seed(separams.seed);
  separams.seed = ran_int.seed();
  ran_int.set_range(0, RAND_MAX);
  ran_dbl.set_seed(ran_int.rnum());
  ran_dbl.set_range(0.0, 1.0);
}

void SEModel::add_possible(const int gene) {
	if (! possible[gene]) {
		possible[gene] = true;
		npossible++;
	}
}

void SEModel::remove_possible(const int gene) {
	if (possible[gene]) {
		possible[gene] = false;
		npossible--;
	}
}

void SEModel::clear_all_possible() {
	for(int g = 0; g < ngenes; g++) {
		possible[g] = false;
	}
	npossible = 0;
}

bool SEModel::is_possible(const int gene) const {
	return possible[gene];
}

int SEModel::possible_size() const {
	return npossible;
}

int SEModel::possible_positions() const {
	return sites.positions_available(possible);
}

void SEModel::clear_sites() {
	sites.clear_sites();
	select_sites.clear_sites();
}

int SEModel::size() const {
	return sites.seqs_with_sites();
}

int SEModel::sites_size() const {
	return sites.number();
}

bool SEModel::is_member(const int gene) const {
	return sites.seq_has_site(gene);
}

void SEModel::genes(int* genes) const {
	int count = 0;
	for (int g = 0; g < ngenes; g++) {
		if (is_member(g)) {
			genes[count] = g;
			count++;
		}
	}
	assert(count == size());
}

void SEModel::seed_random_site() {
	if(npossible < 1) return;
	
	sites.clear_sites();
	int chosen_possible, chosen_seq, chosen_posit;
  bool watson;
	
	chosen_seq = chosen_posit = -1;
	
	/* First choose a sequence */
	ran_int.set_range(0, possible_size() - 1);
	chosen_possible = ran_int.rnum();
	int g;
	for(g = 0; g < ngenes; g++) {
		if(is_possible(g) && chosen_possible == 0) break; 
		if(is_possible(g)) chosen_possible--;
	}
	chosen_seq = g;
	
	/* Now choose a site */
  for(int j = 0; j < 50; j++) {
		ran_int.set_range(0, seqset.len_seq(chosen_seq) - sites.width() - 1);
		double db = ran_dbl.rnum();//random (0,1)
		watson = (db > 0.5);
		chosen_posit = ran_int.rnum();
		if(watson && (chosen_posit > seqset.len_seq(chosen_seq) - sites.width() - 1)) continue;
		if((! watson) && (chosen_posit < sites.width())) continue;
		if(sites.is_open_site(chosen_seq, chosen_posit)) {
      sites.add_site(chosen_seq, chosen_posit, watson);
      break;
    }
  }
}

void SEModel::seed_biased_site(){
  int i,j,k,m;
  double ap=(double)separams.expect/(2.0*sites.positions_available());
  char **ss_seq;
  ss_seq=seqset.seq_ptr();
  double Lw,Lc,Pw,Pc,T,sc;
  int matpos,col;
  for(i=0;i<seqset.num_seqs();i++){
    for(j=0;j<seqset.len_seq(i)-sites.ncols()+1;j++){
      site_bias[i][j]=1.0;
    }
  }
  for(m=0;m<max_motifs;m++){
    sc=get_best_motif(m);
    if(sc<=0.0) break;
    calc_matrix();
    for(i=0;i<seqset.num_seqs();i++){
      for(j=0;j<seqset.len_seq(i)-sites.width()+1;j++){
				Lw=1.0;matpos=0;col=0;
				for(k=0;k<sites.ncols();k++){
					int seq=ss_seq[i][j+col];
					Lw*=score_matrix[matpos+seq];
					col=sites.next_column(col);
					matpos+=sites.depth();
				}
				Lc=1.0;matpos=sites.depth()-1;col=0;
				for(k=0;k<sites.ncols();k++){
					int seq=ss_seq[i][j+sites.width()-1-col];
					Lc*=score_matrix[matpos-seq];
					col=sites.next_column(col);
					matpos+=sites.depth();
				}
				Pw=Lw*ap/(1.0-ap+Lw*ap);
				Pc=Lc*ap/(1.0-ap+Lc*ap);
				//if(Pw<0.8&&Pc<0.8) continue;
				for(k=j;k<=(j+sites.width()-sites.ncols());k++){
					site_bias[i][k]*=(1-Pw)*(1-Pc);
					//site_bias[i][k]=0.0;
				}
      }
    }
  }
  if(m==0){
    seed_random_site();
    return;
  }
	
  for(T=0.0,i=0;i<seqset.num_seqs();i++){
    for(j=0;j<seqset.len_seq(i)-sites.ncols()+1;j++){
      T+=site_bias[i][j];
      //      if(m>3) cerr<<site_bias[i][j]<<'\t'<<T<<'\n';
    }
  }
	
  cerr<<"seeding "<<T<<" biased\n";
	
  sites.clear_sites();
	
  T*=ran_dbl.rnum();
  for(i=0;i<seqset.num_seqs();i++){
    for(j=0;j<seqset.len_seq(i)-sites.ncols()+1;j++){
      T-=site_bias[i][j];
      if(T<=0.0) break;
    }
    if(T<=0.0) break;
  }
	
  double db=ran_dbl.rnum();//random (0,1)
		bool watson=true;
		if(db>0.5) watson=false;
		sites.add_site(i,j,watson);
}

void SEModel::calc_matrix() {
  int d = sites.depth();
  double tot = (double) sites.number() + separams.npseudo;
  sites.calc_freq_matrix(seqset, freq_matrix);
  for(int i = 0; i < d * sites.ncols();i += d){
    score_matrix[i] = score_matrix[i + 5] = 1.0;
    for(int j = 1;j <= 4; j++){
      int x = freq_matrix[i] + freq_matrix[i + 5];
      score_matrix[i + j] = (freq_matrix[i + j] + x * separams.backfreq[j] 
				+ separams.pseudo[j])/(tot * separams.backfreq[j]);
    }
  }
}

void SEModel::single_pass(const double minprob) {
	double ap = (separams.weight * separams.expect * 10
							+ (1 - separams.weight) * sites.number())
							/(2.0 * sites.positions_available(possible));
  calc_matrix();
	sites.remove_all_sites();
  select_sites.remove_all_sites();
  //will only update once per pass
	
  char **ss_seq;
  ss_seq = seqset.seq_ptr();
  double Lw, Lc, Pw, Pc, F;
	int matpos, col;
  int considered = 0;
  int gadd = -1,jadd = -1;
  for(int g = 0; g < seqset.num_seqs(); g++){
		if (! is_possible(g)) continue;
		for(int j = 0; j <= seqset.len_seq(g) - sites.width(); j++){
			Lw = 1.0;
			matpos = 0;
			col = 0;
      for(int k = 0; k < sites.ncols(); k++){
				assert(j + col >= 0);
				assert(j + col <= seqset.len_seq(g));
				int seq = ss_seq[g][j + col];
				Lw *= score_matrix[matpos + seq];
				col = sites.next_column(col);
				matpos += sites.depth();
      }
      Lc = 1.0;
			matpos = sites.depth() - 1;
			col = 0;
      for(int k = 0; k < sites.ncols(); k++){
				assert(j + sites.width() - 1 - col >= 0);
				assert(j + sites.width() - 1 - col <= seqset.len_seq(g));
				int seq = ss_seq[g][j + sites.width() - 1 - col];
				Lc *= score_matrix[matpos - seq];
				col = sites.next_column(col);
				matpos += sites.depth();
      }
      Pw = Lw * ap/(1.0 - ap + Lw * ap);
      Pc = Lc * ap/(1.0 - ap + Lc * ap);
      F = Pw + Pc - Pw * Pc;//probability of either
			if(F > (minprob * 0.025/separams.select)) {
				select_sites.add_site(g, j, true);
			}
			//strand irrelevant for select_sites
			if(g == gadd && j < jadd + sites.width()) continue;
			if(F < minprob) continue;
			considered++;
			Pw = F * Pw / (Pw + Pc);
			Pc = F - Pw;
			double r = ran_dbl.rnum();
			if (r > F) continue;
			else if (r < Pw) {
				assert(j >= 0);
				assert(j <= seqset.len_seq(g) - sites.width());
				sites.add_site(g, j, true);
				gadd = g;
				jadd = j;
			} else {
				assert(j >= 0);
				assert(j <= seqset.len_seq(g) - sites.width());
				sites.add_site(g, j, false);
				gadd = g;
				jadd = j;
			}
    }
  }
  if(verbose)cerr<<"sampling "<<considered<<"/"<<select_sites.number()<<"\n";
}

void SEModel::single_pass_select(const double minprob){
	int i,j,k,n;
  double ap= (separams.weight * separams.expect + (1-separams.weight) * sites.number())
							/ (2.0 * sites.positions_available(possible));
  calc_matrix();
	calc_mean();
	sites.remove_all_sites();
  //will only update once per pass
	
  char **ss_seq;
  ss_seq=seqset.seq_ptr();
  double Lw, Lc, Pw, Pc, F;
  int matpos,col;
  int iadd=-1,jadd=-1;
  for(n=0;n<select_sites.number();n++){
    i=select_sites.chrom(n);
		if(! is_possible(i)) continue;
		j=select_sites.posit(n);
    if(i==iadd&&j<jadd+sites.width()) continue;
    if(j<0||j>seqset.len_seq(i)-sites.width()) continue;
    //could be screwed up with column sampling
    Lw=1.0;matpos=0;col=0;
    for(k=0;k<sites.ncols();k++){
      int seq=ss_seq[i][j+col];
      Lw*=score_matrix[matpos+seq];
      col=sites.next_column(col);
      matpos+=sites.depth();
    }
    Lc=1.0;matpos=sites.depth()-1;col=0;
    for(k=0;k<sites.ncols();k++){
      int seq=ss_seq[i][j+sites.width()-1-col];
      Lc*=score_matrix[matpos-seq];
      col=sites.next_column(col);
      matpos+=sites.depth();
    }
		Pw = Lw * ap / (1.0 - ap + Lw * ap);
    Pc = Lc * ap / (1.0 - ap + Lc * ap);
    F = (Pw + Pc - Pw * Pc);//probability of either
		if(F < minprob) continue;
		Pw = F * Pw/(Pw + Pc);
		Pc = F - Pw;
		double r=ran_dbl.rnum();
		if (r>F)
			continue;
    else if (r<Pw) {
      sites.add_site(i,j,true);
      add_possible(i);
			iadd=i;jadd=j;
    } else {
      sites.add_site(i,j,false);
			add_possible(i);
      iadd=i;jadd=j;
    }
  }
}

bool SEModel::column_sample(const int c, const bool sample){
	//cerr << "\t\t\tRunning column sample with c = " << c << "... ";
	//sample default to true, if false then replace with best column
  //just consider throwing out the worst column, sample for replacement, no need to sample both ways, unless column is specified
  int i, j;
  int col = 0, col_worst = 0, col_removed = 0;
  int *freq = new int[sites.depth()];
  double wt, wt_worst = DBL_MAX;
	
  if(c!=-1)
		col_worst = c;   // user chosen, hopefully a real column
  else {
    for(i = 0; i < sites.ncols(); i++){
      sites.column_freq(col, seqset, freq);
      wt = 0.0;
      for(j = 0; j < sites.depth(); j++) {
				wt += gammaln(freq[j] + separams.pseudo[j]);
				wt -= (double)freq[j] * log(separams.backfreq[j]);
      }
      if(wt < wt_worst) {
				col_worst = col;
				wt_worst = wt;
      }
      col = sites.next_column(col);
    }
  }
	
  col_removed = sites.remove_col(col_worst);
  i = select_sites.remove_col(col_worst);
  if(i != col_removed) {
		cerr << i << "  " << col_removed << " wrong assumption in column_sample\n";
		abort();
	}
	
  int max_left, max_right;
  max_left = max_right = (sites.max_width()-sites.width())/2;
  sites.columns_open(max_left,max_right);
  int cs_span=max_left+max_right+sites.width();
  double *wtx=new double[cs_span];
  int x=max_left;
  //wtx[x+c] will refer to the weight of pos c in the usual numbering
  col=0;
  double best_wt=0.0;
  for(i=0;i<cs_span;i++){
    wtx[i]=0.0;
    if((i-x)==col){
      col=sites.next_column(col);
      continue;
    }
    if(sites.column_freq(i-x,seqset,freq)){
      wt=0.0;
      for(j=0;j<sites.depth();j++){
				wt+=gammaln(freq[j]+separams.pseudo[j]);
				wt-=(double)freq[j]*log(separams.backfreq[j]);
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
		int newwidth=sites.width();
		if(i<x) newwidth+=(x-i);
		else if(i>(x+sites.width()-1)) newwidth+=(i-x-sites.width()+1);
		wtx[i]/=bico(newwidth-2,sites.ncols()-2);
		tot2+=wtx[i];
	}
	
	double pick;
	int col_pick;
	double cutoff=.01;
	if(sample){
		if(1-(wtx[x+col_removed]/tot2)<cutoff){
			sites.add_col(col_removed);
			select_sites.add_col(col_removed);
			delete [] wtx;
			delete [] freq;
			return false;
		} 
		pick=ran_dbl.rnum()*tot2;
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
	sites.add_col(col_pick);
	select_sites.add_col(col_pick);
	
	delete [] wtx;
	delete [] freq;
	
	// cerr << "done." <<  endl;
	
	if(col_removed==col_pick) return false;
  else return true;
}

double SEModel::map_score(){
  int i,j,k;
  double ms = 0.0;
  double map_N = sites.positions_available(possible);  
  double w = separams.weight/(1.0-separams.weight);
  double map_alpha = (double) separams.expect * w;
  double map_beta = map_N * w - map_alpha;
  double map_success = (double)sites.number();
  ms += ( gammaln(map_success+map_alpha)+gammaln(map_N-map_success+map_beta) );
  ms -= ( gammaln(map_alpha)+gammaln(map_N + map_beta) );
	
  sites.calc_freq_matrix(seqset,freq_matrix);
  double sc[6]={0.0,0.0,0.0,0.0,0.0,0.0};
  int d=sites.depth();
  for(k=0;k!=d*sites.ncols();k+=d) {
    int x=freq_matrix[k]+freq_matrix[k+5];
    for(j=1;j<=4;j++){
      ms += gammaln((double)freq_matrix[k+j]+x*separams.backfreq[j]+separams.pseudo[j]);
      sc[j]+=freq_matrix[k+j]+x*separams.backfreq[j];
    }
  }
  ms-=sites.ncols()*gammaln((double)sites.number()+separams.npseudo);
  for (k=1;k<=4;k++) {
    ms -= sc[k]*log(separams.backfreq[k]);
  }
  /*This factor arises from a modification of the model of Liu, et al in which the background frequencies of DNA bases are taken to be constant for the organism under consideration*/
  double  vg;
  ms -= lnbico(sites.width()-2, sites.ncols()-2);
  for(vg=0.0,k=1;k<=4;k++) {
    vg += gammaln(separams.pseudo[k]);
  }
  vg-=gammaln((double)(separams.npseudo));
  ms-=((double)sites.ncols()*vg);
  return ms;
}

double SEModel::get_best_motif(int i){
  return archive.return_best(sites,i);
}

void SEModel::optimize_columns(){
  while(column_sample(false)){}
  //replaces worst column with best column until these are the same
}

void SEModel::optimize_sites(){
  int poss_sites = sites.positions_available()+1;
  int *pos = new int[poss_sites];
  int *chr = new int[poss_sites];
  bool *str = new bool[poss_sites];
  Heap hp(poss_sites,3);
  double ap = (separams.weight * separams.expect
	            + (1 - separams.weight) * sites.number())
							/ (2.0 * sites.positions_available(possible));
	
  calc_matrix();
	calc_mean();
  char **ss_seq;
  ss_seq = seqset.seq_ptr();
  double Lw, Lc, Pw, Pc;
  int matpos, col;
  int h = 1;
  double cutoff = separams.sitecut[3];
  for(int i = 0; i < seqset.num_seqs(); i++) {
    if(! is_possible(i)) continue;
		for(int j = 0; j < seqset.len_seq(i) - sites.width() + 1; j++) {
      Lw = 1.0;
			matpos = 0;
			col = 0;
      for(int k = 0;k < sites.ncols(); k++){
				int seq = ss_seq[i][j + col];
				Lw *= score_matrix[matpos + seq];
				col = sites.next_column(col);
				matpos += sites.depth();
      }
      Lc = 1.0;
			matpos = sites.depth()-1;
			col = 0;
      for(int k = 0; k < sites.ncols(); k++){
				int seq = ss_seq[i][j + sites.width()-1-col];
				Lc *= score_matrix[matpos - seq];
				col = sites.next_column(col);
				matpos += sites.depth();
      }
      Pw = Lw * ap / (1.0 - ap + Lw * ap);
      Pc = Lc * ap / (1.0 - ap + Lc * ap);
      //F=(Pw+Pc-Pw*Pc);//probability of either
      if(Pw > cutoff * 0.025 || Pc > cutoff * 0.025) { 
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
  
  sites.remove_all_sites();
  Sites best = sites;
  double ms, ms_best = -DBL_MAX;
	for(;;){
    h = hp.del_min();
		if(h == -1) break;
    if(sites.is_open_site(chr[h], pos[h])) {
      sites.add_site(chr[h], pos[h], str[h]);
      if(sites.number() <= 3) continue;
      if(hp.value(h) < -0.7) continue;
      ms = map_score();
      if(ms > ms_best){
				ms_best = ms;
				best = sites;
      }
    }
  }
	sites = best;
  delete [] pos;
  delete [] chr;
  delete [] str;
}

void SEModel::orient_motif(){
  double *info = new double[6];
  double *freq = new double[6];
  for(int i = 0; i < 6; i++) info[i] = 0.0;
  int d = sites.depth();
	
  double tot = (double)sites.number() + separams.npseudo;
  sites.calc_freq_matrix(seqset, freq_matrix);
  for(int i = 0; i < d * sites.ncols(); i += d){
    double ii = 0.0;
    for(int j = 1; j <= 4; j++) {
      int x = freq_matrix[i] + freq_matrix[i+5];
      freq[j] = (freq_matrix[i + j] + x * separams.backfreq[j] + separams.pseudo[j]) / tot;
      ii += freq[j]*log(freq[j]);
    }
    ii = 2 + ii;
    for(int j = 1; j <= 4; j++) info[j] += freq[j] * ii;
  }
	
  double flip = 1.5 * info[3] + 1.0 * info[1] - 1.0 * info[4] - 1.5 * info[2];
  //for(i=1;i<5;i++) cerr<<info[i]<<'\t';
  //cerr<<flip<<'\n';
  if(flip < 0.0) sites.flip_sites();
  delete [] info;
  delete [] freq;
}

void SEModel::orient_print_motif(){
  double *info = new double[6];
  double *freq = new double[6];
  for(int i = 0; i < 6; i++) info[i] = 0.0;
  int d = print_sites.depth();
	
  double tot = (double) print_sites.number() + separams.npseudo;
  print_sites.calc_freq_matrix(seqset, freq_matrix);
  for(int i = 0; i < d * print_sites.ncols(); i += d){
    double ii = 0.0;
    for(int j = 1; j <= 4; j++) {
      int x = freq_matrix[i] + freq_matrix[i+5];
      freq[j] = (freq_matrix[i + j] + x * separams.backfreq[j] + separams.pseudo[j]) / tot;
      ii += freq[j]*log(freq[j]);
    }
    ii = 2 + ii;
    for(int j = 1; j <= 4; j++) info[j] += freq[j] * ii;
  }
	
  double flip = 1.5 * info[3] + 1.0 * info[1] - 1.0 * info[4] - 1.5 * info[2];
  //for(i=1;i<5;i++) cerr<<info[i]<<'\t';
  //cerr<<flip<<'\n';
  if(flip < 0.0) print_sites.flip_sites();
  delete [] info;
  delete [] freq;
}

string SEModel::consensus() {
	map<char,char> nt;
  nt[0]=nt[5]='N';
  nt[1]='A';nt[2]='C';nt[3]='G';nt[4]='T';
	char** ss_seq=seqset.seq_ptr();
	
	int numsites = sites.number();
	if(numsites < 1) return "";
	
	// cerr << "Computing consensus with " << numsites << " sites" << endl;
	vector<string> hits(numsites);
	for(int i = 0; i < numsites; i++){
		int c = sites.chrom(i);
    int p = sites.posit(i);
    bool s = sites.strand(i);
    for(int j = 0; j < sites.width(); j++){
      if(s) {
				if(p + j >= 0 && p+j < seqset.len_seq(c))
					 hits[i] += nt[ss_seq[c][p+j]];
				else hits[i] += ' ';
      } else {
				if(p + sites.width() - 1 - j >= 0 && p+sites.width() - 1 - j < seqset.len_seq(c))
					hits[i] += nt[sites.depth()-1-ss_seq[c][p+sites.width()-1-j]];
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

void SEModel::output(ostream &fout){
  map<char,char> nt;
  nt[0] = nt[5] = 'N';
  nt[1] = 'A';
	nt[2] = 'C';
	nt[3] = 'G';
	nt[4] = 'T';
  char** ss_seq = seqset.seq_ptr();
  int x = separams.flanking;
  for(int i = 0; i < print_sites.number(); i++){
    int c = print_sites.chrom(i);
    int p = print_sites.posit(i);
    bool s = print_sites.strand(i);
    for(int j = -x; j < print_sites.width() + x; j++){
      if(s) {
				if(p + j >= 0 && p + j < seqset.len_seq(c))
					fout << nt[ss_seq[c][p + j]];
				else fout << ' ';
      }
      else {
				if(p + print_sites.width() - 1 - j >= 0 && p + print_sites.width()-1-j < seqset.len_seq(c))
					fout << nt[print_sites.depth() - 1 - ss_seq[c][p + print_sites.width() - 1 - j]];
				else fout << ' ';
      }
    }
    fout << '\t' << c << '\t' << p << '\t' << s << '\n';
  }
  for(int i = 0; i < x; i++) fout << ' ';
  int j = 0;
	for(int i = 0;;){
    j = print_sites.next_column(i);
    fout << '*';
    if(i == print_sites.width() - 1) break;
    for(int k = 0; k < (j - i - 1); k++) fout << ' ';
    i = j;
  }
  fout << "\n";
}

void SEModel::full_output(ostream &fout){
  for(int j = 0; j < max_motifs; j++){
    double sc = archive.return_best(print_sites, j);
    orient_print_motif();
    if(sc > 0.0){
      fout << "Motif " << j + 1 << '\n';
      output(fout);
      fout << "MAP Score: " << sc << "\n\n";
    }
    else break;
  }
}

void SEModel::full_output(char *name){
  ofstream fout(name);
  full_output(fout);
}

void SEModel::output_params(ostream &fout){
  fout<<" expect =      \t"<<separams.expect<<'\n';
  fout<<" gcback =      \t"<<separams.gcback<<'\n';
  fout<<" minpass =     \t"<<separams.minpass[0]<<'\n';
  fout<<" seed =        \t"<<separams.seed<<'\n';
  fout<<" numcols =     \t"<<sites.ncols()<<'\n';
  fout<<" undersample = \t"<<separams.undersample<<'\n';
  fout<<" oversample = \t"<<separams.oversample<<'\n';
}

void SEModel::modify_params(int argc, char *argv[]){
  GetArg2(argc,argv,"-expect",separams.expect);
  GetArg2(argc,argv,"-gcback",separams.gcback);
  GetArg2(argc,argv,"-minpass",separams.minpass[0]);
  GetArg2(argc,argv,"-seed",separams.seed);
  GetArg2(argc,argv,"-undersample",separams.undersample);
  GetArg2(argc,argv,"-oversample",separams.oversample);
}

float* SEModel::get_mean() {
	return mean;
}

void SEModel::calc_mean() {
	int i, j;
	for (i = 0; i < npoints; i++) {
		mean[i] = 0;
		for (j = 0; j < ngenes; j++)
			if(is_member(j))
				mean[i] += expr[j][i];
		mean[i] /= size();
	}
}

void SEModel::calc_stdev() {
	calc_mean();
	for(int i = 0; i < npoints; i++) {
		stdev[i] = 0;
		for(int j = 0; j < ngenes; j++)
			if(is_member(j))
				stdev[i] += (expr[j][i] - mean[i]) * (expr[j][i] - mean[i]);
		stdev[i] /= size() - 1;
		stdev[i] = sqrt(stdev[i]);
	} 
}

float SEModel::get_pcorr(const int g1, const int g2) {
	if(pcorr[g1][g2] == -2)
		pcorr[g1][g2] = pcorr[g2][g1] = corr(expr[g1], expr[g2], npoints);
	return pcorr[g1][g2];
}

float SEModel::get_corr_with_mean(const float* pattern) const {
	return corr(mean, pattern, npoints);
}

float SEModel::prob_gene_given_model(int g) const {
	return 1.0;
	/*
	float p = 1;
	float c = 0.39894228; // = 1/sqrt(2Pi)
	for(int i = 0; i < npoints; i++)
		p *= c * exp(-((expr[g][i] - mean[i]) * (expr[g][i] - mean[i]))/(stdev[i] * stdev[i]))/stdev[i];
	return p;
	*/
}

void SEModel::expand_search_around_mean(const double corr_cutoff) {
	calc_mean();
	clear_all_possible();
	for(int g = 0; g < ngenes; g++)
		if(get_corr_with_mean(expr[g]) > corr_cutoff) add_possible(g);
}

void SEModel::expand_search_min_pcorr(const double corr_cutoff) {
 	// First add all genes to list of candidates
 	list<int> candidates;
 	for(int g = 0; g < ngenes; g++)
 		candidates.push_back(g);
 	
 	// Now run through list of genes with sites, and only keep genes within minimum corr_cutoff
 	for(int g = 0; g < ngenes; g++) {
 		if(! is_member(g)) continue;
 		list<int> survivors;
 		for(list<int>::iterator iter = candidates.begin(); iter != candidates.end(); iter++)
 			if(get_pcorr(g, *iter) > corr_cutoff) survivors.push_back(*iter);
 		candidates.assign(survivors.begin(), survivors.end());
 	}
 	
 	// Start with clean slate
 	clear_all_possible();
 	
 	// Add genes which already have sites to the search space
 	for(int g = 0; g < ngenes; g++)
 		if(is_member(g)) add_possible(g);
	for(list<int>::iterator iter = candidates.begin(); iter != candidates.end(); iter++)
		add_possible(*iter);
}

void SEModel::expand_search_avg_pcorr(const double corr_cutoff) {
	if(size() > 0) {
		float avg_pcorr;
		clear_all_possible();
		for(int g1 = 0; g1 < ngenes; g1++) {
			if(is_member(g1)) continue;
			avg_pcorr = 0;
			for(int g2 = 0; g2 < ngenes; g2++) {
				if(! is_member(g2)) continue;
				avg_pcorr += get_pcorr(g1, g2);
			}
			avg_pcorr /= size();
			if(avg_pcorr > corr_cutoff) add_possible(g1);
		}
		for(int g1 = 0; g1 < ngenes; g1++)
			if(is_member(g1)) add_possible(g1);
	}
}

void SEModel::search_for_motif(const double minsize, const double mincorr) {
	double corr_cutoff[4];
	corr_cutoff[0] = 0.70;
	corr_cutoff[3] = mincorr;
	corr_cutoff[2] = sqrt(corr_cutoff[0] * corr_cutoff[3]);
  corr_cutoff[1] = sqrt(corr_cutoff[0] * corr_cutoff[2]);
	double sc, cmp, sc_best_i;
  int i_worse = 0;
	sc_best_i = map_cutoff;
	i_worse = 0;
	int phase = 0;
	int old_phase = 0;
	
	sites.clear_sites();
	select_sites.clear_sites();
	Sites best_sites = sites;
	for(int g = 0; g < ngenes; g++)
		add_possible(g);
	seed_random_site();
	expand_search_avg_pcorr(corr_cutoff[0]);
	while(possible_size() < 2 && phase < 3) {
		phase++;
		old_phase = phase;
		print_status(cerr, 0, old_phase, corr_cutoff[old_phase], sc);
		expand_search_avg_pcorr(corr_cutoff[phase]);
	}
	if(possible_size() < 2) {
		cerr << "Bad search start -- no genes within " << mincorr << endl;
		return;
	}
	set_cutoffs();
	
	for(int i = 1; i <= separams.npass; i++){
		if(old_phase < phase) {
			print_status(cerr, i, old_phase, corr_cutoff[old_phase], sc);
			expand_search_avg_pcorr(corr_cutoff[phase]);
			if(possible_size() < 2 && phase < 3) { 
				phase++; 
				continue;
			}
			set_cutoffs();
			old_phase = phase;
		}
		if(phase == 3) {
			if(size() < minsize/2) {
				cerr << "\t\t\tReached phase " << phase << " with less than " << minsize/2 << " sequences with sites. Restarting..." << endl;
				break;
			}
			double sc1 = map_score();
			sc = 0.0;
			for(int z = 0; sc < sc1 && z < 5; z++){
				optimize_columns();
				optimize_sites();
				sc = map_score();
				print_status(cerr, i, phase, corr_cutoff[phase], sc);
			}
			if(sc < sc1) {
				sites = best_sites;
				sc = sc1;
			}
			sc = map_score();
			print_status(cerr, i, phase, corr_cutoff[phase], sc);
			if(size() < minsize) {
				cerr << "\t\t\tCompleted phase " << phase << " with less than " << minsize << " sequences with sites. Restarting..." << endl;
				break;
			}
			archive.consider_motif(sites, sc);
			cerr << "\t\t\tCompleted phase " << phase << "! Restarting..." << endl;
			break;
		}
		if(i_worse == 0)
			single_pass(separams.sitecut[phase]);
		else 
			single_pass_select(separams.sitecut[phase]);
		if(sites_size() == 0) {
			if(sc_best_i == map_cutoff) {
				cerr << "\t\t\tNo sites and best score matched cutoff! Restarting..." << endl;
				break;
			}
			//if(best_sites.number()<4) break;
			sites = best_sites;
			select_sites = best_sites;
			phase++;
			i_worse = 0;
			continue;
		}
		if(i<=3) continue;
		if(phase < 3 && size() > minsize) {
			if(column_sample(0)) {}
			if(column_sample(sites.width()-1)) {}
			for(int m = 0; m < 3; m++) {
				if(!(column_sample())) break;
			}
		}
		sc = map_score();
		if(sc - sc_best_i > 1e-3){
			i_worse=0;
			if(size() > minsize * 2) {
				cmp = archive.check_motif(sites, sc);
				if(cmp > sim_cutoff) {
					print_status(cerr, i, phase, corr_cutoff[phase], sc);
					cerr <<"\t\t\tToo similar! Restarting..." << endl;
					break;
				}
			}
			sc_best_i = sc;
			best_sites = sites;
		}
		else i_worse++;
		if(i_worse > separams.minpass[phase]){
			if(sc_best_i == map_cutoff) {
				print_status(cerr, i, phase, corr_cutoff[phase], sc);
				cerr << "\t\t\ti_worse is greater than cutoff and best score at cutoff! Restarting..." << endl;
				break;
			}
			if(best_sites.number() < 2) {
				print_status(cerr, i, phase, corr_cutoff[phase], sc);
				cerr << "\t\t\ti_worse is greater than cutoff and only 1 site! Restarting..." << endl;
				break;
			}
			sites = best_sites;
			select_sites = best_sites;
			phase++;
			i_worse = 0;
		}
		
		if(i % 50 == 0) print_status(cerr, i, phase, corr_cutoff[phase], sc);
	}
}

int SEModel::search_for_motif_near_seed(const double minsize, const double mincorr) {
	double corr_cutoff[4];
	corr_cutoff[0] = 0.70;
	corr_cutoff[3] = mincorr;
	corr_cutoff[2] = sqrt(corr_cutoff[0] * corr_cutoff[3]);
  corr_cutoff[1] = sqrt(corr_cutoff[0] * corr_cutoff[2]);
	double sc, cmp, sc_best_i;
  int i_worse = 0;
	sc_best_i = map_cutoff;
	i_worse = 0;
	int phase = 0;
	int old_phase = 0;
	
	Sites best_sites = sites;
	expand_search_avg_pcorr(corr_cutoff[0]);
	while(possible_size() < 2 && phase < 3) {
		phase++;
		old_phase = phase;
		print_status(cerr, 0, old_phase, corr_cutoff[old_phase], sc);
		expand_search_avg_pcorr(corr_cutoff[phase]);
	}
	if(possible_size() < 2) {
		cerr << "Bad search start -- no genes within " << mincorr << endl;
		return BAD_SEARCH_SPACE;
	}
	set_cutoffs();
	
	for(int i = 1; i <= separams.npass; i++){
		if(old_phase < phase) {
			print_status(cerr, i, old_phase, corr_cutoff[old_phase], sc);
			expand_search_avg_pcorr(corr_cutoff[phase]);
			if(possible_size() < 2 && phase < 3) { 
				phase++; 
				continue;
			}
			set_cutoffs();
			old_phase = phase;
		}
		if(phase == 3) {
			if(size() < minsize/2) {
				cerr << "\t\t\tReached phase " << phase << " with less than " << minsize/2 << " sequences with sites. Restarting..." << endl;
				return TOO_FEW_SITES;
			}
			double sc1 = map_score();
			sc = 0.0;
			for(int z = 0; sc < sc1 && z < 5; z++){
				optimize_columns();
				optimize_sites();
				sc = map_score();
				print_status(cerr, i, phase, corr_cutoff[phase], sc);
			}
			if(sc < sc1) {
				sites = best_sites;
				sc = sc1;
			}
			sc = map_score();
			print_status(cerr, i, phase, corr_cutoff[phase], sc);
			if(size() < minsize) {
				cerr << "\t\t\tCompleted phase " << phase << " with less than " << minsize << " sequences with sites. Restarting..." << endl;
				return TOO_FEW_SITES;
			}
			archive.consider_motif(sites, sc);
			cerr << "\t\t\tCompleted phase " << phase << "! Restarting..." << endl;
			return 0;
		}
		if(i_worse == 0)
			single_pass(separams.sitecut[phase]);
		else 
			single_pass_select(separams.sitecut[phase]);
		if(sites_size() == 0) {
			if(sc_best_i == map_cutoff) {
				cerr << "\t\t\tNo sites and best score matched cutoff! Restarting..." << endl;
				return BAD_SEED;
			}
			//if(best_sites.number()<4) break;
			sites = best_sites;
			select_sites = best_sites;
			phase++;
			i_worse = 0;
			continue;
		}
		if(i<=3) continue;
		if(phase < 3 && size() > minsize) {
			if(column_sample(0)) {}
			if(column_sample(sites.width()-1)) {}
			for(int m = 0; m < 3; m++) {
				if(!(column_sample())) break;
			}
		}
		sc = map_score();
		if(sc - sc_best_i > 1e-3){
			i_worse=0;
			if(size() > minsize * 2) {
				cmp = archive.check_motif(sites, sc);
				if(cmp > sim_cutoff) {
					print_status(cerr, i, phase, corr_cutoff[phase], sc);
					cerr <<"\t\t\tToo similar! Restarting..." << endl;
					return TOO_SIMILAR;
				}
			}
			sc_best_i = sc;
			best_sites = sites;
		}
		else i_worse++;
		if(i_worse > separams.minpass[phase]){
			if(sc_best_i == map_cutoff) {
				print_status(cerr, i, phase, corr_cutoff[phase], sc);
				cerr << "\t\t\ti_worse is greater than cutoff and best score at cutoff! Restarting..." << endl;
				return BAD_SEED;
			}
			if(best_sites.number() < 2) {
				print_status(cerr, i, phase, corr_cutoff[phase], sc);
				cerr << "\t\t\ti_worse is greater than cutoff and only 1 site! Restarting..." << endl;
				return BAD_SEED;
			}
			sites = best_sites;
			select_sites = best_sites;
			phase++;
			i_worse = 0;
		}
		
		if(i % 50 == 0) print_status(cerr, i, phase, corr_cutoff[phase], sc);
	}
	
	return 0;
}

void SEModel::print_status(ostream& out, const int i, const int phase, const double cutoff, const double sc) {
	if(size() > 0) {
		out << "\t\t\t" << setw(5) << i;
		out << setw(3) << phase;
		int prec = cerr.precision(2);
		out << setw(5) << setprecision(2) << cutoff;
		out << setw(10) << setprecision(2) << separams.sitecut[phase];
		cerr.precision(prec);
		out << setw(5) << sites_size();
		out << setw(5) << size();
		out << setw(5) << possible_size();
		out << setw(40) << consensus();
		out << setw(10) << sc;
		out << endl;
	} else {
		out << "\t\t\tNo sites!" << endl;
	}
}

void SEModel::set_cutoffs() {
	double ap = (double) separams.expect/(2.0 * sites.positions_available(possible));
	separams.sitecut[0] = 5.0 * ap;//10.0*ap;
	separams.sitecut[3] = 0.2;
	separams.sitecut[2] = 0.2;
	separams.sitecut[1] = sqrt(separams.sitecut[0] * separams.sitecut[2]);
}

void SEModel::print_possible(ostream& out) {
	for(int g = 0; g < ngenes; g++)
		if(is_possible(g)) out << nameset[g] << endl;
}

