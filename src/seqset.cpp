//Copyright 1999 President and Fellows of Harvard University
//seqset.cpp

#include "seqset.h"

Seqset::Seqset() : 
ss_num_seqs(0),
total_seq_len(0),
ss_seq(ss_num_seqs),
gc_genome(0),
gc(ss_num_seqs),
bgmodel0(4),
bgmodel1(16),
bgmodel2(64),
bgmodel3(256),
wbgscores(ss_num_seqs),
cbgscores(ss_num_seqs) {
}

Seqset::Seqset(const vector<string>& v) :
ss_num_seqs(v.size()),
total_seq_len(0),
ss_seq(ss_num_seqs),
gc_genome(0),
gc(ss_num_seqs),
bgmodel0(4),
bgmodel1(16),
bgmodel2(64),
bgmodel3(256),
wbgscores(ss_num_seqs),
cbgscores(ss_num_seqs) {
	map<char, int> nt;
	nt['A'] = nt['a'] = 0;
	nt['C'] = nt['c'] = 1;
	nt['G'] = nt['g'] = 2;
	nt['T'] = nt['t'] = 3;
  int len;
	for(int i = 0; i < ss_num_seqs; i++) {
    gc[i] = 0.0;
		len = v[i].length();
		ss_seq[i].reserve(len);
    for(int j = 0; j < len; j++) {
			ss_seq[i].push_back(nt[v[i][j]]);
      if(v[i][j] == 'C' || v[i][j] == 'G') {
				gc[i]++;
			}
    }
		gc_genome += gc[i];
		gc[i] /= len;
		total_seq_len += len;
  }
	gc_genome /= total_seq_len;
	train_background3();
	train_background2();
	train_background1();
	train_background0();
	calc_bg_scores();
}

Seqset::~Seqset(){
}

float Seqset::bgscore(const int g, const int p, const bool s) {
	if(s) {
		return wbgscores[g][p];
	} else {
		return cbgscores[g][p];
	}
}

void Seqset::train_background3() {
	for(int i = 0; i < 256; i++) {
		bgmodel3[i] = 0;
	}
	
	// Add pseudocounts
	for(int i = 0; i < 256; i++) {
		float atpseudo = 10 * (1 - gc_genome)/2;
		float gcpseudo = 10 * gc_genome/2;
		switch(i % 4) {
			case 1:
			case 4:
				bgmodel3[i] += atpseudo;
				break;
			case 2:
			case 3:
				bgmodel3[i] += gcpseudo;
				break;
		}
	}
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = ss_seq[i].size();
		// Compute counts for forward strand
		for(int j = 3; j < len; j++) {
			bgmodel3[ss_seq[i][j - 3] * 64
								+ ss_seq[i][j - 2] * 16 
								+ ss_seq[i][j - 1] * 4
								+ ss_seq[i][j]]++;
		}
		// Compute counts for reverse strand
		for(int j = len - 4; j >= 0; j--) {
			bgmodel3[(3 - ss_seq[i][j + 3]) * 64
								+ (3 - ss_seq[i][j + 2]) * 16
								+ (3 - ss_seq[i][j + 1]) * 4
								+ (3 - ss_seq[i][j])]++;
		}

	}
	
	// Normalize
	float total;
	for(int i = 0; i < 64; i++) {
		total = bgmodel3[4 * i] + bgmodel3[4 * i + 1] + bgmodel3[4 * i + 2] + bgmodel3[4 * i + 3];
		bgmodel3[4 * i] /= total;
		bgmodel3[4 * i + 1] /= total;
		bgmodel3[4 * i + 2] /= total;
		bgmodel3[4 * i + 3] /= total;
	}
}

void Seqset::train_background2() {
	for(int i = 0; i < 64; i++) {
		bgmodel2[i] = 0;
	}
	
	// Add pseudocounts
	for(int i = 0; i < 64; i++) {
		float atpseudo = 10 * (1 - gc_genome)/2;
		float gcpseudo = 10 * gc_genome/2;
		switch(i % 4) {
			case 1:
			case 4:
				bgmodel2[i] += atpseudo;
				break;
			case 2:
			case 3:
				bgmodel2[i] += gcpseudo;
				break;
		}
	}
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = ss_seq[i].size();
		// Compute counts for forward strand
		for(int j = 2; j < len; j++) {
			bgmodel2[ss_seq[i][j - 2] * 16 
								+ ss_seq[i][j - 1] * 4
								+ ss_seq[i][j]]++;
		}
		// Compute counts for reverse strand
		for(int j = len - 3; j >= 0; j--) {
			bgmodel2[(3 - ss_seq[i][j + 2]) * 16
								+ (3 - ss_seq[i][j + 1]) * 4
								+ (3 - ss_seq[i][j])]++;
		}

	}
	
	// Normalize
	float total;
	for(int i = 0; i < 16; i++) {
		total = bgmodel2[4 * i] + bgmodel2[4 * i + 1] + bgmodel2[4 * i + 2] + bgmodel2[4 * i + 3];
		bgmodel2[4 * i] /= total;
		bgmodel2[4 * i + 1] /= total;
		bgmodel2[4 * i + 2] /= total;
		bgmodel2[4 * i + 3] /= total;
	}
}

void Seqset::train_background1() {
	for(int i = 0; i < 16; i++) {
		bgmodel1[i] = 0;
	}
	
	// Add pseudocounts
	for(int i = 0; i < 16; i++) {
		float atpseudo = 10 * (1 - gc_genome)/2;
		float gcpseudo = 10 * gc_genome/2;
		switch(i % 4) {
			case 1:
			case 4:
				bgmodel1[i] += atpseudo;
				break;
			case 2:
			case 3:
				bgmodel1[i] += gcpseudo;
				break;
		}
	}
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = ss_seq[i].size();
		// Compute counts for forward strand
		for(int j = 1; j < len; j++) {
			bgmodel1[ss_seq[i][j - 1] * 4
								+ ss_seq[i][j]]++;
		}
		// Compute counts for reverse strand
		for(int j = len - 2; j >= 0; j--) {
			bgmodel1[(3 - ss_seq[i][j + 1]) * 4
								+ (3 - ss_seq[i][j])]++;
		}

	}
	
	// Normalize
	float total;
	for(int i = 0; i < 4; i++) {
		total = bgmodel1[4 * i] + bgmodel1[4 * i + 1] + bgmodel1[4 * i + 2] + bgmodel1[4 * i + 3];
		bgmodel1[4 * i] /= total;
		bgmodel1[4 * i + 1] /= total;
		bgmodel1[4 * i + 2] /= total;
		bgmodel1[4 * i + 3] /= total;
	}
}

void Seqset::train_background0() {
	for(int i = 0; i < 4; i++) {
		bgmodel0[i] = 0;
	}
	
	// Just use genome-wide GC content
	bgmodel0[0] = (1 - gc_genome)/2;
	bgmodel0[1] = gc_genome/2;
	bgmodel0[2] = bgmodel0[1];
	bgmodel0[3] = bgmodel0[0];
}

void Seqset::calc_bg_scores() {
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = ss_seq[i].size();
		wbgscores[i].reserve(len);
		
		// Use lower order models for first few Watson bases
		wbgscores[i].push_back(bgmodel0[ss_seq[i][0]]);
		wbgscores[i].push_back(bgmodel1[ss_seq[i][0] * 4 
															+ ss_seq[i][1]]);
		wbgscores[i].push_back(bgmodel2[ss_seq[i][0] * 16 
															+ ss_seq[i][1] * 4 
															+ ss_seq[i][2]]);
				
		// Use third-order model for most bases
		for(int j = 3; j < len; j++) {
			wbgscores[i].push_back(bgmodel3[ss_seq[i][j - 3] * 64
																+ ss_seq[i][j - 2] * 16 
																+ ss_seq[i][j - 1] * 4
																+ ss_seq[i][j]]);
		}
		cbgscores[i].reserve(len);
		for(int j = len - 4; j >= 0; j--) {
			cbgscores[i].push_back(bgmodel3[(3 - ss_seq[i][j + 3]) * 64
																+ (3 - ss_seq[i][j + 2]) * 16
																+ (3 - ss_seq[i][j + 1]) * 4
																+ (3 - ss_seq[i][j])]);
		}
		
		// Use lower order models for first few Watson bases
		cbgscores[i].push_back(bgmodel2[(3 - ss_seq[i][len - 1]) * 16 
															+ (3 - ss_seq[len - 2][1]) * 4 
															+ (3 - ss_seq[len - 1][2])]);
		cbgscores[i].push_back(bgmodel2[(3 - ss_seq[len - 2][1]) * 4 
															+ (3 - ss_seq[len - 1][2])]);
		cbgscores[i].push_back(bgmodel0[3 - ss_seq[len - 1][2]]);
		
	}
}
