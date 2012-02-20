#include "bgmodel.h"

BGModel::BGModel(const Seqset& s, const int ord) :
seqset(s),
total_seq_len(0),
order(ord),
gc_genome(0),
gc(seqset.num_seqs()),
model0(4),
model1(16),
model2(64),
model3(256),
model4(1024),
model5(4096),
wbgscores(seqset.num_seqs()),
cbgscores(seqset.num_seqs()),
cumulscores(seqset.num_seqs()) {
	const vector<vector <char> >& seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		gc[i] = 0.0;
		len = seqset.len_seq(i);
		for(int j = 0; j < len; j++) {
			if(seq[i][j] == 1 || seq[i][j] == 2) {
				gc[i]++;
			}
		}
		gc_genome += gc[i];
		gc[i] /= len;
		total_seq_len += len;
	}
	gc_genome /= total_seq_len;
	
	train_background[0] = &BGModel::train_background_0;
	train_background[1] = &BGModel::train_background_1;
	train_background[2] = &BGModel::train_background_2;
	train_background[3] = &BGModel::train_background_3;
	train_background[4] = &BGModel::train_background_4;
	train_background[5] = &BGModel::train_background_5;
	
	calc_bg_scores[0] = &BGModel::calc_bg_scores_0;
	calc_bg_scores[1] = &BGModel::calc_bg_scores_1;
	calc_bg_scores[2] = &BGModel::calc_bg_scores_2;
	calc_bg_scores[3] = &BGModel::calc_bg_scores_3;
	calc_bg_scores[4] = &BGModel::calc_bg_scores_4;
	calc_bg_scores[5] = &BGModel::calc_bg_scores_5;
	
	for(int i = 0; i <= order; i++) {
		(*this.*train_background[i])();
	}
	
	(*this.*calc_bg_scores[order])();
	
	float last_cumul_score = 0.0;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = seqset.len_seq(i);
		cumulscores[i].push_back(last_cumul_score - wbgscores[i][0] - cbgscores[i][0]);
		for(int j = 1; j < len; j++)
			cumulscores[i].push_back(cumulscores[i][j-1] - wbgscores[i][j] - cbgscores[i][j]);
		assert(cumulscores[i].size() == (unsigned int) len);
		last_cumul_score = cumulscores[i][len - 1];
	}
}

double BGModel::score_site(vector<int>::const_iterator first_col, vector<int>::const_iterator last_col, const int width,
													 const int c, const int p, const bool s) const {
	double L = 0.0;
	int matpos, len;
	vector<int>::const_iterator col_iter = first_col;
	if(s) {
		matpos = 0;
		len = seqset.len_seq(c);
		for(; col_iter != last_col; ++col_iter) {
			assert(p + *col_iter >= 0);
			assert(p + *col_iter < len);
			L += wbgscores[c][p + *col_iter];
			matpos += 4;
		}
	} else {
		matpos = 0;
		len = seqset.len_seq(c);
		for(; col_iter != last_col; ++col_iter) {
			assert(p + width - 1 - *col_iter >= 0);
			assert(p + width - 1 - *col_iter < len);
			L += cbgscores[c][p + width - 1 - *col_iter];
			matpos += 4;
		}
	}
	return L;
}

void BGModel::train_background_5() {
	const vector<vector <char> >& seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	for(int i = 0; i < 4096; i++) {
		model5[i] = 0;
	}
	
	// Add pseudocounts
	float atpseudo = 10 * (1 - gc_genome)/2;
	float gcpseudo = 10 * gc_genome/2;
	for(int i = 0; i < 4096; i++) {
		switch(i % 4) {
			case 1:
			case 4:
				model5[i] += atpseudo;
				break;
			case 2:
			case 3:
				model5[i] += gcpseudo;
				break;
		}
	}
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = seq[i].size();
		// Compute counts for forward strand
		for(int j = 5; j < len; j++) {
			model5[seq[i][j - 5] * 1024
				+ seq[i][j - 4] * 256
				+ seq[i][j - 3] * 64
				+ seq[i][j - 2] * 16 
				+ seq[i][j - 1] * 4
				+ seq[i][j]]++;
		}
		// Compute counts for reverse strand
		for(int j = len - 6; j >= 0; j--) {
			model5[(3 - seq[i][j + 5]) * 1024
				+ (3 - seq[i][j + 4]) * 256
				+ (3 - seq[i][j + 3]) * 64
				+ (3 - seq[i][j + 2]) * 16
				+ (3 - seq[i][j + 1]) * 4
				+ (3 - seq[i][j])]++;
		}

	}
	
	// Normalize
	float total;
	for(int i = 0; i < 1024; i++) {
		total = model5[4 * i] + model5[4 * i + 1] + model5[4 * i + 2] + model5[4 * i + 3];
		model5[4 * i] /= total;
		model5[4 * i + 1] /= total;
		model5[4 * i + 2] /= total;
		model5[4 * i + 3] /= total;
	}
}

void BGModel::train_background_4() {
	const vector<vector <char> >& seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	for(int i = 0; i < 1024; i++) {
		model4[i] = 0;
	}
	
	// Add pseudocounts
	float atpseudo = 10 * (1 - gc_genome)/2;
	float gcpseudo = 10 * gc_genome/2;
	for(int i = 0; i < 1024; i++) {
		switch(i % 4) {
			case 1:
			case 4:
				model4[i] += atpseudo;
				break;
			case 2:
			case 3:
				model4[i] += gcpseudo;
				break;
		}
	}
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = seq[i].size();
		// Compute counts for forward strand
		for(int j = 4; j < len; j++) {
			model4[seq[i][j - 4] * 256
				+ seq[i][j - 3] * 64
				+ seq[i][j - 2] * 16 
				+ seq[i][j - 1] * 4
				+ seq[i][j]]++;
		}
		// Compute counts for reverse strand
		for(int j = len - 5; j >= 0; j--) {
			model4[(3 - seq[i][j + 4]) * 256
				+ (3 - seq[i][j + 3]) * 64
				+ (3 - seq[i][j + 2]) * 16
				+ (3 - seq[i][j + 1]) * 4
				+ (3 - seq[i][j])]++;
		}

	}
	
	// Normalize
	float total;
	for(int i = 0; i < 256; i++) {
		total = model4[4 * i] + model4[4 * i + 1] + model4[4 * i + 2] + model4[4 * i + 3];
		model4[4 * i] /= total;
		model4[4 * i + 1] /= total;
		model4[4 * i + 2] /= total;
		model4[4 * i + 3] /= total;
	}
}

void BGModel::train_background_3() {
	const vector<vector <char> >& seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	for(int i = 0; i < 256; i++) {
		model3[i] = 0;
	}
	
	// Add pseudocounts
	float atpseudo = 10 * (1 - gc_genome)/2;
	float gcpseudo = 10 * gc_genome/2;
	for(int i = 0; i < 256; i++) {
		switch(i % 4) {
			case 1:
			case 4:
				model3[i] += atpseudo;
				break;
			case 2:
			case 3:
				model3[i] += gcpseudo;
				break;
		}
	}
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = seq[i].size();
		// Compute counts for forward strand
		for(int j = 3; j < len; j++) {
			model3[seq[i][j - 3] * 64
								+ seq[i][j - 2] * 16 
								+ seq[i][j - 1] * 4
								+ seq[i][j]]++;
		}
		// Compute counts for reverse strand
		for(int j = len - 4; j >= 0; j--) {
			model3[(3 - seq[i][j + 3]) * 64
								+ (3 - seq[i][j + 2]) * 16
								+ (3 - seq[i][j + 1]) * 4
								+ (3 - seq[i][j])]++;
		}

	}
	
	// Normalize
	float total;
	for(int i = 0; i < 64; i++) {
		total = model3[4 * i] + model3[4 * i + 1] + model3[4 * i + 2] + model3[4 * i + 3];
		model3[4 * i] /= total;
		model3[4 * i + 1] /= total;
		model3[4 * i + 2] /= total;
		model3[4 * i + 3] /= total;
	}
}

void BGModel::train_background_2() {
	const vector<vector <char> >& seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	for(int i = 0; i < 64; i++) {
		model2[i] = 0;
	}
	
	// Add pseudocounts
	float atpseudo = 10 * (1 - gc_genome)/2;
	float gcpseudo = 10 * gc_genome/2;
	for(int i = 0; i < 64; i++) {
		switch(i % 4) {
			case 1:
			case 4:
				model2[i] += atpseudo;
				break;
			case 2:
			case 3:
				model2[i] += gcpseudo;
				break;
		}
	}
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = seq[i].size();
		// Compute counts for forward strand
		for(int j = 2; j < len; j++) {
			model2[seq[i][j - 2] * 16 
								+ seq[i][j - 1] * 4
								+ seq[i][j]]++;
		}
		// Compute counts for reverse strand
		for(int j = len - 3; j >= 0; j--) {
			model2[(3 - seq[i][j + 2]) * 16
								+ (3 - seq[i][j + 1]) * 4
								+ (3 - seq[i][j])]++;
		}

	}
	
	// Normalize
	float total;
	for(int i = 0; i < 16; i++) {
		total = model2[4 * i] + model2[4 * i + 1] + model2[4 * i + 2] + model2[4 * i + 3];
		model2[4 * i] /= total;
		model2[4 * i + 1] /= total;
		model2[4 * i + 2] /= total;
		model2[4 * i + 3] /= total;
	}
}

void BGModel::train_background_1() {
	const vector<vector <char> >& seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	for(int i = 0; i < 16; i++) {
		model1[i] = 0;
	}
	
	// Add pseudocounts
	float atpseudo = 10 * (1 - gc_genome)/2;
	float gcpseudo = 10 * gc_genome/2;
	for(int i = 0; i < 16; i++) {
		switch(i % 4) {
			case 1:
			case 4:
				model1[i] += atpseudo;
				break;
			case 2:
			case 3:
				model1[i] += gcpseudo;
				break;
		}
	}
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = seq[i].size();
		// Compute counts for forward strand
		for(int j = 1; j < len; j++) {
			model1[seq[i][j - 1] * 4
								+ seq[i][j]]++;
		}
		// Compute counts for reverse strand
		for(int j = len - 2; j >= 0; j--) {
			model1[(3 - seq[i][j + 1]) * 4
								+ (3 - seq[i][j])]++;
		}

	}
	
	// Normalize
	float total;
	for(int i = 0; i < 4; i++) {
		total = model1[4 * i] + model1[4 * i + 1] + model1[4 * i + 2] + model1[4 * i + 3];
		model1[4 * i] /= total;
		model1[4 * i + 1] /= total;
		model1[4 * i + 2] /= total;
		model1[4 * i + 3] /= total;
	}
}

void BGModel::train_background_0() {
	for(int i = 0; i < 4; i++) {
		model0[i] = 0;
	}
	
	// Just use genome-wide GC content
	model0[0] = (1 - gc_genome)/2;
	model0[1] = gc_genome/2;
	model0[2] = model0[1];
	model0[3] = model0[0];
}

void BGModel::calc_bg_scores_5() {
	const vector<vector <char> >& seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = seq[i].size();
		wbgscores[i].clear();
		cbgscores[i].clear();
		wbgscores[i].reserve(len);
		cbgscores[i].reserve(len);
		
		// Use lower order models for first five Watson bases
		wbgscores[i].push_back(log(model0[seq[i][0]]));
		wbgscores[i].push_back(log(model1[seq[i][0] * 4 
															+ seq[i][1]]));
		wbgscores[i].push_back(log(model2[seq[i][0] * 16 
															+ seq[i][1] * 4 
															+ seq[i][2]]));
		if(len > 3) {
			wbgscores[i].push_back(log(model3[seq[i][0] * 64 
																+ seq[i][1] * 16 
																+ seq[i][2] * 4
																+ seq[i][3]]));
		}
		if(len > 4) {
			wbgscores[i].push_back(log(model4[seq[i][0] * 256
																+ seq[i][1] * 64 
																+ seq[i][2] * 16 
																+ seq[i][3] * 4
																+ seq[i][4]]));
		}
		
		// Use fifth-order model for most bases
		for(int j = 5; j < len; j++) {
			wbgscores[i].push_back(log(model5[seq[i][j - 5] * 1024
																+ seq[i][j - 4] * 256
			                          + seq[i][j - 3] * 64
																+ seq[i][j - 2] * 16 
																+ seq[i][j - 1] * 4
																+ seq[i][j]]));
		}
		for(int j = 0; j < len - 5; j++) {
			cbgscores[i].push_back(log(model5[(3 - seq[i][j + 5]) * 1024
																+ (3 - seq[i][j + 4]) * 256
			                          + (3 - seq[i][j + 3]) * 64
																+ (3 - seq[i][j + 2]) * 16
																+ (3 - seq[i][j + 1]) * 4
																+ (3 - seq[i][j])]));
		}
		
		// Use lower order models for last four Crick bases
		if(len > 4) {
			cbgscores[i].push_back(log(model4[(3 - seq[i][len - 1]) * 256
																+ (3 - seq[i][len - 2]) * 64
																+ (3 - seq[i][len - 3]) * 16
																+ (3 - seq[i][len - 4]) * 4
																+ (3 - seq[i][len - 5])]));
		}
		if(len > 3) {
		cbgscores[i].push_back(log(model3[(3 - seq[i][len - 1]) * 64
															+ (3 - seq[i][len - 2]) * 16
															+ (3 - seq[i][len - 3]) * 4
															+ (3 - seq[i][len - 4])]));
		}
		cbgscores[i].push_back(log(model2[(3 - seq[i][len - 1]) * 16 
															+ (3 - seq[i][len - 2]) * 4 
															+ (3 - seq[i][len - 3])]));
		cbgscores[i].push_back(log(model1[(3 - seq[i][len - 1]) * 4 
															+ (3 - seq[i][len - 2])]));
		cbgscores[i].push_back(log(model0[3 - seq[i][len - 1]]));
	
		assert(wbgscores[i].size() == (unsigned int) len);
		assert(cbgscores[i].size() == (unsigned int) len);
	}
}

void BGModel::calc_bg_scores_4() {
	const vector<vector <char> >& seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = seq[i].size();
		wbgscores[i].clear();
		cbgscores[i].clear();
		wbgscores[i].reserve(len);
		cbgscores[i].reserve(len);
		
		// Use lower order models for first four Watson bases
		wbgscores[i].push_back(log(model0[seq[i][0]]));
		wbgscores[i].push_back(log(model1[seq[i][0] * 4 
															+ seq[i][1]]));
		wbgscores[i].push_back(log(model2[seq[i][0] * 16 
															+ seq[i][1] * 4 
															+ seq[i][2]]));
		if(len > 4) {
			wbgscores[i].push_back(log(model3[seq[i][0] * 64 
																+ seq[i][1] * 16 
																+ seq[i][2] * 4
																+ seq[i][3]]));
		}
				
		// Use third-order model for most bases
		for(int j = 4; j < len; j++) {
			wbgscores[i].push_back(log(model4[seq[i][j - 4] * 256
			                          + seq[i][j - 3] * 64
																+ seq[i][j - 2] * 16 
																+ seq[i][j - 1] * 4
																+ seq[i][j]]));
		}
		for(int j = 0; j < len - 4; j++) {
			cbgscores[i].push_back(log(model4[(3 - seq[i][j + 4]) * 256
			                          + (3 - seq[i][j + 3]) * 64
																+ (3 - seq[i][j + 2]) * 16
																+ (3 - seq[i][j + 1]) * 4
																+ (3 - seq[i][j])]));
		}
		
		// Use lower order models for last four Crick bases
		if(len > 4) {
			cbgscores[i].push_back(log(model3[(3 - seq[i][len - 1]) * 64
																+ (3 - seq[i][len - 2]) * 16
																+ (3 - seq[i][len - 3]) * 4
																+ (3 - seq[i][len - 4])]));
		}
		cbgscores[i].push_back(log(model2[(3 - seq[i][len - 1]) * 16 
															+ (3 - seq[i][len - 2]) * 4 
															+ (3 - seq[i][len - 3])]));
		cbgscores[i].push_back(log(model1[(3 - seq[i][len - 1]) * 4 
															+ (3 - seq[i][len - 2])]));
		cbgscores[i].push_back(log(model0[3 - seq[i][len - 1]]));
	
		assert(wbgscores[i].size() == (unsigned int) len);
		assert(cbgscores[i].size() == (unsigned int) len);
	}
}

void BGModel::calc_bg_scores_3() {
	const vector<vector <char> >& seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = seq[i].size();
		wbgscores[i].clear();
		cbgscores[i].clear();
		wbgscores[i].reserve(len);
		cbgscores[i].reserve(len);
		
		// Use lower order models for first three Watson bases
		wbgscores[i].push_back(log(model0[seq[i][0]]));
		wbgscores[i].push_back(log(model1[seq[i][0] * 4 
															+ seq[i][1]]));
		wbgscores[i].push_back(log(model2[seq[i][0] * 16 
															+ seq[i][1] * 4 
															+ seq[i][2]]));
				
		// Use third-order model for most bases
		for(int j = 3; j < len; j++) {
			wbgscores[i].push_back(log(model3[seq[i][j - 3] * 64
																+ seq[i][j - 2] * 16 
																+ seq[i][j - 1] * 4
																+ seq[i][j]]));
		}
		for(int j = 0; j < len - 3; j++) {
			cbgscores[i].push_back(log(model3[(3 - seq[i][j + 3]) * 64
																+ (3 - seq[i][j + 2]) * 16
																+ (3 - seq[i][j + 1]) * 4
																+ (3 - seq[i][j])]));
		}
		
		// Use lower order models for last three Crick bases
		cbgscores[i].push_back(log(model2[(3 - seq[i][len - 1]) * 16 
															+ (3 - seq[i][len - 2]) * 4 
															+ (3 - seq[i][len - 3])]));
		cbgscores[i].push_back(log(model1[(3 - seq[i][len - 1]) * 4 
															+ (3 - seq[i][len - 2])]));
		cbgscores[i].push_back(log(model0[3 - seq[i][len - 1]]));
	
		assert(wbgscores[i].size() == (unsigned int) len);
		assert(cbgscores[i].size() == (unsigned int) len);
	}
}

void BGModel::calc_bg_scores_2() {
	const vector<vector <char> >& seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = seq[i].size();
		wbgscores[i].clear();
		cbgscores[i].clear();
		wbgscores[i].reserve(len);
		cbgscores[i].reserve(len);
		
		// Use lower order models for first two Watson bases
		wbgscores[i].push_back(log(model0[seq[i][0]]));
		wbgscores[i].push_back(log(model1[seq[i][0] * 4 
															+ seq[i][1]]));
				
		// Use 2nd order model for most bases
		for(int j = 2; j < len; j++) {
			wbgscores[i].push_back(log(model2[seq[i][j - 2] * 16 
																+ seq[i][j - 1] * 4
																+ seq[i][j]]));
		}
		for(int j = 0; j < len - 2; j++) {
			cbgscores[i].push_back(log(model2[(3 - seq[i][j + 2]) * 16
																+ (3 - seq[i][j + 1]) * 4
																+ (3 - seq[i][j])]));
		}
		
		// Use lower order models for last two Crick bases
		cbgscores[i].push_back(log(model1[(3 - seq[i][len - 1]) * 4 
															+ (3 - seq[i][len - 2])]));
		cbgscores[i].push_back(log(model0[3 - seq[i][len - 1]]));
	
		assert(wbgscores[i].size() == (unsigned int) len);
		assert(cbgscores[i].size() == (unsigned int) len);
	}
}

void BGModel::calc_bg_scores_1() {
	const vector<vector <char> >& seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = seq[i].size();
		wbgscores[i].clear();
		cbgscores[i].clear();
		wbgscores[i].reserve(len);
		cbgscores[i].reserve(len);
		
		// Use 0th order model for first Watson base
		wbgscores[i].push_back(log(model0[seq[i][0]]));
				
		// Use 1st-order model for most bases
		for(int j = 1; j < len; j++) {
			wbgscores[i].push_back(log(model1[seq[i][j - 1] * 4
																+ seq[i][j]]));
		}
		for(int j = 0; j < len - 1; j++) {
			cbgscores[i].push_back(log(model1[(3 - seq[i][j + 1]) * 4
																+ (3 - seq[i][j])]));
		}
		
		// Use 0th order model for last Crick base
		cbgscores[i].push_back(log(model0[3 - seq[i][len - 1]]));
	
		assert(wbgscores[i].size() == (unsigned int) len);
		assert(cbgscores[i].size() == (unsigned int) len);
	}
}

void BGModel::calc_bg_scores_0() {
	const vector<vector <char> >& seq = seqset.seq();
	int ss_num_seqs = seqset.num_seqs();
	
	int len;
	for(int i = 0; i < ss_num_seqs; i++) {
		len = seq[i].size();
		wbgscores[i].clear();
		cbgscores[i].clear();
		wbgscores[i].reserve(len);
		cbgscores[i].reserve(len);
		
		// Use 0th order model for all bases
		for(int j = 0; j < len; j++) {
			wbgscores[i].push_back(log(model0[seq[i][j]]));
		}
		for(int j = 0; j < len; j++) {
			cbgscores[i].push_back(log(model0[(3 - seq[i][j])]));
		}
		
		assert(wbgscores[i].size() == (unsigned int) len);
		assert(cbgscores[i].size() == (unsigned int) len);
	}
}
