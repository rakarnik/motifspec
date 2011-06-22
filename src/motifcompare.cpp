#include "motifcompare.h"

MotifCompare::MotifCompare(const Seqset& s, const BGModel& bgm) 
: seqset(s), bgmodel(bgm)
{
}

float MotifCompare::compare(const Motif& m1, const Motif& m2) const {
	vector<float> scores1;
	vector<float> scores2;
	int c, p, len, window_size;
	int olap = 2; // overlap required between motifs
	bool s;
	float sc, max_sc_1, max_sc_2;
	
	// Set up score matrices
	double* sm1 = new double[4 * m1.ncols()];
	m1.calc_score_matrix(sm1);
	double* sm2 = new double[4 * m2.ncols()];
	m2.calc_score_matrix(sm2);
	
	int w1 = m1.get_width();
	int w2 = m2.get_width();
	
	// Score sites for this motif with both PWMs and background model
	vector<Site>::const_iterator siter1 = m1.sites().begin();
	for(; siter1 != m1.sites().end(); ++siter1) {
		c = siter1->chrom();
		p = siter1->posit();
		s = siter1->strand();
		len = seqset.len_seq(c);
		
		max_sc_1 = -1000;
		window_size = w2 - olap;
		for(int i = max(0, p - window_size); i < min(len - w1 - 1, p + window_size); i++) {
			sc = - m1.score_site(sm1, c, i, s) + bgmodel.score_site(m1.first_column(), m1.last_column(), w1, c, i, s);
			if(sc > max_sc_1)
				max_sc_1 = sc;
		}
		scores1.push_back(max_sc_1);
		
		max_sc_2 = -1000;
		window_size = w1 - olap;
		for(int i = max(0, p - window_size); i < min(len - w2 - 1, p + window_size); i++) {
			sc = - m2.score_site(sm2, c, i, s) + bgmodel.score_site(m2.first_column(), m2.last_column(), w2, c, i, s);
			if(sc > max_sc_2)
				max_sc_2 = sc;
		}
		scores2.push_back(max_sc_2);
	}
	
	// Score sites for the other motif with both PWMs
	vector<Site>::const_iterator siter2 = m2.sites().begin();
	for(; siter2 != m2.sites().end(); ++siter2) {
		c = siter2->chrom();
		p = siter2->posit();
		s = siter2->strand();
		len = seqset.len_seq(c);
		
		max_sc_1 = -1000;
		window_size = w2 - olap;
		for(int i = max(0, p - window_size); i < min(len - w1 - 1, p + window_size); i++) {
			sc = - m1.score_site(sm1, c, i, s) + bgmodel.score_site(m1.first_column(), m1.last_column(), w1, c, i, s);
			if(sc > max_sc_1)
				max_sc_1 = sc;
		}
		scores1.push_back(max_sc_1);
		
		max_sc_2 = -1000;
		window_size = w1 - olap;
		for(int i = max(0, p - window_size); i < min(len - w2 - 1, p + window_size); i++) {
			sc = - m2.score_site(sm2, c, i, s) + bgmodel.score_site(m2.first_column(), m2.last_column(), w2, c, i, s);
			if(sc > max_sc_2)
				max_sc_2 = sc;
		}
		scores2.push_back(max_sc_2);
	}
	delete [] sm1;
	delete [] sm2;
	
	assert(scores1.size() == (unsigned int) m1.number() + m2.number());
	assert(scores2.size() == (unsigned int) m1.number() + m2.number());
	
	float mcorr = corr(scores1, scores2);
	int df = scores1.size() - 2;
	float t = 0.0;
	if(mcorr > 0.999) {
		t = 1000;
	} else {
		t = mcorr * sqrt((df - 2)/(1 - mcorr * mcorr));
	}
	return t;
}

