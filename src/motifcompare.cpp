#include "motifcompare.h"

MotifCompare::MotifCompare() {
}
	
float MotifCompare::compare(const Motif& m1, const Motif& m2) const {
	int cols = 6;
	
	int fmsize1 = m1.get_width() + 2 * m1.ncols();
	vector<float> fm1(fmsize1 * 4);
	m1.freq_matrix_extended(fm1);
	
	int fmsize2 = m2.get_width() + 2 * m2.ncols();
	vector<float> fm2(fmsize2 * 4);
	m2.freq_matrix_extended(fm2);
	
	// Find columns with the highest information content
	vector<struct idscore> csc(fmsize1);
	float ent;
	for(int i = 0; i < fmsize1; i++) {
		ent = 0.0;
		for(int j = i * 4; j < (i + 1) * 4; j++) {
			ent += fm1[j] > 0? -fm1[j] * log(fm1[j]) : 0;
		}
		csc[i].id = i;
		csc[i].score = ent;
	}
	sort(csc.begin(), csc.end(), isc);
	vector<int> chosencols1;
	for(int i = 0; i < cols; i++)
		chosencols1.push_back(csc[i].id);
	sort(chosencols1.begin(), chosencols1.end());
	
	vector<float> chosenfm1;
	copy_subfreq(fm1, chosencols1, chosenfm1);
	
	vector<int>::iterator coliter = chosencols1.begin();
	int offset = chosencols1[0];
	for(; coliter != chosencols1.end(); ++coliter)
		*coliter -= offset;
	
	double bestc = -1.1;
	double c = 0.0;
	for(int i = 0; i < fmsize2; i++) {
		vector<int> chosencols2;
		coliter = chosencols1.begin();
		for(; coliter != chosencols1.end(); ++coliter)
			chosencols2.push_back(i + *coliter);
		vector<float> chosenfm2;
		copy_subfreq(fm2, chosencols2, chosenfm2);
		c = corr(chosenfm1, chosenfm2);
		if(c > bestc)
			bestc = c;
	}
	
	return bestc;
}

void MotifCompare::copy_subfreq(const vector<float>& fm, const vector<int>& cols, vector<float>& subfm) const {
	vector<int>::const_iterator col_iter = cols.begin();
	for(; col_iter != cols.end(); ++col_iter)
		for(int j = 0; j < 4; j++)
			subfm.push_back(fm[*col_iter * 4 + j]);
}

