#ifndef _bgmodel
#define _bgmodel

#include "seqset.h"

class BGModel {
	const Seqset& seqset;
	int total_seq_len;
	int order;
	float gc_genome;
	vector<float> gc;
	vector<float> model0;
	vector<float> model1;
	vector<float> model2;
	vector<float> model3;
	vector<float> model4;
	vector<float> model5;
	vector<vector <float> > wbgscores;
	vector<vector <float> > cbgscores;
	vector<vector <float> > cumulscores;
	void (BGModel::*train_background[6])();
	void (BGModel::*calc_bg_scores[6])();

	void train_background_5();                                     // Train 5th order background model
	void train_background_4();                                     // Train 4th order background model
	void train_background_3();                                     // Train 3rd order background model
	void train_background_2();                                     // Train 2nd order background model
	void train_background_1();                                     // Train 1st order background model
	void train_background_0();                                     // Train 0th order background model
	void calc_bg_scores_5();                                       // Calculate scores using 5th order background model
	void calc_bg_scores_4();                                       // Calculate scores using 4th order background model
	void calc_bg_scores_3();                                       // Calculate scores using 3rd order background model
	void calc_bg_scores_2();                                       // Calculate scores using 2nd order background model
	void calc_bg_scores_1();                                       // Calculate scores using 1st order background model
	void calc_bg_scores_0();                                       // Calculate scores using 0th order background model

public:
	BGModel(const Seqset& s, const int ord = 3);
	float tot_seq_len() const { return total_seq_len; }           // Return total length of all sequences
	float gcgenome() const { return gc_genome; }                  // Return overall GC content
	float gccontent(const int i) const { return gc[i]; }          // Return GC content of a specified sequence
	double score_site(vector<int>::const_iterator first_col, vector<int>::const_iterator last_col, const int width, const int c, const int p, const bool s) const;
	vector<vector <float> > const& get_cumulscores() const { return cumulscores; }
};


#endif
