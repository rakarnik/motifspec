#ifndef _motifcompare
#define _motifcompare

#include "motif.h"

class MotifCompare {
	const Seqset& seqset;
	const BGModel& bgmodel;
	
	struct idscore {
		int id;
		float score;
	};
	
	struct iscomp {
		bool operator() (struct idscore is1, struct idscore is2) { return (is1.score < is2.score); }
	} isc;
	
	void copy_subfreq(const vector<float>& fm, const vector<int>& cols, vector<float>& subfm) const;
	
public:
	MotifCompare(const Seqset& s, const BGModel& bgm);
	float compare(const Motif& m1, const Motif& m2) const;
};

#endif

