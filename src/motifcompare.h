#ifndef _motifcompare
#define _motifcompare

#include "motif.h"

class MotifCompare {
	const Seqset& seqset;
	const BGModel& bgmodel;
	
public:
	MotifCompare(const Seqset& s, const BGModel& bgm);
	float compare(const Motif& m1, const Motif& m2) const;
};

#endif

