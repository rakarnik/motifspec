#ifndef _motifcompare
#define _motifcompare

#include "motif.h"

class MotifCompare {
	const Seqset& seqset;
	
public:
	MotifCompare(const Seqset& s);
	bool compare(const Motif& m1, const Motif& m2, const BGModel& bgm) const;
};

#endif

