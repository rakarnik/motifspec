#include "cluster.h"

Cluster::~Cluster() {
	delete [] membership;
	delete [] mean;
}

Cluster::Cluster(const Cluster& c) {
	init(c.expr, c.npoints, c.nameset);
	nmembers = c.nmembers;
	for(int i = 0; i < nameset.size(); i++) {
		membership[i] = c.membership[i];
	}
	dirty = true;
}

const Cluster& Cluster::operator= (const Cluster& c) {
	delete [] membership;
	delete [] mean;
	init(c.expr, c.npoints, c.nameset);
	nmembers = c.nmembers;
	for(int i = 0; i < nameset.size(); i++) {
		membership[i] = c.membership[i];
	}
	dirty = true;
	return *this;
}

void Cluster::init (float** exprtab, int numexpr, const vector<string>& names) {
	expr = exprtab;
	nmembers = 0;
	ngenes = names.size();
	npoints = numexpr;
	membership = new int[ngenes];
	for(int g = 0; g < ngenes; g++) {
		membership[g] = 0;
	}
	mean = new float[npoints];
	for(int p = 0; p < npoints; p++) {
		mean[p] = 0;
	}
	nameset = names;
	dirty = false;
}

void Cluster::add_gene(const int gene) {
	if (membership[gene] == 0) {
		membership[gene] = 1;
		nmembers++;
	}
	dirty = true;
}

void Cluster::remove_gene(const int gene) {
	if (membership[gene] == 1) {
		membership[gene] = 0;
		nmembers--;
	}
	dirty = true;
}

bool Cluster::is_member(const int gene) {
	return (membership[gene] == 1);
}

float* Cluster::get_mean() {
	if (dirty) calc_mean();
	return mean;
}

int Cluster::size() const {
	return nmembers;
}

void Cluster::calc_mean() {
	if (dirty) {
		int i, j;
		if (npoints > 1) {
			for (i = 0; i < npoints; i++) {
				mean[i] = 0;
				for (j = 0; j < nameset.size(); j++) {
					mean[i] += membership[j] * expr[j][i];
				}
				mean[i] /= nmembers;
			}
		}
		dirty = false;
	}
}

float Cluster::corr (const float* pattern) const {
	int i;
	float u1, u2;
	float s1, s2;
	float c;
	u1 = u2 = s1 = s2 = c = 0.0;
	for (i = 0; i < npoints; i++) {
		u1 += mean[i];
		u2 += pattern[i];
		s1 += mean[i] * mean[i];
		s2 += pattern[i] * pattern[i];
	}
	u1 /= npoints;
	u2 /= npoints;
	for (i = 0; i < npoints; i++) {
		c += (mean[i] - u1) * (pattern[i] - u2);
	}
	s1 = s1/npoints - u1 * u1;
	s2 = s2/npoints - u2 * u2;
	if (s1 > 0 && s2 > 0) {
		c /= sqrt(s1 * s2);
		c /= npoints;
	} else {
		c = 0;
	}
	return c;
}

void Cluster::genes(vector<int>& genes) const {
	for (int g = 0; g < nameset.size(); g++) {
		if (membership[g] == 1) {
			genes.push_back(g);
		}
	}
}