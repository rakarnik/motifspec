#include "cluster.h"

Cluster::Cluster (const vector<vector<float> >& exprtab, const vector<string>& names):
expr(exprtab),
ngenes(0),
membership(exprtab.size(), 0),
mean(exprtab[0].size()),
nameset(names),
dirty(false) {
}

void Cluster::add_gene(const int gene) {
	if (membership[gene] == 0) {
		membership[gene] = 1;
		ngenes++;
	}
	dirty = true;
}

void Cluster::remove_gene(const int gene) {
	if (membership[gene] == 1) {
		membership[gene] = 0;
		ngenes--;
	}
	dirty = true;
}

void Cluster::set_mean(const vector<float>& new_mean) {
	mean.assign(new_mean.begin(), new_mean.end());
}

vector<float>& Cluster::get_mean() {
	return mean;
}

int Cluster::size() const {
	return ngenes;
}

void Cluster::calc_mean() {
	if (dirty) {
		int i, j;
		if (expr.size() > 1) {
			for (i = 0; i < expr[0].size(); i++) {
				mean[i] = 0;
				for (j = 0; j < expr.size(); j++) {
					mean[i] += membership[j] * expr[j][i];
				}
				mean[i] /= ngenes;
			}
		}
		dirty = false;
	}
}

float Cluster::corr (const vector<float>& pattern) const {
	int i;
	float u1, u2;
	float s1, s2;
	float c;
	for (i = 0; i < mean.size(); i++) {
		u1 += mean[i];
		u2 += pattern[i];
		s1 += mean[i] * mean[i];
		s2 += pattern[i] * pattern[i];
	}
	u1 /= mean.size();
	u2 /= pattern.size();
	for (i = 0; i < mean.size(); i++) {
		c += (mean[i] - u1) * (pattern[i] - u2);
	}
	s1 = s1/mean.size() - u1 * u1;
	s2 = s2/pattern.size() - u2 * u2;
	if (s1 > 0 && s2 > 0) {
		c /= sqrt(s1 * s2);
		c /= mean.size();
	} else {
		c = 0;
	}
	return c;
}

void Cluster::genes(vector<int>& genes) const {
	for (int g = 0; g < membership.size(); g++) {
		if (membership[g] == 1) {
			genes.push_back(g);
		}
	}
}