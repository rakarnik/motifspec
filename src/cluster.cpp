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

void Cluster::add_genes(const int* genes, const int count) {
	for(int g = 0; g < count; g++) {
		add_gene(genes[g]);
	}
}

void Cluster::remove_gene(const int gene) {
	if (membership[gene] == 1) {
		membership[gene] = 0;
		nmembers--;
	}
	dirty = true;
}

void Cluster::remove_all_genes() {
	for(int g = 0; g < ngenes; g++) {
		remove_gene(g);
	}
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
		for (i = 0; i < npoints; i++) {
			mean[i] = 0;
			for (j = 0; j < ngenes; j++) {
				mean[i] += membership[j] * expr[j][i];
			}
			mean[i] /= nmembers;
		}
		dirty = false;
	}
}

float Cluster::corrmean(const float* pattern) const {
	return corr(mean, pattern, npoints);
}

void Cluster::genes(int* genes) const {
	int count = 0;
	for (int g = 0; g < ngenes; g++) {
		if (membership[g] == 1) {
			genes[count] = g;
			count++;
		}
	}
	assert(count == nmembers);
}

