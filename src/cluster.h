/*
 *  cluster.h
 *  Represents a k-means cluster
 */

#ifndef _cluster
#define _cluster
#include "standard.h"


class Cluster {
	float** expr;
	vector<string> nameset;
	int nmembers;
	int ngenes;
	int npoints;
	int* membership;
	float* mean;
	bool dirty;
	
public:
	Cluster() {};
	Cluster(const Cluster& c);
	~Cluster();
	const Cluster& operator= (const Cluster& c);
	void init (float** exprtab, int numexpr, const vector<string>& names);
	void add_gene(const int gene);                          // Add gene to this cluster
	void add_genes(const int* genes, const int count);      // Add several genes to this cluster
	void remove_gene(const int gene);                       // Remove gene from this cluster
	void remove_all_genes();                                // Remove all genes from this cluster
	bool is_member(const int gene);                         // Return whether the gene is a member of this cluster
	int size() const;                                       // Return size of this cluster
	void calc_mean();                                       // Calculate the mean for this cluster
	float* get_mean();                                      // Get the mean for this cluster
	float corr (const float* pattern) const;                // Calculate the correlation between the cluster mean and 'pattern'
	void genes(int* genes) const;                           // Return the genes that are assigned to this cluster
};

#endif

