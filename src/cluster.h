/*
 *  cluster.h
 *  Represents a k-means cluster
 */

#ifndef _cluster
#define _cluster
#include "standard.h"


class Cluster {
	vector<string> nameset;
	vector<vector<float> > expr;
	int ngenes;
	vector<int> membership;
	vector<float> mean;
	bool dirty;
	
public:
	Cluster (const vector<vector<float> >& exprtab, const vector<string>& names);
	void add_gene (const int gene);                         // Add gene to this cluster
	void remove_gene (const int gene);                      // Remove gene from this cluster
	int size() const;                                       // Return size of this cluster
	void calc_mean();                                       // Calculate the mean for this cluster
	void set_mean (const vector<float>& new_mean);          // Set the new mean for this cluster
	vector<float>& get_mean();                              // Get the mean for this cluster
	float corr (const vector<float>& pattern) const;        // Calculate the correlation between the cluster mean and 'pattern'
	void genes(vector<int>& genes) const;                   // Return the genes that are assigned to this cluster             
};

#endif