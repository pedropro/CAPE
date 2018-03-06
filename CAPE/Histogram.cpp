/*
 * Copyright 2018 Pedro Proenza <p.proenca@surrey.ac.uk> (University of Surrey)
 *
 */

#include "Histogram.h"

Histogram::Histogram(int nr_bins_per_coord)
{
	this->nr_bins_per_coord = nr_bins_per_coord;
	this->nr_bins= nr_bins_per_coord*nr_bins_per_coord;
	this->H.assign(nr_bins,0);
}

void Histogram::initHistogram(Eigen::MatrixXd & P, vector<bool> & Flags){

	int nr_points = P.rows();
	this->nr_points = nr_points;
	this->B.assign(nr_points,-1); 

	// Set limits
	// Polar angle [0 pi]
	double min_X(0), max_X(3.14);
	// Azimuth angle [-pi pi]
	double min_Y(-3.14), max_Y(3.14);

	// Fill structures
	int X_q, Y_q;
	for (int i=0; i<nr_points; i++){
		if (Flags[i]){
			X_q = (nr_bins_per_coord-1)*(P(i,0)-min_X)/(max_X-min_X);
			// Dealing with degeneracy
			if (X_q>0){
				Y_q = (nr_bins_per_coord-1)*(P(i,1)-min_Y)/(max_Y-min_Y);
			}else{
				Y_q = 0;
			}
			int bin  = Y_q*nr_bins_per_coord + X_q;
			B[i] = bin;
			H[bin]++;
		}
	}
}

vector<int> Histogram::getPointsFromMostFrequentBin(){

	vector<int> point_ids;

	int most_frequent_bin = -1;
	int max_nr_occurrences = 0;
	for (int i=0; i<nr_bins; i++){
		if (H[i]>max_nr_occurrences){
			most_frequent_bin = i;
			max_nr_occurrences = H[i];
		}
	}
	if(max_nr_occurrences>0){
		for (int i=0; i<nr_points; i++){
			if (B[i]==most_frequent_bin){
				point_ids.push_back(i);
			}
		}
	}
	return point_ids;
}

void Histogram::removePoint(int point_id){
	H[B[point_id]]--;
	B[point_id] = -1;
}

Histogram::~Histogram(void)
{
}
