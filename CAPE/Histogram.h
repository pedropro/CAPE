/*
 * Copyright 2018 Pedro Proenza <p.proenca@surrey.ac.uk> (University of Surrey)
 *
 */

#pragma once
#include <Eigen/Dense>
#include <vector>
#include <iostream>

using namespace std;
#define DMAX std::numeric_limits<float>::max()
#define DMIN std::numeric_limits<float>::min()

class Histogram
{
private:
	vector<int> H;
	vector<int> B;
	int nr_bins_per_coord;
	int nr_bins;
	int nr_points;
public:
	Histogram(int nr_bins_per_coord);
	void initHistogram(Eigen::MatrixXd & Points, vector<bool> & Flags);
	vector<int> getPointsFromMostFrequentBin();
	void removePoint(int point_id);
	~Histogram(void);
};

