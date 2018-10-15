/*
 * Copyright 2018 Pedro Proenza <p.proenca@surrey.ac.uk> (University of Surrey)
 *
 */
#pragma once
#include <math.h>
#include <iostream>
#include <iterator>
#include <algorithm>
#include "PlaneSeg.h"
#include "CylinderSeg.h"
#include "Histogram.h"
#include <opencv2/opencv.hpp>
#include <Eigen/Dense>

typedef Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic>  MatrixXb;

class CAPE
{
private:
	int cell_width;
	int cell_height;
	int depth_height;
	int depth_width;
	float max_merge_dist;
	float min_cos_angle_4_merge;
	bool cylinder_detection;
	std::vector<cv::Vec3b> color_code;
	std::vector<PlaneSeg*> Grid;
	cv::Mat_<int> grid_plane_seg_map;
	cv::Mat_<uchar> grid_plane_seg_map_eroded;
	cv::Mat_<int> grid_cylinder_seg_map;
	cv::Mat_<uchar> grid_cylinder_seg_map_eroded;
	cv::Mat mask;
	cv::Mat mask_eroded;
	cv::Mat mask_square_eroded;
	cv::Mat mask_dilated;
	cv::Mat mask_diff;
	cv::Mat mask_square_kernel;
	cv::Mat mask_cross_kernel;
	float * distances_stacked;
	Eigen::ArrayXf distances_cell_stacked;
	unsigned char * seg_map_stacked;
	bool* activation_map;
	bool* unassigned_mask;
public:
	CAPE(int depth_height, int depth_width, int cell_width, int cell_height, bool cylinder_detection, float min_cos_angle_4_merge = 0.97814, float max_merge_dist = 900);
    void process(Eigen::MatrixXf & cloud_array, int & nr_planes, int & nr_cylinders, cv::Mat & seg_output, vector<PlaneSeg> & plane_segments_final, vector<CylinderSeg> & cylinder_segments_final);
	void RegionGrowing(unsigned short width, unsigned short height, bool* input, bool* output, vector<PlaneSeg*> & Grid, vector<float> & cell_dist_tols, unsigned short x, unsigned short y, double * normal, double d);
	void getConnectedComponents(cv::Mat & segment_map, MatrixXb & planes_association_matrix);
	~CAPE(void);
};

