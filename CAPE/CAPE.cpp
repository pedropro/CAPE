/*
 * Copyright 2018 Pedro Proenza <p.proenca@surrey.ac.uk> (University of Surrey)
 *
 */

#include "CAPE.h"

CAPE::CAPE(int depth_height, int depth_width, int cell_width, int cell_height, bool cylinder_detection, float min_cos_angle_4_merge, float max_merge_dist)
{
	this->depth_height = depth_height;
	this->depth_width = depth_width;
	this->cell_width = cell_width;
	this->cell_height = cell_height;
	this->max_merge_dist = max_merge_dist;
	this->min_cos_angle_4_merge = min_cos_angle_4_merge;
	this->cylinder_detection = cylinder_detection;

	int nr_horizontal_cells = depth_width/cell_width;
	int nr_vertical_cells = depth_height/cell_height;
	grid_plane_seg_map = cv::Mat_<int>(nr_vertical_cells,nr_horizontal_cells,0);
	grid_plane_seg_map_eroded = cv::Mat_<uchar>(nr_vertical_cells,nr_horizontal_cells,uchar(0));
	grid_cylinder_seg_map = cv::Mat_<int>(nr_vertical_cells,nr_horizontal_cells,0);
	grid_cylinder_seg_map_eroded = cv::Mat_<uchar>(nr_vertical_cells,nr_horizontal_cells,uchar(0));
	mask = cv::Mat(nr_vertical_cells,nr_horizontal_cells,CV_8U);
	mask_eroded = cv::Mat(nr_vertical_cells,nr_horizontal_cells,CV_8U);
	mask_dilated = cv::Mat(nr_vertical_cells,nr_horizontal_cells,CV_8U);
	mask_diff = cv::Mat(nr_vertical_cells,nr_horizontal_cells,CV_8U);
	
	int nr_total_cells = nr_vertical_cells*nr_horizontal_cells;
	int nr_pts_per_cell = cell_width*cell_height;
	
	for(int i=0;i<nr_total_cells;i++){
		Grid.push_back(NULL);
	}

	distances_stacked = (float*)malloc(depth_height*depth_width*sizeof(float));
	distances_cell_stacked = Eigen::ArrayXf::Zero(nr_pts_per_cell,1);
	seg_map_stacked = (unsigned char*)malloc(depth_height*depth_width*sizeof(unsigned char));
	this->activation_map  = (bool*)malloc(nr_total_cells*sizeof(bool));
	this->unassigned_mask  = (bool*)malloc(nr_total_cells*sizeof(bool));
	mask_square_kernel = cv::Mat::ones(3,3,CV_8U);
	mask_cross_kernel = cv::Mat::ones(3,3,CV_8U);
	mask_cross_kernel.at<uchar>(0,0) = 0; mask_cross_kernel.at<uchar>(2,2) = 0;
	mask_cross_kernel.at<uchar>(0,2) = 0; mask_cross_kernel.at<uchar>(2,0) = 0;
}

void CAPE::process(Eigen::MatrixXf & cloud_array, int & nr_planes_final, int & nr_cylinders_final,  cv::Mat & seg_out, vector<PlaneSeg> & plane_segments_final, vector<CylinderSeg> & cylinder_segments_final){

	int nr_horizontal_cells = depth_width/cell_width;
	int nr_vertical_cells = depth_height/cell_height;
	int nr_total_cells = nr_vertical_cells*nr_horizontal_cells;
	int nr_pts_per_cell = cell_width*cell_height;
	int cylinder_code_offset = 50;

	grid_plane_seg_map = 0;
	grid_plane_seg_map_eroded = 0;
	grid_cylinder_seg_map = 0;
	grid_cylinder_seg_map_eroded = 0;
	std::memset(seg_map_stacked, (uchar)0, depth_height*depth_width*sizeof(unsigned char));
	std::memset(distances_stacked, 100, depth_height*depth_width*sizeof(float)); /* = to really high float*/

	/*--------------------------------- Planar cell fitting  ---------------------------------*/
	std::vector<float> cell_distance_tols(nr_total_cells,0);
	int stacked_cell_id = 0;
	float cell_diameter;
	float sin_cos_angle_4_merge = sqrt(1-pow(min_cos_angle_4_merge,2));
	for (int cell_r=0; cell_r<nr_vertical_cells; cell_r++){
		for (int cell_c=0; cell_c<nr_horizontal_cells; cell_c++){
			Grid[stacked_cell_id] = new PlaneSeg(cloud_array, stacked_cell_id, nr_pts_per_cell,cell_width);
			if (Grid[stacked_cell_id]->planar){
				cell_diameter = (cloud_array.block(stacked_cell_id*nr_pts_per_cell+nr_pts_per_cell-1,0,1,3)-cloud_array.block(stacked_cell_id*nr_pts_per_cell,0,1,3)).norm();
				// Add truncated distance
				cell_distance_tols[stacked_cell_id]=pow(min(max(cell_diameter*sin_cos_angle_4_merge,20.0f),max_merge_dist),2);
			}

			stacked_cell_id++;
		}
	}
	/*------------------------------- Initialize histogram -----------------------------------*/
	//double t3 = cv::getTickCount();
	// Spherical coordinates
	Eigen::MatrixXd C(nr_total_cells,2);
	vector<bool> planar_flags(nr_total_cells,false);
	vector<float> scores_stacked(nr_total_cells,0.0);
	int nr_remaining_planar_cells = 0;
	double nx,ny,nz;
	for (int cell_id=0; cell_id<nr_total_cells; cell_id++){
		if (Grid[cell_id]->planar){
			nx = Grid[cell_id]->normal[0];
			ny = Grid[cell_id]->normal[1];
			nz = Grid[cell_id]->normal[2];
			double n_proj_norm =  sqrt(nx*nx+ny*ny);
			C(cell_id,0) = acos(-nz);
			C(cell_id,1) = atan2(nx/n_proj_norm,ny/n_proj_norm); 
			planar_flags[cell_id] = true;
			scores_stacked[cell_id] = Grid[cell_id]->score;
			nr_remaining_planar_cells++;
		}
	}
	Histogram H(20);
	H.initHistogram(C,planar_flags);

	// Initialization for cell-wise region growing and model fitting
	vector<PlaneSeg> plane_segments;
	vector<CylinderSeg> cylinder_segments;
	bool stop = false;
	int nr_cylinders = 0;
	vector<pair<int,int> > cylinder2region_map;

	for (int cell_id=0; cell_id<nr_total_cells; cell_id++){
		unassigned_mask[cell_id] = planar_flags[cell_id];
	}
	/*------------ Cell-wise region growing with embedded model fitting ---------------------*/
	while(nr_remaining_planar_cells>0){  
		// 1. Seeding
		// Pick up seed candidates
		vector<int> seed_candidates = H.getPointsFromMostFrequentBin();

		// Checkpoint 1
		if (seed_candidates.size()<5)
			break;

		// Select seed based on score
		int seed_id;
		float min_MSE = INT_MAX;
		for(int i=0; i<seed_candidates.size();i++){
			int seed_candidate = seed_candidates[i];
			if(Grid[seed_candidate]->MSE<min_MSE){
				seed_id =  seed_candidate;
				min_MSE = Grid[i]->MSE;
			}
		}

		PlaneSeg new_ps = *Grid[seed_id];

		// 2. Actual seed cell growing
		int y = seed_id/nr_horizontal_cells;
		int x = seed_id%nr_horizontal_cells;
		double * seed_n = new_ps.normal;
		double seed_d = new_ps.d;
		memset(activation_map, false, sizeof(bool)*nr_total_cells);
		RegionGrowing(nr_horizontal_cells, nr_vertical_cells, unassigned_mask, activation_map, Grid, cell_distance_tols, x, y, seed_n, seed_d);

		// 3. Merge activated cells & remove them from histogram and list of remaining cells
		int nr_cells_activated = 0;
		for(int i=0; i<nr_total_cells;i++){
			if (activation_map[i]){
				new_ps.expandSegment(Grid[i]);
				nr_cells_activated++;
				H.removePoint(i);
				unassigned_mask[i] = false;
				nr_remaining_planar_cells--;
			}
		}

		// Checkpoint 2
		if (nr_cells_activated<4)
			continue;

		new_ps.fitPlane();

		// 4. Model fitting
		if(new_ps.score>100){
			// It is a plane
			plane_segments.push_back(new_ps);
			int nr_curr_planes = plane_segments.size();
			// Mark cells
			int i=0;
			int *row;
			for(int r=0; r<nr_vertical_cells; r++){
				row = grid_plane_seg_map.ptr<int>(r);
				for(int c=0; c<nr_horizontal_cells; c++){
					if (activation_map[i]){
						row[c] = nr_curr_planes;
					}
					i++;
				}
			}
		}else{
			if(cylinder_detection && nr_cells_activated>5){
                // It is an extrusion
				CylinderSeg cy(Grid, activation_map, nr_cells_activated);
				cylinder_segments.push_back(cy);
				// Fit planes to subsegments
				for(int seg_id=0; seg_id<cy.nr_segments; seg_id++){
					new_ps.clearPoints();
					for(int c=0;c<nr_cells_activated;c++){
						if (cy.inliers[seg_id](c)){
							new_ps.expandSegment(Grid[cy.local2global_map[c]]);
						}
					}
					new_ps.fitPlane();
					// Model selection based on MSE
					if(new_ps.MSE<cy.MSEs[seg_id]){
						plane_segments.push_back(new_ps);
						int nr_curr_planes = plane_segments.size();
						for(int c=0; c<nr_cells_activated; c++){
							if (cy.inliers[seg_id](c)){
								int cell_id = cy.local2global_map[c];
								grid_plane_seg_map.at<int>(cell_id/nr_horizontal_cells,cell_id%nr_horizontal_cells) = nr_curr_planes;
							}
						}
						cy.cylindrical_mask[seg_id] = false;
					}else{
						nr_cylinders++;
						cylinder2region_map.push_back(make_pair(cylinder_segments.size()-1,seg_id));
						for(int c=0; c<nr_cells_activated; c++){
							if (cy.inliers[seg_id](c)){
								int cell_id = cy.local2global_map[c];
								grid_cylinder_seg_map.at<int>(cell_id/nr_horizontal_cells,cell_id%nr_horizontal_cells) = nr_cylinders;
							}
						}
						cy.cylindrical_mask[seg_id] = true;
					}
				}
			}
		}
	}

	/*------------------------------------ Plane merging ------------------------------------*/
	int nr_planes =  plane_segments.size();
	MatrixXb planes_association_matrix= MatrixXb::Zero(nr_planes,nr_planes);
	getConnectedComponents(grid_plane_seg_map, planes_association_matrix);

	vector<int> plane_merge_labels;
	for(int i=0; i<nr_planes; i++)
		plane_merge_labels.push_back(i);

	// Connect compatible planes
	for(int r=0; r<planes_association_matrix.rows(); r++){
		int plane_id = plane_merge_labels[r];
		bool plane_expanded = false;
		for(int c=r+1; c<planes_association_matrix.cols(); c++){
			if(planes_association_matrix(r,c)){
				double cos_angle = plane_segments[plane_id].normal[0]*plane_segments[c].normal[0]
				+ plane_segments[plane_id].normal[1]*plane_segments[c].normal[1]
				+ plane_segments[plane_id].normal[2]*plane_segments[c].normal[2];
				double distance = pow(plane_segments[r].normal[0]*plane_segments[c].mean[0]
				+ plane_segments[plane_id].normal[1]*plane_segments[c].mean[1]  
				+ plane_segments[plane_id].normal[2]*plane_segments[c].mean[2] + plane_segments[plane_id].d,2);
				if (cos_angle > min_cos_angle_4_merge && distance<max_merge_dist){
					plane_segments[plane_id].expandSegment(&plane_segments[c]);
					plane_merge_labels[c] = plane_id;
					plane_expanded = true;
				}else{
					planes_association_matrix(r,c) = false;
				}
			}
		}
		if(plane_expanded)
			plane_segments[plane_id].fitPlane();
	}

	/*------------------------------- Refine plane boundaries -------------------------------*/
    //vector<PlaneSeg> plane_segments_joint;
	for(int i=0; i<nr_planes;i++){

		if(i!=plane_merge_labels[i])
			continue;

		// Build mask w/ merged segments
		mask = cv::Scalar(0);
		for(int j=i; j<nr_planes;j++){
			if(plane_merge_labels[j]==plane_merge_labels[i]){
				mask.setTo(1, grid_plane_seg_map==j+1);
			}
		}

		// Erode with cross to obtain borders and check support
		cv::erode(mask,mask_eroded,mask_cross_kernel);
		double min, max;
		cv::minMaxLoc(mask_eroded, &min, &max);

		// If completely eroded ignore plane
		if (max==0){
			continue;
		}

        plane_segments_final.push_back(plane_segments[i]);

		// Dilate to obtain borders
		cv::dilate(mask,mask_dilated,mask_square_kernel);
		mask_diff = mask_dilated-mask_eroded;

		int stacked_cell_id = 0;
        uchar plane_nr = (unsigned char)plane_segments_final.size();
		float nx = (float)plane_segments[i].normal[0]; float ny = (float)plane_segments[i].normal[1];
		float nz = (float)plane_segments[i].normal[2]; float d = (float)plane_segments[i].d;
		unsigned char * row_ptr;

		grid_plane_seg_map_eroded.setTo(plane_nr, mask_eroded > 0);

		// Cell refinement
		for (int cell_r=0; cell_r<nr_vertical_cells; cell_r++){
			row_ptr = mask_diff.ptr<uchar>(cell_r);
			for (int cell_c=0; cell_c<nr_horizontal_cells; cell_c++){
				int offset = stacked_cell_id*nr_pts_per_cell;
				int next_offset = offset+nr_pts_per_cell;
				if(row_ptr[cell_c]>0){
					float max_dist = 9*plane_segments[i].MSE;
					// Compute distance block
					distances_cell_stacked = cloud_array.block(offset,0,nr_pts_per_cell,1).array()*nx
						+ cloud_array.block(offset,1,nr_pts_per_cell,1).array()*ny
						+ cloud_array.block(offset,2,nr_pts_per_cell,1).array()*nz
						+ d;
					// Assign pixels
					int j=0;
					float dist;
					for(int pt = offset; pt<next_offset;j++,pt++){
						dist = pow(distances_cell_stacked(j),2);
						if(dist<max_dist && dist<distances_stacked[pt]){ 
							distances_stacked[pt] = dist;
							seg_map_stacked[pt] = plane_nr;
						}
					}
				}
				stacked_cell_id++;
			}
		}
	}
    nr_planes_final = plane_segments_final.size();

	/*------------------------------ Refine cylinder boundaries -----------------------------*/
	nr_cylinders_final = 0;
	if (cylinder_detection){
		Eigen::Vector3f point;
		float dist;
		for(int i=0; i<nr_cylinders;i++){

			int cylinder_nr = i+1;
			int reg_id = cylinder2region_map[i].first;
			int sub_reg_id = cylinder2region_map[i].second;

			// Build mask
			mask = cv::Scalar(0);
			mask.setTo(1, grid_cylinder_seg_map==cylinder_nr);

			// Erode to obtain borders
			cv::erode(mask,mask_eroded,mask_cross_kernel);
			double min, max;
			cv::minMaxLoc(mask_eroded, &min, &max);

			// If completely eroded ignore cylinder
			if (max==0){
				continue;
			}

			nr_cylinders_final++;

			// Dilate to obtain borders
			cv::dilate(mask,mask_dilated,mask_square_kernel);
			mask_diff = mask_dilated-mask_eroded;

			int stacked_cell_id = 0;
			unsigned char * row_ptr;

			grid_cylinder_seg_map_eroded.setTo((unsigned char)cylinder_code_offset+nr_cylinders_final, mask_eroded > 0);

			// Get variables needed for point-surface distance computation
			Eigen::Vector3f P2 = cylinder_segments[reg_id].P2[sub_reg_id];
			Eigen::Vector3f P1P2 = P2 - cylinder_segments[reg_id].P1[sub_reg_id];
			double P1P2_norm = cylinder_segments[reg_id].P1P2_norm[sub_reg_id];
			double radius = cylinder_segments[reg_id].radii[sub_reg_id];

			// Cell refinement
			for (int cell_r=0; cell_r<nr_vertical_cells; cell_r++){
				row_ptr = mask_diff.ptr<uchar>(cell_r);
				for (int cell_c=0; cell_c<nr_horizontal_cells; cell_c++){
					int offset = stacked_cell_id*nr_pts_per_cell;
					int next_offset = offset+nr_pts_per_cell;
					if(row_ptr[cell_c]>0){
						float max_dist = 9*cylinder_segments[reg_id].MSEs[sub_reg_id];
						// Update cells
						int j=0;
						for(int pt = offset; pt<next_offset;j++,pt++){
							point(2) = cloud_array(pt,2);
							if(point(2)>0){
								point(0) = cloud_array(pt,0);
								point(1) = cloud_array(pt,1);
								dist = P1P2.cross(point-P2).norm()/P1P2_norm - radius;
								dist *= dist;
								if(dist<max_dist && dist<distances_stacked[pt]){ 
									distances_stacked[pt] = dist;
									seg_map_stacked[pt] = cylinder_code_offset + nr_cylinders_final;
								}
							}
						}
					}
					stacked_cell_id++;
				}
			}
		}
	}

	/*------------------------- Copying and rearranging segment data ------------------------*/

	uchar * row_ptr, *stack_ptr;
	uchar * grid_plane_eroded_row_ptr, * grid_cylinder_eroded_row_ptr;
	// Copy inlier list to matrix form
	for (int cell_r=0; cell_r<nr_vertical_cells; cell_r++){
		row_ptr = seg_out.ptr<uchar>(cell_r);
		grid_plane_eroded_row_ptr = grid_plane_seg_map_eroded.ptr<uchar>(cell_r);
		grid_cylinder_eroded_row_ptr = grid_cylinder_seg_map_eroded.ptr<uchar>(cell_r);
		int r_offset = cell_r*cell_height;
		int r_limit = r_offset+cell_height;
		for (int cell_c=0; cell_c<nr_horizontal_cells; cell_c++){
			int c_offset = cell_c*cell_width;
			int c_limit = c_offset+cell_width;

			if (grid_plane_eroded_row_ptr[cell_c]>0){
				// Set rectangle equal to assigned cell
				seg_out(cv::Rect(c_offset,r_offset,cell_width,cell_height)).setTo(grid_plane_eroded_row_ptr[cell_c]);
			}else{
				if(grid_cylinder_eroded_row_ptr[cell_c]>0){
					// Set rectangle equal to assigned cell
					seg_out(cv::Rect(c_offset,r_offset,cell_width,cell_height)).setTo(grid_cylinder_eroded_row_ptr[cell_c]);
				}else{
					// Set cell pixels one by one
					stack_ptr = &seg_map_stacked[nr_pts_per_cell*cell_r*nr_horizontal_cells+nr_pts_per_cell*cell_c];
					for(int r=r_offset;r<r_limit;r++){
						row_ptr = seg_out.ptr<uchar>(r);
						for(int c=c_offset;c<c_limit;c++){
							if(*stack_ptr>0){
								row_ptr[c] = *stack_ptr;
							}
							stack_ptr++;
						}
					}
				}
			}
		}
	}

    for(int i=0;i<nr_cylinders;i++){
        int reg_id = cylinder2region_map[i].first;
        if (reg_id>-1){
            int sub_reg_id = cylinder2region_map[i].second;
            CylinderSeg cy;
            cy.radii.push_back(cylinder_segments[reg_id].radii[sub_reg_id]);
            cy.centers.push_back(cylinder_segments[reg_id].centers[sub_reg_id]);
            copy(cylinder_segments[reg_id].axis, cylinder_segments[reg_id].axis+3, cy.axis);
            cout<<cylinder_segments[reg_id].axis[2]<<endl;
            cylinder_segments_final.push_back(cy);
        }
    }


	// Cleaning data
	stacked_cell_id=0;
	for (int cell_r=0; cell_r<nr_vertical_cells; cell_r++){
		for (int cell_c=0; cell_c<nr_horizontal_cells; cell_c++){
			delete Grid[stacked_cell_id];
			stacked_cell_id++;
		}
	}
	//Grid.clear();
}

void CAPE::getConnectedComponents(cv::Mat & segment_map, MatrixXb & planes_association_matrix){

	int *row, *row_below;
	int nr_rows_2_scan = segment_map.rows-1;
	int nr_cols_2_scan = segment_map.cols-1;
	int px_val;
	for(int r=0; r<nr_rows_2_scan; r++){
		row = segment_map.ptr<int>(r);
		row_below = segment_map.ptr<int>(r+1);
		for(int c=0; c<nr_cols_2_scan; c++){
			px_val = row[c];
			if (px_val>0){
				if (row[c+1]>0 && px_val != row[c+1]){planes_association_matrix(px_val-1,row[c+1]-1) = true;}
				if (row_below[c]>0 && px_val != row_below[c]) planes_association_matrix(px_val-1,row_below[c]-1) = true;
			}
		}
	}
	for(int r=0; r<planes_association_matrix.rows(); r++){
		for(int c=r+1; c<planes_association_matrix.cols(); c++){
			planes_association_matrix(r,c) = planes_association_matrix(r,c) || planes_association_matrix(c,r);
		}
	}
}

// Recursive implementation
// TODO: Test instead iterative implementation
void CAPE::RegionGrowing(unsigned short width, unsigned short height, bool* input, bool* output, vector<PlaneSeg*> & Grid, vector<float> & cell_dist_tols, unsigned short x, unsigned short y, double * normal_1, double d)
{
	int index = x + width*y;
	// If pixel is not part of a component or has already been labelled
	if (input[index]==false || output[index]==true)
		return;

	double * normal_2 = Grid[index]->normal;
	double * m = Grid[index]->mean;
	double d_2 = Grid[index]->d;

	if (normal_1[0]*normal_2[0]+normal_1[1]*normal_2[1]+normal_1[2]*normal_2[2]<min_cos_angle_4_merge || 
		pow(normal_1[0]*m[0] + normal_1[1]*m[1] + normal_1[2]*m[2] + d,2)>cell_dist_tols[index])//max_merge_dist
		return;
	output[index] = true;

	// Now label the 4 neighbours:
	if (x > 0)        RegionGrowing(width, height, input, output, Grid, cell_dist_tols, x-1, y, normal_2, d_2);   // left  pixel
	if (x < width-1)  RegionGrowing(width, height, input, output, Grid, cell_dist_tols, x+1, y, normal_2, d_2);  // right pixel
	if (y > 0)        RegionGrowing(width, height, input, output, Grid, cell_dist_tols, x, y-1, normal_2, d_2);   // upper pixel 
	if (y < height-1) RegionGrowing(width, height, input, output, Grid, cell_dist_tols, x, y+1, normal_2, d_2);   // lower pixel
}



CAPE::~CAPE(void)
{
	Grid.clear();
}
