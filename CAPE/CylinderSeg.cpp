/*
 * Copyright 2018 Pedro Proenza <p.proenca@surrey.ac.uk> (University of Surrey)
 *
 */
#include "CylinderSeg.h"

CylinderSeg::CylinderSeg(vector<PlaneSeg*> & Grid, bool * activated_mask,int nr_cells_activated)
{

	int nr_samples = Grid.size();
	int m = nr_cells_activated;
	nr_segments = 0;
	local2global_map = new int[nr_cells_activated];

	Eigen::MatrixXd N(3,2*m);
	Eigen::MatrixXd P(3,m);

	// Init. P and N
	int j =0;
	for(int i=0; i<nr_samples;i++){
		if (activated_mask[i]){
			N(0,j) = Grid[i]->normal[0];
			N(1,j) = Grid[i]->normal[1];
			N(2,j) = Grid[i]->normal[2];
			P(0,j) = Grid[i]->mean[0];
			P(1,j) = Grid[i]->mean[1];
			P(2,j) = Grid[i]->mean[2];
			local2global_map[j] = i;
			j++;
		}
	}

	// Concatenate [N -N]
	for(int i=0; i<nr_samples;i++){
		if (activated_mask[i]){
			N(0,j) = -Grid[i]->normal[0];
			N(1,j) = -Grid[i]->normal[1];
			N(2,j) = -Grid[i]->normal[2];
			j++;
		}
	}

	// Compute covariance
	Eigen::MatrixXd cov = (N*N.adjoint()) / double(N.cols() - 1);

	// PCA using QR decomposition for symmetric matrices
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
	Eigen::Vector3d S = es.eigenvalues();
	double score = S(2)/S(0);

	// Checkpoint 1
	if(score<cylinder_score_min){
		return;
	}

	Eigen::Vector3d vec = es.eigenvectors().col(0);
	axis[0] = vec(0); axis[1] = vec(1); axis[2] = vec(2);

	Eigen::MatrixXd N_cpy = N.block(0,0,3,m);
	N = N_cpy; /* This avoids memory issues */
	Eigen::MatrixXd P_proj(3,m);

	// Projection to plane: P' = P-theta*<P.theta>
	Eigen::MatrixXd P_dot_theta = vec.transpose()*P;
	P_proj.row(0) = P.row(0)-P_dot_theta*vec(0);
	P_proj.row(1) = P.row(1)-P_dot_theta*vec(1);
	P_proj.row(2) = P.row(2)-P_dot_theta*vec(2);
	Eigen::MatrixXd N_dot_theta = vec.transpose()*N;
	N.row(0) -= N_dot_theta*vec(0);
	N.row(1) -= N_dot_theta*vec(1);
	N.row(2) -= N_dot_theta*vec(2);

	// Normalize projected normals
	Eigen::MatrixXd Normals_norm = N.colwise().norm();
	N.row(0) = N.row(0).array()/Normals_norm.array();
	N.row(1) = N.row(1).array()/Normals_norm.array();
	N.row(2) = N.row(2).array()/Normals_norm.array();

	// Ransac params
	float p_success = 0.8;
	float w = 0.33;
	float K = log(1-p_success)/log(1-pow(w,3));

	int m_left = m;
	vector<int> ids_left;
	MatrixXb ids_left_mask(1,m);
	for(int i=0; i<m; i++){
		ids_left.push_back(i);
		ids_left_mask(i) = true;
	}
	// Sequential RANSAC main loop
	while(m_left>5 && m_left>0.1*m){

		int id_1, id_2, id_3;
		Eigen::MatrixXd A, B, center, e1, e2;
		Eigen::MatrixXd D(1,m);
		MatrixXb I(1,m);
		MatrixXb I_final(1,m);
		double a, b, r;
		double dist = 0;
		double min_hypothesis_dist = cylinder_RANSAC_sqr_max_dist*m_left;
		int nr_inliers_accepted = 0.9*m_left;
		int nr_inliers(0), max_nr_inliers(0);

		// RANSAC loop
		int k = 0;
		while(k<K){

			// Random triplet
			id_1 = ids_left[rand()%m_left];
			id_2 = ids_left[rand()%m_left];
			id_3 = ids_left[rand()%m_left];

			e1 = (N.col(id_1)+ N.col(id_2)+ N.col(id_3));
			e2 = (P_proj.col(id_1)+ P_proj.col(id_2)+ P_proj.col(id_3));

			// LLS solution for triplets
			A = e1.transpose()*e1;
			a = 1-A(0)/9;
			b = (N.col(id_1).array()*P_proj.col(id_1).array() + N.col(id_2).array()*P_proj.col(id_2).array() + N.col(id_3).array()*P_proj.col(id_3).array()).sum()/3 
				- (e1.transpose()*e2)(0)/9;
			r = b/a;

			center = (e2- r*e1)/3;

			// Unnecessary calculations here
			// Normal dist
			D = ((P_proj-r*N).colwise()-center.col(0)).colwise().squaredNorm()/(r*r);

			// Rectify radius if concave
			if (r<0)
				r = -r;

			// Inliers
			I = D.array()<cylinder_RANSAC_sqr_max_dist;

			//MSAC truncated distance
			dist = 0;
			nr_inliers = 0;
			for(int i=0;i<m;i++){
				if(ids_left_mask(i)){
					if(I(i)){
						nr_inliers++;
						dist += D(i);
					}else{
						dist += cylinder_RANSAC_sqr_max_dist;
					}
				}
			}

			if(dist<min_hypothesis_dist){
				min_hypothesis_dist = dist;
				max_nr_inliers = nr_inliers;
				for(int i=0;i<m;i++){
					if(ids_left_mask(i)){
						I_final(i) = I(i);
					}else{
						I_final(i) = false;
					}
				}
				if(nr_inliers>nr_inliers_accepted)
					break;
			}
			k++;
		}

		// Checkpoint 2
		if(max_nr_inliers<6)
			break;

		// Increase prob. of finding inlier for next RANSAC runs
		K = log(1-p_success)/log(1-pow(0.5,3));

		// Remove cells from list of remaining cells
		ids_left.clear();
		for(int i=0; i<m; i++){
			if(I_final(i)){
				ids_left_mask(i) = false;
				m_left--;
			}else{
				if(ids_left_mask(i))
					ids_left.push_back(i);
			}
		}

		// LLS solution using all inliers
		e1.setZero();
		e2.setZero();
		b = 0;

		for(int i=0; i<m;i++){
			if(I_final(i)){
				e1 += N.col(i);
				e2 += P_proj.col(i);
				b += (N.col(i).array()*P_proj.col(i).array()).sum();
			}
		}

		A = e1.transpose()*e1;

		a = 1-A(0)/(max_nr_inliers*max_nr_inliers);
		b /= max_nr_inliers;
		b -= (e1.transpose()*e2)(0)/(max_nr_inliers*max_nr_inliers);
		r = b/a;
		center = (e2- r*e1)/max_nr_inliers;

		// Rectify radius if concave
		if (r<0)
			r = -r;

		// Add cylinder
		nr_segments++;
		radii.push_back((float)r);
		centers.push_back(center);
		inliers.push_back(I_final);

		// Save points on axis
		Eigen::Vector3d P1d = center;
		Eigen::Vector3d P2d= center+vec;
		double P1P2d = (P2d-P1d).norm();
		// Use point-to-line distances (same as cylinder distances)
		Eigen::Vector3d P3;
		for(int i=0; i<m;i++){
			if(I_final(i)){
				P3 = P.block<3,1>(0,i);
				//D(i) = (P3-center).norm()-r;
				D(i) = ((P2d-P1d).cross(P3-P2d)).norm()/P1P2d-r;
			}
		}
		D = D.array().square();

		double MSE = 0; 
		for(int i=0; i<m;i++){
			if(I_final(i))
				MSE += D(i);
		}
		MSE = MSE/max_nr_inliers;
		MSEs.push_back(MSE);

		// Save points on axis, useful for computing distances afterwards
		P1.push_back(P1d.cast<float>()); 
		P2.push_back(P2d.cast<float>());
		P1P2_norm.push_back((float)P1P2d);
		cylindrical_mask.push_back(true);
	}
}

float CylinderSeg::distance(Eigen::Vector3f & P3, int id){
	return ((P2[id]-P1[id]).cross(P3-P2[id])).norm()/P1P2_norm[id]-radii[id];

}

CylinderSeg::~CylinderSeg(void)
{
}
