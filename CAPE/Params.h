/*
 * Copyright 2018 Pedro Proenza <p.proenca@surrey.ac.uk> (University of Surrey)
 *
 */
// Kinect 1 sensor model uncertainty according to Khoshelham and Elberink:
const double DEPTH_SIGMA_COEFF =  0.000001425;
const double DEPTH_SIGMA_MARGIN =  10;
const double cylinder_score_min = 100;
const double cylinder_RANSAC_sqr_max_dist = 0.0225; /* square of 15 %*/
