#include "../psr.hpp"
#include "../kdt.hpp"

#include <cmath>
#include <iostream>
#include <vector>
#include <bits/stdc++.h>

void Transform2D::transform(std::vector<Point2D>& set){
    for (int i = 0; i < set.size(); i++){
        //Applies a single matrix multiplication with the provided x, y vector from the set
        double temp_x = rot_mat[0] * set[i].x + rot_mat[1] * set[i].y + trans_vec[0];
        double temp_y = rot_mat[2] * set[i].x + rot_mat[3] * set[i].y + trans_vec[1];
        set[i].x = temp_x;
        set[i].y = temp_y;
    }
}
void Transform2D::print_transform(){
    std::cout << "⎡" << rot_mat[0] << "\t " << rot_mat[1] << "\t " << trans_vec[0] << "⎤\n" \
              << "|" << rot_mat[2] << "\t " << rot_mat[3] << "\t " << trans_vec[1] << "|\n" \
              << "⎣0\t 0\t 1⎦\n";
}

void Transform2D::add_transform(Transform2D& added_transform){
    Transform2D temp_transform(0, 0, 0);
    temp_transform.rot_mat[0] = rot_mat[0] * added_transform.rot_mat[0] + rot_mat[1] * added_transform.rot_mat[2];
    temp_transform.rot_mat[1] = rot_mat[0] * added_transform.rot_mat[1] + rot_mat[1] * added_transform.rot_mat[3];
    temp_transform.rot_mat[2] = rot_mat[2] * added_transform.rot_mat[0] + rot_mat[3] * added_transform.rot_mat[2];
    temp_transform.rot_mat[3] = rot_mat[2] * added_transform.rot_mat[1] + rot_mat[3] * added_transform.rot_mat[3];
    temp_transform.trans_vec[0] = rot_mat[0] * added_transform.trans_vec[0] + rot_mat[1] * added_transform.trans_vec[1] + trans_vec[0];
    temp_transform.trans_vec[1] = rot_mat[2] * added_transform.trans_vec[0] + rot_mat[3] * added_transform.trans_vec[1] + trans_vec[1];
    this->rot_mat[0] = temp_transform.rot_mat[0];
    this->rot_mat[1] = temp_transform.rot_mat[1];
    this->rot_mat[2] = temp_transform.rot_mat[2];
    this->rot_mat[3] = temp_transform.rot_mat[3];
    this->trans_vec[0] = temp_transform.trans_vec[0];
    this->trans_vec[1] = temp_transform.trans_vec[1];
}

Transform2D::Transform2D(double x_trans, double y_trans, double z_rot){
    //Constructor
    //Sets the two elements of the Transform2Ds as per the desired movement
    rot_mat[0] = cos(z_rot);
    rot_mat[1] = -sin(z_rot);
    rot_mat[2] = -rot_mat[1];
    rot_mat[3] = rot_mat[0];
    trans_vec[0] = x_trans;
    trans_vec[1] = y_trans;
    

    this->z_rot = z_rot;
}
std::vector<Point2D> map_scan_points(Transform2D& transform, std::vector<double> las_vec, double scan_period){
    std::vector<Point2D> mapped_pts_vec;
    int j = 0;
    bool skip = true;
    for (int i = 0; i < las_vec.size(); i++){
        Point2D new_point(0,0);
        if (las_vec[i] < 5){
            double temp_scan_x = las_vec[i] * cos(scan_period * i);
            double temp_scan_y = las_vec[i] * sin(scan_period * i);
            new_point.x = transform.trans_vec[0] + temp_scan_x * cos(transform.z_rot) - temp_scan_y * sin(transform.z_rot);
            new_point.y = transform.trans_vec[1] + temp_scan_y * cos(transform.z_rot) + temp_scan_x * sin(transform.z_rot);
            j++;
            mapped_pts_vec.push_back(new_point);
        }
    }
    return mapped_pts_vec;
}

void Transform2D::update_transform(std::vector<double>& vector){
    rot_mat[0] = vector[0];
    rot_mat[1] = vector[1];
    rot_mat[2] = vector[2];
    rot_mat[3] = vector[3];
    trans_vec[0] = vector[4];
    trans_vec[1] = vector[5];

    z_rot = std::atan2(rot_mat[2],rot_mat[0]);
}

bool Transform2D::is_significant(double threshold){
    std::vector<double> comparison = {0,0,0,0,0,0};
    comparison[0] = std::abs(1-std::abs(rot_mat[0]));
    comparison[1] = std::abs(rot_mat[1]);
    comparison[2] = std::abs(rot_mat[2]);
    comparison[3] = std::abs(1-std::abs(rot_mat[3]));
   
    comparison[4] = std::abs(trans_vec[0]);
    comparison[5] = std::abs(trans_vec[1]);

    if (*max_element(comparison.begin(), comparison.end()) < threshold){
        return false;
    }
    else{
        return true;
    }

}

double Transform2D::compare_transform(std::vector<double> vector){
    double transform_diff = 0;
    transform_diff += (this->rot_mat[0] - vector[0]) * (this->rot_mat[0] - vector[0]);
    transform_diff += (this->rot_mat[1] - vector[1]) * (this->rot_mat[1] - vector[1]);
    transform_diff += (this->rot_mat[2] - vector[2]) * (this->rot_mat[2] - vector[2]);
    transform_diff += (this->rot_mat[3] - vector[3]) * (this->rot_mat[3] - vector[3]);
    transform_diff += (this->trans_vec[0] - vector[4]) * (this->trans_vec[0] - vector[4]);
    transform_diff += (this->trans_vec[1] - vector[5]) * (this->trans_vec[1] - vector[5]);
    return transform_diff/6;
}

Correlation::Correlation(KdTree& map_kdt, Point2D& scan_point, char corr_type, Transform2D& g_transf){
    scan = scan_point;
    std::vector<int> ids = map_kdt.find_closest_point(scan_point, 0, corr_type);
    mcor1.x = map_kdt.kd_node_array[ids[0]].x;
    mcor1.y = map_kdt.kd_node_array[ids[0]].y;
    mcor2.x = map_kdt.kd_node_array[ids[1]].x;
    mcor2.y = map_kdt.kd_node_array[ids[1]].y;

    norm[0] = mcor1.y - mcor2.y;
    norm[1] = -(mcor1.x - mcor2.x);
    double norm_size = std::sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
    norm[0] = norm[0]/norm_size;
    norm[1] = norm[1]/norm_size;

    corrected_value = norm[0]*(g_transf.rot_mat[0]*scan_point.x + g_transf.rot_mat[1]*scan_point.y + g_transf.trans_vec[0] - mcor1.x) + \
                                norm[1]*(g_transf.rot_mat[2]*scan_point.x + g_transf.rot_mat[3]*scan_point.y + g_transf.trans_vec[1] - mcor1.y);
}

std::vector<Point2D> make_map(std::vector<double> map_vec, int grid_size_x, int grid_size_y, double pix_res){
    std::vector<Point2D> map;
    if (map_vec.size() != grid_size_x*grid_size_y){
        std::cerr << "Map vector size does not match given grid dimensions! Check the function inputs!\n";
        return map;
    }
    for (int i = 0; i < map_vec.size(); i++){
        if (map_vec[i] > 0){
            double point_x = pix_res * (i%grid_size_y);
            double point_y = pix_res * (int)(i/grid_size_x);
            Point2D map_point(point_x, point_y);
            map.push_back(map_point);
        }
    }
    return map;
}

std::vector<Point2D> place_map_points(std::vector<int> map_vec){
    std::vector<Point2D> map;

    return map;
}