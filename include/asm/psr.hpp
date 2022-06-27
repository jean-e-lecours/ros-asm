#ifndef PSR_HPP
#define PSR_HPP

#include <iostream>
#include <vector>

#include "kdt.hpp"

enum corr_type {SINGLE = 0, DOUBLE = 1};

class Transform2D{
    public:
        double rot_mat[4];
        double trans_vec[2];
        double z_rot;

        void print_transform();

        void transform(std::vector<Point2D>& set);

        void add_transform(Transform2D& added_transform);

        void update_transform(std::vector<double>& vector);

        bool is_significant(double threshold);

        Transform2D(double x_trans, double y_trans, double z_rot);
};

class Correlation{
    public:
        Point2D scan;
        Point2D mcor1;
        Point2D mcor2;

        std::vector<double> norm = {0,0};
        double corrected_value;

        double get_distance();

        std::vector<double> get_trans();

        Correlation(KdTree& map_kdt, Point2D& scan_point, char corr_type, Transform2D& g_transf);
};

class LaserData{
    
    public:
        std::vector<double> las_vec;
        int use_ranges_size;
        double angle_incr;

        std::vector<Point2D> map_scan_points(Transform2D& transform);
};

std::vector<Point2D> map_scan_points(Transform2D& transform, std::vector<double> las_vec, double scan_period);

std::vector<Point2D> make_map(std::vector<double> map_vec, int grid_size_x, int grid_size_y, double pix_res);
#endif