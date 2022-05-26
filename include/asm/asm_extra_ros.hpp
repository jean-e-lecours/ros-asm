#ifndef ASM_H
#define ASM_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "kdtree.hpp"

class Transform2D{
    public:
        double rot_mat[4];
        double trans_vec[2];
        double z_rot;

        double compare_transform(double* vector);

        void update_transform(double* vector);

        void print_transform();

        void transform_set(Set* set);

        Point transform_referenced_point(double x, double y);

        void add_transform(Transform2D* added_transform);

        void set_transform(Transform2D* transform);

        Transform2D(double x_trans, double y_trans, double z_rot);
};

class Correlation{
    
    Point* model_corr_2;
    Node** temp_node;
    Point** temp_point;
    double norm_size;

    public:
        Point current_point;
        Point* model_corr_1;
        double norm_x;
        double norm_y;
        double corrected_value;

        static double std_dev;
        static double mean;

        Correlation(KdTree* model_tree, Set* data_set, int index, Transform2D* transform);
        Correlation();
};

class LaserData{
    
    public:
        std::vector<double> las_vec;
        int use_ranges_size;
        double angle_incr;

        Set map_scan_points(Transform2D* transform);
};

class TextMapData{
    bool done = false;

    public:
        Set* map_set;
        int number_of_points;

        int read_from_file(std::string filename, int size_x, int size_y, double pix_res);
};

double* make_g_vector(Correlation* corrs, Set* data_set, double std_dev, double correntropy_factor, int data_size, bool print);
double* make_G_matrix(Correlation* corrs, Set* data_set, double std_dev, double correntropy_factor, int data_size, bool print);
double* solve_system(double* g, double* G, bool print);

#endif