#ifndef KDT_HPP
#define KDT_HPP

#include <cmath>
#include <iostream>
#include <vector>

enum {APPROXIMATE = 0, BALANCED = 1, THOROUGH = 2};

class Point2D{
    public:
        double x;
        double y;

        double val(bool dim);

        Point2D(double x, double y);
        Point2D(std::vector<double>);
        Point2D();
};

class Node2D{
    public:
        int dim;
        int layer;
        int sl_index = -1;
        int bl_index = -1;
        int p_index = -1;

        double x;
        double y;
        double val(bool dim);

        Node2D(Point2D point){
            this->x = point.x;
            this->y = point.y;
        }
};

class KdTree{
    
    int node_index = 0;

    Node2D populate_node(bool dim, int start, int points_left, int this_index, int parent_index);

    public:
        std::vector<Node2D> kd_node_array;

        int max_layer = 0;
        std::vector<Point2D> points;
        int dim = 0;

        void print_kd_tree();

        void print_fancy_kd_tree();

        std::vector<int> find_closest_point(Point2D& target_point, int dim, char type);

        void make_kd_tree(int dim, int start, int points_left, int parent_index);

        KdTree(std::vector<Point2D> parent_points);
};

#endif