#ifndef KDTREE_H
#define KDTREE_H

#include <cmath>
#include <iostream>

#include "data.hpp"

class Node: public Point{
    public:
        int dim;
        int layer;
        Node* small_link;
        Node* big_link;
        Node* parent;
};

class KdTree{
    Node* kd_node_array;
    int node_index = 0;

    void populate_node(int dim, int start, int points_left, Node* this_node, Node* starting_node);

    int cycle_dimensions(int dim, bool dir);

    public:
        int data_size;
        int max_layer = 0;
        Point* points;
        static int dims;
        int dim = 0;

        void print_kd_tree();

        Node* find_approximate_closest_point(double* target_point, int dim);

        Node** find_closest_point(Point* target_point, int dim, bool find_second_closest);

        void make_kd_tree(int dim, int start, int points_left, Node* starting_node);

        KdTree(Set* parent_set);
};

#endif