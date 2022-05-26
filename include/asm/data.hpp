#ifndef DATA_H
#define DATA_H

#include <iostream>
#include <vector>

class Point{
    public:
        static int dims;
        double* val;

        explicit Point(int set_dims);
        explicit Point();
};

class Set{
    public:
        static int dims;
        int data_size;
        int max_layer = 0;
        Point* points;

        Set* copy_set();

        void print_set();

    Set(int data_size, double** point_vals);

    Set(int data_size, Point* points);

    Set(std::vector<Point> points);
};

#endif