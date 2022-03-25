#include "../data.hpp"
#include <iostream>
#include <vector>

Point::Point(){
    this->val = new double[dims];
}
Point::Point(int set_dims){
    this->val = new double[set_dims];
}


Set* Set::copy_set(){
    Set* new_set = new Set(data_size, points);
    return new_set;
}
void Set::print_set(){
    for (int i = 0; i < data_size; i++){
        for (int d = 0; d < dims; d++){
            std::cout << points[i].val[d];
            if (d+1 == dims){
                std::cout << '\n';
            }
            else{
                std::cout << ", ";
            }
        }
    }
}
Set::Set(int data_size, double** point_vals){
    //Constructor
    this->data_size = data_size;
    this->points = new Point[data_size];

    for (int i = 0; i < data_size; i++){
        for (int d = 0; d < dims; d++){
            this->points[i].val[d] = point_vals[d][i];
        }
    }
}
Set::Set(int data_size, Point* points){
    //Constuctor for copy only
    this->data_size = data_size;
    this->points = new Point[data_size];

    for (int i = 0; i < data_size; i++){
        for (int d = 0; d < dims; d++){
            this->points[i].val[d] = points[i].val[d];
        }
    }
}
Set::Set(std::vector<Point> points){
    this->data_size = points.size();
    this->points = new Point[this->data_size];

    for (int i = 0; i < this->data_size; i++){
        for (int d = 0; d < dims; d++){
            this->points[i].val[d] = points[i].val[d];
        }
    }
}