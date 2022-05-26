#include "../kdt.hpp"
#include <cmath>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

double Point2D::val(bool dim){
    if (!dim){
        return x;
    }
    else{
        return y;
    }
}
Point2D::Point2D(double x_val, double y_val){
    x = x_val;
    y = y_val;
}
Point2D::Point2D(std::vector<double> val_vec){
    int init_size = val_vec.size();
    if (init_size > 2){
        std::cerr << "Bad number of input values for 2D vector, tried to use " << init_size << " values\n";
    }
    else {
        x = val_vec[0];
        y = val_vec[1];
    }
}
Point2D::Point2D(){
    x = 0;
    y = 0;
}

double Node2D::val(bool dim){
    if (!dim){
        return x;
    }
    else{
        return y;
    }
}

Node2D KdTree::populate_node(bool dim, int start, int points_left, int this_index, int parent_index){
    int point_index = start + std::ceil(points_left/2.0)- 1;
    Node2D new_node(points[point_index]);
    new_node.dim = dim;
    new_node.x = points[point_index].x;
    new_node.y = points[point_index].y;

    new_node.p_index = parent_index;
    bool size = false;
    if (parent_index >= 0){
        new_node.layer = kd_node_array[parent_index].layer + 1;
        //Node linking, not needed for the first node (has a NULL starting_node)
        dim = !dim;
        if (new_node.val(dim) < kd_node_array[parent_index].val(dim)){
            kd_node_array[parent_index].sl_index = this_index;
        }
        else{
            kd_node_array[parent_index].bl_index = this_index;
        }
    }
    else{
        new_node.layer = 0;
    }
    return new_node;
}

void KdTree::make_kd_tree(int dim, int start, int points_left, int parent_index){
    for (int i = start; i < (points_left + start); i++){   
        for (int j = i; j < (points_left + start); j++){
            if (points[i].val(dim) > points[j].val(dim)){
                Point2D temp(points[i].x, points[i].y);
                points[i] = points[j];
                points[j] = temp;
            }
        }
    }
    kd_node_array.push_back(populate_node(dim, start, points_left, node_index, parent_index));

    if (kd_node_array[node_index].layer > max_layer){
        max_layer = kd_node_array[node_index].layer;
    }
    int this_node_index = node_index;
    node_index++;
    if (points_left <= 1){
        //If we hit a dead end we return and leave the links NULL
        //std::cout << kd_node_array[node_index].x_val << '\n';
        return;
    }
    else if (points_left <= 2) {
        dim = !dim;
        make_kd_tree(dim, start + 1, 1, this_node_index);
    }
    else{
        dim = !dim;
        make_kd_tree(dim, start, std::ceil(points_left/2.0 - 1), this_node_index);
        make_kd_tree(dim, start + std::ceil(points_left/2.0), std::floor(points_left/2.0), this_node_index);
        return;
    }
}

std::vector<int> KdTree::find_closest_point(Point2D& target_point, int dim, char type){
    std::vector<int> possible_node_vec;
    std::vector<int> best_node_ind{-1,-1};
    int possible_node_index = 0;
    int node_index = 0;
    double distance = 0;
    double border_distance = 0;
    double closest_distance = INFINITY;
    double second_closest_distance = INFINITY;

    bool border_cross = false;
    if (type == APPROXIMATE){
        border_cross = true;
    }

    int counter = 0;

    while (true){

        while (node_index >= 0){
            //We only go forwards and end when the index dies, unless there is
            //another branch to inspect.
            distance = sqrt((target_point.x - kd_node_array[node_index].x) * (target_point.x - kd_node_array[node_index].x) + (target_point.y - kd_node_array[node_index].y) * (target_point.y - kd_node_array[node_index].y));
            if (distance < closest_distance){
                best_node_ind[1] = best_node_ind[0];
                second_closest_distance = closest_distance;

                best_node_ind[0] = node_index;
                closest_distance = distance;
            }
            else if (distance < second_closest_distance){
                best_node_ind[1] = node_index;
                second_closest_distance = distance;
            }

            dim = kd_node_array[node_index].dim;
            border_distance = target_point.val(dim) - kd_node_array[node_index].val(dim);
            //If the distance to the node is larger than the distance to the border,
            //there might be a closer point on the other side.
            if (border_distance < 0){
                //The point is to the left of the node, need to refer to small link
                if (closest_distance > std::abs(border_distance) && !border_cross){
                    //Point is close to the border, will keep the big link in mind for later as well
                    possible_node_vec.push_back(kd_node_array[node_index].bl_index);
                    possible_node_index++;
                }
                node_index = kd_node_array[node_index].sl_index;
            }
            else{
                //The point is equal to or to the right of the node, need to refer to big link
                if (closest_distance > std::abs(border_distance) && !border_cross){
                    //Point is close to the border, will keep the small link in mind for later as well
                    possible_node_vec.push_back(kd_node_array[node_index].sl_index);
                    possible_node_index++;
                }
                node_index = kd_node_array[node_index].bl_index;
            }
        }

        if (possible_node_index > 0){
            //The possible node pointer will always be 1 ahead of where it should be.
            //If it is 0, then we would want to look at -1 and that is no good.
            if (type == BALANCED){
                border_cross = true;
            }
            node_index = possible_node_vec[possible_node_index - 1];
            possible_node_vec.pop_back();
            possible_node_index--;
        }
        else{
            //This is the exit condition for the loop, better than a do while in my
            //opinion (saves a branch).
            break;
        }
        counter ++;
    }
    return best_node_ind;
}

void KdTree::print_fancy_kd_tree(){
    if (max_layer > 10){
        std::cerr << "Kd Tree is too big! Printing it like this will not look good!\n";
        return;
    }
    std::vector<int> indices{0};
    int pos = 0;
    for (int l = 0; l < max_layer; l++){
        for (int p = pos; p < (pos + pow(2,l)); p++){
            int curr_index = indices[p];
            if (curr_index < 0){
                indices.push_back(-1);
                indices.push_back(-1);
            }
            else{
                indices.push_back(kd_node_array[indices[p]].sl_index);
                indices.push_back(kd_node_array[indices[p]].bl_index);
            }  
        }
        pos += pow(2,l);
    }
    int total_spaces = 17;
    for (int i = 1; i <= max_layer; i++){
        total_spaces = total_spaces*2+1;
    }
    int num_index = 0;
    for (int l = 0; l <= max_layer; l++){
        for (int s = 0; s <= pow(2,l); s++){
            int spaces = 0;
            bool print_number = true;
            spaces = (total_spaces + 1 - (11 * pow(2, l)))/(pow(2,l));
            if (s == 0 || s == pow(2,l)){
                spaces = (spaces - 1)/2;
                if (s > 0){
                    print_number = false;
                }
            }
            for (int i = 0; i < spaces; i++){
                std::cout << ' ';
            }
            if (print_number){
                if (indices[num_index] < 0){
                    std::cout << "<         >";
                }
                else{
                    double tx = kd_node_array[indices[num_index]].x;
                    tx = trunc(tx*1000)/1000;
                    double ty = kd_node_array[indices[num_index]].y;
                    ty = trunc(ty*1000)/1000;
                    std::cout << tx << ',' << ty;
                }
                num_index++;
            }
            else {
                std::cout << '\n';
            }
        }
    }

    return;
}

void KdTree::print_kd_tree(){
    for (int i = 0; i < kd_node_array.size(); i++) {
        std::cout << kd_node_array[i].x << ", " << kd_node_array[i].y;
        std::cout << "[Layer " << kd_node_array[i].layer << "] ";
        if (kd_node_array[i].bl_index >= 0){
            std::cout << " =+=> " << kd_node_array[kd_node_array[i].bl_index].x << ", " << kd_node_array[kd_node_array[i].bl_index].y;
        }
        if (kd_node_array[i].sl_index >= 0){
            std::cout << " =-=> " << kd_node_array[kd_node_array[i].sl_index].x << ", " << kd_node_array[kd_node_array[i].sl_index].y;
        }
        std::cout << '\n';
    }
}

KdTree::KdTree(std::vector<Point2D> parent_set){
    //Constructor
    this->points = parent_set;

    //Calls a recursive function to make create the object and all its properties
    std::cout << "making tree\n";
    make_kd_tree(false, 0, points.size(), -1);
}