#include <cmath>
#include <iostream>

#include "../data.hpp"
#include "../kdtree.hpp"

void KdTree::populate_node(int dim, int start, int points_left, Node* this_node, Node* starting_node){
    int index = start + std::ceil(points_left/2.0)- 1;
    this_node->dim = dim;
    for (int d = 0; d < dims; d++){
        this_node->val[d] = points[index].val[d];
    }
    this_node->parent = starting_node;
    bool size = false;
    if (starting_node != NULL){
        this_node->layer = starting_node->layer + 1;
        //Node linking, not needed for the first node (has a NULL starting_node)
        dim = cycle_dimensions(dim, false);
        if (this_node->val[dim] < starting_node->val[dim]){
            starting_node->small_link = this_node;
        }
        else{
            starting_node->big_link = this_node;
        }
    }
    else{
        this_node->layer = 0;
    }
}
int KdTree::cycle_dimensions(int dim, bool dir){
    if (dir){
        dim++;
        if (dim >= dims){
            dim = 0;
        }
    }
    else{
        dim--;
        if (dim < 0){
            dim = dims - 1;
        }
    }
    return dim;
}
void KdTree::print_kd_tree(){
    for (int i = 0; i < data_size; i++) {
        std::cout << kd_node_array[i].val[0] << ", " << kd_node_array[i].val[1];
        std::cout << "[Layer " << kd_node_array[i].layer << "] ";
        if (kd_node_array[i].big_link != NULL){
            std::cout << " =+=> " << kd_node_array[i].big_link->val[0] << ", " << kd_node_array[i].big_link->val[1];
        }
        if (kd_node_array[i].small_link != NULL){
            std::cout << " =-=> " << kd_node_array[i].small_link->val[0] << ", " << kd_node_array[i].small_link->val[1];
        }
        std::cout << '\n';
    }
}
Node* KdTree::find_approximate_closest_point(double* target_point, int dim){
    Node* node_ptr = NULL;
    Node* next_node_ptr = &kd_node_array[0];

    while (next_node_ptr != NULL){
        node_ptr = next_node_ptr;
        if (target_point[0] < node_ptr->val[dim]){
            next_node_ptr = node_ptr->small_link;
        }
        else{
            next_node_ptr = node_ptr->big_link;
        }
        dim = cycle_dimensions(dim, true);
    }
    return node_ptr;
}
Node** KdTree::find_closest_point(Point* target_point, int dim, bool find_second_closest){
    Node** possible_node_ptr = new Node*[(max_layer * (max_layer + 1))/2];
    Node** best_node_ptr = new Node*[2];
    best_node_ptr[0] = NULL;
    best_node_ptr[1] = NULL;
    int possible_node_index = 0;
    Node* node_ptr = &kd_node_array[0];
    double distance = 0;
    double border_distance = 0;
    double closest_distance = INFINITY;
    double second_closest_distance = INFINITY;

    while (true){
        while (node_ptr != NULL){
            //We only go forwards and end when the pointer dies, unless there is
            //another branch to inspect.
            distance = sqrt((target_point->val[0] - node_ptr->val[0]) * (target_point->val[0] - node_ptr->val[0]) + (target_point->val[1] - node_ptr->val[1]) * (target_point->val[1] - node_ptr->val[1]));
            if (distance < closest_distance){
                if (find_second_closest){
                    best_node_ptr[1] = best_node_ptr[0];
                    second_closest_distance = closest_distance;
                }
                best_node_ptr[0] = node_ptr;
                closest_distance = distance;
            }
            else if (find_second_closest && distance < second_closest_distance){
                best_node_ptr[1] = node_ptr;
                second_closest_distance = distance;
            }

            dim = node_ptr->dim;
            border_distance = target_point->val[dim] - node_ptr->val[dim];
            //If the distance to the node is larger than the distance to the border,
            //there might be a closer point on the other side.
            if (border_distance < 0){
                //The point is to the left of the node, need to refer to small link
                if (distance > std::abs(border_distance)){
                    //Point is close to the border, will keep the big link in mind for later as well
                    possible_node_ptr[possible_node_index] = node_ptr->big_link;
                    possible_node_index++;
                }
                node_ptr = node_ptr->small_link;
            }
            else{
                //The point is equal to or to the right of the node, need to refer to big link
                if (distance > std::abs(border_distance)){
                    //Point is close to the border, will keep the small link in mind for later as well
                    possible_node_ptr[possible_node_index] = node_ptr->small_link;
                    possible_node_index++;
                }
                node_ptr = node_ptr->big_link;
            }
        }

        if (possible_node_index != 0){
            //The possible node pointer will always be 1 ahead of where it should be.
            //If it is 0, then we would want to look at -1 and that is no good.
            node_ptr = possible_node_ptr[possible_node_index - 1];
            possible_node_index--;
        }
        else{
            //This is the exit condition for the loop, better than a do while in my
            //opinion (saves a branch).
            break;
        }
    }
    if (find_second_closest){
        return best_node_ptr;
    }
    else{
        return &best_node_ptr[0];
    }  
}
void KdTree::make_kd_tree(int dim, int start, int points_left, Node *starting_node){
    for (int i = start; i < (points_left + start); i++){   
        for (int j = i; j < (points_left + start); j++){
            if (points[i].val[dim] > points[j].val[dim]){
                double* temp = new double[dims];
                for (int d = 0; d < dims; d++){
                    temp[d] = points[i].val[d];
                }
                for (int d = 0; d < dims; d++){
                    points[i].val[d] = points[j].val[d];
                }
                for (int d = 0; d < dims; d++){
                    points[j].val[d] = temp[d];
                }
            }
        }
    }
    populate_node(dim, start, points_left, &kd_node_array[node_index], starting_node);
    Node* node = &kd_node_array[node_index];
    if (node->layer > max_layer){
        max_layer = node->layer;
    }
    node_index++;
    if (points_left <= 1){
        //If we hit a dead end we return and leave the links NULL
        //std::cout << kd_node_array[node_index].x_val << '\n';
        return;
    }
    else if (points_left <= 2) {
        dim = cycle_dimensions(dim, true);
        make_kd_tree(dim, start + 1, 1, node);
    }
    else{
        dim = cycle_dimensions(dim, true);
        make_kd_tree(dim, start, std::ceil(points_left/2.0 - 1), node);
        make_kd_tree(dim, start + std::ceil(points_left/2.0), std::floor(points_left/2.0), node);
        return;
    }
}
KdTree::KdTree(Set* parent_set){
    //Constructor
    this->data_size = parent_set->data_size;
    this->dim = 0;
    this->points = new Point[data_size];

    for (int i = 0; i < data_size; i++){
        for (int d = 0; d < dims; d++){
            this->points[i].val[d] = parent_set->points[i].val[d];
        }
    }

    this->kd_node_array = new Node[data_size];

    //Calls a recursive function to make create the object and all its properties
    std::cout << "making tree\n";
    make_kd_tree(false, 0, data_size, NULL);
}