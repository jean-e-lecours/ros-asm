#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include "../include/gnuplot.h"

#include "ros/ros.h"
#include "std_msgs/String.h"
#include "sensor_msgs/LaserScan.h"
#include "geometry_msgs/Pose.h"
#include "geometry_msgs/TransformStamped.h"
#include "tf2_ros/transform_listener.h"
#include "tf2_geometry_msgs/tf2_geometry_msgs.h"

struct PARAMS{
    int dimensions;
    int max_guesses;
    double corr_factor;
    double transf_tresh;
    bool init = true;
    bool pub_ready = false;
    double start_x;
    double start_y;
    double start_t;
};

PARAMS params;

class Point{
    public:
        double* val;

        explicit Point(int dims){
            this->val = new double[dims];
        }
        explicit Point(){
            this->val = new double[params.dimensions];
        }
};

class Node: public Point{
    public:
        int dim;
        int layer;
        Node* small_link;
        Node* big_link;
        Node* parent;
};

class Set{
    public:
        int data_size;
        int max_layer = 0;
        Point* points;

        Set* copy_set(){
            Set* new_set = new Set(data_size, points);
            return new_set;
        }

        void print_data(){
            std::cout << data_size;
            for (int i = 0; i < data_size; i++){
                for (int d = 0; d < params.dimensions; d++){
                    std::cout << points[i].val[d];
                    if (d+1 == params.dimensions){
                        std::cout << '\n';
                    }
                    else{
                        std::cout << ", ";
                    }
                }
            }
        }

    Set(int data_size, double** point_vals){
        //Constructor
        this->data_size = data_size;
        this->points = new Point[data_size];

        for (int i = 0; i < data_size; i++){
            for (int d = 0; d < params.dimensions; d++){
                this->points[i].val[d] = point_vals[d][i];
            }
        }
    }
    Set(int data_size, Point* points){
        //Constuctor for copy only
        this->data_size = data_size;
        this->points = new Point[data_size];

        for (int i = 0; i < data_size; i++){
            for (int d = 0; d < params.dimensions; d++){
                this->points[i].val[d] = points[i].val[d];
            }
        }
    }
    ~Set(){
        //Deconstructor
    }
};

class KdTree{
    Node* kd_node_array;
    int node_index = 0;
    int dims;

    void populate_node(int dim, int start, int points_left, Node* this_node, Node* starting_node){
        int index = start + std::ceil(points_left/2.0)- 1;
        this_node->dim = dim;
        for (int d = 0; d < params.dimensions; d++){
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

    int cycle_dimensions(int dim, bool dir){
        if (dir){
            dim++;
            if (dim >= params.dimensions){
                dim = 0;
            }
        }
        else{
            dim--;
            if (dim < 0){
                dim = params.dimensions - 1;
            }
        }
        return dim;
    }

    public:
        int data_size;
        int max_layer = 0;
        Point* points;
        int dim = 0;

        void print_kd_tree(){
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

        Node* find_approximate_closest_point(double* target_point, int dim){
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

        Node** find_closest_point(Point* target_point, int dim, bool find_second_closest){
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
                        if (distance > abs(border_distance)){
                            //Point is close to the border, will keep the big link in mind for later as well
                            possible_node_ptr[possible_node_index] = node_ptr->big_link;
                            possible_node_index++;
                        }
                        node_ptr = node_ptr->small_link;
                    }
                    else{
                        //The point is equal to or to the right of the node, need to refer to big link
                        if (distance > abs(border_distance)){
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

        void make_kd_tree(int dim, int start, int points_left, Node* starting_node){
            for (int i = start; i < (points_left + start); i++){
                
                for (int j = i; j < (points_left + start); j++){
                    if (points[i].val[dim] > points[j].val[dim]){
                        double* temp = new double[params.dimensions];
                        for (int d = 0; d < params.dimensions; d++){
                            temp[d] = points[i].val[d];
                        }
                        for (int d = 0; d < params.dimensions; d++){
                            points[i].val[d] = points[j].val[d];
                        }
                        for (int d = 0; d < params.dimensions; d++){
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

    KdTree(Set* parent_set){
        //Constructor
        this->data_size = parent_set->data_size;
        this->dims = params.dimensions-1;
        this->dim = 0;
        this->points = new Point[data_size];

        for (int i = 0; i < data_size; i++){
            for (int d = 0; d < params.dimensions; d++){
                this->points[i].val[d] = parent_set->points[i].val[d];
            }
        }

        this->kd_node_array = new Node[data_size];

        //Calls a recursive function to make create the object and all its properties
        //std::cout << "making tree\n";
        make_kd_tree(false, 0, data_size, NULL);
    }
    ~KdTree(){
        //Destructor
    }
};

class Transform2D{
    public:
        double rot_mat[4];
        double trans_vec[2];
        double z_rot;

        double compare_transform(double* vector){
            double transform_diff = 0;
            transform_diff += (this->rot_mat[0] - vector[0]) * (this->rot_mat[0] - vector[0]);
            transform_diff += (this->rot_mat[1] - vector[1]) * (this->rot_mat[1] - vector[1]);
            transform_diff += (this->rot_mat[2] - vector[2]) * (this->rot_mat[2] - vector[2]);
            transform_diff += (this->rot_mat[3] - vector[3]) * (this->rot_mat[3] - vector[3]);
            transform_diff += (this->trans_vec[0] - vector[4]) * (this->trans_vec[0] - vector[4]);
            transform_diff += (this->trans_vec[1] - vector[5]) * (this->trans_vec[1] - vector[5]);
            return transform_diff/6;
        }

        void update_transform(double* vector){
            this->rot_mat[0] = vector[0];
            this->rot_mat[1] = vector[1];
            this->rot_mat[2] = vector[2];
            this->rot_mat[3] = vector[3];
            this->trans_vec[0] = vector[4];
            this->trans_vec[1] = vector[5];
        }

        void print_transform(){
            std::cout << "⎡" << rot_mat[0] << "\t " << rot_mat[1] << "\t " << trans_vec[0] << "⎤\n" \
                      << "|" << rot_mat[2] << "\t " << rot_mat[3] << "\t " << trans_vec[1] << "|\n" \
                      << "⎣0\t 0\t 1⎦\n";
        }

        void transform_set(Set* set){
            for (int i = 0; i < set->data_size; i++){
                //Applies a single matrix multiplication with the provided x, y vector from the set
                double temp_x = rot_mat[0] * set->points[i].val[0] + rot_mat[1] * set->points[i].val[1] + trans_vec[0];
                double temp_y = rot_mat[2] * set->points[i].val[0] + rot_mat[3] * set->points[i].val[1] + trans_vec[1];
                set->points[i].val[0] = temp_x;
                set->points[i].val[1] = temp_y;
            }
        }
        Point transform_referenced_point(double x, double y){
            double temp_x = rot_mat[0] * x + rot_mat[1] * y + trans_vec[0];
            double temp_y = rot_mat[2] * x + rot_mat[3] * y + trans_vec[1];
            Point point(2);
            point.val[0] = temp_x;
            point.val[1] = temp_y;
            return point;
        }
        void add_transform(Transform2D* added_transform){
            Transform2D temp_transform(0, 0, 0);
            temp_transform.rot_mat[0] = rot_mat[0] * added_transform->rot_mat[0] + rot_mat[1] * added_transform->rot_mat[2];
            temp_transform.rot_mat[1] = rot_mat[0] * added_transform->rot_mat[1] + rot_mat[1] * added_transform->rot_mat[3];
            temp_transform.rot_mat[2] = rot_mat[2] * added_transform->rot_mat[0] + rot_mat[3] * added_transform->rot_mat[2];
            temp_transform.rot_mat[3] = rot_mat[2] * added_transform->rot_mat[1] + rot_mat[3] * added_transform->rot_mat[3];
            temp_transform.trans_vec[0] = rot_mat[0] * added_transform->trans_vec[0] + rot_mat[1] * added_transform->trans_vec[1] + trans_vec[0];
            temp_transform.trans_vec[1] = rot_mat[2] * added_transform->trans_vec[0] + rot_mat[3] * added_transform->trans_vec[1] + trans_vec[1];
            this->rot_mat[0] = temp_transform.rot_mat[0];
            this->rot_mat[1] = temp_transform.rot_mat[1];
            this->rot_mat[2] = temp_transform.rot_mat[2];
            this->rot_mat[3] = temp_transform.rot_mat[3];
            this->trans_vec[0] = temp_transform.trans_vec[0];
            this->trans_vec[1] = temp_transform.trans_vec[1];
        }

    Transform2D(double x_trans, double y_trans, double z_rot){
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
    ~Transform2D(){
        //Deconstructor
    }
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

        Correlation(KdTree* model_tree, Set* data_set, int index, Transform2D* transform){
            //Making copies of individual point
            this->current_point.val[0] = data_set->points[index].val[0];
            this->current_point.val[1] = data_set->points[index].val[1];
            //Doing some weird pointer thing. Looks like this is working.
            this->temp_node = model_tree->find_closest_point(&current_point, false, true);
            model_corr_1 = new Point();
            model_corr_2 = new Point();
            this->model_corr_1->val[0] = temp_node[0]->val[0]; //mi1
            this->model_corr_1->val[1] = temp_node[0]->val[1]; //mi1
            this->model_corr_2->val[0] = temp_node[1]->val[0]; //mi2
            this->model_corr_2->val[1] = temp_node[1]->val[1]; //mi2
            
            //Sign of the norm does not matter, using cross product with [0;0;1]
            this->norm_x = model_corr_1->val[1] - model_corr_2->val[1];
            this->norm_y = - (model_corr_1->val[0] - model_corr_2->val[0]);
            this->norm_size = sqrt(norm_x*norm_x + norm_y*norm_y);
            this->norm_x = norm_x/norm_size;
            this->norm_y = norm_y/norm_size;

            //CHECK: I assume this is correct for now, to be investigated further if stuff doesn't work as expected
            this->corrected_value = norm_x*(transform->rot_mat[0]*current_point.val[0] + transform->rot_mat[1]*current_point.val[1] + transform->trans_vec[0] - model_corr_1->val[0]) + \
                                norm_y*(transform->rot_mat[2]*current_point.val[0] + transform->rot_mat[3]*current_point.val[1] + transform->trans_vec[1] - model_corr_1->val[1]);
        }
        Correlation(){
        }

};

double* make_g_vector(Correlation* corrs, Set* data_set, double std_dev, double correntropy_factor, int data_size, bool print){
    double* g = new double[6];
    double delta = 0;
    double g_delta = 0;
    for (int i = 0; i < data_size; i++){
            delta = exp(- corrs[i].corrected_value * corrs[i].corrected_value / (2 * correntropy_factor * std_dev));
            g_delta = (corrs[i].model_corr_1->val[0] * corrs[i].norm_x + corrs[i].model_corr_1->val[1] * corrs[i].norm_y) * delta;
            //g = [nx⋅px⋅(mx⋅nx + my⋅ny)  nx⋅py⋅(mx⋅nx + my⋅ny)  ny⋅px⋅(mx⋅nx + my⋅ny)  ny⋅py⋅(mx⋅nx + my⋅ny)  nx⋅(mx⋅nx + my⋅ny)  ny⋅(mx⋅nx+ my⋅ny)]
            g[0] += corrs[i].norm_x * data_set->points[i].val[0] * g_delta;
            g[1] += corrs[i].norm_x * data_set->points[i].val[1] * g_delta;
            g[2] += corrs[i].norm_y * data_set->points[i].val[0] * g_delta;
            g[3] += corrs[i].norm_y * data_set->points[i].val[1] * g_delta;
            g[4] += corrs[i].norm_x * g_delta;
            g[5] += corrs[i].norm_y * g_delta;
    }
    if (print){
        //Print vector so that it can easily be imported into matlab if need be
        std::cout << "g=[" << g[0] << ',' << g[1] << ',' << g[2] << ',' << g[3] << ',' << g[4] << ',' << g[5] << "]\n\n";
    }
    return g;
}

double* make_G_matrix(Correlation* corrs, Set* data_set, double std_dev, double correntropy_factor, int data_size, bool print){
    double* G = new double[36];
    double delta = 0;
    for (int i = 0; i < data_size; i++){
        delta = exp(- corrs[i].corrected_value * corrs[i].corrected_value / (2 * correntropy_factor * std_dev));

    /*    ⎡     2   2        2                   2                       2                    ⎤
        ⎢ 0 nx ⋅px     1 nx ⋅px⋅py   2 nx⋅ny⋅px   3 nx⋅ny⋅px⋅py    4 nx ⋅px       5 nx⋅ny⋅px⎥
        ⎢                                                                                   ⎥
        ⎢     2              2   2                            2        2                    ⎥
        ⎢ 6 nx ⋅px⋅py    7 nx ⋅py    8 nx⋅ny⋅px⋅py  9 nx⋅ny⋅py     10 nx ⋅py     11 nx⋅ny⋅py⎥
        ⎢                                                                                   ⎥
        ⎢           2                        2   2         2                           2    ⎥
        ⎢12 nx⋅ny⋅px   13 nx⋅ny⋅px⋅py   14 ny ⋅px     15 ny ⋅px⋅py  16 nx⋅ny⋅px   17 ny ⋅px ⎥
        ⎢                                                                                   ⎥
        ⎢                           2        2               2   2                     2    ⎥
        ⎢18 nx⋅ny⋅px⋅py  19 nx⋅ny⋅py    20 ny ⋅px⋅py    21 ny ⋅py    22 nx⋅ny⋅py  23 ny ⋅py ⎥
        ⎢                                                                                   ⎥
        ⎢     2               2                                            2                ⎥
        ⎢24 nx ⋅px       25 nx ⋅py     26 nx⋅ny⋅px    27 nx⋅ny⋅py     28 nx       29 nx⋅ny  ⎥
        ⎢                                                                                   ⎥
        ⎢                                    2              2                           2   ⎥
        ⎣30 nx⋅ny⋅px    31 nx⋅ny⋅py     32 ny ⋅px      33 ny ⋅py     34 nx⋅ny     35 ny     ⎦ */

        G[0] += corrs[i].norm_x * corrs[i].norm_x * data_set->points[i].val[0] * data_set->points[i].val[0] * delta;
        G[7] += corrs[i].norm_x * corrs[i].norm_x * data_set->points[i].val[1] * data_set->points[i].val[1] * delta;
        G[14] += corrs[i].norm_y * corrs[i].norm_y * data_set->points[i].val[0] * data_set->points[i].val[0] * delta;
        G[21] += corrs[i].norm_y * corrs[i].norm_y * data_set->points[i].val[1] * data_set->points[i].val[1] * delta;
        G[28] += corrs[i].norm_x * corrs[i].norm_x * delta;
        G[35] += corrs[i].norm_y * corrs[i].norm_y * delta;

        G[1] += corrs[i].norm_x * corrs[i].norm_x * data_set->points[i].val[0] * data_set->points[i].val[1] * delta; G[6] = G[1];//xx xy
        G[2] += corrs[i].norm_x * corrs[i].norm_y * data_set->points[i].val[0] * data_set->points[i].val[0] * delta; G[12] = G[2];//xy xx
        G[3] += corrs[i].norm_x * corrs[i].norm_y * data_set->points[i].val[0] * data_set->points[i].val[1] * delta; G[18] = G[3];//xy xy
        G[4] += corrs[i].norm_x * corrs[i].norm_x * data_set->points[i].val[0] * delta;                              G[24] = G[4];//xx x
        G[5] += corrs[i].norm_x * corrs[i].norm_y * data_set->points[i].val[0] * delta;                              G[30] = G[5];//xy x

        //G[8]           -- This is equal to G[3] so we don't compute it here --                                      //xy xy
        G[9] += corrs[i].norm_x * corrs[i].norm_y * data_set->points[i].val[1] * data_set->points[i].val[1] * delta; G[19] = G[9];//xy yy
        G[10] += corrs[i].norm_x * corrs[i].norm_x * data_set->points[i].val[1] * delta;                             G[25] = G[10];//xx y
        G[11] += corrs[i].norm_x * corrs[i].norm_y * data_set->points[i].val[1] * delta;                             G[31] = G[11];//xy y

        G[15] += corrs[i].norm_y * corrs[i].norm_y * data_set->points[i].val[0] * data_set->points[i].val[1] * delta; G[20] = G[15];//yy xy
        //G[16]          -- This is equal to G[5] so we don't compute it here --                         //xy x
        G[17] += corrs[i].norm_y * corrs[i].norm_y * data_set->points[i].val[0] * delta;                             G[32] = G[17];//yy x
        
        //G[22]          -- This is equal to G[11] so we don't compute it here --                         //xy y
        G[23] += corrs[i].norm_y * corrs[i].norm_y * data_set->points[i].val[1] * delta;                             G[33] = G[23];//yy y
        G[29] += corrs[i].norm_x * corrs[i].norm_y * delta;                                                          G[34] = G[29];//xy

    }
    //Set values that reappear
    G[8] = G[3]; G[13] = G[3]; G[18] = G[3];
    G[16] = G[5]; G[26] = G[5]; G[30] = G[5];
    G[22] = G[11]; G[27] = G[11]; G[31] = G[11];

    //Print matrix so it can be easily imported into matlab if need be;
    if (print){
        std::cout << std::setprecision(30);
        std::cout << "G=[";
        for (int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                std::cout << G[6*(i) + j];
                if (j+1 < 6){
                    std::cout << ',';
                }
            }
            if (i+1 < 6){
                std::cout << ';';
            }
        }
        std::cout << "]\n\n";
    }

    return G;
}

double* solve_system(double* g, double* G, bool print){
    double* x = new double[6];
    double mult = 0;
    //This will make the G matrix triangular using Gaussian Elimination
    //TODO: swap lines in case of bad pivot point
    for (int j = 0; j < 5; j++){
        for (int i = j+1; i < 6; i++){
            mult = G[6*i+j]/G[6*j+j];
            for (int k = 0; k < 6; k++){
                G[6*i+k] = G[6*i+k] - mult * G[6*j+k];
            }
            g[i] = g[i] - (mult * g[j]);
        }
    }
    //Solve the now triangular system
    double sum = 0;
    x[5] = g[5]/G[35];
    for (int i = 4; i >= 0; i--){
        sum = 0;
        for (int j = 5; j >= i+1; j--){
            sum += G[6*i+j] * x[j];
        }
        x[i] = (g[i] - sum)/G[6*i+i];
    }
    //Print the solution vector to double check if need be
    if (print){
        std::cout << "Solution vector is: ";
        for (int i = 0; i < 6; i++){
            std::cout << x[i] << '\t';
        }
        std::cout << "\n\n";
    }
    return x;
}

class LaserData{
    
    public:
        std::vector<double> las_vec;
        int use_ranges_size;
        double angle_incr;

        Set map_scan_points(Transform2D* transform){
        double* las_x = new double[use_ranges_size];
        double* las_y = new double[use_ranges_size];
        int j = 0;
        int las_size = las_vec.size();
        for (int i = 0; i < las_size; i++){
            if (las_vec[i] < 5){
            double temp_scan_x = las_vec[i] * cos(angle_incr * i);
            double temp_scan_y = las_vec[i] * sin(angle_incr * i);
            las_x[j] = transform->trans_vec[0] + temp_scan_x * cos(transform->z_rot) - temp_scan_y * sin(transform->z_rot);
            las_y[j] = transform->trans_vec[1] + temp_scan_y * cos(transform->z_rot) + temp_scan_x * sin(transform->z_rot);
            j++;
            }
        }
        double* las_ar[2] = {las_x, las_y};
        Set las_set(use_ranges_size, las_ar);
        return las_set;
        }
};

struct DATA{
    Set* old_set;
    int old_data_size;
    Transform2D* main_transf;
    geometry_msgs::Pose msg;
};

DATA perm_data;

int func(LaserData* las_data){
    int guesses = 0;
    double transf_diff = INFINITY;
    int data_size = las_data->use_ranges_size;

    Transform2D base_transform(0, 0, 0);

    if (params.init){
        Set new_set = las_data->map_scan_points(&base_transform);
        perm_data.old_set = new_set.copy_set();
        //plot.add_data(perm_data.old_set);
        params.init = false;
        //plot.plot_data();
        return 0;
    }
    Set new_set = las_data->map_scan_points(&base_transform); //TODO: consider odometry, also maybe this is jank

    KdTree model_tree(perm_data.old_set);

    perm_data.old_set = new_set.copy_set();
 
    Transform2D guess_transf(0,0,0);
    Transform2D total_transf(0,0,0);
    

    while (transf_diff > params.transf_tresh && guesses < params.max_guesses){
        //Apply guess transform on the set and add it to the total transform
        //The guess transform is with respect to the previous guess
        guess_transf.transform_set(&new_set);
        total_transf.add_transform(&guess_transf);

        Correlation* corrs = new Correlation[data_size];
        //Get correlation between points and average for the standard deviation
        double corr_mean = 0;
        //We need to go through the old data as it is the one being compared to the kd tree!
        for (int i = 0; i < data_size; i++){
            corrs[i] = Correlation(&model_tree, &new_set, i, &guess_transf);
            corr_mean += corrs[i].corrected_value;
        }
        corr_mean = corr_mean/data_size;
        /* if (perm_data.main_transf->trans_vec[0] > 1){
            std::cout << "adding more\n";
            plot.add_data(&new_set);
            plot.plot_data();
            return 0;
        } */
    
        //Calculate the standard deviation of the corrected value estimate
        double std_dev = 0;
        for (int i = 0; i < data_size; i++){
            std_dev += abs(corrs[i].corrected_value - corr_mean);
        }
        std_dev = sqrt(std_dev/data_size);

        //For two dimensions, g is a 1x6 vector and G is a 6x6 matrix
        double* g = make_g_vector(corrs, &new_set, std_dev, params.corr_factor, data_size, false);
        double* G = make_G_matrix(corrs, &new_set, std_dev, params.corr_factor, data_size, false);
        //Solve the system using gaussian elimination
        //TODO: skip or/and error message on receiving unsolvable system
        double* x = solve_system(g, G, false);
        if (x[0] != x[0]){
            std::cout << "System has no solution! Assuming the robot is stationary..\n";
            x[0] = 1; x[1] = 0; x[2] = 0; x[3] = 1; x[4] = 0; x[5] = 0;
        }

        //Evaluate difference and update guess Transform2D and update
        transf_diff = guess_transf.compare_transform(x);
        guess_transf.update_transform(x);

        guesses++;
    }
    
    //guess_transf.add_transform(perm_data.trans_tranf);
    guess_transf.transform_set(&new_set);
    total_transf.add_transform(&guess_transf);
    
    perm_data.main_transf->add_transform(&total_transf);
    
    data_size = data_size;
    /* std::cout << "New transform:\n";
    total_transf.print_transform();
    std::cout << "Total transform is now:\n";
    perm_data.main_transf->print_transform(); */

    //Go from Affine to Rotation matrix and publish the rigid transform
    double angle = std::atan2(perm_data.main_transf->rot_mat[2],perm_data.main_transf->rot_mat[0]);

    Transform2D rigid_transf(perm_data.main_transf->trans_vec[0], perm_data.main_transf->trans_vec[1], angle);

    perm_data.msg.position.x = rigid_transf.trans_vec[0];
    perm_data.msg.position.y = rigid_transf.trans_vec[1];
    perm_data.msg.position.z = 0;

    perm_data.msg.orientation.w = 0.5 * sqrt(2 + rigid_transf.rot_mat[0] + rigid_transf.rot_mat[3]);
    perm_data.msg.orientation.x = 0;
    perm_data.msg.orientation.y = 0;
    perm_data.msg.orientation.z = 1/(4.0 * perm_data.msg.orientation.w) * (rigid_transf.rot_mat[2] - rigid_transf.rot_mat[1]);

    params.pub_ready = true;
    
    return 0;
}

class LaserSubscriber{
  ros::Subscriber sub;

  public:
    static bool ready;
    static LaserData data;

    static void laser_callback(const sensor_msgs::LaserScan::ConstPtr &msg){
        if (ready){
            //std::cout << "receiving data\n";
            data.angle_incr = msg->angle_increment;
            data.las_vec.clear();
            int true_ranges_size = msg->ranges.size();
            int use_ranges_size = 0;

            for (int i = 0; i < true_ranges_size; i++){
                data.las_vec.push_back(msg->ranges[i]);
                /* if (i < true_ranges_size-1){
                    if (true){
                        std::cout << msg->ranges[i] << ", ";
                    }
                }
                else{
                    if (true){
                        std::cout << msg->ranges[i] << "\n";
                    }
                } */
                if (msg->ranges[i] < INFINITY){
                    use_ranges_size++; 
                }
            }
            ready = false;
            data.use_ranges_size = use_ranges_size;
            //std::cout << "doing func\n";
            func(&data);
        }
    }
    LaserSubscriber(ros::NodeHandle n){
      //Constructor
      this->sub = n.subscribe("scan", 1000, laser_callback);
    }
};

bool LaserSubscriber::ready = true;
LaserData LaserSubscriber::data;

int main(int argc, char **argv)
{
  ros::init(argc, argv, "listener");
  ros::NodeHandle n;
  ros::Publisher scan_pos = n.advertise<geometry_msgs::Pose>("scan_pos", 1000);
  ros::Rate rate(1);  
  params.start_x = n.param<double>("starting_x_pos", 0.25);
  params.start_y = n.param<double>("starting_y_pos", 0.25);
  params.start_t = n.param<double>("starting_rad_or", 0);
  params.dimensions = n.param<int>("dimensions", 2);
  params.max_guesses = n.param<int>("max_guesses", 10);
  params.corr_factor = n.param<double>("correntropy_factor", 0.1);
  params.transf_tresh = n.param<double>("transform_match_treshold", 0.0001);

  Transform2D starting_transf(params.start_x, params.start_y, params.start_t);
  perm_data.main_transf = &starting_transf;

    //The following technically implements a tf listener but it's not actually needed...

  /* geometry_msgs::TransformStamped transform_stamped;
  tf2_ros::Buffer tfBuffer;
  tf2_ros::TransformListener tfListener(tfBuffer);

  while(true){
      try{
      transform_stamped = tfBuffer.lookupTransform("base_scan", "base_link",ros::Time(0));
    }
    catch (tf2::TransformException &ex) {
      ROS_WARN("%s",ex.what());
      ros::Duration(1.0).sleep();
      continue;
    }
    ROS_INFO("Found frame");
    break;
  }

  Transform2D scan_transf(transform_stamped.transform.translation.x, transform_stamped.transform.translation.y, 
    atan2(2*(transform_stamped.transform.rotation.w * transform_stamped.transform.rotation.z + transform_stamped.transform.rotation.x * transform_stamped.transform.rotation.y),
        transform_stamped.transform.rotation.w*transform_stamped.transform.rotation.w + transform_stamped.transform.rotation.x*transform_stamped.transform.rotation.x 
        - transform_stamped.transform.rotation.y*transform_stamped.transform.rotation.y - transform_stamped.transform.rotation.z*transform_stamped.transform.rotation.z));
    
    perm_data.main_transf->add_transform(&scan_transf); */

  LaserSubscriber laser_sub(n);

  while(ros::ok()){
      laser_sub.ready = true;
      ros::spinOnce();
      if (params.pub_ready){
          scan_pos.publish(perm_data.msg);
          params.pub_ready = false;
      }
      rate.sleep();
      }
    return 0;
}