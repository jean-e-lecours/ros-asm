#include <iomanip>

#include "../data.hpp"
#include "../asm_extra_ros.hpp"


double Transform2D::compare_transform(double* vector){
    double transform_diff = 0;
    transform_diff += (this->rot_mat[0] - vector[0]) * (this->rot_mat[0] - vector[0]);
    transform_diff += (this->rot_mat[1] - vector[1]) * (this->rot_mat[1] - vector[1]);
    transform_diff += (this->rot_mat[2] - vector[2]) * (this->rot_mat[2] - vector[2]);
    transform_diff += (this->rot_mat[3] - vector[3]) * (this->rot_mat[3] - vector[3]);
    transform_diff += (this->trans_vec[0] - vector[4]) * (this->trans_vec[0] - vector[4]);
    transform_diff += (this->trans_vec[1] - vector[5]) * (this->trans_vec[1] - vector[5]);
    return transform_diff/6;
}
void Transform2D::update_transform(double *vector){
    this->rot_mat[0] = vector[0];
    this->rot_mat[1] = vector[1];
    this->rot_mat[2] = vector[2];
    this->rot_mat[3] = vector[3];
    this->trans_vec[0] = vector[4];
    this->trans_vec[1] = vector[5];
}
void Transform2D::print_transform(){
    std::cout << "⎡" << rot_mat[0] << "\t " << rot_mat[1] << "\t " << trans_vec[0] << "⎤\n" \
              << "|" << rot_mat[2] << "\t " << rot_mat[3] << "\t " << trans_vec[1] << "|\n" \
              << "⎣0\t 0\t 1⎦\n";
}
void Transform2D::transform_set(Set *set){
    for (int i = 0; i < set->data_size; i++){
        //Applies a single matrix multiplication with the provided x, y vector from the set
        double temp_x = rot_mat[0] * set->points[i].val[0] + rot_mat[1] * set->points[i].val[1] + trans_vec[0];
        double temp_y = rot_mat[2] * set->points[i].val[0] + rot_mat[3] * set->points[i].val[1] + trans_vec[1];
        set->points[i].val[0] = temp_x;
        set->points[i].val[1] = temp_y;
    }
}
Point Transform2D::transform_referenced_point(double x, double y){
    double temp_x = rot_mat[0] * x + rot_mat[1] * y + trans_vec[0];
    double temp_y = rot_mat[2] * x + rot_mat[3] * y + trans_vec[1];
    Point point(2);
    point.val[0] = temp_x;
    point.val[1] = temp_y;
    return point;
}
void Transform2D::add_transform(Transform2D *added_transform){
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

Correlation::Correlation(KdTree* model_tree, Set* data_set, int index, Transform2D* transform){
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
Correlation::Correlation(){
}

Set LaserData::map_scan_points(Transform2D *transform){
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

int TextMapData::read_from_file(std::string filename, int size_x, int size_y, double pix_res){
    std::vector<Point> map_points;
    std::ifstream imapf{filename};
    if (!imapf){
        std::cerr << "Could not read map data!\n";
        return 0;
    }
    else{
        std::cout << "Reading map data\n";
    }
    int i = 0;
    while (!done){
        std::string map_read;
        imapf >> map_read;
        int map_length = map_read.length()-1;
        if (map_read[map_length] == ','){
            map_read.erase(map_read.end()-1);
        }
        else{
            done = true;
        }
        int map_val = std::stoi(map_read);
        if (map_val > 0){
            Point curr_point(2);
            curr_point.val[0] = i%size_y * pix_res;
            curr_point.val[1] = trunc(float(i)/size_x) * pix_res;
            map_points.push_back(curr_point);
        }
        i++;
    }
    std::cout << "done reading map\n";
    this->map_set = new Set(map_points);
    imapf.close();
    return 1;
}


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