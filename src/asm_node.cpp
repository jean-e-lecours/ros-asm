#include <iostream>
#include <memory>
#include <ros/node_handle.h>
#include <ros/publisher.h>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>

#include "ros/ros.h"
#include "std_msgs/String.h"
#include "sensor_msgs/LaserScan.h"
#include "geometry_msgs/PoseStamped.h"
#include "geometry_msgs/TransformStamped.h"
#include "tf2_ros/transform_listener.h"
#include "tf2_geometry_msgs/tf2_geometry_msgs.h"
#include "nav_msgs/OccupancyGrid.h"
#include "nav_msgs/Odometry.h"

#include "../include/asm/asm.hpp"
#include "../include/asm/kdt.hpp"
#include "../include/asm/psr.hpp"
#include "../include/asm/test.hpp"

//This structs holds all ROS parameters and is made global to be accessed everywhere in the node
struct PARAMS{
    int max_guesses;
    double corr_factor;
    double transf_tresh;
    bool pub_ready = false;
    double start_x;
    double start_y;
    double start_t;
    bool gen_plots;
};

//Various global quantities

PARAMS params;

double map_res = 0;

std::vector<Point2D> map_set;

KdTree* map_tree;

Transform2D scan_transf(0,0,0);
Transform2D scan_transf_inv(0,0,0);
Transform2D odom_transf(0,0,0);
Transform2D main_transf(0,0,0);
Transform2D las_pos_transf(0,0,0);

geometry_msgs::PoseStamped msg;

bool computing = false;

//This is the main function which applies the ICP
int func(std::vector<Point2D> las_data){
    
    //std::cout << "Doing func\n";
    int guesses = 0;

    Transform2D base_transf(0,0,0);
 
    Transform2D guess_transf(0,0,0);
    Transform2D total_transf(0,0,0);

    //Correlations vector is made global to the function, it is used after the loop
    std::vector<Correlation> correlations;
    
    do{
        correlations.clear();
        double corr_mean = 0;
        double corr_stdev = 0;

        //Get correlation between points and average for the standard deviation
        for (int i = 0; i < las_data.size(); i++){
            //The data set is only modified when calculating this using the current total transform.
            //The output of every step is thus an *incremental* transform that should be summed up
            Correlation corr(*map_tree, las_data[i], THOROUGH, total_transf);
            //Note: DO NOT use anything but thorough search, the algorithm will not converge!
            corr_mean += corr.corrected_value;
            correlations.push_back(corr);
        }
        corr_mean = corr_mean/correlations.size();

        //Calculate the standard deviation of the corrected value estimate
        for (int i = 0; i < correlations.size(); i++){ 
            corr_stdev += std::abs(correlations[i].corrected_value - corr_mean);
        }
        corr_stdev = sqrt(corr_stdev/correlations.size());

        //Solve the system of linear equation using Gaussian
        std::vector<double> g = make_g_vector(correlations, corr_stdev, params.corr_factor, false);
        std::vector<double> G = make_G_matrix(correlations, corr_stdev, params.corr_factor, false);
        std::vector<double> x = solve_system(g, G, false);
        if (x[0] != x[0]){
            //If there is no solution. This should not happen, but in case it does this branch provides added robustness
            std::cout << "System has no solution! Assuming the robot is stationary..\n";
            x[0] = 1; x[1] = 0; x[2] = 0; x[3] = 1; x[4] = 0; x[5] = 0;
        }
        guess_transf.update_transform(x);

        //Extract rotation angle from affine transform to project into SO2. This makes the angle right every time.
        double angle = std::atan2(guess_transf.rot_mat[2],guess_transf.rot_mat[0]);
        Transform2D rigid_guess(guess_transf.trans_vec[0], guess_transf.trans_vec[1], angle);

        //The guess transform is with respect to the previous guess. Must be added to total transform
        guess_transf = rigid_guess;
        total_transf.add_transform(guess_transf);

        guesses++;
    } while (guess_transf.is_significant(params.transf_tresh) && guesses < params.max_guesses);

    //Not much need to find new correlations because the last guess transform is not significant <-- AT LEAST NOT IN THEORY

    //Evaluate the distance between data and map points, remember the ones that are "converged"
    std::vector<int> sig_ind;
    for (int i = 0; i < correlations.size(); i++){

        if (correlations[i].get_distance() < map_res){
            sig_ind.push_back(i);
        }
    }
    //Project affine transform into SO2
    double angle = std::atan2(total_transf.rot_mat[2],total_transf.rot_mat[0]);
    Transform2D rigid_transf(total_transf.trans_vec[0], total_transf.trans_vec[1], angle);

    

    //Make new correlations based on the rigid total transform
    correlations.clear();
    for (int i = 0; i < las_data.size(); i++){
        Correlation corr(*map_tree, las_data[i], THOROUGH, rigid_transf);
        correlations.push_back(corr);
    }

    //If data points that were converged are now away from the corresponding map point, consider moving
    std::vector<double> sig_vec = {0,0};
    int sig_count = 0;
    double dist_score = 0;
    for (int i = 0; i < sig_ind.size(); i++){
        if (correlations[i].get_distance() > map_res){
            //Using an average right now, may need a max value in the future
            sig_vec[0] += correlations[i].get_trans()[0];
            sig_vec[1] += correlations[i].get_trans()[1];
            sig_count++;
        }
    }
    sig_vec[0] = sig_vec[0]/sig_count;
    sig_vec[1] = sig_vec[1]/sig_count;
    std::cout << "Sig vec is: [" << sig_vec[0] << ", " << sig_vec[1] << "]\n";

    //Use obtained sig vec as a maximum value for a bisection search
    double bisec = sqrt(sig_vec[0]*sig_vec[0] + sig_vec[1]*sig_vec[1]);
    std::cout << bisec << '\n';
    Transform2D base_t = rigid_transf;
    Transform2D far_t = rigid_transf;
    std::vector<double> dist_vec = sig_vec;
    far_t.trans_vec[0] -= dist_vec[0];
    far_t.trans_vec[1] -= dist_vec[1];

    //Evaluate both sides
    double dist_score_base = 0;
    double dist_score_far = 0;
    for (int i = 0; i < sig_ind.size(); i++){
        Correlation corr_base(*map_tree, las_data[sig_ind[i]], THOROUGH, base_t);
        dist_score_base += corr_base.get_distance();
        Correlation corr_far(*map_tree, las_data[sig_ind[i]], THOROUGH, far_t);
        dist_score_far += corr_far.get_distance();
    }

    Transform2D mid_t = rigid_transf;
    std::cout << "Entering bisection...\n";
    while(bisec > map_res/2){
        mid_t.trans_vec[0] = (base_t.trans_vec[0] + far_t.trans_vec[0])/2;
        mid_t.trans_vec[1] = (base_t.trans_vec[1] + far_t.trans_vec[1])/2;
        //Evaluate the middle
        double dist_score_mid = 0;
        for (int i = 0; i < sig_ind.size(); i++){
            Correlation corr_mid(*map_tree, las_data[sig_ind[i]], THOROUGH, mid_t);
            dist_score_mid += corr_mid.get_distance();
        }
        if (dist_score_mid < dist_score_base && dist_score_mid > dist_score_far){
            //Here the answer lies between mid and far, re-assign mid to base and try again
            base_t = mid_t;
            bisec = bisec/2;
        }
        else if (dist_score_mid > dist_score_base && dist_score_mid < dist_score_far){
            //Here the answer lies between base and mid, re-assign far to mid and try again
            far_t = mid_t;
            bisec = bisec/2;
        }
        else if ((dist_score_mid < dist_score_base && dist_score_mid < dist_score_far) ||\
         (dist_score_mid > dist_score_base && dist_score_mid > dist_score_far)){
            //There are multiple minima in the range, need to pick the best one
            if (dist_score_base < dist_score_far){
                far_t = mid_t;
                bisec = bisec/2;
            }
            else{
                base_t = mid_t;
                bisec = bisec/2;
            }
            break;
        }
    }

    /* std::vector<std::vector<double>> sig_vec;
    for (int i = 0; i < correlations.size(); i++){
        sig_vec.push_back(correlations[i].get_trans());
    }

    std::vector<double> sig_trans = {0,0};
    double sign_x;
    double sign_y;
    for (int i = 0; i < sig_vec.size(); i++){
        if (std::abs(sig_vec[i][0]) > 3*map_res){
            sig_trans[0] += sig_vec[i][0];
            sign_x ++;
        }
        if (std::abs(sig_vec[i][1]) > 3*map_res){
            sig_trans[1] += sig_vec[i][1];
            sign_y ++;
        }
    }
    if (sign_x > 50){
        sig_trans[0] = sig_trans[0]/sign_x;
        std::cout << "Moving in X, " << sig_trans[0] << "\n"; 
    }
    else{
        sig_trans[0] = 0;
    }
    if (sign_y > 50){
        sig_trans[1] = sig_trans[1]/sign_y;
        std::cout << "Moving in Y, " << sig_trans[1] << "\n";
    }
    else{
        sig_trans[1] = 0;
    }
    
    //Shape the translation correction transform based on the previous evaluation and add it
    Transform2D add_transf(-sig_trans[0], -sig_trans[1], 0); */
   /*  if (params.gen_plots){
        Plot2D plot(false, "data");
        plot.add_data(map_set);
        plot.add_data(las_data);
        total_transf.transform(las_data);
        plot.add_data(las_data);
        add_transf.transform(las_data);
        plot.add_data(las_data);
        plot.plot_data();
    } */
    total_transf = mid_t;
    
    //add_transf.print_transform();

    //Before this total_transf is just adjustment from odometry, odometry is added here to get global pose
    total_transf.add_transform(odom_transf);
    //scan_transf_inv.print_transform();
    main_transf = total_transf;

    std::cout << "Total transform is now:\n";
    ROS_INFO("NEW ASM ESTIMATE AVAILABLE");
    
    main_transf.print_transform();
    //rigid_transf.print_transform();

    msg.pose.position.x = main_transf.trans_vec[0];
    msg.pose.position.y = main_transf.trans_vec[1];
    msg.pose.position.z = 0;

    //Converting z axis rotation to quaternion because this is the standard for ROS
    msg.pose.orientation.w = 0.5 * sqrt(2 + main_transf.rot_mat[0] + main_transf.rot_mat[3]);
    msg.pose.orientation.x = 0;
    msg.pose.orientation.y = 0;
    msg.pose.orientation.z = 1/(4.0 * msg.pose.orientation.w) * (main_transf.rot_mat[2] - main_transf.rot_mat[1]);

    //We have a new estimate with the correct format, time to publish
    params.pub_ready = true;
    
    return 0;
}

class LaserSubscriber{
  ros::Subscriber sub;
    public:
    static void laser_callback(const sensor_msgs::LaserScan::ConstPtr &msg){
        if (!computing){
            //Receive data from the message and run the algorithm, then freeze the subscriber
            //Laser scan subscriber is unfrozen when new estimate is published
            std::vector<double> las_vec;
            for (int i = 0; i < msg->ranges.size(); i++){
                las_vec.push_back(msg->ranges[i]);
            }
            computing = true;
            //TODO: add components of previous guesses to initial position of the data set
            std::vector<Point2D> read_vec = map_scan_points(las_pos_transf, las_vec, msg->angle_increment);
            func(read_vec);
        }
    }
    LaserSubscriber(ros::NodeHandle n){
      //Constructor
      this->sub = n.subscribe("scan", 1000, laser_callback);
    }
};

class MapSubscriber{
    ros::Subscriber sub;
    static std::vector<Point2D> map_points;
    public:
        int number_of_points;
        static bool got_data;
        static void map_callback(const nav_msgs::OccupancyGrid::ConstPtr &msg){
            if (!got_data){
                std::cout << "Listening for map...\n";
                for (int i = 0; i < msg->info.width*msg->info.height; i++){
                    if (msg->data[i] > 0){
                        map_res = msg->info.resolution;
                        Point2D curr_point(0,0);
                        //Convert position in the array to coordinates on a 2D grid
                        curr_point.x = i%msg->info.width * msg->info.resolution + msg->info.origin.position.x;
                        curr_point.y = trunc(float(i)/msg->info.width) * msg->info.resolution + msg->info.origin.position.y;
                        map_points.push_back(curr_point);
                    }
                }
                /* for (int i = 0; i < msg->info.width*msg->info.height; i++){
                    if (msg->data[i] > 0){
                        map_res = msg->info.resolution;
                        Point2D curr_point(0,0);
                        if (i < msg->info.width && msg->data[i+1] <= 0){
                            //Free space detected in the next positive x position
                            curr_point.x = i%msg->info.width * msg->info.resolution + msg->info.origin.position.x + map_res/2;
                            curr_point.y = trunc(float(i)/msg->info.width) * msg->info.resolution + msg->info.origin.position.y;
                        }
                        else if (i > 0 && msg->data[i-1] <= 0){
                            //Free space detected in the previous negative x position
                            curr_point.x = i%msg->info.width * msg->info.resolution + msg->info.origin.position.x;
                            curr_point.y = trunc(float(i)/msg->info.width) * msg->info.resolution;
                        }
                        else if (i > msg->info.width && msg->data[i-msg->info.width]){
                            //Free space detected in the previous negative y direction
                            curr_point.x = i%msg->info.width * msg->info.resolution + msg->info.origin.position.x;
                            curr_point.y = trunc(float(i)/msg->info.width) * msg->info.resolution + msg->info.origin.position.y - map_res/2;
                        }
                        else if (i < msg->info.width * (msg->info.height)-1 && msg->data[i+msg->info.width]){
                            //Free space detected in the next positive y direction
                            curr_point.x = i%msg->info.width * msg->info.resolution + msg->info.origin.position.x;
                            curr_point.y = trunc(float(i)/msg->info.width) * msg->info.resolution + msg->info.origin.position.y + map_res/2;
                        }
                        map_points.push_back(curr_point);
                    }
                } */
                map_set = map_points;
                got_data = true;
                map_tree = new KdTree(map_set);
                Plot2D map_plot(false, "map");
                map_plot.add_data(map_points);
                map_plot.plot_data();
                std::cout << "Made tree\n";
            }
        }
        MapSubscriber(ros::NodeHandle n){
            //Constructor
            this->sub = n.subscribe("map", 1000, map_callback);
        }
};

class OdomSubscriber{
    ros::Subscriber sub;
    public:
        static void odom_callback(const nav_msgs::Odometry::ConstPtr &msg){
            if (!computing){
                //Odom subscriber freezes with laser scan subscriber to prevent time shift
                //This is unfrozen when the new estimate is published

                //Calculate angle from quaternion and put it in a transform
                double odom_angle = atan2(2*(msg->pose.pose.orientation.w * msg->pose.pose.orientation.z + msg->pose.pose.orientation.x * msg->pose.pose.orientation.y),
                                msg->pose.pose.orientation.w*msg->pose.pose.orientation.w + msg->pose.pose.orientation.x*msg->pose.pose.orientation.x 
                                - msg->pose.pose.orientation.y*msg->pose.pose.orientation.y - msg->pose.pose.orientation.z*msg->pose.pose.orientation.z);
                Transform2D temp_transf(msg->pose.pose.position.x, msg->pose.pose.position.y, odom_angle);
                //Get odometry transform considering only base robot frame
                odom_transf = temp_transf;
                //Get odometry transform, adding in the transform from base to scan frame
                temp_transf.add_transform(scan_transf);
                las_pos_transf = temp_transf;
            }
        }
        OdomSubscriber(ros::NodeHandle n){
            //Constructor
            this->sub = n.subscribe("odom", 1000, odom_callback);
        }
};

bool MapSubscriber::got_data = false;
std::vector<Point2D> MapSubscriber::map_points;

int main(int argc, char **argv){
    //Starting set-up
    std::cout << "asm_node started, reading parameters...";
    ros::init(argc, argv, "listener");
    ros::NodeHandle n;
    ros::Publisher scan_pos = n.advertise<geometry_msgs::PoseStamped>("asm_pos", 1000);
    ros::Rate rate(1);

    //Listening to and setting parameters
    //For now, these first three parameters are not really necessary, may come in useful for algorithm fusion
    params.start_x = n.param<double>("starting_x_pos", 0.25);
    params.start_y = n.param<double>("starting_y_pos", 0.25);
    params.start_t = n.param<double>("starting_rad_or", 0);

    //ICP parameters
    params.max_guesses = n.param<int>("max_guesses", 20);
    params.corr_factor = n.param<double>("correntropy_factor", 0.1);
    params.transf_tresh = n.param<double>("transform_match_treshold", 0.0001);

    //Other parameters
    params.gen_plots = n.param<bool>("generate_plots", false);

    //The following technically implements a tf listener but it's not actually needed, the transform holds
    //a static value which may or may not be reliable on a real platform

    geometry_msgs::TransformStamped transform_stamped;
    tf2_ros::Buffer tfBuffer;
    tf2_ros::TransformListener tfListener(tfBuffer);

    std::cout << "Looking for transforms...\n";

    //ROS rate of 1s is on the slow side but accounts for less powerful hardware
    ros::Duration(1.0).sleep();

    while(true){
        try{
            transform_stamped = tfBuffer.lookupTransform("base_footprint", "base_scan",ros::Time(0));
        }
        catch (tf2::TransformException &ex) {
            ROS_WARN("%s",ex.what());
            ros::Duration(1.0).sleep();
            continue;
        }
        ROS_INFO("Found frame");
        break;
    }
    //Converting the tf2 transform into my Transform2D object
    double scan_t_x = transform_stamped.transform.translation.x;
    double scan_t_y = transform_stamped.transform.translation.y;
    double scan_t_a = atan2(2*(transform_stamped.transform.rotation.w * transform_stamped.transform.rotation.z + transform_stamped.transform.rotation.x * transform_stamped.transform.rotation.y),
                            transform_stamped.transform.rotation.w*transform_stamped.transform.rotation.w + transform_stamped.transform.rotation.x*transform_stamped.transform.rotation.x 
                            - transform_stamped.transform.rotation.y*transform_stamped.transform.rotation.y - transform_stamped.transform.rotation.z*transform_stamped.transform.rotation.z);
    scan_transf = Transform2D(scan_t_x, scan_t_y, scan_t_a);
    scan_transf.print_transform();

    std::cout << "Starting subscribers\n";
    //The map subscriber reads the map server and generates a KD-Tree
    MapSubscriber map_sub(n);
    std::cout << "made map_sub\n";
    //The laser scan subscriber starts looking for messages
    LaserSubscriber laser_sub(n);
    std::cout << "made laser_sub\n";
    //The odometry subscriber starts looking for messages
    OdomSubscriber odom_sub(n);
    std::cout << "made odom_sub\n";

    //Main ROS loop
    while(ros::ok()){
        computing = false;
        ros::spinOnce();
        //Runs the algorithm if it is not being ran
        if (params.pub_ready){
            //Applies the PoseStamped message to the map (absolute) frame
            msg.header.frame_id = "map";
            //Publish a message on the asm_pos topic
            scan_pos.publish(msg);
            //Resets the publish flag
            params.pub_ready = false;
        }
        rate.sleep();
    }
    return 0;
}