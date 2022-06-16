#include <boost/function/function_base.hpp>
#include <boost/mpl/count_fwd.hpp>
#include <iostream>
#include <memory>
#include <ros/node_handle.h>
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

struct PARAMS{
    int max_guesses;
    double corr_factor;
    double transf_tresh;
    bool pub_ready = false;
    double start_x;
    double start_y;
    double start_t;
};

double map_res = 0;

PARAMS params;

std::vector<Point2D> map_set;

KdTree* map_tree;

Transform2D scan_transf(0,0,0);
Transform2D scan_transf_inv(0,0,0);
Transform2D odom_transf(0,0,0);
Transform2D main_transf(0,0,0);

geometry_msgs::PoseStamped msg;

bool nope = false;

int func(std::vector<Point2D> las_data){
    if (nope){
        return 0;
    }
    std::cout << "Doing func\n";
    int guesses = 0;
    double transf_diff = INFINITY;

    Transform2D base_transf(0,0,0);

    /* if (!nope && perm_data.odom_transf->trans_vec[1] > 1){

            std::cout << perm_data.odom_transf->z_rot << '\n';
            plot.add_data(perm_data.map_set);
            plot.add_data(&data_set);
            plot.plot_data();
            nope = true;
            return 0;
        } */
 
    Transform2D guess_transf(0,0,0);
    Transform2D total_transf(0,0,0);

    std::vector<Correlation> correlations;
    
    while (transf_diff > params.transf_tresh && guesses < params.max_guesses){
        //Apply guess transform on the set and add it to the total transform
        //The guess transform is with respect to the previous guess
        guess_transf.transform(las_data);
        total_transf.add_transform(guess_transf);
        
        correlations.clear();
        double corr_mean = 0;
        double corr_stdev = 0;

        for (int i = 0; i < las_data.size(); i++){
            Correlation corr(*map_tree, las_data[i], THOROUGH, guess_transf);
            corr_mean += corr.corrected_value;
            correlations.push_back(corr);
        }
        
        //Get correlation between points and average for the standard deviation

        corr_mean = corr_mean/correlations.size();
        /* if (perm_data.main_transf->trans_vec[0] > 1){
            std::cout << "adding more\n";
            plot.add_data(perm_data.map_set);
            plot.add_data(&data_set);
            plot.plot_data();
            return 0;
        } */
    
        //Calculate the standard deviation of the corrected value estimate
        corr_mean = corr_mean/correlations.size();

        for (int i = 0; i < correlations.size(); i++){
            corr_stdev += std::abs(correlations[i].corrected_value - corr_mean);
        }
        corr_stdev = sqrt(corr_stdev/correlations.size());
        std::cout << corr_stdev << '\n';

        std::vector<double> g = make_g_vector(correlations, las_data, corr_stdev, params.corr_factor, false);
        std::vector<double> G = make_G_matrix(correlations, las_data, corr_stdev, params.corr_factor, false);
        std::vector<double> x = solve_system(g, G, false);
        if (x[0] != x[0]){
            std::cout << "System has no solution! Assuming the robot is stationary..\n";
            x[0] = 1; x[1] = 0; x[2] = 0; x[3] = 1; x[4] = 0; x[5] = 0;
        }

        //Evaluate difference and update guess Transform2D and update
        double angle = std::atan2(x[2],x[0]);
        x[0] = cos(angle);
        x[1] = -sin(angle);
        x[2] = sin(angle);
        x[3] = cos(angle);
        transf_diff = guess_transf.compare_transform(x);
        guess_transf.update_transform(x);

        guesses++;
    }
    std::vector<double> sig_vec = {0,0};

    for (int i = 0; i < correlations.size(); i++){
        double dist_x = correlations[i].mcor1.x - correlations[i].scan.x;
        double dist_y = correlations[i].mcor1.y - correlations[i].scan.y;
        if (sqrt(dist_x*dist_x+dist_y*dist_y) > 1.5*map_res){
            sig_vec[0] += dist_x;
            sig_vec[1] += dist_y;
        }
    }
    sig_vec[0] = sig_vec[0]/correlations.size();
    sig_vec[1] = sig_vec[1]/correlations.size();
    if (sig_vec[0] > map_res && sig_vec[1] > map_res){
        Transform2D sig_correct(sig_vec[0],sig_vec[1],0);
        guess_transf.add_transform(sig_correct);
    }
    else if (sig_vec[0] > map_res){
        Transform2D sig_correct(sig_vec[0],0,0);
        guess_transf.add_transform(sig_correct);
    }
    else if (sig_vec[1] > map_res) {
        Transform2D sig_correct(0,sig_vec[1],0);
        guess_transf.add_transform(sig_correct);
    }

    //guess_transf.add_transform(perm_data.trans_tranf);
    guess_transf.transform(las_data);
    total_transf.add_transform(guess_transf);

    

    //NOTE: Doing this now may cause jank
    total_transf.add_transform(odom_transf);
    total_transf.add_transform(scan_transf_inv);
    main_transf = total_transf;

    std::cout << "Total transform is now:\n";
    main_transf.print_transform();

    //Go from Affine to Rotation matrix and publish the rigid transform
    double angle = std::atan2(main_transf.rot_mat[2],main_transf.rot_mat[0]);

    Transform2D rigid_transf(main_transf.trans_vec[0], main_transf.trans_vec[1], angle);

    msg.pose.position.x = main_transf.trans_vec[0];
    msg.pose.position.y = main_transf.trans_vec[1];
    msg.pose.position.z = 0;

    msg.pose.orientation.w = 0.5 * sqrt(2 + main_transf.rot_mat[0] + main_transf.rot_mat[3]);
    msg.pose.orientation.x = 0;
    msg.pose.orientation.y = 0;
    msg.pose.orientation.z = 1/(4.0 * msg.pose.orientation.w) * (main_transf.rot_mat[2] - main_transf.rot_mat[1]);

    params.pub_ready = true;
    
    return 0;
}

class LaserSubscriber{
  ros::Subscriber sub;

  public:
    static bool ready;

    static void laser_callback(const sensor_msgs::LaserScan::ConstPtr &msg){
        if (ready){
            //std::cout << "receiving data\n";
            std::vector<double> las_vec;

            for (int i = 0; i < msg->ranges.size(); i++){
                las_vec.push_back(msg->ranges[i]);
            }
            ready = false;
            //std::cout << "doing func\n";
            //CHECK: Having to use a dereferencing * here is kinda sus...
            std::vector<Point2D> read_vec = map_scan_points(odom_transf, las_vec, msg->angle_increment);
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
                        curr_point.x = i%msg->info.height * msg->info.resolution + msg->info.origin.position.x;
                        curr_point.y = trunc(float(i)/msg->info.width) * msg->info.resolution + msg->info.origin.position.y;
                        map_points.push_back(curr_point);
                    }
                }
                map_set = map_points;
                got_data = true;
                map_tree = new KdTree(map_set);
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
            //std::cout << "Got some odom\n";
            double odom_angle = atan2(2*(msg->pose.pose.orientation.w * msg->pose.pose.orientation.z + msg->pose.pose.orientation.x * msg->pose.pose.orientation.y),
                            msg->pose.pose.orientation.w*msg->pose.pose.orientation.w + msg->pose.pose.orientation.x*msg->pose.pose.orientation.x 
                            - msg->pose.pose.orientation.y*msg->pose.pose.orientation.y - msg->pose.pose.orientation.z*msg->pose.pose.orientation.z);
            Transform2D temp_transf(msg->pose.pose.position.x, msg->pose.pose.position.y, odom_angle);
            //std::cout << odom_angle << '\n';
            temp_transf.add_transform(scan_transf);
            //temp_transf.print_transform();
            odom_transf = temp_transf;
            //std::cout << "Set some odom\n";
        }
        OdomSubscriber(ros::NodeHandle n){
            //Constructor
            this->sub = n.subscribe("odom", 1000, odom_callback);
        }
};

bool LaserSubscriber::ready = true;
bool MapSubscriber::got_data = false;
std::vector<Point2D> MapSubscriber::map_points;

int main(int argc, char **argv){
    std::cout << "asm_node started, reading parameters...";
    ros::init(argc, argv, "listener");
    ros::NodeHandle n;
    ros::Publisher scan_pos = n.advertise<geometry_msgs::PoseStamped>("asm_pos", 1000);
    ros::Rate rate(1);  
    params.start_x = n.param<double>("starting_x_pos", 0.25);
    params.start_y = n.param<double>("starting_y_pos", 0.25);
    params.start_t = n.param<double>("starting_rad_or", 0);
    params.max_guesses = n.param<int>("max_guesses", 10);
    params.corr_factor = n.param<double>("correntropy_factor", 0.1);
    params.transf_tresh = n.param<double>("transform_match_treshold", 0.0001);

    //std::cout << "done It dies here...";
    Transform2D starting_transf(params.start_x, params.start_y, params.start_t);
    
    main_transf = starting_transf;
    //std::cout << "Actually no\n";

    //The following technically implements a tf listener but it's not actually needed...

    geometry_msgs::TransformStamped transform_stamped;
    tf2_ros::Buffer tfBuffer;
    tf2_ros::TransformListener tfListener(tfBuffer);

    std::cout << "Looking for transforms...\n";

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
    double scan_t_x = transform_stamped.transform.translation.x;
    double scan_t_y = transform_stamped.transform.translation.y;
    double scan_t_a = atan2(2*(transform_stamped.transform.rotation.w * transform_stamped.transform.rotation.z + transform_stamped.transform.rotation.x * transform_stamped.transform.rotation.y),
                            transform_stamped.transform.rotation.w*transform_stamped.transform.rotation.w + transform_stamped.transform.rotation.x*transform_stamped.transform.rotation.x 
                            - transform_stamped.transform.rotation.y*transform_stamped.transform.rotation.y - transform_stamped.transform.rotation.z*transform_stamped.transform.rotation.z);
    scan_transf = Transform2D(scan_t_x, scan_t_y, scan_t_a);
    scan_transf_inv = Transform2D(-scan_t_x, -scan_t_y, -scan_t_a);
    scan_transf.print_transform();
    std::cout << "done.\n Starting subscribers";

    LaserSubscriber laser_sub(n);
    std::cout << "made laser_sub\n";
    MapSubscriber map_sub(n);
    std::cout << "made map_sub\n";
    OdomSubscriber odom_sub(n);
    std::cout << "made odom_sub\n";
    while(ros::ok()){
        laser_sub.ready = true;
        ros::spinOnce();
        if (params.pub_ready){
            msg.header.frame_id = "base_footprint";
            scan_pos.publish(msg);
            params.pub_ready = false;
        }
        rate.sleep();
    }
    return 0;
}