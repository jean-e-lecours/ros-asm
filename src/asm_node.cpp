#include <boost/mpl/count_fwd.hpp>
#include <iostream>
#include <ros/node_handle.h>
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

#include "nav_msgs/OccupancyGrid.h"

#include "../include/asm/asm_extra_ros.hpp"
#include "../include/asm/data.hpp"
#include "../include/asm/kdtree.hpp"
#include "../include/asm/mat.hpp"

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

int Point::dims = 2;
int Set::dims = 2;
int KdTree::dims = 2;

PARAMS params;

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

class MapSubscriber{
    ros::Subscriber sub;
    static std::vector<Point> map_points;
    public:
        static Set* map_set;
        int number_of_points;
        static bool got_data;
        static void map_callback(const nav_msgs::OccupancyGrid::ConstPtr &msg){
            if (!got_data){
                std::cout << "Listening for map...\n";
                Point curr_point(2);
                for (int i = 0; i < msg->info.width*msg->info.height; i++){
                    curr_point.val[0] = i%msg->info.height * msg->info.resolution + msg->info.origin.position.x;
                    curr_point.val[1] = trunc(float(i)/msg->info.width) * msg->info.resolution + msg->info.origin.position.y;
                    map_points.push_back(curr_point);
                }
                map_set = new Set(map_points);
                got_data = true;
            }
        }
        MapSubscriber(ros::NodeHandle n){
            //Constructor
            this->sub = n.subscribe("map", 1000, map_callback);
        }

};

bool LaserSubscriber::ready = true;
LaserData LaserSubscriber::data;

int main(int argc, char **argv)
{
  ros::init(argc, argv, "listener");
  ros::NodeHandle n;
  ros::Publisher scan_pos = n.advertise<geometry_msgs::Pose>("scan_pos", 1000);
  ros::Rate rate(5);  
  params.start_x = n.param<double>("starting_x_pos", 0.25);
  params.start_y = n.param<double>("starting_y_pos", 0.25);
  params.start_t = n.param<double>("starting_rad_or", 0);
  params.dimensions = n.param<int>("dimensions", 2);
  Point::dims = params.dimensions;
  Set::dims = params.dimensions;
  KdTree::dims = params.dimensions;
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
    MapSubscriber map_sub(n);


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