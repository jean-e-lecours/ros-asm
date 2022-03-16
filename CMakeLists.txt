cmake_minimum_required(VERSION 3.0.2)
project(asm)

## Find catkin and any catkin packages
find_package(catkin REQUIRED COMPONENTS 
    roscpp 
    rospy 
    std_msgs 
    tf2_ros 
    geometry_msgs
)

## Declare a catkin package
catkin_package()

## Build talker and listener
include_directories(include ${catkin_INCLUDE_DIRS})

add_executable(asm_node src/asm_node.cpp)
target_link_libraries(asm_node ${catkin_LIBRARIES})