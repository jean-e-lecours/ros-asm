cmake_minimum_required(VERSION 3.0.2)
project(asm)

## Find catkin and any catkin packages
find_package(catkin REQUIRED COMPONENTS 
    roscpp 
    rospy 
    std_msgs 
    tf2_ros 
    geometry_msgs
    message_generation
)

add_message_files(
    FILES
    info.msg
)

generate_messages(
    DEPENDENCIES
    std_msgs
    geometry_msgs 
)

## Declare a catkin package
catkin_package(
    CATKIN_DEPENDS message_runtime
)



include_directories(include include/asm ${catkin_INCLUDE_DIRS})

add_executable(asm_node src/asm_node.cpp)

set(ASM_LIBS
    include/asm/src/psr.cpp
    include/asm/src/asm.cpp
    include/asm/src/kdt.cpp
    include/asm/src/test.cpp)


add_library(asm_lib ${ASM_LIBS})
add_library(gnuplotcpp INTERFACE)

target_link_libraries(asm_node asm_lib ${catkin_LIBRARIES})