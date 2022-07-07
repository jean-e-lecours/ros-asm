#include <bits/types/FILE.h>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "../../gnuplot.hpp"
#include "../psr.hpp"
#include "../test.hpp"

void Plot2D::add_data(std::vector<Point2D>& data_points){
    //Add file name to the array and create the file itself
    std::string filename = this->file_header + std::to_string(number_of_plots)+ ".dat";
    filenames.push_back(filename);

    std::ofstream osetf(filename);
    //Go through all the data to find what is the maximum and minimum in x
    //and y so that we can size the plot window properly. This is not
    //necessary but we are going through all the data anyways so why not.
    int data_size = data_points.size();
    for (int i = 0; i < data_size; i++){ 
        if (data_points[i].x < extremas[0]){
            extremas[0] = data_points[i].x; 
        }
        else if (data_points[i].x > extremas[1]){
            extremas[1] = data_points[i].x; 
        }
        if (data_points[i].y < extremas[2]){
            extremas[2] = data_points[i].y; 
        }
        else if (data_points[i].y > extremas[3]){
            extremas[3] = data_points[i].y;
        }
        //This is where we actually place the data into the file
        osetf << data_points[i].x << ' ' << data_points[i].y << '\n';
    }
    number_of_plots++;
    gnup_line += "\"" + filename + "\",";
    osetf.close();
}

void Plot2D::add_vec_data_2d(std::vector<std::vector<double>>& data_points){
    //Add file name to the array and create the file itself
    std::string filename = this->file_header + std::to_string(number_of_plots)+ ".dat";
    filenames.push_back(filename);

    std::ofstream osetf(filename);
    //Go through all the data to find what is the maximum and minimum in x
    //and y so that we can size the plot window properly. This is not
    //necessary but we are going through all the data anyways so why not.
    int data_size = data_points.size();
    for (int i = 0; i < data_size; i++){ 
        if (data_points[i][0] < extremas[0]){
            extremas[0] = data_points[i][0]; 
        }
        else if (data_points[i][0] > extremas[1]){
            extremas[1] = data_points[i][0]; 
        }
        if (data_points[i][1] < extremas[2]){
            extremas[2] = data_points[i][1]; 
        }
        else if (data_points[i][1] > extremas[3]){
            extremas[3] = data_points[i][1];
        }
        //This is where we actually place the data into the file
        osetf << data_points[i][0] << ' ' << data_points[i][1] << '\n';
    }
    number_of_plots++;
    gnup_line += "\"" + filename + "\",";
    osetf.close();
}

void Plot2D::add_line(Point2D& point1, Point2D& point2){
    double x1 = point1.x;
    double x2 = point2.x;
    double y1 = point1.y;
    double y2 = point2.y;

    double a = (y2 - y1)/(x2 - x1);
    double b = y1 - a*x1;

    if (a != INFINITY && b != INFINITY && a == a){
        gnup_line += "[" + std::to_string(x1) + ":" + std::to_string(x2) + "] " + std::to_string(a) + "*x+" + std::to_string(b) + " notitle,";
    }
}

void Plot2D::add_corrs(std::vector<Correlation>& corrs, char corr_type){
    for (int i = 0; i < corrs.size(); i++){
        add_line(corrs[i].scan, corrs[i].mcor1);
        if (corr_type == DOUBLE){
            add_line(corrs[i].scan, corrs[i].mcor2);
        }
    }
}

void Plot2D::plot_data(){
    //The plot is drawn when the pipe is destroyed. It is created here so that it
    //goes out of scope right after the function is done and before the files are
    //deleted when the Plot2D object is destroyed.
    GnuplotPipe gp;
            
    //Set graph bounds from the data gathered previously
    std::string gnup_range = "set xrange[" + std::to_string(extremas[0] - padding) + ":" + std::to_string(extremas[1] + padding) + "]\n" + \
                             "set yrange[" + std::to_string(extremas[2] - padding) + ":" + std::to_string(extremas[3] + padding) + "]\n";
    gp.sendLine(gnup_range);
    gnup_line.erase(gnup_line.end() - 1);

    std::ofstream gnup_file{file_header + "gnup_line"};
    gnup_file << gnup_line << '\n';
    gnup_file.close();

    gp.sendLine(gnup_line);
}

Plot2D::Plot2D(bool keep_data, std::string file_header){
    //Constructor
    this->file_header = file_header;
    this->data_permanence = keep_data;
    this->number_of_plots = 0;
}

Plot2D::~Plot2D(){
    //Deconstructor, removes data files if desired
    if (!data_permanence){
        for (int i = 0; i < filenames.size(); i++){
            remove(filenames[i].c_str());
        }
    }
}

int TextData::read_from_file(std::string filename){
    std::ifstream idatf{filename};
    if (!idatf){
        std::cerr << "Could not read laser scan data!\n";
        return 0;
    }
    else{
        std::cout << "Reading data\n";
    }
    while (!done){
        std::string dat_read;
        idatf >> dat_read;
        int dat_length = dat_read.length()-1;
        if (dat_read[dat_length] == ','){
            dat_read.erase(dat_read.end()-1);
        }
        else{
            done = true;
        }
        double dat_val = 0;
        dat_val = std::stod(dat_read);
        dat_vec.push_back(dat_val);                
    }
    std::cout << "done reading\n";
    idatf.close();
    return 1;
}

void make_las_data(){
    
}