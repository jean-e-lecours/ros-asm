#include <cmath>
#include <string>
#include <fstream>

#include "../data.hpp"
#include "../gnuplot.h"
#include "../cpplot.hpp"

void Plot2D::add_data(Set *set){
    //Add file name to the array and create the file itself
    std::string filename = "dat/temp2d_" + std::to_string(number_of_plots)+ ".dat";
    filenames += filename + "|";

    std::ofstream osetf{filename};
    //Go through all the data to find what is the maximum and minimum in x
    //and y so that we can size the plot window properly. This is not
    //necessary but we are going through all the data anyways so why not.
    for (int i = 0; i < set->data_size; i++){ 
        if (set->points[i].val[0] < extremas[0]){
            extremas[0] = set->points[i].val[0]; 
        }
        else if (set->points[i].val[0] > extremas[1]){
            extremas[1] = set->points[i].val[0];
        }
        if (set->points[i].val[1] < extremas[2]){
            extremas[2] = set->points[i].val[1]; 
        }
        else if (set->points[i].val[1] > extremas[3]){
            extremas[3] = set->points[i].val[1];
        }
        //This is where we actually place the data into the file
        osetf << set->points[i].val[0] << ' ' << set->points[i].val[1] << '\n';
    }
    number_of_plots++;
    gnup_line += "\"" + filename + "\",";
}
void Plot2D::add_line(Point *point1, Point *point2){
    double x1 = point1->val[0];
    double x2 = point2->val[0];
    double y1 = point1->val[1];
    double y2 = point2->val[1];

    double a = (y2 - y1)/(x2 - x1);
    double b = y1 - a*x1;

    if (a != INFINITY && b != INFINITY && a == a){
        gnup_line += "[" + std::to_string(x1) + ":" + std::to_string(x2) + "] " + std::to_string(a) + "*x+" + std::to_string(b) + " notitle,";
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

    std::ofstream gnup_file{"gnup_line"};
    gnup_file << gnup_line << '\n';
    gnup_file.close();

    gp.sendLine(gnup_line);
}
Plot2D::Plot2D(bool keep_data, int max_guesses){
    //Constructor
    this->data_permanence = keep_data;
    this->max_plots = max_guesses;
}
Plot2D::~Plot2D(){
    //Deconstructor, removes data files if desired
    if (!data_permanence){
        int s = 0;
        for (int i = 0; i < max_plots+1; i++){
            std::string filename = "";
            while(filenames[s] != '|'){
                filename += filenames[s];
                s++;
            }
            remove(filename.c_str());
            s++;
        }
    }
}