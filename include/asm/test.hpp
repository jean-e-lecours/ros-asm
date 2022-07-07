#ifndef TEST_HPP
#define TEST_HPP

#include <cmath>
#include <string>
#include <vector>

#include "kdt.hpp"
#include "psr.hpp"

class Plot2D{
    bool data_permanence = false;
    double padding = 1;
    double extremas[4] = {INFINITY, -INFINITY, INFINITY, -INFINITY};

    int number_of_plots;

    std::string gnup_line = "plot ";
    std::string file_header = "";
    std::vector<std::string> filenames;

    public:
        void add_data(std::vector<Point2D>& data_points);

        void add_vec_data_2d(std::vector<std::vector<double>>& data_points);

        void add_line(Point2D& point1, Point2D& point2);

        void add_corrs(std::vector<Correlation>& corr, char corr_type);
        
        void plot_data();

        Plot2D(bool keep_data, std::string file_header);
    
    ~Plot2D();
};



class TextData{
    
    bool done = false;

    public:
        std::vector<double> dat_vec;

        int read_from_file(std::string filename);

        void make_las_data();

        void make_map_data();
};

#endif