#ifndef CPPLOT_H
#define CPPLOT_H

#include <string>
#include <cmath>

#include "data.hpp"

class Plot2D{
    bool data_permanence = false;
    int max_plots = 0;
    int number_of_plots = 0;
    double padding = 1;
    double extremas[4] = {INFINITY, 0, INFINITY, 0};

    std::string gnup_line = "plot ";
    std::string filenames = "";

    public:
        int total_plots;

        void add_data(Set* set);

        void add_line(Point* point1, Point* point2);

        void plot_data();

        Plot2D(bool keep_data, int number_of_plots);

    ~Plot2D();
};

#endif