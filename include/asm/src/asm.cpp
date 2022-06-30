#include "../asm.hpp"
#include <iomanip>
#include <iostream>

std::vector<double> make_g_vector(std::vector<Correlation>& corrs, double stdev, double correntropy_factor, bool print){
    std::vector<double> g = {0,0,0,0,0,0};
    double delta = 0;
    double g_delta = 0;
    int data_size = corrs.size();
    for (int i = 0; i < data_size; i++){
            delta = exp(- corrs[i].corrected_value * corrs[i].corrected_value / (2 * correntropy_factor * stdev));
            g_delta = (corrs[i].mcor1.x * corrs[i].norm[0] + corrs[i].mcor1.y * corrs[i].norm[1]) * delta;
            //g = [nx⋅px⋅(mx⋅nx + my⋅ny)  nx⋅py⋅(mx⋅nx + my⋅ny)  ny⋅px⋅(mx⋅nx + my⋅ny)  ny⋅py⋅(mx⋅nx + my⋅ny)  nx⋅(mx⋅nx + my⋅ny)  ny⋅(mx⋅nx+ my⋅ny)]
            g[0] += corrs[i].norm[0] * corrs[i].scan.x * g_delta;
            g[1] += corrs[i].norm[0] * corrs[i].scan.y * g_delta;
            g[2] += corrs[i].norm[1] * corrs[i].scan.x * g_delta;
            g[3] += corrs[i].norm[1] * corrs[i].scan.y * g_delta;
            g[4] += corrs[i].norm[0] * g_delta;
            g[5] += corrs[i].norm[1] * g_delta;
    }
    if (print){
        //Print vector so that it can easily be imported into matlab if need be
        std::cout << "g=[" << g[0] << ',' << g[1] << ',' << g[2] << ',' << g[3] << ',' << g[4] << ',' << g[5] << "]\n\n";
    }
    return g;
}
std::vector<double> make_G_matrix(std::vector<Correlation>& corrs, double stdev, double correntropy_factor, bool print){
    std::vector<double> G = {0,0,0,0,0,0,
                             0,0,0,0,0,0,
                             0,0,0,0,0,0,
                             0,0,0,0,0,0,
                             0,0,0,0,0,0,
                             0,0,0,0,0,0};
    double delta = 0;
    int data_size = corrs.size();
    for (int i = 0; i < data_size; i++){
        delta = exp(- corrs[i].corrected_value * corrs[i].corrected_value / (2 * correntropy_factor * stdev));

  /*    ⎡     2   2        2                   2                       2                ⎤
        ⎢ 0 nx ⋅px     1 nx ⋅px⋅py   2 nx⋅ny⋅px   3 nx⋅ny⋅px⋅py    4 nx ⋅px       5 nx⋅ny⋅px⎥
        ⎢                                                                               ⎥
        ⎢     2              2   2                            2        2                ⎥
        ⎢ 6 nx ⋅px⋅py    7 nx ⋅py    8 nx⋅ny⋅px⋅py  9 nx⋅ny⋅py     10 nx ⋅py     11 nx⋅ny⋅py⎥
        ⎢                                                                               ⎥
        ⎢           2                      2   2         2                             2⎥
        ⎢12 nx⋅ny⋅px   13 nx⋅ny⋅px⋅py   14 ny ⋅px     15 ny ⋅px⋅py  16 nx⋅ny⋅px   17 ny ⋅px ⎥
        ⎢                                                                               ⎥
        ⎢                         2        2              2   2                     2   ⎥
        ⎢18 nx⋅ny⋅px⋅py  19 nx⋅ny⋅py    20 ny ⋅px⋅py    21 ny ⋅py    22 nx⋅ny⋅py  23 ny ⋅py ⎥
        ⎢                                                                               ⎥
        ⎢     2               2                                          2              ⎥
        ⎢24 nx ⋅px       25 nx ⋅py     26 nx⋅ny⋅px    27 nx⋅ny⋅py     28 nx      29 nx⋅ny  ⎥
        ⎢                                                                               ⎥
        ⎢                                    2              2                       2   ⎥
        ⎣30 nx⋅ny⋅px    31 nx⋅ny⋅py     32 ny ⋅px      33 ny ⋅py     34 nx⋅ny     35 ny    ⎦ */

        G[0] += corrs[i].norm[0] * corrs[i].norm[0] * corrs[i].scan.x * corrs[i].scan.x * delta;
        G[7] += corrs[i].norm[0] * corrs[i].norm[0] * corrs[i].scan.y * corrs[i].scan.y * delta;
        G[14] += corrs[i].norm[1] * corrs[i].norm[1] * corrs[i].scan.x * corrs[i].scan.x * delta;
        G[21] += corrs[i].norm[1] * corrs[i].norm[1] * corrs[i].scan.y * corrs[i].scan.y * delta;
        G[28] += corrs[i].norm[0] * corrs[i].norm[0] * delta;
        G[35] += corrs[i].norm[1] * corrs[i].norm[1] * delta;

        G[1] += corrs[i].norm[0] * corrs[i].norm[0] * corrs[i].scan.x * corrs[i].scan.y * delta; G[6] = G[1];//xx xy
        G[2] += corrs[i].norm[0] * corrs[i].norm[1] * corrs[i].scan.x * corrs[i].scan.x * delta; G[12] = G[2];//xy xx
        G[3] += corrs[i].norm[0] * corrs[i].norm[1] * corrs[i].scan.x * corrs[i].scan.y * delta; G[18] = G[3];//xy xy
        G[4] += corrs[i].norm[0] * corrs[i].norm[0] * corrs[i].scan.x * delta;                              G[24] = G[4];//xx x
        G[5] += corrs[i].norm[0] * corrs[i].norm[1] * corrs[i].scan.x * delta;                              G[30] = G[5];//xy x

        //G[8]           -- This is equal to G[3] so we don't compute it here --                                      //xy xy
        G[9] += corrs[i].norm[0] * corrs[i].norm[1] * corrs[i].scan.y * corrs[i].scan.y * delta; G[19] = G[9];//xy yy
        G[10] += corrs[i].norm[0] * corrs[i].norm[0] * corrs[i].scan.y * delta;                             G[25] = G[10];//xx y
        G[11] += corrs[i].norm[0] * corrs[i].norm[1] * corrs[i].scan.y * delta;                             G[31] = G[11];//xy y

        G[15] += corrs[i].norm[1] * corrs[i].norm[1] * corrs[i].scan.x * corrs[i].scan.y * delta; G[20] = G[15];//yy xy
        //G[16]          -- This is equal to G[5] so we don't compute it here --                         //xy x
        G[17] += corrs[i].norm[1] * corrs[i].norm[1] * corrs[i].scan.x * delta;                             G[32] = G[17];//yy x
        
        //G[22]          -- This is equal to G[11] so we don't compute it here --                         //xy y
        G[23] += corrs[i].norm[1] * corrs[i].norm[1] * corrs[i].scan.y * delta;                             G[33] = G[23];//yy y
        G[29] += corrs[i].norm[0] * corrs[i].norm[1] * delta;                                             G[34] = G[29];//xy

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
std::vector<double> solve_system(std::vector<double>& g, std::vector<double>& G, bool print){
    std::vector<double> x = {0,0,0,0,0,0};
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