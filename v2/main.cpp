#include <iostream>
#include <fstream>
#include<cmath>
#include <itpp/itbase.h>
#include <chrono>
#include <NTL/GF2X.h>
#include "func.hpp"
#include "weilei_lib/weilei_lib.h" 

int main(int argc, char **argv){
    itpp::Parser parser;
    parser.init(argc,argv);
    int w = 3;
    parser.get(w,"w");

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    std::ofstream data;
    data.open ("data.txt");
    

    for(int l = 2*w; l <= 10; l++){ 

        std::cout << "w = " << w << ", l = " << l << "\n";

        itpp::bvec a(l);
        int wt, n;

        a.zeros();
        a(0) = 1;

        set_weight(a, 1, w, l);

    }

    data.close();

    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_span = (end - start)/1000;
    std::cout << "Run-time " << time_span.count() << " seconds.\n";
    return 0;
}