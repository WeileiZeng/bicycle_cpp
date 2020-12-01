#include <iostream>
#include <fstream>
#include<cmath>
#include <itpp/itbase.h>
#include <chrono>
#include <NTL/GF2X.h>
#include "func.hpp" 

int main(int argc, char **argv){
    itpp::Parser parser;
    parser.init(argc,argv);
    int w = 3;
    parser.get(w,"w");

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    
    int max = pow(w,2);
    std::ofstream data;
    data.open ("data.txt");

    for(int l = 2*w; l <= 500; l++){ 

        std::cout << "w = " << w << ", l = " << l << "\n";

        itpp::bvec a(l),vtemp(l);
        itpp::GF2mat Aa(l,l), TAa(l,l), Hx(l,2*l);
        int tot = pow(2,l);
        int wt, cof, k, n;
        NTL::GF2X a_pol, g;
        NTL::GF2X l_pol = NTL::GF2X(l,1) + 1;
        long coff;
        

        for(int i = w; i <= tot; i++){

            a = binary_num(l, i);
            wt = weight(a);

            if(wt == w){
                
                for(int j = 0; j < l; j++){
                    cof = a(j);
                    coff = cof;
                    a_pol += NTL::GF2X(j,coff);
                }

                g = GCD(a_pol, l_pol);

                if(g != 1){
        
                    k = 2 * NTL::deg(g);
                    if(k != 0){

                        Aa = circulant_mat(a);
                        TAa = Aa.transpose();

                        for(int m = 0; m < l; m++){
                            vtemp = Aa.get_col(m);
                            Hx.set_col(m, vtemp);
                        }

                        for(int o = l; o < 2*l; o++){
                            vtemp = TAa.get_col(o-l);
                            Hx.set_col(o, vtemp);
                        }
                        std::cout << "n = " << 2*l <<", k = " << k << "\n"; 
                        std::cout << "a = " << a << "\n";
                    
                        //data << "n = " << 2*l <<", k = " << k << "\n";
                        //data << "Hx = Hz = " << Hx << "\n";
                    }
                }

                
                a_pol = 0;
                      
            }
            
        }
        
    }

    data.close();

    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_span = (end - start)/1000;
    std::cout << "Run-time " << time_span.count() << " seconds.\n";
    return 0;
}