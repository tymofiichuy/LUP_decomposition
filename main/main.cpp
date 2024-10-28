#include "eq_system.hpp"
#include<iostream>

using namespace std;

int main(){
    eq_system sys(6);
    matrix solution;

    sys.coef_mtr.add_row(string(""), 0);
    sys.coef_mtr.add_row(string(""), 1);
    sys.coef_mtr.add_row(string(""), 2);
    sys.coef_mtr.add_row(string(""), 3);
    sys.coef_mtr.add_row(string(""), 4);
    sys.coef_mtr.add_row(string(""), 5);

    sys.res_mtr.add_row(string(""), 0);
    
    sys.print();
    try{
        sys.LUP_decomp();
        solution = sys.solve();
        solution.print();
    }
    catch(invalid_argument error){
        cout << error.what();
    }
    
    return 0;
}