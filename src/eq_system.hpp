#pragma once
#include<string>

class matrix{
    float** mtr = nullptr;
    int rows;
    int columns;    
public:
    matrix operator+(matrix& in);
    matrix operator-(matrix& in);
    matrix operator*(matrix& in);
    void swap_rows(int ind1, int ind2);
    //next method return is an index, not value
    int column_abs_max(int col, int from);
    void print();
    void resize(int new_rows, int new_cols);
    void add_row(std::string& in, int row);

    matrix& operator=(matrix&& other);
    matrix(int rows, int cols);
    matrix(matrix&& other);
    matrix();    
    ~matrix();

    friend class eq_system;
};

class eq_system{
public:
    matrix coef_mtr;
    matrix res_mtr;
    
    void print();
    matrix LUP_decomp();
    matrix solve();

    eq_system(int eq_num);
};