#include "eq_system.hpp"
#include<iostream>
#include<sstream>
#include<vector>

using namespace std;

matrix::matrix(){
    mtr = nullptr;
    rows = 0;
    columns = 0;
}

matrix::matrix(int rows, int cols): rows(rows), columns(cols){
    if(rows<=0 || cols<=0){
        throw invalid_argument("Invalid rows or columns number");
    }
    else{
        mtr = new float*[rows];
        for(int i = 0; i < rows; i++){
            mtr[i] = new float[cols]();
        }
    }        
}

matrix::matrix(matrix&& other): rows(other.rows), columns(other.columns), mtr(other.mtr){
    other.mtr = nullptr;
    other.rows = 0;
    other.columns = 0;
}

matrix& matrix::operator=(matrix&& other){
    if(this != &other){
        if(mtr){
            for(int i = 0; i < rows; i++){
                delete[] mtr[i];
            }    
            delete[] mtr; 
            mtr = nullptr;               
        }

        mtr = other.mtr;
        rows = other.rows;
        columns = other.columns;

        other.mtr = nullptr;
        other.rows = 0;
        other.columns = 0;
    }

    return *this;
}

matrix::~matrix(){
    if(mtr){
        for(int i = 0; i < rows; i++){
            delete[] mtr[i];
        }  
        delete[] mtr;              
    }
}

matrix matrix::operator+(matrix& in){
    if(rows != in.rows || columns != in.columns){
        throw invalid_argument("Matrices must have the same size");
    }
    else{
        matrix out(rows, columns);
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < columns; j++){
                out.mtr[i][j] = mtr[i][j] + in.mtr[i][j];
            }
        }
        return move(out);        
    }
}

matrix matrix::operator-(matrix& in){
    if(rows != in.rows || columns != in.columns){
        throw invalid_argument("Matrices must have the same size");
    }
    else{
        matrix out(rows, columns);
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < columns; j++){
                out.mtr[i][j] = mtr[i][j] - in.mtr[i][j];
            }
        }
        return move(out);        
    }
}

matrix matrix::operator*(matrix& in){
    if (columns != in.rows){
        throw invalid_argument("Invalid matrices sizes");
    }
    else{
        matrix out(rows, in.columns);
        for(int i = 0; i < out.rows; i++){
            for(int j = 0; j < out.columns; j++){
                for(int k = 0; k < columns; k++){
                    out.mtr[i][j] += mtr[i][k]*in.mtr[k][j];
                }
            }
        }
        return move(out);
    }
}

void matrix::resize(int new_rows, int new_cols){
    if(mtr){
        for(int i = 0; i < rows; i++){
            delete[] mtr[i];
        }
        delete[] mtr;
    }
    rows = new_rows;
    columns = new_cols;

    mtr = new float*[new_rows];
    for(int i = 0; i < new_rows; i++){
        mtr[i] = new float[new_cols]();
    }
}

void matrix::add_row(string& in, int row){
    stringstream in_stream (in);
    vector<string> el_vec;
    string el;

    while(getline(in_stream, el, ' ')){
        if(!el.empty()){
            el_vec.push_back(el);                
        }
    }

    float* temp = new float[columns];
    for(int i = 0; i < columns; i++){
        temp[i] = stof(el_vec[i]);
    }

    delete[] mtr[row];
    mtr[row] = temp;
}

void matrix::swap_rows(int ind1, int ind2){
    float* temp = mtr[ind1];
    mtr[ind1] = mtr[ind2];
    mtr[ind2] = temp;
}

int matrix::column_abs_max(int col, int from){
    int max_index = from;
    float max_value = abs(mtr[from][col]);
    float temp = 0;

    for(int i = from+1; i < rows; i++){
        temp = abs(mtr[i][col]);
        if(max_value < temp){
            max_index = i;
            max_value = temp;
        }
    }
    return max_index;
}

void matrix::print(){
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < columns; j++){
            cout << mtr[i][j] << " ";
        }
        cout << "\n";
    }
}

matrix eq_system::LUP_decomp(){
    matrix P(1, res_mtr.columns);
    int col_max = 0;
    float temp_int;
    int temp_ext;

    for(int i = 0; i < coef_mtr.columns; i++){
        P.mtr[0][i] = i;
        col_max = coef_mtr.column_abs_max(i, i);
        if(coef_mtr.mtr[col_max][i] == 0){
            throw invalid_argument("No solution");
        }
        else{
            coef_mtr.swap_rows(col_max, i);
            temp_ext = P.mtr[0][i];
            temp_int = res_mtr.mtr[0][i];
            P.mtr[0][i] = P.mtr[0][col_max];
            res_mtr.mtr[0][i] = res_mtr.mtr[0][col_max];
            P.mtr[0][col_max] = temp_ext;
            res_mtr.mtr[0][col_max] = temp_int;
            for(int j = i+1; j < coef_mtr.columns; j++){
                coef_mtr.mtr[j][i] = coef_mtr.mtr[j][i]/coef_mtr.mtr[i][i];
                for(int k = i+1; k < coef_mtr.columns; k++){
                    coef_mtr.mtr[j][k] = coef_mtr.mtr[j][k] - coef_mtr.mtr[j][i]*coef_mtr.mtr[i][k];
                }
            }
        }
    }
    return move(P);
}

matrix eq_system::solve(){
    float sum;
    matrix Y(1, res_mtr.columns);
    matrix X(1, res_mtr.columns);

    Y.mtr[0][0] = res_mtr.mtr[0][0]; 
    for(int i = 1; i < res_mtr.columns; i++){
        sum = 0;
        for(int j = 0; j < i; j++){
            sum += coef_mtr.mtr[i][j]*Y.mtr[0][j];
        }
        Y.mtr[0][i] = res_mtr.mtr[0][i] - sum;
    }

    X.mtr[0][res_mtr.columns-1] = Y.mtr[0][res_mtr.columns-1];
    for(int i = res_mtr.columns-1; i >= 0; i--){
        sum = 0;
        for(int j = i+1; j < res_mtr.columns; j++){
            sum += coef_mtr.mtr[i][j]*X.mtr[0][j];
        }
        X.mtr[0][i] = (Y.mtr[0][i] - sum)/coef_mtr.mtr[i][i];
    }

    return move(X);
}

void eq_system::print(){
    cout << "Coeficients: \n";
    coef_mtr.print();
    cout << "\nConstants: \n";
    res_mtr.print();
    cout << "\n";
}

eq_system::eq_system(int eq_num){
    coef_mtr.resize(eq_num, eq_num);
    res_mtr.resize(1, eq_num);
}