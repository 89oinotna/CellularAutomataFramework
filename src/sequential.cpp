#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <functional>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstring>
#define cimg_use_png

#include "./cimg/CImg.h"
#include <cstdint>
#include <cassert>
#include "utimer.cpp"

using namespace std;
using namespace cimg_library;

/**
 * Cellular Automata abstract implementation
 * @tparam T state type
 * @tparam C CImg type to represent the image
 */
template <class T, class C>
class CellularAutomata{
        
  
    int _n; //number of rows
    int _m; //number of columns
    vector<vector<T>> matrices; //the two matrices as alternating buffers
    int _nIterations; 
    vector<CImg<C>> images;
    
    
    /**
     * Computes the new state of a cell
     * @param matrix old state of the matrix 1-d
     * @param index index in the matrix 1-d
     * @param n rows number
     * @param m columns number
     * @return the new state
     */
    virtual inline T rule(const vector<T>& matrix, int const&index, int const& n, int const& m)=0;

    /**
     * Computes the representation of the state and inserts it in the image object
     * @param img the image object
     * @param i row index
     * @param j column index
     * @param state value of the cell state
     */
    virtual inline void repr(CImg<C> &img, int const&i, int const &j, T const& state)=0;

    /**
     * Used to initialize the image objects
     * @param n rows
     * @param m columns
     * @return the image object
     */
    virtual inline CImg<C> imgBuilder(int const& n, int const &m)=0;

    public:
    CellularAutomata(vector<T>& initialState, int n, int m, int nIterations) 
        : _n(n), _m(m), _nIterations(nIterations){
        matrices = vector<vector<T>>(2, initialState);
    }

    public:
     void run(){  
        bool index=0; //index used to alternate the matrices
        for(int j=0;j<_nIterations;j++){                    
            for(int k = 0; k < _n*_m; k++){
                #ifdef WIMG
                auto res=rule(matrices[index], k, _n, _m);
                matrices[!index][k]=res; 
                repr(images[j], k/_m, k%_m, res);
                #endif

                #ifndef WIMG
                matrices[!index][k]=rule(matrices[index], k, _n, _m);
                #endif
            }
            #ifdef WIMG
            string filename="./frames/"+to_string(j)+".png";
            char name[filename.size()+1];
            strcpy(name, filename.c_str());
            images[j].save(name);
            #endif
            index=!index;
        }            
    }    

    void init(){
        #ifdef WIMG
        images=vector<CImg<C>>(_nIterations, imgBuilder(_n,_m));
        #endif
    }
};



int random_init(){
    return (std::rand())%2;
}

class MyCa : public CellularAutomata<int, unsigned char> {
    public:
    MyCa(vector<int>& initialState, 
                    int n, int m, 
                    int nIterations)
        :  CellularAutomata(initialState, n, m, nIterations){}

    inline int pmod(int v, int m){
        return v % m < 0 ? v % m + m : v % m;
    }

    int rule(const vector<int>& matrix, int const& index, int const& n, int const& m){
        int row = index/m;
        int col = index%m;
        int rm1 = pmod(row-1, m);
        int cm1 = pmod(col-1, n);
        int rp1 = pmod(row+1, m);
        int cp1 = pmod(col+1, m);
        int sum=matrix[rm1*m + col];    //up
        sum += matrix[rm1*m + cm1];     //up_left   
        sum += matrix[rm1*m + cp1];     //up_right       
        sum += matrix[rp1*m + col];     //down
        sum += matrix[row*m + cm1];     //left
        sum += matrix[row*m + cp1];     //right
        sum += matrix[rp1*m + cm1];     //down_left
        sum += matrix[rp1*m + cp1];     //down_right

        int s=matrix[index];

        if(sum==3){
            return 1;
        }
        if(s==1 && sum==2){
            return s;
        }
        return 0;
    }

    void repr(CImg<unsigned char> &img, int const& i, int const &j, int const& s){
            img(i,j)= s==0 ? 0 : 255;
    }
    
    CImg<unsigned char> imgBuilder(int const& n, int const &m){
        return CImg<unsigned char>(n, m);
    }
};

int main(int argc, char * argv[]) {

    if(argc != 4) {
        std::cout << "Usage is: " << argv[0] << " n m iterations" << std::endl;
        return(-1);
    }
    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    int iter = atoi(argv[3]);
    srand(0);
    vector<int> matrix (n*m);  
    std::generate(matrix.begin(), matrix.end(), random_init); 

    utimer tp("completion time");
    MyCa ca(matrix, n, m, iter);
    ca.init();
    //utimer tp("run time");
    ca.run();
    return(0);
}