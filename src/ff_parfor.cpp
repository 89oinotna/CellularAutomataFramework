#include <ff/ff.hpp>
#include <ff/parallel_for.hpp>
#include <ff/barrier.hpp>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <functional>
#include <vector>
#include <thread>
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
        
    typedef struct {
        int start;
        int end;
    } range;

    ff::Barrier ba; //the FF barrier
    int _n; //number of rows
    int _m; //number of columns
    vector<vector<T>> matrices; //the two matrices as alternating buffers
    int _nIterations;
    int _parallelism;
    int _nworkers;
    vector<range> ranges; //list of ranges to be assigned to each worker
    #ifdef WIMG
    vector<CImg<C>> images;
    #endif
    ff::ParallelFor* pf;

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

    /**
     * Computes and initializes the ranges that will be assigned to workers
     */
    void initRanges(){
        int delta { int(_n)*int(_m) / int(_nworkers) };
        for(int i=0; i<_nworkers; i++) { 
            ranges[i].start = i*delta;
            ranges[i].end   = (i != (_nworkers-1) ? (i+1)*delta : _n*_m); 
        }        
    }   

    public:
    CellularAutomata(vector<T>& initialState, 
                    int n, int m, 
                    int nIterations,  int nworkers){   
        _n=n;
        _m=m;
        matrices = vector<vector<T>>(2, initialState);
        _nIterations=nIterations;
        _nworkers = nworkers;
        ranges= vector<range>(_nworkers);

        ba.barrierSetup(_nworkers);
        
        pf = new ff::ParallelFor(_nworkers);        
        pf->disableScheduler(true);  
    }

    public:
    #ifdef WIMG
    /**
     * Get the iterations that the given worker has to write and writes them
     */
    inline void writeImages(int const& thid){
        int nprint=ceil(double(_nIterations) / double(_nworkers));
        int wstart=thid*nprint;
        int wend= min(int(_nIterations), (int(thid)+1) * nprint);
        for(int k=wstart; k<wend; k++){
            string path="./frames/"+to_string(k)+".png";
            char name[path.size()+1];
            strcpy(name, path.c_str());
            images[k].save(name);
        }
    }
    #endif
    /**
     * Initializes ranges and images
     */
    void init(){
        initRanges();
        #ifdef WIMG
        images = vector<CImg<C>>(_nIterations, imgBuilder(_n,_m));
        #endif
    }
    void run(){ 
        pf->parallel_for_thid(0,ranges.size(),1,0,[&](const long i, const int thid) { 
            bool index=0;
            for(int j=0;j<_nIterations;j++){ 
                for(int k = ranges[i].start; k < ranges[i].end; k++){
                    #ifdef WIMG
                    auto res=rule(matrices[index], k, _n, _m);
                    matrices[!index][k]=res; 
                    repr(images[j], k/_m, k%_m, res);
                    #endif
                    #ifndef WIMG
                    matrices[!index][k]=rule(matrices[index], k, _n, _m);
                    #endif
                }                  
                ba.doBarrier(thid); 
                 
                index=!index; //change the index of the matrix
            }
            #ifdef WIMG
            writeImages(thid);
            #endif
        },_nworkers);
        
    }
};

/**
 * @return random number 0 or 1
 */
int random_init(){
    return (std::rand())%2;
}

class MyCa : public CellularAutomata<int, unsigned char> {
    public:
    MyCa(vector<int>& initialState, 
                    int n, int m, 
                    int nIterations, int nworkers)
        :  CellularAutomata(initialState, n, m, nIterations, nworkers){}

    /**
     * @param v
     * @param m
     * @return positive modulo
     */
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
            img(i,j)= s==0 ? 0 : 255; //black & white
    }
    
    CImg<unsigned char> imgBuilder(int const& n, int const &m){
        return CImg<unsigned char>(n, m);
    }
};

int main(int argc, char* argv[]){
    if(argc != 5) {
        cout << "Usage is: " << argv[0] << " N M number_step number_worker" << endl;
        return(-1);
    }
    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    int iter = atoi(argv[3]);
    int nw = atoi(argv[4]);
    
    std::srand(0);
    vector<int> matrix (n*m);  
    
    std::generate(matrix.begin(), matrix.end(), random_init);

    utimer tp("completion time");
    MyCa ca(matrix,n,m, 
        iter,
        nw
    );
    ca.init();
    //utimer tp("run time");
    ca.run();
    return 0;
}
