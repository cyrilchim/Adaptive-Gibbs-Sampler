// modified code by leemes user at Stack Overflow: http://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf

#ifndef progress_bar_hpp
#define progress_bar_hpp
#define ARMA_DONT_PRINT_ERRORS
// [[Rcpp::depends(RcppArmadillo)]]



#include <RcppArmadillo.h>
#include <cstdlib>


using namespace Rcpp;
using namespace arma;
using namespace std;


int bar_width = 70; // width of the progress bar

void print_bar(int bar_width, int start, int N)// display progress bar
{
    int i;
    double progress_bar = double(start)/double(N); // update progress bar
    Rcout << "[";
    int pos = int(bar_width * progress_bar);
    
    for (i = 0; i < bar_width; ++i) {
        if (i < pos) Rcout << "=";
        else if(i == pos) Rcout << ">";
        else Rcout<<" ";
    }

    Rcout << "] " << int(progress_bar * 100.0) << " %\r";
    Rcout.flush();
}


#endif /* progress_bar_hpp */
