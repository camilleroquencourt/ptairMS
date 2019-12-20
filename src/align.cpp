// vi: fdm=marker ts=4 et cc=80 tw=80
#include <Rcpp.h>

// Sign function
int sign_d(double x) {
    return((x>0)-(x<0));
}

// [[Rcpp::export]]
Rcpp::List findLimDensity(Rcpp::NumericVector dens, int istart, int state) {

    --istart; // convert from R index to C index
    int difference;
    int linflex = istart;
    int rinflex = istart;

    // Loop from istart to end of vector
    for(int i = istart ; i < dens.size() - 1 ; ++i)
    {
        difference = sign_d(dens[i+1]-dens[i])-sign_d(dens[i]-dens[i-1]);
        if(difference>=1)  //inflex point
        {
            if(state==2)
            {
                state=1;
                rinflex=i;
                istart=i;
                break;
            }
            else//First encounter with an inflex
            {
                state=1;
                linflex=i;
            }
        }
        else if(difference == -2 && state == 1)
            state=2;
    }

    Rcpp::List ret = Rcpp::List::create(Rcpp::Named("linflex") = linflex + 1,
                                  Rcpp::Named("rinflex") = rinflex + 1,
                                  Rcpp::Named("state") = state);
    return(ret);
}
