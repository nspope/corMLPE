#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <chrono>
#include "corMLPE.h"

using namespace Rcpp;

arma::mat makeBlocks(arma::vec v, arma::vec& x, arma::uword p){
	/*
	v is 3x1 with coefficients
	x is nx1 with values to multiply
	p is integer, gives number of populations
	*/
  arma::uword i, j, s;

  arma::uvec l(p-1);
	for(i=1;i<p;i++){
		l[i-1] = p-i;
	}
	
  arma::mat blocks(p+1, p+1, arma::fill::zeros);
	// first calculate empty blocks
	s=0;
	for(i=0;i<p-1;i++){	
		for(j=0;j<l[i];j++){
			blocks(i,i) += x[s++]*v[0];	
		}
	}
	// multiply all x by v1
	double vd = v[1] - v[0];
	for(i=0;i<x.n_elem;i++){
		x[i] *= vd;
	}
	
	// get full blocks
	s=0;
	for(i=0;i<p-1;i++){
		blocks(i+1,i) = blocks(i,i);
		for(j=0;j<l[i];j++){
			blocks(i+1,i) += x[s++];
		}
	}
	// get diagonal blocks
	s=0;
	for(i=0;i<p-1;i++){ // loop over horizontal
		for(j=i+2;j<p+1;j++){ // loop over vertical
			blocks(j,i) = blocks(i,i) + x[s++];
		}
	}
	return(blocks);
}

arma::vec multBlocks(arma::mat blocks, arma::uword p, arma::vec x){
	int st, s, z, i, k, j;
  arma::vec o(x.n_elem, arma::fill::zeros);
	// upper triangle first
	z=0;
	st=p-1;
	for(i=0;i<p-1;i++){ // horizontal over blocks
		s = st;
		for(k=i+1;k<p-1;k++){ // vertical over blocks
			for(j=k+2;j<p+1;j++){ // vertical within blocks
				o(s) += blocks(j,i) + x(z);
				s++;
			}
			z++;
		}
		z++;
		st += p-i-2;
	}
	// hopefully that does it.

	
	// lower triangle
	// first block
  arma::vec Sm(p-1, arma::fill::zeros);
	for(i=2;i<p+1;i++){ // loop over rows
		for(j=1;j<i;j++){ // loop over columns
			Sm(i-2) += blocks(i,j);
		}
		for(j=i;j<p-1;j++){
			Sm(i-2) += blocks(j,j);
		}
	}

	// add this to running count
	s=0;
	for( i=0;i<p-2;i++ ){ // loop vertically over blocks
		for(j=i;j<p-1;j++){ // loop vertically within blocks
			o[s++] += Sm(j) + blocks(i+1,i);
			// subtract out
			Sm(j) -= blocks(j+2,i+1);
		}
	}
	o[s++] += blocks(p-1,p-2);	


	return(o);
}

// [[Rcpp::export]]
arma::vec MultLambda(arma::vec x, arma::vec v, arma::uword p){
	double vd = v[2]-v[1];
  arma::vec x2 = vd*x;
  arma::mat bl = makeBlocks(v, x, p);
  arma::vec o = multBlocks(bl, p, x);
	o += x2;
	return(o);
}

// [[Rcpp::export]]
arma::vec MultLambdaTrio(const arma::vec& x, const arma::vec& v){
  arma::vec y(3, arma::fill::zeros);
   y[0] = (x[1]+x[2])*v[1] + x[0]*v[2];
   y[1] = (x[0]+x[2])*v[1] + x[1]*v[2];
   y[2] = (x[0]+x[1])*v[1] + x[2]*v[2];
   return(y);
}

// [[Rcpp::export]]
arma::mat MultLambdaGroups(arma::mat x, arma::mat v, arma::uvec n, arma::uvec p){
  arma::uword s = 0;
  arma::uword i,j;
    for(i=0;i<p.n_elem;++i){
        for(j=0;j<x.n_cols;++j){
            if(p(i) > 3){
                x( arma::span(s,s+n(i)-1), j ) = MultLambda( x( arma::span(s,s+n(i)-1), j ) , v.col(i), p(i) );
            } else {
                x( arma::span(s,s+n(i)-1), j ) = MultLambdaTrio( x( arma::span(s,s+n(i)-1), j ) , v.col(i) );
            }
        }
        s += n(i);
    }    
    return(x);
}

// [[Rcpp::export]]
arma::vec eigenVals(double tau, arma::uword n){
	double stau = tau/0.5;
  arma::vec L(3, arma::fill::zeros);
	    L[0] = 1 - stau;
	    L[1] = 1 + stau * (double(n) - 4)/2;
	    L[2] = 1 + stau * (double(n) - 2);
	return(L);
}

// [[Rcpp::export]]
arma::mat eigenVecs(arma::uword n){
	double nd = double(n);
	double N = nd*(nd-1)/2.0;
  arma::mat V(3,3);
	V(0,0) = 1 / (0.5 * (nd-1) * (nd-2) );
	V(0,1) = -(nd-3) / ( (nd-1) * (nd-2) );
	V(0,2) = (nd-3) / (nd-1);
	V(1,0) = -8 * (N-nd) / ( nd*nd * (nd-2)* (nd-3) );
	V(1,1) = (2*N/nd-3) /  (nd*(nd-2));
	V(1,2) = 2/nd;
	V(2,0) = 2 / ((nd-1)*nd);
	V(2,1) = 2 / ((nd-1)*nd);
	V(2,2) = 2 / ((nd-1)*nd);
	return(V);
}

// [[Rcpp::export]]
arma::uvec eigenCount(arma::uword n){
  arma::uvec lcnt(3);
	lcnt[0] = (n-1)*n/2 - (n-1) - 1;
	lcnt[1] = n-1;
	lcnt[2] = 1;
	return(lcnt);	
}

// [[Rcpp::export]]
arma::uvec matrixCount(arma::uword n){
  arma::uvec lcnt(3);
	lcnt[0] = 0;
	for(int i=n-3;i>0;i--){
		lcnt[0] += i;
	}
	lcnt[1] = (n-1)*n/2 - lcnt[0] - 1;
	lcnt[2] = 1;
	return(lcnt);
}

// [[Rcpp::export]]
arma::vec matrixVals(arma::vec lambda, arma::mat e){
		return( trans(e.row(0))/lambda[0] + trans(e.row(1))/lambda[1] + trans(e.row(2))/lambda[2] );
}

// [[Rcpp::export]]
double likelihood(arma::vec v, arma::uvec cnt, arma::vec l, arma::uvec lcnt, arma::vec x, arma::vec y, arma::uword p, double sigma){
  // x is assumed to be scaled and centered
  int i;
  arma::vec LambdaY = MultLambda(y, v, p);
  arma::vec LambdaX = MultLambda(x, v, p);
  double Int = 0;
  for(i=0;i<3;i++){
    Int += double(cnt[i])*v[i];
  }

  arma::vec YLambdaX(2, arma::fill::zeros);
  YLambdaX[0] = arma::sum(y*Int);
  YLambdaX[1] = arma::dot(y, LambdaX);

  double V = arma::dot(x, LambdaX);
  double M = double(x.n_elem)*Int;

  arma::vec L(3, arma::fill::zeros);
  for(i=0;i<3;i++){
    L[i] = double(lcnt[i])*log(l[i]/sigma);
  }

  return( -0.5 * sigma * ( arma::dot(y, LambdaY) - YLambdaX[0]*YLambdaX[0]/M - YLambdaX[1]*YLambdaX[1]/V) - (0.5*double(x.n_elem)-1)*log(2*arma::datum::pi) - 0.5*log(M/sigma) - 0.5*log(V/sigma) - 0.5*arma::sum(L) );
}

double priorTauSigma(double tau, double sigma){
	return( 2*R::dnorm( sqrt(1/sigma), 0, 10, 1) );
}

double logit(double p, double base=1){ 
		return( log(p) - log(base-p) );
	}
	
double rlogit(double x, double base=1){
	return( base*exp(x) / (1 + exp(x)) );
}

// [[Rcpp::export]]
arma::vec sampleTauSigma(arma::vec lltausig, arma::uword iter, double tune, arma::mat e, arma::uvec cnt, arma::uvec lcnt, arma::vec x, arma::vec y, arma::uword p){
	arma::vec prop(3);
	arma::vec l(3);
	arma::vec v(3);
	double r;
	for(int i=0;i<iter;i++){
		prop[1] = rlogit( R::rnorm( logit(lltausig[1], 0.5), tune) , 0.5);
		prop[2] = exp( R::rnorm( log(lltausig[2]), tune) ); 
		l = eigenVals(prop[1], p);
		v = matrixVals(l, e);
		prop[0] = likelihood(v, cnt, l, lcnt, x, y, p, prop[2]);
		r = prop[0] + priorTauSigma(prop[1], prop[2]) - lltausig[0] - priorTauSigma(lltausig[1], lltausig[2]);
		if( log(R::runif(0,1)) < r ){
			lltausig = prop;
		}
	}
	return(lltausig);
}

// [[Rcpp::export]]
arma::ivec build_index(int nr, int nc){
	int j;
	int stopping = nr*nc-2;
	int i = 1;
  arma::ivec vn(2*nr*nc - nr - nc);
	for(j=0;j<nr-1;++j){
		vn[j] = j;
	}
	j = nr-1;
	vn[j] = 0;
	while(j < vn.n_elem){
		++j;
		vn[j] = i;
		++j;
		vn[j] = i + (nr-1); 
		++i;
		if( (i % nr) == 0 ){
			++j;
			vn[j] = i;
			++i;
		}
	}
	return(vn);
}

// [[Rcpp::export]]
arma::ivec build_pointer(int nr, int nc){
  arma::ivec block(nr*nc+1);
	block[0]=0;
	int j = 1;
	int k, i;
	int a = 1;
	for(k=0;k<nc;k++){
		block[j]=block[j-1]+1-a;
		j++;
		for(i=1;i<nr;i++){
			block[j]=block[j-1]+2-a;
			j++;
		}
		if(k==0){ a=0; }
	}
	return(block);
}

// [[Rcpp::export]]
arma::vec fill_values(int nr, int nc, arma::ivec ind, arma::ivec p, arma::vec R){
	int i, j, pn, np;
  arma::vec x(2*nr*nc - nr - nc);
	int k = 0;
	for(i=0;i<nr*nc;i++){
		pn = p[i];
		np = p[i+1];
		for(j=pn;j<np;j++){			
			if( i != ind[j] ){
				x[k] = -(R[i] + R[ind[j]])/2;
			}		
			k++;
		}
	}
	return(x);
}

// [[Rcpp::export]]
arma::umat all_pairs(arma::uword p){
	int i,j,z;
  arma::umat out(p*(p-1)/2, 2, arma::fill::zeros);

	z=0;
	for(i=0;i<p-1;++i){
		for(j=i+1;j<p;++j){
			out(z,0) = i;
			out(z,1) = j;
			++z;
		}
	}
	return(out);
}

