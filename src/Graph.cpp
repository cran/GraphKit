// [[Rcpp::depends(RcppArmadillo)]]
#include "helpers.h"
#include <vector>
#include <array>
#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

double quantileCpp(arma::vec x, double probs) {
    Rcpp::Environment stats("package:stats");
    Rcpp::Function quantile = stats["quantile"];

    double  ans = Rcpp::as<double>(quantile(x, probs));
    
    return ans;
}

// [[Rcpp::export]]
arma::mat bootstrap (arma::mat thetaHat, arma::mat data, int numB) {
  int d = thetaHat.n_cols;
  int n = data.n_rows;
  int z = ((d-1)*d)/2;
  arma::mat M = data * thetaHat;
  arma::mat mat_norm = arma::mat(numB,n,arma::fill::randn);

  arma::mat BSR = arma::mat(numB, z, arma::fill::zeros);
  arma::mat Mt = M.t();

  //extract xi
  for(int i = 0; i < numB; i++) {
    arma::rowvec xi = mat_norm.row(i);
    double sum = 0;
    int c = 0;

    for(int j = 0; j < d; j++) {sum += xi.at(j);} //sum row

    arma::mat tmp = Mt;
    arma::mat tmp2 = M;


    for(int j = 0; j < n; j++) {
      for(int k = 0; k < d; k++) {
        tmp2.at(j,k) = tmp2.at(j,k) * xi.at(j);
      }
    }

    tmp = tmp * tmp2;
    tmp = tmp - thetaHat * sum;

    //idx process
    for(int j = 0; j < d; j++) {
      for(int k = 0; k < j; k++) {

        BSR.at(i,c) = tmp.at(j,k)/double(n);
        c += 1;
      }
    }
  }
  return BSR;
}

// [[Rcpp::export]]
arma::mat debias (arma::mat sigmaHat, arma::mat thetaHat) {
    arma::mat denom = thetaHat.t() * sigmaHat;
    arma::mat num = denom * thetaHat - thetaHat;
    arma:: vec bottom = arma::diagvec(denom);
    arma::mat tmp = arma::repmat(bottom,1,thetaHat.n_rows);
    arma::mat inner = arma::trans(thetaHat - arma::trans(num.t()/tmp));
    return inner;
}

// [[Rcpp::export]]
bool skipDownConn(arma::mat deOmega, arma::mat bootStats, double alpha, int numConn) {
    int v = deOmega.n_cols;
    UnionFind components = UnionFind(v);
    arma::vec na = arma::zeros<arma::vec>(bootStats.n_rows);
    
    while(components.getCount() > numConn) {
        vector <array<int, 2> > critical;
        
        //build the critical set
        for(int i = 0; i < v; i++) {
            for(int j= i + 1; j < v; j++) {
                if(!components.connected(i,j)) {
                    array<int,2> edge = {i,j};
                    critical.push_back(edge);
                }
            }
        }
        
        arma::uvec cols(bootStats.n_cols,arma::fill::zeros);
        
        bool edited = false;
        
        //build the rejected set
        for(int i = 0; i < critical.size(); i++) {
            int e_i = critical[i][0];
            int e_j = critical[i][1];
            
            cols(edgeNum(e_i,e_j)) = 1;
        }
        
        arma::mat bootFilt = bootStats;
        
        for(int i = 0; i < bootFilt.n_cols; i++) {
            if(!cols(i)) {bootFilt.col(i) = na;}
        }
        
        arma::vec maxes(bootFilt.n_rows);
        
        for(int i = 0; i < maxes.size(); i++) {
            maxes(i) = arma::abs(bootFilt.row(i)).max();
        }
        
        double statsMM = quantileCpp(maxes,1-alpha);
        for(int i = 0; i < critical.size(); i++) {
            int e_i = critical[i][0];
            int e_j = critical[i][1];
            
            if(abs(deOmega(e_i,e_j)) > statsMM && !components.connected(e_i, e_j)) {
                edited = true;
                components.connect(e_i,e_j);
                if(components.getCount() <= numConn) {return true;}
            }
        }
        
        if(!edited) {return components.getCount() <= numConn;}
    }
    
    return true;
}

// [[Rcpp::export]]
arma::vec skipDownConnCI (arma::mat deOmega, arma::mat bootStats, double alpha) {
    int v = deOmega.n_cols;
    UnionFind components = UnionFind(v);
    arma::vec na = arma::zeros<arma::vec>(bootStats.n_rows);
    arma::vec gen = arma::randu<arma::vec>(1);
    bool upper = gen(0) > 1 - alpha/2;

    while(true) {
        vector < array<int, 2> > critical;

        //build the critical set
        for(int i = 0; i < v; i++) {
            for(int j = i + 1; j < v; j++) {
                if(!components.connected(i,j)) {
                    array<int,2> edge = {i,j};
                    critical.push_back(edge);
                }
            }
        }
        
      
        
        arma::uvec cols(bootStats.n_cols,arma::fill::zeros);
        
        bool edited = false;

        //build the rejected set
        for(int i = 0; i < critical.size(); i++) {
            int e_i = critical[i][0];
            int e_j = critical[i][1];
            cols(0) = 1;
            cols(edgeNum(e_i,e_j)) = 1;
        }

        arma::mat bootFilt = bootStats;

        for(int i = 0; i < bootFilt.n_cols; i++) {
            if(!cols(i)) {bootFilt.col(i) = na;}
        }

        arma::vec maxes(bootFilt.n_rows);

        for(int i = 0; i < maxes.size(); i++) {
            maxes(i) = arma::abs(bootFilt.row(i)).max();
        }
        
        double statsMM = quantileCpp(maxes, 1 - alpha/2);
        for(int i = 0; i < critical.size(); i++) {
            int e_i = critical[i][0];
            int e_j = critical[i][1];

            if(abs(deOmega(e_i,e_j)) > statsMM && !components.connected(e_i, e_j)) {
                edited = true;
                components.connect(e_i,e_j);
                }
        }

        if(!edited || components.getCount() == 1) {
            arma::vec CI = arma::zeros<arma::vec>(2);
            CI(0) = components.getCount();
            CI(1) = upper ? CI(0) : 1;
            return CI;
        }
    }

    return arma::zeros<arma::vec>(2);
}

// [[Rcpp::export]]
bool skipDownChain (arma::mat deOmega, arma::mat bootStats, double alpha , int k) {
    return false;
}

// [[Rcpp::export]]
arma::vec skipDownChainCI (arma::mat deOmega, arma::mat bootStats, double alpha) {
    return arma::zeros<arma::vec>(2);
}

// [[Rcpp::export]]
bool skipDownDeg (arma::mat deOmega, arma::mat bootStats, double alpha , int k) {
    int v = deOmega.n_cols;
    int maxDegree = 0;
    arma::vec degrees = arma::zeros<arma::vec>(v);
    arma::mat edges = arma::zeros<arma::mat>(v,v);
    arma::vec na = arma::zeros<arma::vec>(bootStats.n_rows);
    
    
    while(true) {
        vector <array<int, 2> > critical;
        
        //build the critical set
        for(int i = 0; i < v; i++) {
            for(int j= i + 1; j < v; j++) {
                if(!edges(i,j)) {
                    array<int,2> edge = {i,j};
                    critical.push_back(edge);
                }
            }
        }
        
        arma::uvec cols(bootStats.n_cols,arma::fill::zeros);
        
        bool edited = false;
        
        //build the rejected set
        for(int i = 0; i < critical.size(); i++) {
            int e_i = critical[i][0];
            int e_j = critical[i][1];
            
            cols(edgeNum(e_i,e_j)) = 1;
        }
        
        arma::mat bootFilt = bootStats;
        
        for(int i = 0; i < bootFilt.n_cols; i++) {
            if(!cols(i)) {bootFilt.col(i) = na;}
        }
        
        arma::vec maxes(bootFilt.n_rows);
        
        for(int i = 0; i < maxes.size(); i++) {
            maxes(i) = arma::abs(bootFilt.row(i)).max();
        }
        
        double statsMM = quantileCpp(maxes,1-alpha);
        for(int i = 0; i < critical.size(); i++) {
            int e_i = critical[i][0];
            int e_j = critical[i][1];
            
            if(abs(deOmega(e_i,e_j)) > statsMM) {
                edited = true;
                edges(e_i,e_j) = 1;
                degrees(e_i) += 1;
                degrees(e_j) += 1;
                maxDegree = degrees(e_i) > maxDegree ? degrees(e_i):maxDegree;
                maxDegree = degrees(e_j) > maxDegree ? degrees(e_j):maxDegree;
                
                if(maxDegree > k) {return false;}
            }
        }
        
        if(!edited) {return maxDegree <= k;}
    }
    
    return true;
}

// [[Rcpp::export]]
arma::vec skipDownDegCI (arma::mat deOmega, arma::mat bootStats, double alpha) {
    int v = deOmega.n_cols;
    int maxDegree = 0;
    arma::vec na = arma::zeros<arma::vec>(bootStats.n_rows);
    arma::vec gen = arma::randu<arma::vec>(1);
    arma::vec degrees = arma::zeros<arma::vec>(v);
    arma::mat edges = arma::zeros<arma::mat>(v,v);

    bool upper = gen(0) > 1 - alpha/2;
    
    while(true) {
        vector < array<int, 2> > critical;
        
        //build the critical set
        for(int i = 0; i < v; i++) {
            for(int j = i + 1; j < v; j++) {
                if(!edges(i,j)) {
                    array<int,2> edge = {i,j};
                    critical.push_back(edge);
                }
            }
        }
        
        
        
        arma::uvec cols(bootStats.n_cols,arma::fill::zeros);
        
        bool edited = false;
        
        //build the rejected set
        for(int i = 0; i < critical.size(); i++) {
            
            int e_i = critical[i][0];
            int e_j = critical[i][1];
            
            cols(edgeNum(e_i,e_j)) = 1;
        }
        
        arma::mat bootFilt = bootStats;
        
        for(int i = 0; i < bootFilt.n_cols; i++) {
            if(!cols(i)) {bootFilt.col(i) = na;}
        }
        
        arma::vec maxes(bootFilt.n_rows);
        
        for(int i = 0; i < maxes.size(); i++) {
            maxes(i) = arma::abs(bootFilt.row(i)).max();
        }
        
        double statsMM = quantileCpp(maxes, 1 - alpha/2);
        for(int i = 0; i < critical.size(); i++) {
            int e_i = critical[i][0];
            int e_j = critical[i][1];
            
            if(abs(deOmega(e_i,e_j)) > statsMM) {
                edited = true;
                edges(e_i,e_j) = 1;
                degrees(e_i) += 1;
                degrees(e_j) += 1;
                maxDegree = degrees(e_i) > maxDegree ? degrees(e_i):maxDegree;
                maxDegree = degrees(e_j) > maxDegree ? degrees(e_j):maxDegree;
            }
        }
        
        if(!edited || maxDegree == v-1) {
            arma::vec CI = arma::zeros<arma::vec>(2);
            CI(0) = maxDegree;
            CI(1) = upper ? CI(0) : 0;
            return CI;
        }
    }
    
    return arma::zeros<arma::vec>(2);
}

// [[Rcpp::export]]
bool skipDownClique (arma::mat deOmega, arma::mat bootStats, double alpha , int k) {
    return false;
}

// [[Rcpp::export]]
arma::vec skipDownCliqueCI (arma::mat deOmega, arma::mat bootStats, double alpha) {
    return arma::zeros<arma::vec>(2);
}

// [[Rcpp::export]]
bool skipDownChrom (arma::mat deOmega, arma::mat bootStats, double alpha , int k) {
    return false;
}

// [[Rcpp::export]]
arma::vec skipDownChromCI (arma::mat deOmega, arma::mat bootStats, double alpha) {
    return arma::zeros<arma::vec>(2);
}

// [[Rcpp::export]]
bool skipDownSingle (arma::mat deOmega, arma::mat bootStats, double alpha , int k) {
    int v = deOmega.n_cols;
    arma::vec singletons = arma::ones<arma::vec>(v);
    arma::vec na = arma::zeros<arma::vec>(bootStats.n_rows);
    int numSingletons = v;
    
    while(true) {
        vector <array<int, 2> > critical;
        
        //build the critical set
        for(int i = 0; i < v; i++) {
            for(int j= i + 1; j < v; j++) {
                if(singletons(i) || singletons(j)) {
                    array<int,2> edge = {i,j};
                    critical.push_back(edge);
                }
            }
        }
        
        arma::uvec cols(bootStats.n_cols,arma::fill::zeros);
        
        bool edited = false;
        
        //build the rejected set
        for(int i = 0; i < critical.size(); i++) {
            int e_i = critical[i][0];
            int e_j = critical[i][1];
            
            cols(edgeNum(e_i,e_j)) = 1;
        }
        
        arma::mat bootFilt = bootStats;
        
        for(int i = 0; i < bootFilt.n_cols; i++) {
            if(!cols(i)) {bootFilt.col(i) = na;}
        }
        
        arma::vec maxes(bootFilt.n_rows);
        
        for(int i = 0; i < maxes.size(); i++) {
            maxes(i) = arma::abs(bootFilt.row(i)).max();
        }
        
        double statsMM = quantileCpp(maxes,1-alpha);
        for(int i = 0; i < critical.size(); i++) {
            int e_i = critical[i][0];
            int e_j = critical[i][1];
            
            if(abs(deOmega(e_i,e_j)) > statsMM) {
                edited = true;
                if (singletons(e_i)) {
                    singletons(e_i) = 0;
                    numSingletons-= 1;
                }
                
                if (singletons(e_j)) {
                    singletons(e_j) = 0;
                    numSingletons-= 1;
                }
                if(numSingletons <= k) {return true;}
            }
        }
        
        if(!edited) {return arma::sum(singletons) <= k;}
    }
    
    return true;
}

// [[Rcpp::export]]
arma::vec skipDownSingleCI (arma::mat deOmega, arma::mat bootStats, double alpha) {
    int v = deOmega.n_cols;
    arma::vec na = arma::zeros<arma::vec>(bootStats.n_rows);
    arma::vec gen = arma::randu<arma::vec>(1);
    arma::vec singletons = arma::ones<arma::vec>(v);
    int numSingletons = v;
    bool upper = gen(0) > 1 - alpha/2;

    while(true) {
        vector < array<int, 2> > critical;
        
        //build the critical set
        for(int i = 0; i < v; i++) {
            for(int j = i + 1; j < v; j++) {
                if(singletons(i) || singletons(j)) {
                    array<int,2> edge = {i,j};
                    critical.push_back(edge);
                }
            }
        }
        
        
        
        arma::uvec cols(bootStats.n_cols,arma::fill::zeros);
        
        bool edited = false;
        
        //build the rejected set
        for(int i = 0; i < critical.size(); i++) {
            
            int e_i = critical[i][0];
            int e_j = critical[i][1];

            cols(edgeNum(e_i,e_j)) = 1;
        }
        
        arma::mat bootFilt = bootStats;
        
        for(int i = 0; i < bootFilt.n_cols; i++) {
            if(!cols(i)) {bootFilt.col(i) = na;}
        }
        
        arma::vec maxes(bootFilt.n_rows);
        
        for(int i = 0; i < maxes.size(); i++) {
            maxes(i) = arma::abs(bootFilt.row(i)).max();
        }
        
        double statsMM = quantileCpp(maxes, 1 - alpha/2);
        for(int i = 0; i < critical.size(); i++) {
            int e_i = critical[i][0];
            int e_j = critical[i][1];
            
            if(abs(deOmega(e_i,e_j)) > statsMM) {
                edited = true;
                if (singletons(e_i)) {
                    singletons(e_i) = 0;
                    numSingletons-= 1;
                }
                
                if (singletons(e_j)) {
                    singletons(e_j) = 0;
                    numSingletons-= 1;
                }
            }
        }
        
        if(!edited || numSingletons == 0) {
            arma::vec CI = arma::zeros<arma::vec>(2);
            CI(0) = numSingletons;
            CI(1) = upper ? CI(0) : 0;
            return CI;
        }
    }
    
    return arma::zeros<arma::vec>(2);
}

// [[Rcpp::export]]
bool skipDownGirth (arma::mat deOmega, arma::mat bootStats, double alpha , int k) {
    return false;
}

// [[Rcpp::export]]
arma::vec skipDownGirthCI (arma::mat deOmega, arma::mat bootStats, double alpha) {
    return arma::zeros<arma::vec>(2);
}

// [[Rcpp::export]]
bool skipDownMatch (arma::mat deOmega, arma::mat bootStats, double alpha , int k) {
    return false;
}

// [[Rcpp::export]]
bool skipDownPlan (arma::mat deOmega, arma::mat bootStats, double alpha , int k) {
    return false;
}

// [[Rcpp::export]]
bool skipDownBipar (arma::mat deOmega, arma::mat bootStats, double alpha , int k) {
    return false;
}

// [[Rcpp::export]]
bool skipDownCycle (arma::mat deOmega, arma::mat bootStats, double alpha) {
    int v = deOmega.n_cols;
    UnionFind components = UnionFind(v);
    arma::vec na = arma::zeros<arma::vec>(bootStats.n_rows);
    arma::mat edges = arma::zeros<arma::mat>(v,v);
    
    while(true) {
        vector <array<int, 2> > critical;
        
        //build the critical set
        for(int i = 0; i < v; i++) {
            for(int j= i + 1; j < v; j++) {
                if(!edges(i,j)) {
                    array<int,2> edge = {i,j};
                    critical.push_back(edge);
                }
            }
        }
        
        arma::uvec cols(bootStats.n_cols,arma::fill::zeros);
        
        bool edited = false;
        
        //build the rejected set
        for(int i = 0; i < critical.size(); i++) {
            int e_i = critical[i][0];
            int e_j = critical[i][1];
            
            cols(edgeNum(e_i,e_j)) = 1;
        }
        
        arma::mat bootFilt = bootStats;
        
        for(int i = 0; i < bootFilt.n_cols; i++) {
            if(!cols(i)) {bootFilt.col(i) = na;}
        }
        
        arma::vec maxes(bootFilt.n_rows);
        
        for(int i = 0; i < maxes.size(); i++) {
            maxes(i) = arma::abs(bootFilt.row(i)).max();
        }
        
        double statsMM = quantileCpp(maxes,1-alpha);
        for(int i = 0; i < critical.size(); i++) {
            int e_i = critical[i][0];
            int e_j = critical[i][1];
            
            if(abs(deOmega(e_i,e_j)) > statsMM) {
                edited = true;
                if(!components.connected(e_i,e_j)){
                components.connect(e_i,e_j);
                edges(e_i,e_j) = 1;
                }
                else {
                   return false;        
                }
            }
        }
        
        if(!edited) {return true;}
    }
    
    return true;
}

// [[Rcpp::export]]
arma::rowvec stepDown(arma::vec stats, arma::mat bootStats, double alpha) {
    int len = stats.size();
    int sum = len;
    int height = bootStats.n_rows;
    arma::vec Enull = arma::vec(len, arma::fill::ones);
    arma::rowvec maxVec = arma::rowvec(height);

    while(sum) {

        for(int i = 0; i < height; i++) {
            double max = 0;

            for(int j = 0; j < len; j++) {
                if(Enull(j)) {
                    if(std::abs(bootStats(i,j)) > max) {
                        max =  std::abs(bootStats(i,j));
                    }
                }
            }

            maxVec(i) = max;
        }

        double statsMM = quantileCpp(maxVec.t(), 1 - alpha);
        bool edited = false;

        for(int j = 0; j < len; j++) {
            if(Enull(j)) {
                if(std::abs(stats(j)) > statsMM) {
                    Enull(j) = false;
                    sum--;
                    edited = true;
                }
            }
        }

        if(!edited) {break;}
    }

    len = len - sum;
    arma::rowvec diff = arma::rowvec(len);
    int index = 0;
    int j = 0;

    while(index < len) {

        if(!Enull(j)) {
            diff(index) = j + 1;
            index++;
        }

        j++;
    }

    return diff;
}
