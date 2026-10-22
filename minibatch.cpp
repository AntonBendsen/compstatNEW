#include <Rcpp.h>
#include <chrono>
#include <cmath>
using namespace Rcpp ;

/* 8.0: Model. */

// [[Rcpp::export]]
Rcpp::NumericVector dose_response(
    const Rcpp::NumericVector& x,
    const Rcpp::NumericVector& theta
) {
  double alpha = theta[0];
  double beta  = theta[1];
  double gamma = theta[2];
  double delta = theta[3];
  
  int n = x.size();
  Rcpp::NumericVector out(n);

  for (int i = 0; i < n; ++i) {
    double p = R::plogis(
      -beta * (std::log(x[i]) - alpha),
      0.0, 1.0, 1, 0
    );
  
    out[i] = gamma + (delta - gamma) * p;
  }
  
  return out;
}


/* 8.1: Unit loss. */
// [[Rcpp::export]]
double
  loss_i(const Rcpp::NumericVector& theta, double x, double y)
  {
    double residual = y - dose_response(
      Rcpp::NumericVector::create(x), theta)[0] ;
    return 0.5 * residual * residual ;
  }


/* 8.2: Minibatch gradient. */
// [[Rcpp::export]]
Rcpp::NumericVector
  grad_i(const Rcpp::NumericVector theta, double x, double y)
  {
    double alpha = theta[0] ;
    double beta  = theta[1] ;
    double gamma = theta[2] ;
    double delta = theta[3] ;
    
    double lx = std::log(x) ;
    double p = R::plogis(-beta * (lx - alpha), 0, 1, 1, 0) ;
    double q = p * (1 - p) ;
    
    double fit = gamma + (delta - gamma) * p ;
    double residual = y - fit ;
    
    Rcpp::NumericVector df = {(delta - gamma) * beta * q,
                              -(delta - gamma) * (lx - alpha) * q,
                              1 - p,
                              p};
    
    return -residual * df ;
    
  }



// [[Rcpp::export]]
Rcpp::List sgd_batch(
    const Rcpp::NumericVector& x,
    const Rcpp::NumericVector& y,
    Rcpp::NumericVector theta,
    List batches,
    int step,
    double lr0,
    int K,
    double a,
    Rcpp::NumericMatrix trajectory
) {
  const auto start_time = std::chrono::steady_clock::now();
  Rcpp::NumericVector theta_local = Rcpp::clone(theta);
  int n = batches.size();
  Rcpp::NumericVector gradient_norm(n);
  Rcpp::NumericVector gradient_time(n);
  Rcpp::NumericVector elapsed_time(n);
  
  for (int i = 0; i < n; ++i) {
    IntegerVector batch = batches[i];
    step += 1;
    
    double lr = lr0 /
      std::pow(1.0 + static_cast<double>(step) / K, a);
    
    int m = batch.size();
    Rcpp::NumericVector g(theta_local.size());
    
    for (int j = 0; j < m; ++j) {
      g += grad_i(theta_local, x[batch[j] - 1], y[batch[j] - 1]);
    }
    g = g / m;
    gradient_norm[i] = std::sqrt(Rcpp::sum(g * g));
    gradient_time[i] = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start_time
    ).count();

    theta_local += -lr * g;
    trajectory(step, Rcpp::_) = theta_local;
    elapsed_time[i] = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start_time
    ).count();
  }
  
  return List::create(
    Named("theta") = theta_local,
    Named("trajectory") = trajectory,
    Named("step") = step,
    Named("gradient_norm") = gradient_norm,
    Named("gradient_time") = gradient_time,
    Named("elapsed_time") = elapsed_time
  );
}
