#include <cmath>

double evaluateCrystalBallTail(double t, double alpha, double n)
{
   double r = n / alpha;
   double b = r - alpha;
 
   return std::exp(-0.5 * alpha * alpha) * std::pow(r / (b - t), n);
}

double CrystalBall(double* xs, double* pars)
{
   double x = xs[0];
   double x0 = pars[0];
   double sigma = pars[1];
   double alphaL = pars[2];
   double nL = pars[3];
   double alphaR = pars[4];
   double nR = pars[5];
 
   const double t = (x - x0) / sigma;
 
   if (t < -alphaL) {
      return evaluateCrystalBallTail(t, alphaL, nL);
   } else if (t <= alphaR) {
      return std::exp(-0.5 * t * t);
   } else {
      return evaluateCrystalBallTail(-t, alphaR, nR);
   }
}