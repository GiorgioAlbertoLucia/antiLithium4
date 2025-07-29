#pragma once
 
#include <RooAbsPdf.h>
#include <RooRealProxy.h>
 
class RooCustomVoigtian : public RooAbsPdf {
public:
  RooCustomVoigtian() {}
  RooCustomVoigtian(const char *name, const char *title,
              RooAbsReal& _x, RooAbsReal& _mean,
              RooAbsReal& _width, RooAbsReal& _sigma,
              RooAbsReal& _alpha, RooAbsReal& _n,
              bool doFast = false);
  RooCustomVoigtian(const RooCustomVoigtian& other, const char* name=nullptr) ;
  TObject* clone(const char* newname=nullptr) const override { return new RooCustomVoigtian(*this,newname); }
 
  /// Enable the fast evaluation of the complex error function using look-up
  /// tables (default is the "slow" CERNlib algorithm).
  inline void selectFastAlgorithm()    { _doFast = true;  }
 
  /// Disable the fast evaluation of the complex error function using look-up
  /// tables (default is the "slow" CERNlib algorithm).
  inline void selectDefaultAlgorithm() { _doFast = false; }
 
protected:
 
  RooRealProxy x ;
  RooRealProxy mean ;
  RooRealProxy width ;
  RooRealProxy sigma ;  
  RooRealProxy alpha ; // NEW: shape parameter for the Voigtian distribution
  RooRealProxy n ; // NEW: shape parameter for the Voigtian distribution
 
  double evaluate() const override ;
 
private:
 
  bool _doFast = false;
  ClassDefOverride(RooCustomVoigtian,2) // Voigtian PDF (Gauss (x) BreitWigner)
};


///////////////////////// Class implementation /////////////////////////

/** \class RooCustomVoigtian
    \ingroup Roofit

Custom modification of the RooVoigtian class to include the shape parameter.

RooVoigtian is an efficient implementation of the convolution of a
Breit-Wigner with a Gaussian, making use of the complex error function.
RooFitCore provides two algorithms for the evaluation of the complex error
function (the default CERNlib C335 algorithm, and a faster, look-up-table
based method). By default, RooVoigtian employs the default (CERNlib)
algorithm. Select the faster algorithm either in the constructor, or with
the selectFastAlgorithm() method.
 
\note The "width" parameter that determines the Breit-Wigner shape
      represents the **full width at half maximum (FWHM)** of the
      Breit-Wigner (often referred to as \f$\Gamma\f$ or \f$2\gamma\f$).
**/
  
#include <RooMath.h>
 
#include <cmath>
#include <complex>
 
 
////////////////////////////////////////////////////////////////////////////////
/// Construct a RooVoigtian PDF, which represents the convolution of a
/// Breit-Wigner with a Gaussian.
/// \param name Name that identifies the PDF in computations.
/// \param title Title for plotting.
/// \param _x The observable for the PDF.
/// \param _mean The mean of the distribution.
/// \param _width The **full width at half maximum (FWHM)** of the Breit-Wigner
///               (often referred to as \f$\Gamma\f$ or \f$2\gamma\f$).
/// \param _sigma The width of the Gaussian distribution.
/// \param _alpha The shape parameter of the Voigtian distribution.
/// \param _n The shape parameter of the Voigtian distribution (optional, can be
///            used to adjust the shape of the Voigtian).
/// \param doFast Use the faster look-up-table-based method for the evaluation
///               of the complex error function.
 
RooCustomVoigtian::RooCustomVoigtian(const char *name, const char *title,
          RooAbsReal& _x, RooAbsReal& _mean,
          RooAbsReal& _width, RooAbsReal& _sigma,
          RooAbsReal& _alpha, // NEW: shape parameter
          RooAbsReal& _n, // NEW: shape parameter
          bool doFast) :
  RooAbsPdf(name,title),
  x("x","Dependent",this,_x),
  mean("mean","Mean",this,_mean),
  width("width","Breit-Wigner Width",this,_width),
  sigma("sigma","Gauss Width",this,_sigma),
  alpha("alpha","Shape Parameter",this,_alpha),
  n("n","Shape Parameter",this,_n),
  _doFast(doFast)
{
 
}
 
////////////////////////////////////////////////////////////////////////////////
 
RooCustomVoigtian::RooCustomVoigtian(const RooCustomVoigtian& other, const char* name) :
  RooAbsPdf(other,name), x("x",this,other.x), mean("mean",this,other.mean),
  width("width",this,other.width),sigma("sigma",this,other.sigma),
  alpha("alpha",this,other.alpha), // NEW: shape parameter
  n("n",this,other.n), // NEW: shape parameter
  _doFast(other._doFast)
{
 
}
 
////////////////////////////////////////////////////////////////////////////////
 
double RooCustomVoigtian::evaluate() const
{
  double s = std::abs((double)sigma);
  double w = std::abs((double)width);
  double a = (double)alpha;
  double shape_n = (double)n;
  
  double arg = x - mean;
  double s_eff = sigma + alpha * arg;
  double w_eff = width; 

  // Avoid nonphysical sigma
  if (s_eff <= 0.0) return 0.0;  // or clamp to small positive value  

  // Special cases
  if (s_eff == 0. && w_eff == 0.) return 1.0;
  if (s_eff == 0.) return 1.0 / (arg * arg + 0.25 * w_eff * w_eff);
  if (w_eff == 0.) return std::exp(-0.5 * std::pow(arg, 2) / (s_eff * s_eff));  

  // Modify Gaussian tails using n
  s_eff *= std::sqrt(1.0 + shape_n * arg * arg);  

  double c = 1.0 / (std::sqrt(2.0) * s_eff);
  double a_voigt = 0.5 * c * w_eff;
  double u = c * arg; 

  std::complex<double> z(u, a_voigt);
  std::complex<double> v; 

  if (_doFast) {
      v = RooMath::faddeeva_fast(z);
  } else {
      v = RooMath::faddeeva(z);
  } 

  return c * v.real();
}
