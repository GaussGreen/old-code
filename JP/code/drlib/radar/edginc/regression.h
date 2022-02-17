// Author: ???

/****************************************************************************/
/*                                                                          */
/*  The description of the regression objects used in the fitting procedures*/
/*                                                                          */
/****************************************************************************/
/*       basis.H                                                          */
/****************************************************************************/

#ifndef REGRESSION_H
#define REGRESSION_H

#include "edginc/RadarRepUtil.hpp" // TDate
#include <iostream>
#include <fstream>
#include <vector>

DRLIB_BEGIN_NAMESPACE

using namespace std;

// typedef long int TDate;
#define MAX(a,b) ((a) > (b) ? (a) : (b)) 
#define MIN(a,b) ((a) < (b) ? (a) : (b)) 
#define TOL 1.0e-10		// Default value for single precision and variables scaled to order unity.
#define SQR(a) ((a)*(a))
#define SIGN(a,b) ((b) >= 0.? fabs(a): -fabs(a))

//! A class that describes the merit function of a regression
class  RADAR_DLL Regression{
public:
    Regression(){}
//computes the coefficients for the regression using  a container of basis functions
// Basic setup is we have a matrix A and a vector b.  We want to solve for x
// so that Ax = b using regression.
// A(i,*) = row i = EvaluateBasis(x_k), where x_k is fitting variables at a
// point on a path.
// b(i) = instrument value at that point we are trying to match for path i.
// basisValues = A, dealValues = b
// coeffs is result of regression = x (in Ax=b)
    virtual void DoTheRegression(
                             vector<vector<double> > & basisValues, 
                             vector<double> & dealValues, 
                             vector<double> & coeffs) = 0;
    // When quants finish determining choice of basis, regression, etc. for a
    // product, they need to output everything for use in production.  The
    // following write is used to store the type of regression used.
    // TODO:  The type of the regression should really been written one level
    // up so deserialization is cleaner.
    //write the state to a file
    virtual void writeState(std::ostream & ostr) = 0;
    //read the state from a file
    virtual void readState(std::istream & istr) = 0;

    //Get the standard deviation of the error
    double GetChi(){return chisq;}
    //Get the errors of the fitted values
    vector<double> GetErrors(){return errors;}
    virtual ~Regression() {}
    protected:
        //!the standard deviation of the error
        double chisq;
        //!the errors of the fitting for each deal value 
        vector<double> errors;
 };

class  RADAR_DLL PlainRegression: public Regression{
public:
  PlainRegression(){}
  //computes the coefficients for the regression using  a container of basis functions
  void DoTheRegression(vector<vector<double> > & basisValues, vector<double> & dealValues, vector<double> & coeffs);
  //write the state to a file
  virtual void writeState(std::ostream & ostr){
    ostr << "GAUSS_REGRESSION" << endl;
  }
  //read the state from a file
  virtual void readState(std::istream & istr){}
};

class RADAR_DLL  SVDRegression: public Regression{
public:
  SVDRegression(){}
  //computes the coefficients for the regression using  a container of basis functions
  void DoTheRegression(vector<vector<double> > & basisValues, vector<double> & dealValues, vector<double> & coeffs);
  //write the state to a file
  virtual void writeState(std::ostream & ostr){
    ostr << "SVD_REGRESSION" << endl;
  }
  //read the state from a file
  virtual void readState(std::istream & istr){}
};

DRLIB_END_NAMESPACE
#endif
