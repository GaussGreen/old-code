#ifndef RADAR_REP_UTIL_H
#define RADAR_REP_UTIL_H

#include <vector>

DRLIB_BEGIN_NAMESPACE

typedef long TDate;  //FIXME !

typedef std::vector<double> MarketArray;
typedef std::vector<double> FittingArray;
typedef std::vector<double> TransformedArray; // fitting array after fit. var. transformation
typedef std::vector<double> BasisValArray;
//typedef matrix<double, double> BasisDerivativeMatrix;
typedef std::vector<vector<double> > BasisDerivativeMatrix;
typedef std::vector<double> RegressionCoeff;

DRLIB_END_NAMESPACE

#endif
