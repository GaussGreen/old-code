//----------------------------------------------------------------------------
//
//   Group       : Credit QRD
//
//   Filename    : FAMultiQDistribution.hpp
//
//   Description : Helper class 
//                 
//
//   Author      :
//
//   Date        : 
//
//
//----------------------------------------------------------------------------

#ifndef FAMULTIQDISTRIBUTION_HPP
#define FAMULTIQDISTRIBUTION_HPP

#include "edginc/AtomicArray.hpp"
#include "edginc/Object.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/MultiQDistribution.hpp"
#include "edginc/Function.hpp"
#include "edginc/Range.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(FAMultiQDistribution)

/** generic function class*/
template<class T>  class MemberMethod1DDouble: public Function1DDouble
{
public:
  //  MemberMethod1DDouble():Function1DDouble(),r(0,1,0,1){};
  MemberMethod1DDouble(double (T::*psi) (double) const,const T *const fa,Range r)
    : psi(psi),fa(fa),Function1DDouble(r),r(r){};
  //if you do not want to define a range, use this. range will default to ]-inf,+inf[
  MemberMethod1DDouble(double (T::*psi) (double) const,const T *const fa)
    : psi(psi),fa(fa),Function1DDouble(r),r(0,false,0,false){};
  
  
  virtual double operator ()(double x) const;
private:
  double (T::*psi)(double) const;
  const T *const fa;
  Range r;
  
};
typedef MemberMethod1DDouble<FAMultiQDistribution> FAMultiQDistributionFunction;

/**FAMultiQDistribution Class*/
class UTIL_DLL FAMultiQDistribution : public CObject {
public:
    static CClassConstSP const TYPE; // for reflection API
    static void load(CClassSP& clazz);

    virtual ~FAMultiQDistribution();

  /**Constructor */
  FAMultiQDistribution(MultiQDistributionSP mq,
		       double convexityAlpha,
		       double convexityPower,
		       double convexityRecovery,
		       double convexityRiskFreeRate,
		       double convexityTime,
		       double delayAlpha,
		       double delayPower,
		       double delayTime,
		       int RichardsonPolyOrder);

 /**Returns the expected value of the distribution - probably the
       forward price, or spread, depending on what you are using this
       for. */
    double forward() const;
  double norm() const ;
  /**Probability density function of the mapped value: P(y<=Y<=y+dy) = fa.pdf(y) dy.*/
    double pdf(double y) const ;


    /**Computes the price of a vanilla call or put with the M-Q distribution*/
    double vanillaOptionPrice (
        bool isCall, /**<true if call, else put*/
        double         strike  /**<Option strike price */
    ) const;

  // double payoffPrice (Function1DDouble payoff) const;
		   
protected:
  
/**********************************************************************
*************  Fields  *****************************************************************/
  MultiQDistributionSP mq;
  /**Convexity adjustment Fields*/
  double convexityAlpha;
  double convexityPower;
  double convexityRecovery;
  double convexityRiskFreeRate;
  double convexityTime;
  /**Delay Adjustment Fields*/
  double delayAlpha;
  double delayPower;
  double delayTime;
  /** Richardson Poly Order for integration:*/
  int RichardsonPolyOrder;
/**Default Constructor for Reflection interface only*/
  FAMultiQDistribution();

private:
   /**Calibrated Flag*/
  mutable bool calibrated;
  mutable double Norm;
  mutable double lb;
  mutable double ub;
  
  double xpdf(double y) const ;
  static IObject* defaultFAMultiQDistribution();

  friend class FAMultiQDistributionAddin2;
  friend class Q3MQQuasiPricer;
};

typedef smartPtr<FAMultiQDistribution> FAMultiQDistributionSP;

DRLIB_END_NAMESPACE

#endif
