//----------------------------------------------------------------------------
//
//   Group       : Credit QRD
//
//   Filename    : FAMultiQDistribution.cpp
//
//   Description : Helper class 
//                
//   Author      : 
//
//   Date        : 
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_FAMULTIQDISTRIBUTION_CPP

#include "edginc/FAMultiQDistribution.hpp"
#include "edginc/Addin.hpp"  
#include "edginc/mathlib.hpp"
#include "edginc/Black.hpp"
#include "edginc/RootFinderND.hpp"
#include "edginc/Integrator.hpp"

DRLIB_BEGIN_NAMESPACE

#define FAMULTIQ_BDRY 1E-2
#define FAMULTIQ_TINY 1E-6
#define FAMULTIQ_LARGE_YIELD 0.5
#define FAMULTIQ_MAX_STEPS 500

/******************************************************
********* template class MemberMehtod1DDouble ****************
*******************************************************/

/*======================================================
  Implementing operator()
  =======================================================*/
template<class T> double MemberMethod1DDouble<T>::operator() (double x) const
{
   static const string method = "MemberMethod1DDouble::operator()";

   //take care of cases of null pointers.
  //  if (!(fa->*psi)) throw ModelException(method,
// 				      "No function defined");
  return (fa->*psi)(x);

}



/*********************************************************
********* Class FAMultiQDistribution **************************
**********************************************************/

/*=======================================================
  Constructor 
  =======================================================*/
  FAMultiQDistribution::FAMultiQDistribution(MultiQDistributionSP mq,
					     double convexityAlpha,
					     double convexityPower,
					     double convexityRecovery,
					     double convexityRiskFreeRate,
					     double convexityTime,
					     double delayAlpha,
					     double delayPower,
					     double delayTime,
					     int RichardsonPolyOrder) :
    CObject(TYPE), mq(mq),convexityAlpha(convexityAlpha),
    convexityPower(convexityPower),convexityRecovery(convexityRecovery),
    convexityRiskFreeRate(convexityRiskFreeRate),convexityTime(convexityTime),
    delayAlpha(delayAlpha),delayPower(delayPower),delayTime(delayTime),calibrated(false),
    RichardsonPolyOrder(RichardsonPolyOrder) 
	{
	lb = mq->qMap(-5.0 * mq->getNormalVol());
	ub= mq->qMap(5.0 * mq->getNormalVol());
	}



/*=======================================================================
  Probability density function of the mapped value: P(y<=Y<=y+dy) = fa.pdf(y) dy.
  Pdf(y) is not normalized yet.
  ========================================================================*/
double FAMultiQDistribution::pdf(double y) const{

    static const char* method = "FAMutiQDistribution::pdf";
    double yCutoff = FAMULTIQ_BDRY * mq->forward();
    double yCut = (y > Maths::max(FAMULTIQ_TINY,yCutoff) ? y : yCutoff);
    
    //Calculate delay adjustment
    double delayLambda;
    double delayAdjustment;
    if (y < FAMULTIQ_TINY)
      {
	delayAdjustment = 1.0;
      }else{
     delayLambda = delayAlpha * pow(yCut,delayPower);
     delayAdjustment = pow(1+delayLambda,-delayTime);
     double  p = Maths::min(y/yCutoff,1.0);
     delayAdjustment = p * delayAdjustment + (1-p);
    }

    //Calculate convexity adjustment
  
    double convexityLambda = convexityAlpha * pow(yCut,convexityPower);
    double convexityAdjustment = (1-convexityRecovery) / yCut
      * convexityLambda /(convexityLambda + convexityRiskFreeRate)
      * (1- pow((1+convexityLambda)*(1+convexityRiskFreeRate),-convexityTime));

    //Return adjusted distribution. Not normalized.
    return mq->pdf(y)* delayAdjustment / convexityAdjustment;
}

/*=================================================================
  Function y * pdf(y)
  ===================================================================*/
double FAMultiQDistribution::xpdf(double y) const {
  static const char* method = "FAMutiQDistribution::pdf";
  return y * pdf(y);
}

/*==============================================================
  Normalize the pdf.
  Need to divide by this number each time you integrate.
  ===============================================================*/
double FAMultiQDistribution::norm() const {
   static const string method = "FAMultiQDistribution::norm";
   double sum;
  OpenRomberg1D integ(FAMULTIQ_TINY,FAMULTIQ_MAX_STEPS,RichardsonPolyOrder);
  //  Infinity inf(Infinity::Sign::Plus);
OpenBoundary inf(ub);
  OpenBoundary lowerBound(lb);
  Range r(lowerBound,inf);

  FAMultiQDistributionFunction const integrand(&FAMultiQDistribution::pdf,this,r) ;

    try{
      sum = integ.integrate(integrand);
    }catch(exception& e) {
      throw ModelException(&e,method);
    }
    
  calibrated = true;
  Norm = sum;
  if (Norm <= FAMULTIQ_TINY) throw ModelException(method,
						  "Norm of FA distribution is negative or too small");

  return sum;
}

/*====================================================================
  Returns the expected value of the distribution - probably the
  forward price or spread 
  =====================================================================*/
double FAMultiQDistribution::forward() const {
   static const string method = "FAMultiQDistribution::forward";
   double sum; 
  OpenRomberg1D integ(FAMULTIQ_TINY,FAMULTIQ_MAX_STEPS,RichardsonPolyOrder);
  //  Infinity inf(Infinity::Sign::Plus);
  OpenBoundary inf(ub);
  OpenBoundary lowerBound(lb);
  Range r(lowerBound,inf);

  FAMultiQDistributionFunction const integrand(&FAMultiQDistribution::xpdf,this,r) ;
    try{
 sum = integ.integrate(integrand);
    }catch(exception& e) {
      throw ModelException(&e,method);
    }
    return sum / (calibrated? Norm : norm());
}

/*====================================================================
  Calculates option price according to new normalized distribution
  ======================================================================*/
double FAMultiQDistribution::vanillaOptionPrice (
						 bool isCall, /**<true if call, else put*/
						 double         strike  /**<Option strike price */
						 ) const {

   static const string method = "FAMultiQDistribution::vanillaOptionPrice";
   double sum1;
   double sum2;
  OpenRomberg1D integ(FAMULTIQ_TINY,FAMULTIQ_MAX_STEPS,RichardsonPolyOrder);
 
  //Infinity infBound(Infinity::Sign::Plus);
  OpenBoundary infBound(ub);
  OpenBoundary zeroBound(lb);
  OpenBoundary strikeBound(strike);
 
   Range callRange(strikeBound,infBound);
   Range putRange(zeroBound,strikeBound);
   Range r(putRange);

   if(isCall)  r = callRange;

   //if the range is empty, the option price is nil.
   if (!r.isNonEmpty()) return 0.0;

   try{
     FAMultiQDistributionFunction const integrand1(&FAMultiQDistribution::pdf,this,r) ;
   sum1 = integ.integrate(integrand1);
  
   FAMultiQDistributionFunction const integrand2(&FAMultiQDistribution::xpdf,this,r) ;
    sum2 = integ.integrate(integrand2);
    }catch(exception& e) {
      throw ModelException(&e,method);
    }
    
   return (sum2 - strike * sum1)* (isCall? 1 : -1 ) / (calibrated? Norm : norm() );
  
  
    }

/*======================================================================
  Purpose is: price a generic payoff expressed as a function of y.
  Not yet implemented
  ==================================================================*/
// double FAMultiQDistribution::payoffPrice (Function1DDouble payoff)const {
//   return 1;
// }


/**Destructor*/
FAMultiQDistribution::~FAMultiQDistribution(){};


/****************************************************************
********* Analytics End Here. ************************************
*******************************************************************
******************************************************************/



/*=============================================================================
 * Reflection static load method
 *===========================================================================*/
void FAMultiQDistribution::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(FAMultiQDistribution, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultFAMultiQDistribution);

    /* INTERFACE FIELDS */
    FIELD(mq,   "MultiQ Distribution: expresses distribution under annuity measure");
    FIELD(convexityAlpha,    " ");
    FIELD(convexityPower,        " ");
    FIELD(convexityRecovery,      " ");
    FIELD(convexityRiskFreeRate,       " ");
    FIELD(convexityTime,        " ");
    FIELD(delayAlpha," ");
    FIELD(delayPower," ");
    FIELD(delayTime," ");
    FIELD(RichardsonPolyOrder,"Parameter for the Romberg numerical integrator");

	/*Transient fields*/
	FIELD(ub,"upper bound for integration");
	FIELD(lb,"lower bound for integration");
	FIELD(calibrated,"boolean whether the norm has been calculated");
	FIELD(Norm,"Norm of the distribution");
	FIELD_MAKE_TRANSIENT(ub);
	FIELD_MAKE_TRANSIENT(lb);
	FIELD_MAKE_TRANSIENT(Norm);
	FIELD_MAKE_TRANSIENT(calibrated);
};

/*=============================================================================
 * Register class loader
 *===========================================================================*/
CClassConstSP const FAMultiQDistribution::TYPE = 
    CClass::registerClassLoadMethod("FAMultiQDistribution", typeid(FAMultiQDistribution), load);

// bool FAMultiQDistributionLoad(){
//   return (FAMultiQDistribution::TYPE !=0); 
// }

/**Default Constructor*/
FAMultiQDistribution::FAMultiQDistribution() : CObject(TYPE), calibrated(false) {};



/*=============================================================================
 * Reflection callback for default constructor.
 *===========================================================================*/
IObject* FAMultiQDistribution::defaultFAMultiQDistribution() {return new FAMultiQDistribution();};




/*=============================================================================
 * Class to provide add-in functions - values and option prices
 *===========================================================================*/
class FAMultiQDistributionAddin2 : public CObject {
public:
    static CClassConstSP const TYPE;

    // addin parameters
    FAMultiQDistributionSP faMq;
    double value;
    string callOrPut;

  

    double pdf() {
        return faMq->pdf(value);
    };
  
  double norm() {
    return faMq->norm();
  };

  double forward() {
    return faMq->forward();
  }

    double vanillaPrice() {
        bool isCall = false;
        if (callOrPut.size()<1) 
            throw ModelException("FAMultiQDistributionAddin2::vanillaPrice",
            "callOrPut may not be blank and must begin with 'C' for call or 'P' for put.");
        char cp = toupper(callOrPut[0]);
        if (cp!='C' && cp!='P')
            throw ModelException("FAMultiQDistributionAddin2::vanillaPrice",
            "callOrPut must begin with 'C' for call or 'P' for put.");

        if (cp=='C') isCall = true;

        return faMq->vanillaOptionPrice(isCall, value);
    };

   //  double binaryPrice() {
//         bool isCall = false;
//         if (callOrPut.size()<1) 
//             throw ModelException("MultiQDistributionAddin::binaryPrice",
//             "callOrPut may not be blank and must begin with 'C' for call or 'P' for put.");
//         char cp = toupper(callOrPut[0]);
//         if (cp!='C' && cp!='P')
//             throw ModelException("MultiQDistributionAddin::binaryPrice",
//             "callOrPut must begin with 'C' for call or 'P' for put.");

//         if (cp=='C') isCall = true;

//         return mqDist->binaryOptionPrice(isCall, value);
//     }
    
    FAMultiQDistributionAddin2() : CObject(TYPE), value(0), callOrPut("") {};

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(FAMultiQDistributionAddin2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultFAMultiQDistributionAddin2);
        FIELD(faMq, "Forward Adjusted MultiQDistribution object.");
        FIELD(value,"Input value or strike for addin functions");
        FIELD(callOrPut,"Should begin with 'C' for call or 'P' for put.");
        FIELD_MAKE_OPTIONAL(callOrPut);

       

        Addin::registerDoubleMethod("FAMULTIQ_VANILLA_PRICE",
            Addin::UTILITIES,
            "Prices a vanilla call or put with given strike under the forward adjusted measure.",
            &FAMultiQDistributionAddin2::vanillaPrice);

        // Addin::registerDoubleMethod("MULTIQ_BINARY_PRICE",
//             Addin::UTILITIES,
//             "Prices a binary call or put with given strike.",
//             &MultiQDistributionAddin2::binaryPrice);

//         Addin::registerDoubleMethod("MULTIQ_CDF",
//             Addin::UTILITIES,
//             "Returns the cumulative probability density of a given value.",
//             &MultiQDistributionAddin2::cdf);

        Addin::registerDoubleMethod("FAMULTIQ_PDF",
            Addin::UTILITIES,
            "Returns the non-normalized probability density function at a given value.",
            &FAMultiQDistributionAddin2::pdf);

	  Addin::registerDoubleMethod("FAMULTIQ_NORM",
            Addin::UTILITIES,
            "Returns the norm of the (not normalized) probability density function",
            &FAMultiQDistributionAddin2::norm);
	  Addin::registerDoubleMethod("FAMULTIQ_FORWARD",
				      Addin::UTILITIES,
				      "Returns the expectation of the probability density function",
				      &FAMultiQDistributionAddin2::forward);
    }

    static IObject* defaultFAMultiQDistributionAddin2() {
        return new FAMultiQDistributionAddin2();
    }
};

CClassConstSP const FAMultiQDistributionAddin2::TYPE = CClass::registerClassLoadMethod(
    "FAMultiQDistributionAddin2", typeid(FAMultiQDistributionAddin2), FAMultiQDistributionAddin2::load);

bool FAMultiQDistributionLinkIn() {
    return FAMultiQDistribution::TYPE!=0 && FAMultiQDistributionAddin2::TYPE!=0;
}



DRLIB_END_NAMESPACE
