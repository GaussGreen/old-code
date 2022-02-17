/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanillaarg.h
 *
 *  \brief this files enables to convert 
 *      security of the arm kernel to gp priceables
 *      choice at this stage is to plug this directly into
 *      the pricing model pricing function. For a more general 
 *      approach, we may want to create from a arm kernel security
 *      a deal description!
 *
 *
 *	\author  E.M Ezzine, E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */


#ifndef _INGPCALIB_VANILLAARG_H
#define _INGPCALIB_VANILLAARG_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"

#include "gpcalib/typedef.h"

/// gpinfra
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/gplinalgtypedef.h"

#include <string>
CC_USING_NS(std,string)

/// forward declaration
/// kernel object in the global namespace!
class ARM_Security;

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_PricingModel;
struct ARM_ConverterFromKernel;

///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaArg
/// \brief base class for vanilla instrument
///////////////////////////////////////////////////////////////
struct ARM_VanillaArg: public ARM_RootObject
{
	/// use for efficiency reason instread of a dynamic cast!
	enum VanillaArgType
	{
		VANILLA_CAP,
		VANILLA_CORRIDOR,
		VANILLA_DIGITAL,
		VANILLA_EQUITYOPTION,
		VANILLA_FXOPTION,
		VANILLA_SWAPTION,
		VANILLA_IRLEG,
		VANILLA_IRSWAP,
		VANILLA_IRFIXEDLEG,
		VANILLA_INFLEG,
		VANILLA_INFCAP,
		VANILLA_PORTFOLIO,
		VANILLA_SPREADOPTION,
		VANILLA_SUMOPT,
		VANILLA_MEPI,
		VANILLA_IRFXSWAPTION,
		VANILLA_STRIP
	};


	ARM_VanillaArg( 
		const string& curveName = string(), 
		double evalTime = 0.0, 
		int CallPut		= K_CALL,
        double expiry	= 0.0,
		double tenor	= 0.0,
		double spread	= 0.0,
		int decompFreq	= 1)
	:	
		ARM_RootObject(),
		itsCurveName(curveName),
		itsEvalTime(evalTime),
		itsCallPut(CallPut),
		itsMktPrice(-1111111111),
		itsIndex(    1111111111),
        itsExpiry(expiry),
		itsTenor(tenor),
		itsSpread(spread),
		itsDecompFreq(decompFreq)
	{}

    ARM_VanillaArg(const ARM_VanillaArg& arg);
    ARM_VanillaArg& operator=(const ARM_VanillaArg& rhs);
	virtual ~ARM_VanillaArg();
	virtual ARM_CLASS_NAME GetRootName() { return ARM_GP_VANILLAARG; }


	/// the important method!
	virtual double Price(ARM_PricingModel* model) const = 0;
    virtual double ImpliedVol(ARM_PricingModel* model) const = 0;
	virtual double Derivative(ARM_PricingModel* model, 
        const ARM_ModelParam& modelParam, 
        size_t xIndex, 
		size_t yIndex = 0,
		size_t zIndex = 0,
		int factorNb =  1,
        ARM_MktTargetType targetType = ARM_CalibrationTarget::PriceTarget) const;
	
private:
// FIXMEFRED: mig.vc8 (22/05/2007 18:04:38):missing return type
	void CopyNoCleanUp(const ARM_VanillaArg& rhs);
    double NumericDerivative(ARM_PricingModel* model, 
        ARM_ModelParam* modelParam, 
        size_t xIndex, 
		size_t yIndex = 0,
		size_t zIndex = 0,
		int factorNb =  1,
        ARM_MktTargetType targetType = ARM_CalibrationTarget::PriceTarget) const;

    int     itsCallPut;		///capfloor for cap/floor and RecPay for swaption.
    string  itsCurveName;
    double  itsEvalTime;
    double  itsMktPrice;    /// used only to calibrate a portfolio
    int     itsIndex;       /// used only to calibrate a portfolio
    double	itsExpiry;		/// used only to calibrate a portfolio
    double	itsTenor;		/// used only to calibrate a portfolio
	double	itsSpread;
	int		itsDecompFreq;

public:
	/// Set accessors
	inline void SetMktPrice( double mktPrice ) { itsMktPrice = mktPrice; }
	inline void SetIndex( int index ) { itsIndex = index; }
    inline void SetExpiry( double expiry) { itsExpiry = expiry; }
    inline void SetTenor( double tenor) { itsTenor= tenor; }
	inline void SetSpread( double spread) { itsSpread= spread; }
	inline void SetCurveName( const string& curveName ) { itsCurveName=curveName; }
	
	/// Get accessors
	inline string GetCurveName() const { return itsCurveName; }
	inline double GetEvalTime( ) const { return itsEvalTime; }
	inline int	  GetCallPut() const { return itsCallPut; }
	inline int    GetIndex( ) const { return itsIndex; }
    inline double GetExpiry( ) const { return itsExpiry; }
    inline double GetTenor( ) const { return itsTenor; }
    virtual double GetMktPrice( ) const { return itsMktPrice; }
	inline double GetSpread( ) const { return itsSpread; }

	/// method for validation
	void ThrowErrorOnNullObject( const string& objName, ARM_Object* obj ) const;
	void ThrowErrorOnEmptyVector( const string& vectorName, const ARM_VectorVector& vec ) const;

	/// use for efficiency reason instread of a dynamic cast!
	virtual VanillaArgType GetType() const = 0;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
