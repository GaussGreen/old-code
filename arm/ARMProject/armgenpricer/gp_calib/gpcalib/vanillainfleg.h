/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file vanillairleg.h
 *
 *  \brief vanilla interest rates leg
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#ifndef _INGPCALIB_VANILLAINFLEG_H
#define _INGPCALIB_VANILLAINFLEG_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/datestrip.h"

#include "vanillaarg.h"


CC_BEGIN_NAMESPACE( ARM )


///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaInfSwaplegArg
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaInfSwaplegArg:public ARM_VanillaArg
{
    ARM_VanillaInfSwaplegArg(
		const string& curveName,
		const string& infCurveName,
		double evalTime,
		int swapType,
		int CallPut,
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		double expiry	= 0.0,
		double tenor	= 0.0,
		double Spread = 0.,
		int DecompFreq = 1);
		
    /// Constructors/Destructors
	ARM_VanillaInfSwaplegArg(const ARM_VanillaInfSwaplegArg& arg);
    ARM_VanillaInfSwaplegArg& operator=(const ARM_VanillaInfSwaplegArg& rhs);
    virtual ~ARM_VanillaInfSwaplegArg();

	/// Functions for pricing
    virtual double Price(ARM_PricingModel* model) const;
	double ImpliedVol( ARM_PricingModel * model ) const { return 0;}
	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_INFLEG; }

	/// Standard ARM stuff
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;

	/// Accessors (const and non const version)
	inline int getSwapType() const { return itsSwapType; }
	inline void setSwapType( int swapType ) { itsSwapType = swapType; }
	inline ARM_DateStripPtr getNumDateStrip() const { return itsNumDateStrip; }
	inline ARM_DateStripPtr getDenomDateStrip() const { return itsDenomDateStrip; }
	inline double getSpread() const { return itsSpread; }
	inline void setSpread( double spread ) { itsSpread = spread; }
	inline void setNumDateStrip( const ARM_DateStripPtr& numDateStrip ) { itsNumDateStrip = numDateStrip; }
	inline void setDenomDateStrip( const ARM_DateStripPtr& denomDateStrip ) { itsDenomDateStrip = denomDateStrip; }
	inline void setIRCurveName( const string& irCurveName ) { itsIRCurveName = irCurveName; }
	inline void setInfCurveName( const string& infCurveName ) { itsInfCurveName= infCurveName; }
	inline string getIRCurveName() const { return itsIRCurveName; }
	inline string getInfCurveName() const { return itsInfCurveName; }

private :
	void CopyNoCleanUp(const ARM_VanillaInfSwaplegArg& rhs);
	void CleanUp();

	string itsIRCurveName;
	string itsInfCurveName;

	ARM_DateStripPtr itsNumDateStrip;
	ARM_DateStripPtr itsDenomDateStrip;
	
	double itsSpread;
	int itsSwapType;
};

///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaInfCapArg
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaInfCapArg:public ARM_VanillaInfSwaplegArg
{
    ARM_VanillaInfCapArg(
		const string& curveName,
		const string& infCurveName,
		double evalTime,
		int swapType,
		int CallPut,
		double Strike,
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		double expiry	= 0.0,
		double tenor	= 0.0,
		double Spread = 0.,
		int DecompFreq = 1);

	/// Constructors/Destructors
	ARM_VanillaInfCapArg(const ARM_VanillaInfCapArg& arg);
    ARM_VanillaInfCapArg& operator=(const ARM_VanillaInfCapArg& rhs);
    virtual ~ARM_VanillaInfCapArg();

	/// Functions for pricing
	virtual double Price(ARM_PricingModel* model) const;
	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_INFCAP; }

	/// Standard ARM stuff
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;

	/// Accessors
	inline double getStrike() { return itsStrike; }
	inline void setStrike( double strike ) { itsStrike = strike; }

private:
	double itsStrike;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
