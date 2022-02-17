/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file genpricer.h
 *
 *  \brief the generic pricer object
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_GENPRICER_H
#define _INGPINFRA_GENPRICER_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/typedef.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration (class)
class ARM_GenSecurity;
class ARM_PricingModel;

/// struct object
struct ARM_PricerInfo;


/////////////////////////////////////////////////////////
/// \class ARM_GenPricer
/// the ARM_GenPricer contains already
///	- a pointor to the security
/// - a pointor to the model
/// beware that it does not clone the security and model
///		neither delete them!
///
/// I REPEAT
/// BEWARE THAT IT DOES NOT CLONE THE SECURITY AND MODEL
///		NEITHER DELETE THEM!
///
/////////////////////////////////////////////////////////

struct ARM_CVInfo: public ARM_RootObject
{
private:
	ARM_StringVector itsCVColumns;
	string itsRefPriceColumn;
	ARM_GP_Vector itsCVPrices;
	ARM_GP_Vector itsBetas;

public:
	ARM_CVInfo( const ARM_StringVector& cvColumns, const string& refPriceColumn, const ARM_GP_Vector& cvPrices, const ARM_GP_Vector& betas)
	:	itsCVColumns( cvColumns ), itsRefPriceColumn( refPriceColumn ), itsCVPrices( cvPrices ), itsBetas( betas )
	{}
	
	ARM_StringVector GetCVColumns() const { return itsCVColumns; }
	string GetRefPriceColumn() const { return itsRefPriceColumn; }
	ARM_GP_Vector GetCVPrices() const { return itsCVPrices; }
	ARM_GP_Vector GetBetas() const { return itsBetas; }
	void SetBetas (ARM_GP_Vector betas) {itsBetas = betas;}

	/// standard ARM Object support
	virtual ARM_Object* Clone() const { return new ARM_CVInfo(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const {return "ARM_CVInfo";}
};


class ARM_GenPricer : public ARM_RootObject
{
private:
	ARM_GenSecurity* itsGenSecurity;
	ARM_PricingModel* itsPricingMod;
	ARM_PricerInfo* itsPricerInfo;
	ARM_CVInfo* itsCVInfo;

	bool itsDetailMode;
	void CheckColumnNameExists( const string& columnName, const ARM_StringVector& secColumnNames ) const;
	void ComputeCV();

public :
	ARM_GenPricer( ARM_GenSecurity* sec, ARM_PricingModel* mod, const ARM_StringVector& cvsColumns = ARM_StringVector(), const ARM_GP_Vector& cvsPrices=ARM_GP_Vector(), const string& refPriceColumn = string(), const ARM_GP_Vector& betas = ARM_GP_Vector() );
	ARM_GenPricer(const ARM_GenPricer &pricer);
	ARM_GenPricer&operator = (const ARM_GenPricer&pricer);
	~ARM_GenPricer();

	/// pricing functions
	virtual double Price();
	ARM_PricerInfo* GetPricerInfo() const {return itsPricerInfo; }

	/// standard ARM Object support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;
	inline const ARM_CVInfo* GetCVInfo() { return itsCVInfo;}
	/// for easy debugging 
	inline bool SetDetailMode(bool value) {return itsDetailMode=value;}
	static const string CVTag;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/


