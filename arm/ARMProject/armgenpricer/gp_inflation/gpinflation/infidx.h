/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *
 *
 */


/*----------------------------------------------------------------------------*/

/*! \file inflidx.h
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \brief Header for the ARM_InfIdx class
 * To keep the ARM spirit, we built an inflation index
 */

/*----------------------------------------------------------------------------*/

  
#ifndef _INGPINFLATION_INFIDX_H
#define _INGPINFLATION_INFIDX_H

#include <inst/irindex.h>					//ARM_IrIndex
#include <ccy/currency.h>					//ARM_Currency			
#include <string>
#include "gpbase/port.h"
#include <gpbase/countedptr.h>				//ARM_CountedPtr
#include "gpinflation/infcurvmodel.h"
#include "gpinflation/infcurv.h"




CC_BEGIN_NAMESPACE( ARM )

class ARM_InfIdx : public ARM_IRIndex
{
protected:
	std::string itsIndexName;
	ARM_InfCurv* itsInfFwdCurve ;
	std::string itsResetLag;
	std::string itsDCFLag;

public:
	/// constructor
	ARM_InfIdx( ARM_Currency* ccy,
				ARM_InfCurv* infCurve,
				const std::string& IndexName,
				const std::string& ResetLag = GETDEFAULTVALUESTR, 
				const std::string& DCFLag   = GETDEFAULTVALUESTR);
	
	// for backward compatibility
	ARM_InfIdx( const std::string& IndexName);


	ARM_InfIdx( const ARM_InfIdx& Index );
	ARM_InfIdx& operator=( const ARM_InfIdx& Index );
	virtual ~ARM_InfIdx()
	{
		if (itsInfFwdCurve )
			delete itsInfFwdCurve ;
	};
	
	/// Standard ARM_Object support
	virtual ARM_Object* Clone(void);
	virtual void View(char* id = NULL, FILE* ficOut= NULL);

	/// Accessors
	std::string GetResetLag() const { return itsResetLag; }
	std::string GetDCFLag() const { return itsDCFLag; }
	std::string GetIndexName() const { return itsIndexName; }
	ARM_InfCurv* GetInfFwdCurve() const {return itsInfFwdCurve ;}
	
	double FwdCPI( const ARM_Date& fixingDate );
};

CC_END_NAMESPACE()


#endif
