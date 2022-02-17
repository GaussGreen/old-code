/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *
 * $Log: implicitcorrel.h,v $
 * Revision 1.7  2003/12/18 14:37:58  nbelgrade
 * *** empty log message ***
 *
 * Revision 1.3  2003/11/19 17:18:02  nbelgrade
 * *** empty log message ***
 *
 * Revision 1.2  2003/11/17 15:00:04  rguillemot
 * bug Unix
 *
 * Revision 1.1  2003/11/17 11:41:46  nbelgrade
 * Initial revision
 *
 *
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file implicitcorrel.h
 *
 *  \brief computation of boundaries for correlation
 *
 *	\author  Nabyle Belgrade
 *	\version 1.0
 *	\date Novmber 2003
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPINFLATION_IMPLICITCORREL_H
#define _INGPINFLATION_IMPLICITCORREL_H

#include "gpbase/port.h"
#include "gpbase/gplinalgtypedef.h"

#include <string>
CC_USING_NS( std, string )

/// ARM namespace
CC_BEGIN_NAMESPACE( ARM )

class ARM_MarketDataValidator
{
private:
	enum tool
	{
		Case1, 
		Case2, 
		Case3
	};

	void Convert( const string& x, tool& y );
	void Interpol( const double& a, const double& b, const double& c, const double& f_a, double& f_b, const double& f_c);
	void HmgCov( ARM_GP_Vector* pZCVol,  ARM_GP_Vector* pMaturity, const int length, ARM_GP_Vector& Cov );
	void HmgVol( ARM_GP_Vector* pYtYVol, ARM_GP_Vector* pMaturity, ARM_GP_Vector* pCor, const int length,	ARM_GP_Vector& f );

public:
	void Vol_to_Cor(
		const ARM_GP_Vector& ZCVol, 
		const ARM_GP_Vector& YtYVol, 
		const ARM_GP_Vector& Maturity,
		ARM_GP_Vector& Cor );

	void YtYCor_to_ZC( 
		const ARM_GP_Vector& YtYVol, 
		const ARM_GP_Vector& Cor, 
		const ARM_GP_Vector& Maturity,
		ARM_GP_Vector& ZCVol );

	void ZCCor_to_YtY( 
		const ARM_GP_Vector& ZCVol, 
		const ARM_GP_Vector& Cor, 
		const ARM_GP_Vector& Maturity,
		ARM_GP_Vector& YtYVol );

	void Bounds( 
		ARM_GP_Vector* pZCVol,
		ARM_GP_Vector* pYtYVol,
		ARM_GP_Vector* pCor,
		ARM_GP_Vector* pMaturity,
		ARM_GP_Vector& UBound,
		ARM_GP_Vector& LBound,
		const string& choice,
		string& TBound );

	void HmgYtYCor_to_ZC( 
		ARM_GP_Vector* pYtYVol, 
		ARM_GP_Vector* pCor, 
		ARM_GP_Vector* pMaturity,
		const int length,
		ARM_GP_Vector& ZCVol );

	void HmgVol_to_Cor(
		ARM_GP_Vector* pZCVol, 
		ARM_GP_Vector* pYtYVol, 
		ARM_GP_Vector* pMaturity,
		const int length,
		ARM_GP_Vector& Cor );

	void HmgZCCor_to_YtY( 
		ARM_GP_Vector* pZCVol, 
		ARM_GP_Vector* pCor, 
		ARM_GP_Vector* pMaturity,
		const int length,
		ARM_GP_Vector& YtYVol );
	
};


class ARM_Convert_Vol
{
public:
	void VolYoY_to_VolSwp( 
		ARM_GP_Vector* pDFactor, 
		ARM_GP_Vector* pFwdCPI, 
		ARM_GP_Vector* pVol_DF,
		ARM_GP_Vector* pVol_YoY, 
		ARM_GP_Vector* pAvgCor, 
		ARM_GP_Vector* pDates, 
		ARM_GP_Vector* pTenors, 
		const double SwpRate,
		double &Vol_Swp);

};


CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
