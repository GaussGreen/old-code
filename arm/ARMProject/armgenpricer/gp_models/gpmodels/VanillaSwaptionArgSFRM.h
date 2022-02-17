/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: VanillaCapArgSFRM.h,v $
 * Revision 1.1  2004/03/31 07:51:19  ebenhamou
 * Initial revision
 *
 */

/*! \file VanillaSwaptionArgSFRM.h
 *
 *  \brief this files enables to store various data to fast calibrate
 *      swaption for the SFRM model
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#ifndef _INGPMODEL_VANILLASWAPTIONARGSFRM_H
#define _INGPMODEL_VANILLASWAPTIONARGSFRM_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/port.h"
#include "gpcalib/vanillaswaption.h"
#include "gpinfra/typedef.h"
#include <vector>
CC_USING_NS(std,vector)


CC_BEGIN_NAMESPACE( ARM )

class ARM_SFRM;


///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaSwaptionArgSFRM
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaSwaptionArgSFRM: public ARM_VanillaSwaptionArg
{
	ARM_VanillaSwaptionArgSFRM(
		double resetTime,
		double startTime,
		double endTime,
		double averageShift,
		std::vector<double>& fixPayTimes,
		std::vector<double>& fixPayPeriods,
		std::vector<double>& floatResetTimes,
		std::vector<double>& floatStartTimes,
		std::vector<double>& floatEndTimes,
		std::vector<double>& floatIntTerms,
		const ARM_VectorPtr& mu,
		const ARM_VectorPtr& fixAnnuity,
		const ARM_VectorPtr& swapFwd,
        int fixFrequency,
		const ARM_VectorPtr& fixAnnuityWithNominal = ARM_VectorPtr(NULL),
		const ARM_VectorPtr& swapFwdWithNominal = ARM_VectorPtr(NULL),
		std::vector<double>& fixNominal = NULL,
		std::vector<double>& floatNominal = NULL );

	ARM_VanillaSwaptionArgSFRM( const ARM_VanillaSwaptionArg& rhs );
	ARM_VanillaSwaptionArgSFRM( const ARM_VanillaSwaptionArgSFRM& rhs );
    ARM_VanillaSwaptionArgSFRM& operator=(const ARM_VanillaSwaptionArgSFRM& rhs);
// FIXMEFRED: mig.vc8 (30/05/2007 16:10:18): reference missing
	void SetSwapFwd(ARM_VectorPtr &swapFwd) {itsSwapFwd = swapFwd; }
	virtual ~ARM_VanillaSwaptionArgSFRM();
	virtual ARM_Object* Clone();

	/// accessor
	inline void SetAverageShift( double averageShift )          { itsAverageShift= averageShift;    }
	inline void SetMu( const ARM_VectorPtr& mu )                { itsMu = mu;                       }
	inline void SetFixAnnuity( const ARM_VectorPtr& fixAnnuity) { itsFixAnnuity = fixAnnuity;       }
	inline void SetSwapFwd( const ARM_VectorPtr& swapFwd )      { itsSwapFwd= swapFwd;              }
	inline void SetFixAnnuityWithNominal( const ARM_VectorPtr& fixAnnuityWithNominal)  { itsFixAnnuityWithNominal = fixAnnuityWithNominal;       }
	inline void SetSwapFwdWithNominal ( const ARM_VectorPtr& swapFwdWithNominal )      { itsSwapFwdWithNominal	  = swapFwdWithNominal;              }


	inline double  GetAverageShift( ) const     { return itsAverageShift;   }
	inline ARM_VectorPtr GetMu() const          { return itsMu;             }
	inline ARM_VectorPtr GetFixAnnuity( ) const { return itsFixAnnuity;     }
	inline ARM_VectorPtr GetSwapFwd( ) const    { return itsSwapFwd;        }
	inline ARM_VectorPtr GetFixAnnuityWithNominal( ) const { return itsFixAnnuityWithNominal;     }
	inline ARM_VectorPtr GetSwapFwdWithNominal   ( ) const { return itsSwapFwdWithNominal;        }


private :
    double itsAverageShift;
	ARM_VectorPtr itsMu;
	ARM_VectorPtr itsFixAnnuity;
	ARM_VectorPtr itsFixAnnuityWithNominal;  /// takes into account variable nominal
	ARM_VectorPtr itsSwapFwdWithNominal;  /// The Swap Fwd and NOT the SwapRate
	ARM_VectorPtr itsSwapFwd;
// FIXMEFRED: mig.vc8 (25/05/2007 15:44:34):missing return type
	void CopyNoCleanUp(const ARM_VanillaSwaptionArgSFRM& rhs);
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
