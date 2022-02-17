/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: VanillaCapArgSFRM.h,v $
 * Revision 1.1  2004/03/31 07:51:19  ebenhamou
 * Initial revision
 *
 */

/*! \file VanillaSwaptionArgSmiledFRM.h
 *
 *  \brief this files enables to store various data to fast calibrate
 *      swaption for the SFRM model
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#ifndef _INGPMODEL_VanillaSwaptionArgSmiledFRM_H
#define _INGPMODEL_VanillaSwaptionArgSmiledFRM_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/port.h"
#include "gpcalib/vanillaswaption.h"
#include "gpinfra/typedef.h"
#include <vector>
CC_USING_NS(std,vector)


CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaSwaptionArgSmiledFRM
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
class ARM_VanillaSwaptionArgSmiledFRM: public ARM_RootObject
{
public:
	ARM_VanillaSwaptionArgSmiledFRM();
	ARM_VanillaSwaptionArgSmiledFRM( double averageShift,const ARM_VectorPtr& mu,const ARM_VectorPtr& fixAnnuity,const ARM_VectorPtr& swapFwd) ;
	ARM_VanillaSwaptionArgSmiledFRM( const ARM_VanillaSwaptionArgSmiledFRM& rhs );
    virtual ~ARM_VanillaSwaptionArgSmiledFRM();

	inline void SetAverageShift( double averageShift )          { itsAverageShift= averageShift;    }
	inline void SetMu( const ARM_VectorPtr& mu )                { itsMu = mu;                       }
	inline void SetFixAnnuity( const ARM_VectorPtr& fixAnnuity) { itsFixAnnuity = fixAnnuity;       }
	inline void SetSwapFwd( const ARM_VectorPtr& swapFwd )      { itsSwapFwd= swapFwd;              }
	
	inline double  GetAverageShift( ) const				{ return itsAverageShift;   }
	inline const ARM_VectorPtr& GetMu() const			{ return itsMu;             }
	inline const ARM_VectorPtr& GetFixAnnuity( ) const	{ return itsFixAnnuity;     }
	inline const ARM_VectorPtr& GetSwapFwd( ) const		{ return itsSwapFwd;        }

	   /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_VanillaSwaptionArgSmiledFRM(*this); };
	virtual string toString(const string& indent="",const string& nextIndent="") const;

	
private :
    double itsAverageShift;
	ARM_VectorPtr itsMu;
	ARM_VectorPtr itsFixAnnuity;
	ARM_VectorPtr itsSwapFwd;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
