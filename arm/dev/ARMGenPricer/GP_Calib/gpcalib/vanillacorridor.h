/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanillacorridor.h
 *
 *  \brief vanilla corridor
 *
 *	\author  E.M Ezzine, E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */


#ifndef _INGPCALIB_VANILLACORRIDOR_H
#define _INGPCALIB_VANILLACORRIDOR_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "vanillaarg.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaCorridorLegArg
/// \brief CorridorLeg argument simple struct
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

struct ARM_VanillaCorridorLegArg: public ARM_VanillaArg
{
public:
    ARM_VanillaCorridorLegArg(
        const string& curveName             = string(), 
		double evalTime				        = 0.0, 
		int CallPut					        = K_CALL, 
	    ARM_GP_Vector* ResetTimes           = NULL, 
	    ARM_GP_Vector* StartTimes           = NULL,
        ARM_GP_Vector* EndTimes             = NULL,
        ARM_GP_Vector* PayTimes             = NULL, 
        int indexPaymentType                = IDXFIXED,
        ARM_GP_Vector* fwdPaymentPeriod     = NULL,
        ARM_VectorVector RefIdxResettimes   = ARM_VectorVector(0),
        ARM_VectorVector RefIdxStarttimes   = ARM_VectorVector(0),
        ARM_VectorVector RefIdxEndtimes     = ARM_VectorVector(0),
        ARM_VectorVector RefFwdPeriods      = ARM_VectorVector(0),
        ARM_VectorVector RefIndexWeight     = ARM_VectorVector(0),
        ARM_VectorVector DownBarrier        = ARM_VectorVector(0),
        ARM_VectorVector UpBarrier          = ARM_VectorVector(0),
        ARM_GP_Vector*  CouponMargin        = NULL,
        ARM_GP_Vector*  Nominals            = NULL)   
	:
        ARM_VanillaArg      (curveName, evalTime, CallPut ),
        itsResetTimes       (ResetTimes					   ),
        itsStartTimes       (StartTimes                    ),
        itsEndTimes         (EndTimes                      ),
        itsPayTimes         (PayTimes                      ), 
        itsIndexPaymentType (indexPaymentType              ),
        itsFwdPaymentPeriod (fwdPaymentPeriod              ),
        itsRefIdxResettimes (RefIdxResettimes              ),
        itsRefIdxStarttimes (RefIdxStarttimes              ),
        itsRefIdxEndtimes   (RefIdxEndtimes                ),
        itsRefFwdPeriods    (RefFwdPeriods                 ),
        itsRefIndexWeight   (RefIndexWeight                ),
        itsDownBarrier      (DownBarrier                   ),
        itsUpBarrier        (UpBarrier                     ),
        itsCouponMargin     (CouponMargin                  ),
        itsNominals         (Nominals                      )
        
    {}

    ARM_VanillaCorridorLegArg(const ARM_VanillaCorridorLegArg& arg);
    ARM_VanillaCorridorLegArg& operator=(const ARM_VanillaCorridorLegArg& rhs);
    virtual ~ARM_VanillaCorridorLegArg();

    virtual double Price(ARM_PricingModel* model) const;
    virtual double ImpliedVol(ARM_PricingModel* model) const;
	virtual ARM_Object* Clone() const { return new ARM_VanillaCorridorLegArg(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const;

private :
    void CopyNoCleanUp(const ARM_VanillaCorridorLegArg& rhs);
	void CleanUp();

    ARM_GP_Vector* itsResetTimes; 
	ARM_GP_Vector* itsStartTimes;
    ARM_GP_Vector* itsEndTimes;
    ARM_GP_Vector* itsPayTimes; 
    int itsIndexPaymentType;
    ARM_GP_Vector* itsFwdPaymentPeriod;

    ARM_VectorVector itsRefIdxResettimes;
    ARM_VectorVector itsRefIdxStarttimes;
    ARM_VectorVector itsRefIdxEndtimes;
    ARM_VectorVector itsRefFwdPeriods;
    ARM_VectorVector itsRefIndexWeight;

    ARM_VectorVector  itsDownBarrier;
    ARM_VectorVector  itsUpBarrier;
    ARM_GP_Vector*  itsCouponMargin;
    ARM_GP_Vector*  itsNominals;

public:
	/// accessor part
    inline const ARM_GP_Vector* GetResetTimes() const           {return itsResetTimes;       }
    inline const ARM_GP_Vector* GetStartTimes() const           {return itsStartTimes;       }
    inline const ARM_GP_Vector* GetEndTimes() const             {return itsEndTimes;         }
    inline const ARM_GP_Vector* GetPayTimes() const             {return itsPayTimes;         }
    inline const int  GetIndexPaymentType() const			    {return itsIndexPaymentType; }
    inline const ARM_GP_Vector* GetFwdPaymentPeriod() const		{return itsFwdPaymentPeriod; }

    inline const ARM_VectorVector GetRefIdxResettimes() const   {return itsRefIdxResettimes;}
    inline const ARM_VectorVector GetRefIdxStarttimes() const   {return itsRefIdxStarttimes;}
    inline const ARM_VectorVector GetRefIdxEndtimes() const     {return itsRefIdxEndtimes;  }
    inline const ARM_VectorVector GetRefFwdPeriods() const      {return itsRefFwdPeriods;   }
    inline const ARM_VectorVector GetRefIndexWeight() const     {return itsRefIndexWeight;  }

    inline const ARM_VectorVector GetDownBarrier() const		{return itsDownBarrier;		}
    inline const ARM_VectorVector GetUpBarrier() const			{return itsUpBarrier;       }
    inline const ARM_GP_Vector* GetCouponMargin() const			{return itsCouponMargin;    }
    inline const ARM_GP_Vector* GetNominals() const				{return itsNominals;        }

	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_CORRIDOR; }
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
