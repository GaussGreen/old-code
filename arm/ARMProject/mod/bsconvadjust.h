/*
 *
 * Copyright (c) CDC IXIS CM September 2003 Paris
 *
 *
 * $Log: bsconvadjust.h,v $
 * Revision 1.13  2004/03/24 13:44:20  rguillemot
 * SUMMIT Payment Lag
 *
 * Revision 1.12  2004/02/09 08:52:40  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.11  2004/01/20 13:50:10  rguillemot
 * Replication Conbvexity Adjustment
 *
 * Revision 1.10  2003/12/29 10:38:06  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.9  2003/12/29 07:38:39  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.8  2003/12/23 16:55:10  rguillemot
 * Conv Adj Merge
 *
 * Revision 1.7  2003/12/18 08:44:39  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.6  2003/12/11 14:58:38  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.5  2003/12/09 08:29:46  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.4  2003/12/04 14:41:51  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.3  2003/11/28 07:03:53  rguillemot
 * Transfert ajustement de convexite inflation
 *
 * Revision 1.2  2003/11/03 17:27:24  emezzine
 * conv2unix
 *
 * Revision 1.2  2003/10/22 13:07:47  emezzine
 * First version
 *
 */
/*----------------------------------------------------------------------------*/
/*! \file bsconvadjust.h
 *
 *  \brief B&S Convexity adjustment 
 *
 *	\author  El Mostafa EZZINE
 *	\version 1.0
 *	\date October 2003
 */
/*----------------------------------------------------------------------------*/

#ifndef _ARMBSCONVADJUST_H
#define _ARMBSCONVADJUST_H

/*! \class   ARM_BSConvAdjust
 *	\brief  object that manage convexity and payment lag for alla model
 *	
 *	\author  El Mostafa EZZINE
 *	\version 1.0
 */


#include "convadjustmanager.h"
#include "linalg.h"
#include "model.h"

using std::string;



class ARM_BSConvAdjust : public ARM_ConvAdjustManager  
{
    private:

	    ARM_Model* itsUsedModel;

	    int itsSUMMITFormulaeUsed;
	    int itsUseSabrCMS;

    public:
	
        // Init Method
	    void Init(void);

	    ARM_BSConvAdjust(int itsSMUMITFormulaeUsed = 0, int itsUseSabrCMS = 0);// Mode: K_CONV_ADJ_EXP (ARM formulae)
        ARM_BSConvAdjust(const ARM_BSConvAdjust& rhs);
        ARM_BSConvAdjust& operator= (const ARM_BSConvAdjust& rhs);


virtual ~ARM_BSConvAdjust(void);
        
	// return the convexity adjustment of Libor paid in Arrears or in advance
virtual double ConvAdjustLibor(ARM_Model* Model, 
                               const AdjustLiborData& Input);

	// compute the fwd volatility used in the libor adjustment
virtual double ComputeFwdVolatility(ARM_Model* Model, 
                                    const AdjustLiborData& Input);
   
	// return the natural convexity adjustment of CMS 
virtual double NaturalAdjstCMS(ARM_Model* Model, 
                               const NaturalAdjustData& Input, 
                               StoreFwdRateInfo* StoreInfo = NULL);

	// compute the cms volatility used in the cms adjustment
virtual double ComputeCMSVolatility(ARM_Model* Model, const NaturalAdjustData& Input, 
                                    StoreFwdRateInfo* StoreInfo = NULL);
	
	// return the convexity adjustment of CMS paid in Arrears or in advance
virtual double TimeLagAdjustCMS(ARM_Model* Model, 
                                const TimeLagAdjustData& Input, 
                                StoreFwdRateInfo* StoreInfo = NULL);
	
	// return the forward CPI ratio convexity adjustment
virtual double FwdCPIRatioAdjust(ARM_Model* Model, const ARM_Vector* Input);

	
	// compute the cms volatility implied for an CMS oplet Price
	// virtual double ComputeCMSOpletVolatility(ARM_Model* Model, const NaturalAdjustData& Input, StoreFwdRateAndCapletInfo* StoreInfo = NULL);

	// return UseOrNot SabrCMS
virtual int GetUseSabrCMS(void)
        {
		    return(itsUseSabrCMS);          
        }

virtual int GetSUMMITFormulaUsed(void)
{
	return itsSUMMITFormulaeUsed;
}


	// Standard ARM object support
    ARM_Object* Clone(void);        
    
    void View(char* id = NULL, FILE* ficOut = NULL);
};

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/