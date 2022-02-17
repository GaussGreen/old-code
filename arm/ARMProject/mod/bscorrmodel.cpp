/*
 * $Log: bscorrmodel.cpp,v $
 * Revision 1.7  2004/01/14 09:34:03  ebenhamou
 * added more validation for negative fwd in 2LOG
 *
 * Revision 1.4  2003/11/26 16:46:20  sgasquet
 * Ajout exception constructeur
 *
 * Revision 1.3  2003/02/11 16:47:46  mab
 * Added : #include "foncstd.h"
 *
 * Revision 1.2  2002/11/13 10:35:56  mab
 * Formatting
 *
 * Revision 1.1  2002/05/23 15:53:22  sgasquet
 * Initial revision
 *
 */


/*----------------------------------------------------------------------------*
     bscorrmodel.cpp
 
    This file implements the ARM_BSCorrModel class

    (The Black-Sholes Model)

*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "armglob.h"
#include "bscorrmodel.h"
#include "volflat.h"
#include "zeroflat.h"
#include "foncstd.h"




/*---------------------------------------------------------------------------*/
/*  Set                                                                      */
/*---------------------------------------------------------------------------*/

void ARM_BSCorrModel::Set(ARM_VolCurve* capCashVol, // CASH Vol or 
                                                    // Gaussian Spread Vol
                            ARM_VolCurve* SpreadVol,// Spread Vol
                            ARM_VolCurve* Correlation, // Correlations
                            int ModelType,
                            int VolType)
{
/*    if (itsCashVolMatrix)
       delete itsCashVolMatrix;
 
    if (capCashVol)
    {
       itsCashVolMatrix = (ARM_VolCurve *) capCashVol->Clone();
    }
    else 
       itsCashVolMatrix = NULL;
*/
    SetCashVolMatrix(capCashVol); 
  
	SetCorrelationMatrix(Correlation);

	SetModelType(ModelType);

	SetSpreadVolCurve(SpreadVol);

	SetSpreadVolType(VolType);
	

/*    if (itsSpreadVolCurve)
       delete itsSpreadVolCurve;
 
    if (SpreadVol)
    {
       itsSpreadVolCurve = (ARM_VolCurve *) SpreadVol->Clone();
    }
    else 
       itsSpreadVolCurve = NULL;

    if (itsCorrelationMatrix)
       delete itsCorrelationMatrix;
 
    if (Correlation)
    {
       itsCorrelationMatrix = (ARM_VolCurve *) Correlation->Clone();
    }
    else 
       itsCorrelationMatrix = NULL;

    itsModelType = ModelType;
    itsVolType = VolType;*/
}
        


// ARM_BSCORRMODEL
/*---------------------------------------------------------------------------*/
/*   Constructor                                                             */
/*---------------------------------------------------------------------------*/

ARM_BSCorrModel::ARM_BSCorrModel(ARM_Date& startDate,
                                 ARM_ZeroCurve* zeroCurve, 
                                 ARM_VolCurve* spreadLock, // SpreadLock Vol
                                 ARM_VolCurve* capIRGVol,  // IRG Vol
                                 ARM_VolCurve* capCashVol, // CASH Vol
                                 ARM_VolCurve* indexVAdjol,// Adjustment Vol
                                 ARM_VolCurve* SpreadVol,  // Spread Vol
                                 ARM_VolCurve* Correlation,// Correlations
                                 int ModelType,
                                 int VolType )
                                :ARM_BSModel(startDate, zeroCurve,
                                  spreadLock, capIRGVol, indexVAdjol)
{
    Init();

    SetName(ARM_BSMODEL);

	if (capIRGVol->GetName() != ARM_VOL_CUBE)
	{
		if (capIRGVol->GetVolType() != K_ATMF_VOL)
		  throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
             "Vol Cube or ATM Curve expected");                
	}

	if (capCashVol->GetName() != ARM_VOL_CUBE)
	{
		if (capCashVol->GetVolType() != K_ATMF_VOL)
		  throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
             "Vol Cube or ATM Curve expected");                
	}

    Set(capCashVol, SpreadVol, Correlation, ModelType, VolType);
}



double ARM_BSCorrModel::SpreadOptionPrice(double fwd1, double fwd2, 
     double vol1, 
     double vol2, double Correl, 
     double strike, 
     double optMat, int optType, int ComputedFormula /* 0 formule d'Olivier, 1 Formule Classique*/)
{
	double SpreadVol=-1;
	if(GetSpreadVolType() == K_INPUTED )
		SpreadVol = GetSpreadVolCurve()->ComputeVolatility(optMat,strike)*0.01;
    return ARM_BSCorrModel::SpreadOptionFormula(fwd1,fwd2,vol1,vol2,Correl,strike,optMat,optType,GetModelType(),SpreadVol, ComputedFormula);
}

/// static method to allow function call directly from the interface
// fwd and strike are in % (vol and correl have been converted into real values)
double ARM_BSCorrModel::SpreadOptionFormula(double fwd1, double fwd2, 
     double vol1, 
     double vol2, double Correl, 
     double strike, 
     double optMat, int optType, int modelType, double spreadVol, int ComputedFormula /* 0 formule d'Olivier, 1 Formule Classique*/)
{
    double price;
	
	/// compare to one day!
	if( optMat< 0.0001) // one day 1.0/365.0 = 0.0027393
		throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
			"Trying to compute an option with a time shorter than one day!");

	// validation on the lognormal model!
	if(modelType == K_LOG )
	{
		if((fwd2-fwd1)*0.01 <= K_DOUBLE_TOL)
			throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
				"Negative spread: cannot be lognormal");
		if ( strike*0.01<=K_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
				"Negative strike makes no sense with lognormal model!");                
	}

	if ( modelType == K_2LOG )        // 2 Log diffusions for each fwd
	{
		if(fwd1*0.01 <= K_DOUBLE_TOL)
			throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
				"fwd1 cannot be nagative under lognormal assumtion");

		if(fwd2*0.01 <= K_DOUBLE_TOL)
			throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
				"fwd2 cannot be nagative under lognormal assumtion");

		price = 100.0 * SpreadOption(fwd1*0.01, fwd2*0.01, vol1, 
			vol2, 0.0, 0.0, Correl, 0.0, strike*0.01, optMat, optType, ComputedFormula); 
    }
    else
    {
		// we look if we need to compute implied vol for spread
		if( spreadVol == -1 )
		{
			spreadVol = sqrt(pow(vol1*fwd1*0.01, 2)
				+pow(vol2*fwd2*0.01, 2)
				-2.0*vol1*vol2*fwd1*0.01*fwd2*0.01*Correl);

			if( modelType == K_LOG )
				spreadVol /= (fwd2-fwd1)*0.01;

			if(spreadVol<K_DOUBLE_TOL)
				throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
					"Implied vol from lognormal vols is negative");
        }
		
        if ( modelType == K_NOR ) // 1 gaussian diffusions for the spread
        {
			price = 100.0 * GaussianSpreadOption(fwd1*0.01, fwd2*0.01, 
				spreadVol, 0.0, 0.0, Correl, 0.0, strike*0.01, optMat, optType); 
        }
        else if ( modelType == K_LOG ) // 1 Log diffusions for the spread
        {
			price = 100.0 * ::bsOption((fwd2-fwd1)*0.01, strike*0.01, 
				spreadVol, 0.0, 0.0, optMat, optType); 
				
        }
    }

    return(price);
}



 
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
