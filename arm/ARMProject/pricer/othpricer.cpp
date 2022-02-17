/*
 * $Log: othpricer.cpp,v $
 * Revision 1.3  2003/12/02 19:18:43  mab
 * SetModel replaced by : SetModelVariable(..)
 *
 * Revision 1.2  2003/08/20 16:21:26  sgasquet
 * Forcage prise en compte du modele grace a SetModel(NULL) dans la fonction Price
 *
 */


/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : othpricer.cpp                                                */
/*                                                                            */
/* DESCRIPTION : classe pricers pour compatibilite avec existant              */
/*                                                                            */
/* DATE        : Fri Sep 11 1998                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/



#include "othpricer.h"

#include "bsmodel.h"







double ARM_OtherPricer::Price(void)
{
    itsSecurity->SetModelVariable(NULL);

    if ( itsModel->GetName() == ARM_BSPRICINGMODEL )
    {
       ARM_Model* theBSPricingModel = ((ARM_BSPricingModel *)
                                      itsModel)->GetGeneratedRealBSPricingModel(itsSecurity);

       itsSecurity->SetModel(theBSPricingModel);
    }
    else
    {
	   itsSecurity->SetModel(itsModel);
    }

    double thePrice = itsSecurity->ComputePrice();

	return(thePrice);
}

