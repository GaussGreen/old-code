/*
 * $Log: ipricer.cpp,v $
 * Revision 1.18  2003/11/17 17:05:43  mcampet
 * MC quanto dualcap integrated
 *
 * Revision 1.17  2003/07/01 19:35:49  jpriaudel
 * suppression of include vcmarray
 *
 * Revision 1.16  2003/05/15 09:22:28  ykhlif
 *  Add computePrices for pricing option and underlyings
 *
 * Revision 1.15  2003/03/28 14:36:58  mab
 * Added : Treatement of ARM_POWERREVERSE specific case
 *
 * Revision 1.14  2003/03/10 12:15:40  mab
 * Added : Markov tree pricer
 *
 * Revision 1.13  2001/10/10 16:31:08  abizid
 * *** empty log message ***
 *
 * Revision 1.11  2001/10/02 07:56:57  abizid
 * Check cas quanto
 *
 * Revision 1.10  2001/01/31 16:22:55  smysona
 * Supreesion templatage des objets FRM
 *
 * Revision 1.9  2001/01/30 10:27:16  smysona
 * Rajout du MTreePricer
 *
 * Revision 1.8  2000/11/10 15:50:30  smysona
 * Error correction
 *
 * Revision 1.7  2000/11/09 20:55:46  smysona
 * New case SCRMONTECARLO
 *
 * Revision 1.6  2000/09/28 08:51:48  mab
 * Rajout two plus pricer
 *
 */



/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : ipricer.cpp                                                  */
/*                                                                            */
/* DESCRIPTION : Interface pricer classe                                      */
/*                                                                            */
/* DATE        : Fri Sep 11 1998                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#include "ipricer.h"
#include "othpricer.h"
#include "mcfnpricer.h"
#include "mcrnpricer.h"
#include "2dtreepricer.h"
#include "twopluspricer.h"
#include "lsmcpricer.h"
#include "mpricer.h"
#include "markovtreepricer.h"




ARM_IFPricer::ARM_IFPricer(ARM_Security *sec, ARM_Model *mod)
{

//  Test the Currencies of the Model and the Instrument to check compatibility
//  prohibit symetric models
    if (( sec->GetName() != ARM_POWERREVERSE ) && ( sec->GetName() != ARM_DUALCAP ))
    {
       if (mod->IsQuanto())
       {
          if ((strcmp(sec->GetCurrencyUnit()->GetCcyName(),
                      ARM_DEFAULT_COUNTRY)) 
              ||
             (strcmp(mod->GetZeroCurve()->GetCurrencyUnit()->GetCcyName(),
                     ARM_DEFAULT_COUNTRY)))
          {
             if (strcmp(sec->GetCurrencyUnit()->GetCcyName(),
                        mod->GetZeroCurve()->GetCurrencyUnit()->GetCcyName()))
             {
                throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
                        "Invalid Quanto Model to Price this Security!");
             }
          }
       }
    }

    switch (mod->PricerType())
    {
        case PT_NONE : 
        {
            itsPricer = new ARM_OtherPricer(sec, mod);
        };
        break;

        case PT_FNMONTECARLO :
        {
            itsPricer = new ARM_FNMonteCarloPricer(sec, mod);
        };
        break;

        case PT_RNMONTECARLO :
        {
            itsPricer = new ARM_RNMonteCarloPricer(sec, mod);

        };
        break;

        case PT_HWTREEQUANTO :
        {
            itsPricer = new ARM_2DTreePricer(sec, mod);
        };
        break;
    
        case PT_TWOPLUSTREE :
        {
            itsPricer = new ARM_TwoPlusPricer(sec, mod);
        };
        break;

        case PT_SCRMONTECARLO :
        {
             itsPricer = new ARM_LSMonteCarloPricer(sec, mod);
        };
        break;

        case PT_MTree :
        {
            switch (mod->GetDim())
            {
                case 1 :
                {
                    itsPricer = new ARM_MTreePricer1(sec, mod);
                };
                break;

                default :
                {
                   throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                            "MPricer not yet ready for this dimension !");
                }
            }
        }
        break;

        case PT_MARKOVTREE :
        {
            itsPricer = new ARM_MarkovTreePricer(sec, mod);
        };
        break;

        default :
        {
            throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Invalid Pricer Type <Constructor ARM_IFPricer> ");
        };
        break;
    }
}


double ARM_IFPricer::Price(void)
{
    return (itsPricer->Price());
}

ARM_Vector* ARM_IFPricer::ComputePrices(void)
{
    return (itsPricer->ComputePrices());
}
