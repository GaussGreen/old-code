/*
 * $Log: pricer.h,v $
 * Revision 1.10  2003/06/20 11:40:56  ebenhamou
 * remove unused var
 *
 * Revision 1.9  2003/05/15 09:21:38  ykhlif
 *  Add computePrices for pricing option and underlyings
 *
 * Revision 1.8  2003/03/07 15:50:19  mab
 * Added : SetPricingType, GetPricingType
 *
 * Revision 1.7  2001/12/12 12:01:24  arm
 * remplacement de cerr parr printf
 *
 */


/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : pricer.h                                                     */
/*                                                                            */
/* DESCRIPTION : classe mere de tout les pricers                              */
/*                                                                            */
/* DATE        : Fri Sep 11 1998                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#ifndef _PRICER_H
#define _PRICER_H

#include "security.h"
#include "model.h"
#include "time.h"


class ARM_Pricer : public ARM_Object
{
     protected :

         ARM_Security *itsSecurity;

         ARM_Model    *itsModel;

         time_t       itsBegining;


     public :

         ARM_Pricer(void)
         {
             itsSecurity = NULL;

             itsModel = NULL;

             time(&itsBegining);
         }

         ARM_Pricer(ARM_Security *sec, ARM_Model *mod);

         virtual double Price(void)
         {
             throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                  "Unimplemented <Price> method");
         }

         virtual ARM_Vector* ComputePrices(void)
         {
             throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                  "Unimplemented <Price> method");
         }


        long SummitPrice(double& pricevalue)
        {
            try
            {
                 pricevalue = Price();

                 return(RET_OK);
            }

            catch(Exception&)
            {
                printf("\nSomething caught in ARM pricing !\n");

                return(RET_KO);
            }

            catch(...)
            {
                printf("\nSomething caught in ARM pricing !\n");

                return(RET_KO);
            }

        }


        // Services

        ARM_Pricer(const ARM_Pricer &pricer)
        {
            this->BitwiseCopy(&pricer);
        }

        ARM_Pricer &operator = (const ARM_Pricer &pricer)
        {
            Copy(&pricer);

            return(*this);
        }

    
        void BitwiseCopy(const ARM_Object* opricer)
        {
             ARM_Pricer *pricer = (ARM_Pricer *)opricer;

             itsSecurity = pricer->itsSecurity;

             itsModel    = pricer->itsModel;
        }

        void Copy(const ARM_Object* src)
        {                
            ARM_Object* obj = (ARM_Pricer *) src;

            BitwiseCopy(obj);
        }

        time_t GetTime(void)
        {
            return itsBegining; 
        }

        void SetTime(void)
        {
            time(&itsBegining);
        }

        void SetTime(time_t x)
        {
            itsBegining = x;
        }

        virtual void StartPricing(void);

        virtual void EndPricing(void) {}

        virtual void SetPricingType(int ptype)
        { 
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                            "This Function has not been implemanted yet ");
        }

        virtual int GetPricingType(void) 
        { 
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                            "This Function has not been implemanted yet ");
        }
};


#endif
/*-----------------------------------------------------------------*/
/*---- End Of File ----*/
