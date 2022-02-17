#ifndef _CURRENCY_H
#define _CURRENCY_H



#include <stdio.h>
#include <string.h>

//#include "zerocurv.h"
#include "armcalypso.h"
#include "armglob.h"

namespace ARM {
#define CREDIT_LAG 2;

class ARM_ZeroCurve;
class ARM_Currency : public ARM_Object
{
    private:

    // the following variables are initialized when initializing the object
    // and are not to be modified afterwards
    
        char   itsName[4];		// Code ISO 

        double itsCrossValue;   // against security currency price ???

        int itsSpotDays;        // Nb of days after (>0) or prior (<0) 
                                // payment will be made

		int itsCreditStartDateLag;        // Nb of days after (>0) or prior (<0) 
										 // start date of the Credit product

        int itsFxSpotDays;        // Nb of days after (>0) or prior (<0) 
                                // payment will be made for Fx

        int itsLiborIndexDayCount;        // day count basis for libor

        int itsMMDayCount;        // day count basis for MM

        int itsLiborTerm;        // Term of Libor rate

        int itsFwdRule;         // adjust payment date to the previous (-1) 
                                // or next (+1) business day

        int itsFixedPayFreq;     // Payment frequency of the fixed leg 
                                 // of plain vanilla swap

        int itsCashPayFreq;     // Payment frequency for bond

        int itsFixedDayCount;   // day count basis for fixed leg 

        int itsCashDayCount;    // day count basis for bond


        ARM_ZeroCurve* itsZeroCurve; // Risk Free Zero Curve

    public:

        ARM_Currency(void);

        ARM_Currency(char* name);

        ARM_Currency(const char* name);

        ARM_Currency(char* name, ARM_ZeroCurve* zeroCurve);

        ARM_Currency(char* name, int spotDays, int dayCount, int fwdRule, 
                     ARM_ZeroCurve* zeroCurve);

        ARM_Currency& operator = (const ARM_Currency& currency);

        ARM_Currency(const ARM_Currency& currency);

       ~ARM_Currency(void)
        {
            itsZeroCurve = (ARM_ZeroCurve *) NULL;
        }

        int operator == (ARM_Currency& ccy) 
        { 
           if ( strcmp(itsName, ccy.itsName) == 0 )
           {
              return(1);
           }
           else
           {
              return(0);
           }
        }

        void View(char* id, FILE* ficOut);

        void Init(void)
        {
            strcpy(itsName , ARM_DEFAULT_COUNTRY);

            itsSpotDays = 2; 

			itsCreditStartDateLag = CREDIT_LAG; // figé à 1 jour sur toutes les devises CF juristes

            itsFxSpotDays = 2; // For the currencies used so far 

            itsCrossValue = 1.0;

            itsLiborIndexDayCount = KACTUAL_360;

            itsMMDayCount = KACTUAL_360;

            itsLiborTerm = K_SEMIANNUAL;

            itsFwdRule = K_MOD_FOLLOWING; 
	
            itsFixedDayCount = K30_360;

            itsCashDayCount = K30_360;

            itsFixedPayFreq = K_SEMIANNUAL;

            itsCashPayFreq = K_SEMIANNUAL;

            itsZeroCurve = NULL;
        }

        void BitwiseCopy(const ARM_Object* srcCcy)
        {
            ARM_Currency* src = (ARM_Currency *) srcCcy;


            strcpy(itsName, src->itsName);

            itsCrossValue = src->itsCrossValue;

            itsSpotDays = src->itsSpotDays;

			itsCreditStartDateLag = src->itsCreditStartDateLag;

            itsFxSpotDays = src->itsFxSpotDays;

            itsLiborIndexDayCount = src->itsLiborIndexDayCount;

            itsMMDayCount = src->itsMMDayCount;

            itsLiborTerm = src->itsLiborTerm;

            itsFwdRule = src->itsFwdRule;

            itsFixedPayFreq = src->itsFixedPayFreq;

            itsCashPayFreq = src->itsCashPayFreq;

            itsFixedDayCount = src->itsFixedDayCount;

            itsCashDayCount = src->itsCashDayCount;

            itsZeroCurve = src->itsZeroCurve;
        }
 
        void Copy(const ARM_Object* srcCcy)
        {
            ARM_Object::Copy(srcCcy);

            this->BitwiseCopy(srcCcy);
        }

        ARM_Object* Clone(void)
        {
            ARM_Currency* theClone = new ARM_Currency();
 

            theClone->Copy(this);
  
            return(theClone);
        }

        ARM_CLASS_NAME GetRootName(void)
        {
            return(ARM_CURRENCY);
        }

        void Set(char *name, int spotDays, int dayCount, int fwdRule, 
                 ARM_ZeroCurve* zeroCurve);

        ARM_INDEX_TYPE GetVanillaIndexType(void);

        double GetCrossValue(void) const 
        { 
            return(itsCrossValue); 
        }

        void SetCrossValue(double xvalue)
        {
            itsCrossValue = xvalue;
        }

        int GetSpotDays(void) const 
        { 
            return(itsSpotDays); 
        }

        void SetSpotDays(int spotDays)
        { 
            itsSpotDays = spotDays; 
        }

		int GetCreditStartDateLag(void) const 
        { 
            return(itsCreditStartDateLag); 
        }

        void SetCreditStartDateLag(int creditLag)
        { 
            itsCreditStartDateLag = creditLag; 
        }

        int GetFxSpotDays(void) const 
        { 
            return(itsFxSpotDays); 
        }

        void SetFxSpotDays(int spotDays)
        { 
            itsFxSpotDays = spotDays; 
        }

        int GetLiborIndexDayCount(void) const 
        { 
            return(itsLiborIndexDayCount); 
        }

        void SetLiborIndexDayCount(int dayCount)
        { 
            itsLiborIndexDayCount = dayCount; 
        }

        int GetMMDayCount(void) const 
        { 
            return(itsMMDayCount); 
        }

        void SetMMDayCount(int dayCount)
        { 
            itsMMDayCount = dayCount; 
        }

		int GetFixedDayCountFromIndex(ARM_INDEX_TYPE pIndexType) const
        {
			int fixedDayCount = GetFixedDayCount(); // get infos from default index
            
			switch (pIndexType) // specify different value for other non vanilla index
            {
                case EONIA:
                case EURIBOR1M:
                case EUR1M:
                {
					fixedDayCount = KACTUAL_360;
                    return(fixedDayCount);
                };
                break;

                default:
                {
                    return(fixedDayCount);
                };
                break;
            }
        }

        int GetFixedDayCount(void) const 
        { 
            return(itsFixedDayCount); 
        }

        void SetFixedDayCount(int dayCount)
        { 
            itsFixedDayCount = dayCount; 
        }

        int GetFwdRule(void) const 
        { 
            return(itsFwdRule); 
        }

        void SetFwdRule(int rule)
        { 
            itsFwdRule = rule; 
        }

        int GetFixedPayFreq(void) const 
        { 
            return(itsFixedPayFreq); 
        }

        void SetFixedPayFreq(int freq)
        { 
            itsFixedPayFreq = freq; 
        }

        void SetZeroCurve(ARM_ZeroCurve* zeroCurve)
        { 
            itsZeroCurve = zeroCurve; 
        }

        ARM_ZeroCurve* GetZeroCurve(void) const 
        { 
            return(itsZeroCurve);
        }

        void SetCcyName(char* cname)
        { 
        	strcpy(itsName , cname);
        }

        char* GetCcyName(void) const 
        { 
            return((char *) itsName);
        }

        void CalcResetCal(char* resetCal)
        {
            char* resCal = GetResetCalName(GetVanillaIndexType());

            strcpy(resetCal, resCal);

            delete [] resCal;
        }

        char* GetResetCalName(ARM_INDEX_TYPE pLiborType) const
        {
            switch (pLiborType)
            {
                case LIBOR1M:
                case LIBOR2M:
                case LIBOR3M:
                case LIBOR6M:
                case LIBOR1Y:
                {
                    char* newCal = new char [4];

                    if ( strcmp(itsName, "AUD") == 0 )
                    {
                       strcpy(newCal, "AUD"); 
                    }
                    else if ( strcmp(itsName, "SKK") == 0 )
                    {
                        strcpy(newCal, "SKK");
                    }
                    else if ( strcmp(itsName, "CNY") == 0 )
                    {
                       strcpy(newCal, "CNY");
                    }
                    else if ( strcmp(itsName, "KRW") == 0 )
                    {
                       strcpy(newCal, "KRW");
                    }                    else if ( strcmp(itsName, "INR") == 0 )
                    {
                       strcpy(newCal, "INR");
                    }
                    else if ( strcmp(itsName, "TWD") == 0 )
                    {
                       strcpy(newCal, "TWD");
                    }
                    else if ( strcmp(itsName, "HKD") == 0 )
                    {
                       strcpy(newCal, "HKD");
                    }
                    else if ( strcmp(itsName, "SGD") == 0 )
                    {
                       strcpy(newCal, "SGD");
                    }
                    else
                    {
                       strcpy(newCal, "GBP"); 
                    }

                    return(newCal);
                };
                break;

                default:
                {
                    char* newCal = new char [strlen(itsName)+1];

                    strcpy(newCal, itsName);

                    return(newCal);
                };
                break;
            }
        }

        void CalcFloatPayCal(char* payCal)
        {
            char* resPayCal = GetPayCalName(GetVanillaIndexType());

            strcpy(payCal, resPayCal);

            delete [] resPayCal;
        }

        void CalcFixPayCal(char* payCal)
        {
            if (( strcmp(itsName, "AUD") == 0 )
                ||
                ( strcmp(itsName, "GBP") == 0 )
               )
            {
               char buf[100];

               strcpy(buf, itsName);

               strcpy(payCal, buf); 
            }
            else
            {
               CalcFloatPayCal(payCal);
            }
        }

        char* GetPayCalName(ARM_INDEX_TYPE pLiborType) const
        {
            switch (pLiborType)
            {
                case LIBOR1M:
                case LIBOR2M:
                case LIBOR3M:
                case LIBOR6M:
                case LIBOR1Y:
                {
                    if ( strcmp(itsName, "CNY") == 0 ) // In fact: CNY/USD
                    {
                       char* newCal = new char [strlen(itsName)+1];

                       strcpy(newCal, "CZU");

                       return(newCal);
                    }
                    else if ( strcmp(itsName, "KRW") == 0 )
                    {
                       char* newCal = new char [strlen(itsName)+1];

                       strcpy(newCal, "KZU");

                       return(newCal);
                    }

                    if ((strcmp(itsName, "JPY") != 0)
                        && (strcmp(itsName, "USD") != 0))
                    {
                       char* newCal = new char [strlen(itsName)+1];

                       strcpy(newCal, itsName);

                       return(newCal);
                    }
                    else
                    {
						char* newCal = new char [10];
	
						if (GetCALYPSOVersion() == 0)
						{
							char itsPayCalName[4];
							itsPayCalName[0] = 'Z';
							itsPayCalName[1] = 'G';
							itsPayCalName[2] = itsName[0];
							itsPayCalName[3] = '\0';
             
							strcpy(newCal, itsPayCalName);
						}
						else
						{
							strcpy(newCal, "GBP");
							strcat(newCal, itsName);
						}

						return(newCal);
                    };
                    break;
                }

                default:
                {
                    char* newCal = new char [strlen(itsName)+1];

                    strcpy(newCal, itsName);

                    return(newCal);
                };
                break;
            }
        }

        

		void SetCurrencyProperties(void);

        inline int GetLiborTerm(void)
        {
            return itsLiborTerm;
        }


};
}
#endif
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
