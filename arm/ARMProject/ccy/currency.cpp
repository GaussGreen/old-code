#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "currency.h"
#include "paramview.h"

namespace ARM{

ARM_Currency* ARM_DEFAULT_CURRENCY = new ARM_Currency(ARM_DEFAULT_COUNTRY);

ARM_Currency* ARM_FRF_CURRENCY    = new ARM_Currency("FRF");


void CUSTOMISE_CURRENCY(void)
{
    CUSTOMISE_COUNTRY();

    if (ARM_DEFAULT_CURRENCY)
    {
       ARM_Currency ccy(ARM_DEFAULT_COUNTRY);

       *ARM_DEFAULT_CURRENCY = ccy; 
    }
}


void CUSTOMISE_ENV(void)
{
    CUSTOMISE_CURRENCY();
}


ARM_Currency::ARM_Currency(void)
{
    Init();

    SetName(ARM_CURRENCY);

    SetCurrencyProperties();
}



ARM_Currency::ARM_Currency(char* cname)
{ 
    Init();

    SetName(ARM_CURRENCY);

    strcpy(itsName, cname);

    SetCurrencyProperties();
}



ARM_Currency::ARM_Currency(const char* cname)
{ 
    Init();

    SetName(ARM_CURRENCY);

    strcpy(itsName, cname);

    SetCurrencyProperties();
}



ARM_Currency::ARM_Currency(char* cname, ARM_ZeroCurve* zeroCurve)
{
    Init();

    SetName(ARM_CURRENCY);

    strcpy(itsName , cname);

    SetCurrencyProperties();

    itsZeroCurve = zeroCurve;
}



void ARM_Currency::Set(char *cname, int spotDays, 
                       int dayCount, int fwdRule, 
                       ARM_ZeroCurve* zeroCurve)
{
    strcpy(itsName , cname);

    SetCurrencyProperties();

    itsSpotDays = spotDays; 

    itsFxSpotDays = spotDays; 

    itsLiborIndexDayCount = dayCount; 

    itsFwdRule = fwdRule;

    itsZeroCurve = zeroCurve;
}


/* fonction a mettre apres le Init et le set correct de "itsName" */
void ARM_Currency::SetCurrencyProperties(void)
{
   ARM_COUNTRY_SYMBOL ctry = ARM_GetCountrySymbol(itsName);

   switch(ctry)
   {
      case ARM_ZAR :
      {
         itsSpotDays = 0;
         
         itsFixedDayCount = KACTUAL_365; 
         itsLiborIndexDayCount = KACTUAL_365; 
         itsFixedPayFreq = K_QUARTERLY;

         itsLiborTerm = K_QUARTERLY;
      };
      break;

      case ARM_HUF :
      {
         itsFixedDayCount = KACTUAL_360; 
         itsLiborIndexDayCount = KACTUAL_360; 

         itsFixedPayFreq = K_QUARTERLY;
         
         itsLiborTerm = K_QUARTERLY;
      };
      break;

      case ARM_PLN :
      {
         itsFixedDayCount = KACTUAL_365; 
         itsLiborIndexDayCount = KACTUAL_365; 

         itsFixedPayFreq = K_ANNUAL;
      };
      break;

      case ARM_BEF:
      {
         itsCashPayFreq = K_ANNUAL;
         itsFixedDayCount = KACTUAL_365; 
         itsFixedPayFreq = K_ANNUAL;
         itsLiborIndexDayCount = KACTUAL_365; 
      };
      break;

      case ARM_DKK:
      {
         itsSpotDays = 0; // MA Le 9 Oct 2001
                          // Avant itsSpotDays = 1
         itsLiborIndexDayCount = KACTUAL_360; 

         itsFixedPayFreq = K_ANNUAL;
      };
      break;

      case ARM_GRD:
      {
         itsFixedDayCount      = KACTUAL_365;
         itsLiborIndexDayCount = KACTUAL_365;
      };
      break;

      case ARM_ATS:
      case ARM_CHF:
      case ARM_DEM:
      case ARM_NLG:
      case ARM_SEK:
      case ARM_XEU:
      case ARM_NOK:
      {
         itsFixedPayFreq = K_ANNUAL;
         itsCashPayFreq = K_ANNUAL;
      };
      break;

      case ARM_ESP:
      {
         itsFixedPayFreq = K_ANNUAL;
      };
      break;

      case ARM_FRF:
      {
         itsFixedPayFreq = K_ANNUAL;
         itsCashPayFreq = K_ANNUAL;
         itsSpotDays = 1;
         itsCashDayCount = KACTUAL_FEB29; 
         itsLiborTerm = K_QUARTERLY;
      };
      break;

      case ARM_GBP:
      {
         itsSpotDays = 0;
         itsMMDayCount = KACTUAL_365;
         itsCashDayCount = KACTUAL_365; 
         itsFixedDayCount = KACTUAL_365; 
         itsLiborIndexDayCount = KACTUAL_365;
      };
      break;

      case ARM_ITL:
      {
         itsFixedPayFreq = K_ANNUAL;
         itsCashDayCount = K30_360; 
      };
      break;

      case ARM_JPY:
      {
//         itsFixedPayFreq = K_ANNUAL;
         itsFixedDayCount = KACTUAL_365; 
      };
      break;

      case ARM_PTE:
      {
         itsFixedPayFreq = K_ANNUAL;
         itsLiborIndexDayCount = KACTUAL_365; 
      };
      break;

      case ARM_USD:
      {
         itsFxSpotDays = 2; // JMP 7 IV 2005 avant = 0;
         itsLiborTerm = K_QUARTERLY;
      };
      break;

      case ARM_CAD:
      {
         itsSpotDays = 1;

		 itsMMDayCount = KACTUAL_365;

         itsLiborIndexDayCount = KACTUAL_360;

         itsFixedPayFreq = K_SEMIANNUAL;

         itsFixedDayCount = KACTUAL_365;
      };
      break;

      case ARM_CZK:
      {
         itsLiborIndexDayCount = KACTUAL_360;

         itsFixedPayFreq = K_ANNUAL;

         itsFixedDayCount = KACTUAL_360;
      };
      break;

      case ARM_AUD:
      {
         itsFxSpotDays = 2; // 1 before!

         itsMMDayCount = KACTUAL_365;

         itsLiborTerm = K_QUARTERLY;

         itsLiborIndexDayCount = KACTUAL_365; 
         itsFixedDayCount      = KACTUAL_365;


         itsFixedPayFreq = K_SEMIANNUAL;

         itsCashPayFreq  = K_QUARTERLY;

         itsCashPayFreq  = K_QUARTERLY;

         itsCashDayCount = KACTUAL_365;

         itsSpotDays = 1;
      };
      break;

      case ARM_EUR:
      {
         itsFixedPayFreq = K_ANNUAL;
         itsCashPayFreq = K_ANNUAL;
         itsSpotDays = 2;
         itsCashDayCount = KACTUAL_FEB29; 
      }
      break;

      case ARM_SKK:
      {
          itsLiborIndexDayCount = KACTUAL_360;

          itsFixedPayFreq       = K_ANNUAL;

          itsFixedDayCount      = KACTUAL_360;

      };
      break;

      case ARM_HKD:
      {
         itsSpotDays = 1;
         itsLiborIndexDayCount = KACTUAL_365;
         itsMMDayCount    = KACTUAL_365;
         itsLiborTerm     = K_QUARTERLY;
         itsFixedDayCount = KACTUAL_365;
         itsCashDayCount  = KACTUAL_365;
         itsFixedPayFreq  = K_QUARTERLY;
         itsCashPayFreq   = K_QUARTERLY;
      };
      break;

      case ARM_SGD:
      {
         itsLiborIndexDayCount = KACTUAL_365;
         itsMMDayCount    = KACTUAL_365;
         itsLiborTerm     = K_SEMIANNUAL;
         itsFixedDayCount = KACTUAL_365;
         itsCashDayCount  = KACTUAL_365;
         itsFixedPayFreq  = K_SEMIANNUAL;
         itsCashPayFreq   = K_SEMIANNUAL;
      };
      break;

      case ARM_PHP:
      {
         itsSpotDays = 1;
         itsLiborIndexDayCount = KACTUAL_365;
         itsMMDayCount    = KACTUAL_365;
         itsLiborTerm     = K_SEMIANNUAL;
         itsFixedDayCount = KACTUAL_365;
         itsCashDayCount  = KACTUAL_365;
         itsFixedPayFreq  = K_SEMIANNUAL;
      };
      break;

      case ARM_NZD:
      {
         itsLiborIndexDayCount = KACTUAL_365;
         itsMMDayCount = KACTUAL_365;
         itsLiborTerm     = K_SEMIANNUAL;
         itsFixedDayCount = KACTUAL_365;
         itsCashDayCount  = KACTUAL_365;
         itsFixedPayFreq  = K_SEMIANNUAL;
      };
      break;

      case ARM_CNY:
      {
         itsSpotDays = 1;

         itsLiborIndexDayCount = KACTUAL_365;
         itsMMDayCount         = KACTUAL_365;
         itsLiborTerm          = K_QUARTERLY;
         itsFixedDayCount      = KACTUAL_365;
         itsCashDayCount       = KACTUAL_365;
         itsFixedPayFreq       = K_QUARTERLY;
         itsCashPayFreq        = K_QUARTERLY;
      };
      break;

      case ARM_TWD:
      {
         itsLiborIndexDayCount = KACTUAL_365;
         itsMMDayCount         = KACTUAL_365;
         itsLiborTerm          = K_QUARTERLY;
         itsFixedDayCount      = KACTUAL_365;
         itsCashDayCount       = KACTUAL_365;
         itsFixedPayFreq       = K_QUARTERLY;
         itsCashPayFreq        = K_QUARTERLY;
      };
      break;

      case ARM_KRW:
      {
         itsSpotDays = 1;
         itsLiborIndexDayCount = KACTUAL_365;
         itsMMDayCount         = KACTUAL_365;
         itsLiborTerm          = K_QUARTERLY;
         itsFixedDayCount      = KACTUAL_365;
         itsCashDayCount       = KACTUAL_365;
         itsFixedPayFreq       = K_QUARTERLY;
         itsCashPayFreq        = K_QUARTERLY;
      };
      break;

      case ARM_INR:
      {
         itsSpotDays = 1;
         itsLiborIndexDayCount = KACTUAL_365;
         itsMMDayCount         = KACTUAL_365;
         itsLiborTerm          = K_ANNUAL;
         itsFixedDayCount      = KACTUAL_365;
         itsCashDayCount       = KACTUAL_365;
         itsFixedPayFreq       = K_ANNUAL;
         itsCashPayFreq        = K_ANNUAL;
      };
      break;
   }
}


ARM_Currency::ARM_Currency(char* cname, int spotDays, 
                           int dayCount, int fwdRule, 
                           ARM_ZeroCurve* zeroCurve)
{
    Init();

    SetName(ARM_CURRENCY);

    Set(cname, spotDays, dayCount, fwdRule, zeroCurve);
}



ARM_Currency::ARM_Currency(const ARM_Currency& ccy) : ARM_Object(ccy)
{
    Init();

    this->BitwiseCopy(&ccy);    
}



ARM_Currency & ARM_Currency::operator = (const ARM_Currency& ccy)
{
    (*this).ARM_Object::operator = (ccy);

    this->BitwiseCopy(&ccy); 

    return(*this);
}


void ARM_Currency::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[254];

    if ( ficOut == NULL )
    {
        ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);        
        fOut = fopen(fOutName, "w");
    }
    else
    {
        fOut = ficOut;
    }

    try
    {
        fprintf(fOut, "\t ======> Currency Infos <======\n\n");

		fprintf(fOut, "\t Currency Name      : %s \n", itsName);
		fprintf(fOut, "\t Forward Rule       : %s \n", ARM_ParamView::GetMappingName( S_FORWARD_RULES, itsFwdRule ));
		fprintf(fOut, "\t Libor Term         : %s \n", ARM_ParamView::GetMappingName( S_FREQUENCY,  itsLiborTerm ));
		fprintf(fOut, "\t Libor DayCount     : %s \n", ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsLiborIndexDayCount ));
		fprintf(fOut, "\t Spot Days          : %d \n", itsSpotDays ) ;
		fprintf(fOut, "\t MM DayCount        : %s \n", ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsMMDayCount ));
		fprintf(fOut, "\t Fixed PayFrequency : %s \n", ARM_ParamView::GetMappingName( S_FREQUENCY, itsFixedPayFreq ));
		fprintf(fOut, "\t Fixed DayCount     : %s \n", ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsFixedDayCount ));
		fprintf(fOut, "\t Cash PayFrequency  : %s \n", ARM_ParamView::GetMappingName( S_FREQUENCY, itsCashPayFreq ));
		fprintf(fOut, "\t Cash DayCount      : %s \n", ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsCashDayCount ));
        fprintf(fOut, "\n\t ======> End of Currency Infos <======\n");
    
	}catch(.../*Exception&*/ )
    {
        //throw Exception(__LINE__, __FILE__, ERR_CALC_DET_PB, /* the constant ERR_CALC_DET_PB is not used*/
        //                    "At least one of the flags can not be found");
    }

/*
    if (itsZeroCurve)
    {
        itsZeroCurve->View(id, fOut);
    }
*/

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}



ARM_INDEX_TYPE ARM_Currency::GetVanillaIndexType(void)
{
    ARM_COUNTRY_SYMBOL ctry = ARM_GetCountrySymbol(itsName);

    ARM_INDEX_TYPE defaultIndex = (ARM_INDEX_TYPE) LIBOR6M;


    switch(ctry)
    {
        case ARM_ZAR :
        {
            defaultIndex = (ARM_INDEX_TYPE) LIBOR3M;
        };
        break;

        case ARM_HUF :
        {
            defaultIndex = (ARM_INDEX_TYPE) LIBOR3M;
        };
        break;

        case ARM_PLN :
        {
            defaultIndex = (ARM_INDEX_TYPE) LIBOR6M;
        };
        break;

        case ARM_FRF :
        {
            defaultIndex = (ARM_INDEX_TYPE) PIBOR3M;
        };
        break;

        case ARM_USD :
        {
           defaultIndex = (ARM_INDEX_TYPE) LIBOR3M;
        };
        break;

        case ARM_AUD :
        {
           defaultIndex = (ARM_INDEX_TYPE) LIBOR3M;
        };
        break;

        case ARM_SEK :
        {
           defaultIndex = (ARM_INDEX_TYPE) LIBOR3M;
        };
        break;

        case ARM_EUR :
        {
            defaultIndex = (ARM_INDEX_TYPE) EURIBOR6M;
        };
        break;

        case ARM_GBP :
        {
           defaultIndex = (ARM_INDEX_TYPE) LIBOR6M;
        };
        break;

        case ARM_JPY :
        {
           defaultIndex = (ARM_INDEX_TYPE) LIBOR6M;
        };
        break;

        case ARM_DKK :
        {
           defaultIndex = (ARM_INDEX_TYPE) LIBOR6M;
        };
        break;

        case ARM_GRD :
        {
           defaultIndex = (ARM_INDEX_TYPE) LIBOR6M;
        };
        break;

        case ARM_CHF :
        {
           defaultIndex = (ARM_INDEX_TYPE) LIBOR6M;
        };
        break;

        default :
        {
            return(defaultIndex);
        };
        break;
    }

    return(defaultIndex);
}
}


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
