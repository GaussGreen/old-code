/*
 * $Log: armglob.cpp,v $
 * Revision 1.22  2004/06/04 16:21:59  jpriaudel
 * ARM_GetIsoFromCountrySymbol added
 *
 * Revision 1.21  2003/12/10 17:08:47  jpriaudel
 * ajout de ARM_NZD
 *
 * Revision 1.20  2003/11/06 15:02:16  jpriaudel
 * Added : currency : PHP
 *
 * Revision 1.19  2003/09/29 07:40:12  ebenhamou
 * remove second definition
 *
 * Revision 1.17  2003/09/26 13:59:19  mab
 * Changes in "namespace" ARM_Constants
 *
 * Revision 1.16  2003/07/21 12:59:50  mab
 * correction of a typing pb! in:
 * ARM_GetViewFile
 *
 * Revision 1.15  2003/07/01 17:29:59  jpriaudel
 * bug correction
 *
 * Revision 1.14  2003/06/23 07:17:55  ebenhamou
 * remove unused var
 *
 * Revision 1.13  2003/05/21 12:13:35  mab
 * Added : currencies : HKD , SGD
 *
 * Revision 1.12  2003/01/14 13:00:54  mab
 * Ameliorations
 *
 * Revision 1.11  2001/10/09 09:02:43  mab
 * Rajout de : ARM_ZAR, ARM_HUF, ARM_PLN
 *
 * Revision 1.10  2001/06/28 14:16:29  mab
 * Rajout de ifdef unix pour traiter le cas d'un usage
 * sur poste NT
 *
 * Revision 1.9  2001/03/19 08:44:24  mab
 * Rajout de : ARM_AUD ds isoSymbols
 *
 * Revision 1.8  1999/08/13 12:09:18  mab
 * void InitCurveName(void)
 *
 * Revision 1.7  1999/02/04 16:00:36  ypilchen
 * Rajout ARM_IEP, ARM_FIM, ARM_NOK, ARM_GBP
 *
 * Revision 1.6  1998/12/10 15:01:18  nicolasm
 * EURO : Ajout de la mise a jour de la devise a partir de la nouvelle
 * ARM_DEFAULT_COUNTRY dans la fonction ARM_SetDefaultCountry
 *
 * Revision 1.5  1998/12/08 09:24:54  ypilchen
 * Rajout de $Log: armglob.cpp,v $
 * Rajout de Revision 1.22  2004/06/04 16:21:59  jpriaudel
 * Rajout de ARM_GetIsoFromCountrySymbol added
 * Rajout de
 * Rajout de Revision 1.21  2003/12/10 17:08:47  jpriaudel
 * Rajout de ajout de ARM_NZD
 * Rajout de
 * Rajout de Revision 1.20  2003/11/06 15:02:16  jpriaudel
 * Rajout de Added : currency : PHP
 * Rajout de
 * Rajout de Revision 1.19  2003/09/29 07:40:12  ebenhamou
 * Rajout de remove second definition
 * Rajout de
 * Rajout de Revision 1.17  2003/09/26 13:59:19  mab
 * Rajout de Changes in "namespace" ARM_Constants
 * Rajout de
 * Rajout de Revision 1.16  2003/07/21 12:59:50  mab
 * Rajout de correction of a typing pb! in:
 * Rajout de ARM_GetViewFile
 * Rajout de
 * Rajout de Revision 1.15  2003/07/01 17:29:59  jpriaudel
 * Rajout de bug correction
 * Rajout de
 * Rajout de Revision 1.14  2003/06/23 07:17:55  ebenhamou
 * Rajout de remove unused var
 * Rajout de
 * Rajout de Revision 1.13  2003/05/21 12:13:35  mab
 * Rajout de Added : currencies : HKD , SGD
 * Rajout de
 * Rajout de Revision 1.12  2003/01/14 13:00:54  mab
 * Rajout de Ameliorations
 * Rajout de
 * Rajout de Revision 1.11  2001/10/09 09:02:43  mab
 * Rajout de Rajout de : ARM_ZAR, ARM_HUF, ARM_PLN
 * Rajout de
 * Rajout de Revision 1.10  2001/06/28 14:16:29  mab
 * Rajout de Rajout de ifdef unix pour traiter le cas d'un usage
 * Rajout de sur poste NT
 * Rajout de
 * Rajout de Revision 1.9  2001/03/19 08:44:24  mab
 * Rajout de Rajout de : ARM_AUD ds isoSymbols
 * Rajout de
 * Rajout de Revision 1.8  1999/08/13 12:09:18  mab
 * Rajout de void InitCurveName(void)
 * Rajout de
 * Rajout de Revision 1.7  1999/02/04 16:00:36  ypilchen
 * Rajout de Rajout ARM_IEP, ARM_FIM, ARM_NOK, ARM_GBP
 * Rajout de
 * Rajout de Revision 1.6  1998/12/10 15:01:18  nicolasm
 * Rajout de EURO : Ajout de la mise a jour de la devise a partir de la nouvelle
 * Rajout de ARM_DEFAULT_COUNTRY dans la fonction ARM_SetDefaultCountry
 * Rajout de Pour RCS et amelioration de ARM_GetCountrySymbol
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : armglob.cpp                                                  */
/*                                                                            */
/* DESCRIPTION : global utilities                                             */
/*                                                                            */
/* DATE        : Wed Apr 17 1996                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/


 
#include <stdarg.h>

#include <stdio.h>

#include <ctype.h>

#include <string.h>

#include <stdlib.h>

#include "expt.h"

#include "armglob.h"
#include "currency.h"

#ifdef unix
#    include <unistd.h>
#    include <fcntl.h>
#    include <sys/stat.h>
#endif


/*---- Global variables ----*/

int gTrace = 0;


#ifdef __ARM_NO_NAMESPACE

double rateBase   = 100.0;
double volBase    = 100.0;
double correlBase = 100.0;

#else
/// no definition as it is already in the header!

#endif


static ARM_COUNTRY_SYMBOL_ISO
  isoSymbols[100] = {
                       { ARM_ZAR, "ZAR"},
                       { ARM_HUF, "HUF"},
                       { ARM_PLN, "PLN"},
                       { ARM_ATS, "ATS"},
                       { ARM_BEF, "BEF"},
                       { ARM_CAD, "CAD"},
                       { ARM_CHF, "CHF"},
                       { ARM_DEM, "DEM"},
                       { ARM_DKK, "DKK"},
                       { ARM_ESP, "ESP"},
                       { ARM_EUR, "EUR"},
                       { ARM_FRF, "FRF"},
                       { ARM_GBP, "GBP"},
                       { ARM_GRD, "GRD"},
                       { ARM_ITL, "ITL"},
                       { ARM_JPY, "JPY"},
                       { ARM_NLG, "NLG"},
                       { ARM_NOK, "NOK"},
                       { ARM_PTE, "PTE"},
                       { ARM_SEK, "SEK"},
                       { ARM_USD, "USD"},
                       { ARM_XEU, "XEU"},
                       { ARM_IEP, "IEP"},
                       { ARM_FIM, "FIM"},
                       { ARM_AUD, "AUD"},
                       { ARM_CZK, "CZK"},
                       { ARM_HKD, "HKD"},
                       { ARM_SGD, "SGD"},
                       { ARM_PHP, "PHP"},
                       { ARM_NZD, "NZD"},
                       { ARM_SKK, "SKK"},
                       { ARM_CNY, "CNY"},
                       { ARM_TWD, "TWD"},
                       { ARM_KRW, "KRW"},
                       { ARM_INR, "INR"}
                    };




/*---------------------------------------------------------------------------*/


ARM_Country ARM_DEFAULT_COUNTRY = "EUR";




void InitCurveName(char* cvName)
{
    char* SUMMIT_CURVE_NAME = NULL;
 
 
    SUMMIT_CURVE_NAME = getenv("ARM_CV_NAME");
 
    if ( SUMMIT_CURVE_NAME != NULL )
    {
       strcpy(cvName, SUMMIT_CURVE_NAME);
    }
    else
    {
       strcpy(cvName, CV_HISTORY_NAME);
    }
}



char* ARM_GetDefaultCountry(void)
{
    return((char *) ARM_DEFAULT_COUNTRY);
}



char* ARM_SetDefaultCountry(ARM_Country newCountry)
{
    ARM_COUNTRY_SYMBOL isoSymb;
 
 
    try
    {
       isoSymb = ARM_GetCountrySymbol(newCountry);
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
        throw x;
    }

    strcpy(ARM_DEFAULT_COUNTRY, newCountry);

    ARM::ARM_Currency ccy(ARM_DEFAULT_COUNTRY);

    if (ARM_DEFAULT_CURRENCY)
        *ARM_DEFAULT_CURRENCY = ccy;
    else
        ARM_DEFAULT_CURRENCY = new ARM::ARM_Currency(ARM_DEFAULT_COUNTRY);

    return((char *) ARM_DEFAULT_COUNTRY);
}



ARM_COUNTRY_SYMBOL ARM_GetCountrySymbol(char* isoName)
{
    int sz = 0;
    int i;


    if (( isoName == (char *) NULL ) || ( strlen(isoName) < 3 ))
    {
       Exception x(__LINE__, __FILE__, ERR_OBJECT_UNK,
                "Unknown ISO name");
 
       throw x;
    }

    sz = sizeof(isoSymbols)/sizeof(isoSymbols[0]);

   
    i = 0;

    for (i = 0; i < sz; i++)
    {
        if ( strcmp(isoSymbols[i].isoName, isoName) == 0 )
        {
           return(isoSymbols[i].symb);
        }
    } 

    Exception x(__LINE__, __FILE__, ERR_OBJECT_UNK,
             "Unknown ISO name");
 
    throw x;

    return(ARM_FRF);
}


char* ARM_GetIsoFromCountrySymbol(ARM_COUNTRY_SYMBOL isoSymbol)
{
    return isoSymbols[(int)isoSymbol].isoName;
}



void CUSTOMISE_COUNTRY(void)
{
    char* country = NULL;


    country = getenv("ARM_COUNTRY");

    if ( country != NULL )
    {
       strcpy(ARM_DEFAULT_COUNTRY, country);
    }
}



char* ARM_UPPERCASE(char* str)
{
    int i;

    
    if ( str == NULL )
    {
       return(str);
    }

    i = 0;

    while ( str[i] != '\0' )
    {
        str[i] = toupper(str[i]);

        i++;
    } 

    return(str);
}


/*-----------------------------------------------------------------------------*/
/* Create a Lock File                                                          */
/*-----------------------------------------------------------------------------*/
 
void createFLock(char* flagFileName)
{
#ifdef unix
   umask(0);
   close(creat(flagFileName, 0666));
#endif
}
 
 
 

/*

void ARM_ERROR(char* funcName, int traceModule,...)
{
    va_list args;
    char*   format;


    va_start(args, );

    format = va_arg(args, char *);

    if ((gTrace) || (traceModule))
    {
       fprintf(stderr, "\n <#####?? ERR IN : %s : %s \n", 
                       __FILE__, funcName);
        
       (void) vfprintf(stderr, format, args);
    }

    va_end(args);
}


void ARM_PRINTF(char* funcName, int traceModule,...)
{
    va_list args;
    char*   format;
 
 
    va_start(args, );
 
    format = va_arg(args, char *); 
 
    if ((gTrace) || (traceModule))
    {
       fprintf(stderr, "\n <---- IN : %s : %s \n", 
                           __FILE__, 
                           funcName); 
       
       (void) vfprintf(stderr, format, args);
    }
 
    va_end(args);
}


*/




void ARM_GetViewFile(char* fileName, char* id, char* fOutName)
{
    fOutName[0] = '\0';


    char* homeDir;

    homeDir = getenv("ARM_INSTALL_HOME");


#ifdef unix

    if (homeDir)
    {
       strcpy(fOutName, homeDir);
       strcat(fOutName, "/tmp/");
       strcat(fOutName, fileName);
    }
    else
    {
       strcpy(fOutName, "/tmp/");
       strcat(fOutName, fileName);
    }

#else

    if (homeDir)
    {
       strcpy(fOutName, homeDir);
       strcat(fOutName, "\\temp\\");
       strcat(fOutName, fileName);
    }
    else
    {
       strcpy(fOutName, "c:\\temp\\");
       strcat(fOutName, fileName);
    }

#endif


    if ( id != NULL )
    {
       if ( strlen(id) >= 2 )
       {
          strcat(fOutName, id);
       }
    }
}



void ARM_GetTmpAbsFile(char* fileName, char* fOutName)
{
    fOutName[0] = '\0';


    char* homeDir;

    homeDir = getenv("ARM_INSTALL_HOME");


#ifdef unix

    if (homeDir)
    {
       strcpy(fOutName, homeDir);
       strcat(fOutName, "/tmp/");
       strcat(fOutName, fileName);
    }
    else
    {
       strcpy(fOutName, "/tmp/");
       strcat(fOutName, fileName);
    }

#else

    if (homeDir)
    {
       strcpy(fOutName, homeDir);
       strcat(fOutName, "\\temp\\");
       strcat(fOutName, fileName);
    }
    else
    {
       strcpy(fOutName, "c:\\temp\\");
       strcat(fOutName, fileName);
    }

#endif
}



int ARM_GetCalendarFile(char* fileName)
{
    char* calName;
    FILE* calFile = NULL;
 
 
 

#ifdef unix
    calName = getenv("ARM_CALENDAR");
#else
    calName = "C:\\Program Files\\ARM\\CAL_SUMMIT";
#endif
 
    if (calName)
    {
       strcpy(fileName, calName);

       calFile = fopen(fileName, "r");

       if ( calFile == NULL )
       {
          return(RET_KO);
       }
       else
       {
          fclose(calFile);

          return(RET_OK);
       } 
    }
    else
    {
       return(RET_KO);
    }
}



FILE* ARM_GetOutputHomeFile(char* fileName, char* mode)
{
    char* homeDir;
    char  fOutName[200]; 
    FILE* fOut = NULL; 

 
    homeDir = getenv("HOME");
 
    if (homeDir)
    {
       sprintf(fOutName, "%s/%s", homeDir, fileName);

       fOut = fopen(fOutName, mode);

       if (fOut)
          return(fOut);
    }

    throw Exception(__LINE__, __FILE__, ERR_PB_OPEN_FILE,
                    "Openning a file failed");

    return(NULL);
}



double ARM_Object::CalcNumericalObjectSignature(void)
{
    double signature = -1.0; // No signature

    return(signature);
}



/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
