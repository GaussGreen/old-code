#include "firsttoinc.h"


/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE         : volcurv.cpp                                                 */
/*                                                                            */
/* DESCRIPTION  : This file implements the ARM_VolCurve class, a class for    */
/*                managing volatility curves                                  */
/*                                                                            */
/* DATE         : Tue Jan 14 1997                                             */
/*                                                                            */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*---- System Include ----*/

#ifdef unix
#include <sys/types.h>
#include <unistd.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*---- Application Include ----*/

#include "dates.h"
#include "volcurv.h"
#include "expt.h"
#include "interpol.h"
#include "volfxspinterp.h"
#include "fromto.h"

#include "volint.h"
//#include "irindex.h"

/*----------------------------------------------------------------------------*/
/* Constructor                                                                */
/*----------------------------------------------------------------------------*/

ARM_VolCurve::ARM_VolCurve(const ARM_Date& asOf, ARM_Currency* ccy)
{
    Init();

    SetCurrencyUnit(ccy);

    // Set default values of other variables

    itsAsOfDate = asOf;
	itsLastKnownDate = asOf;

    itsStrikeType = K_STK_TYPE_PRICE;
}



ARM_VolCurve::ARM_VolCurve(const ARM_Date& asOf, int KType, int volType,
                           ARM_Currency* ccy)
{
    Init();

    SetCurrencyUnit(ccy);

    // set default values of other variables
 
    itsAsOfDate = asOf;
	itsLastKnownDate = asOf;
    
 
    itsStrikeType = KType;

    itsVolType = volType;
}



ARM_VolCurve::~ARM_VolCurve(void)
{
    if (itsVolatilities)
       delete itsVolatilities;

    if (itsExpiryDates)
       delete itsExpiryDates;
         
    if (itsExpiryTerms)
       delete itsExpiryTerms;

    if (itsCurrency && itsCurrency != ARM_DEFAULT_CURRENCY)
       delete itsCurrency;

    if (itsFxVolSmileInterpolation)
       delete itsFxVolSmileInterpolation;

    // if (itsIndexName)
    //    delete[] itsIndexName;

	if (itsIndex) 
       delete itsIndex ; 
}

void ARM_VolCurve::SetIndex(const ARM_IRIndex&ref)
{
	/*if (itsIndex) delete itsIndex; 
	itsIndex=NULL; 
	itsIndex=(ARM_IRIndex*) const_cast<ARM_IRIndex&>(ref).Clone(); */
}
//	
const ARM_IRIndex&  ARM_VolCurve::GetIndex() const
{
	if (!itsIndex) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "ARM_VolCurve::GetIndex(): No index defined");
	return *itsIndex; 
}



void ARM_VolCurve::BitwiseCopy(const ARM_Object* srcVolCurve)
{
    ARM_VolCurve* vCurve = (ARM_VolCurve *) srcVolCurve;

    itsLastKnownDate = vCurve->itsLastKnownDate;

    itsStrikeType = vCurve->itsStrikeType;

    itsVolType = vCurve->itsVolType;

	itsOptionType = vCurve->itsOptionType;

    // JLA useless itsIndexType = vCurve->itsIndexType;

    itsInterpType = vCurve->itsInterpType;

    SetCurrencyUnit(vCurve->itsCurrency);

    if (itsExpiryDates)
    {
       delete itsExpiryDates;
       itsExpiryDates = NULL;
    }

    if (vCurve->itsExpiryDates)
       itsExpiryDates = (ARM_Vector *) vCurve->itsExpiryDates->Clone();

    if (itsExpiryTerms)
    {
       delete itsExpiryTerms;
       itsExpiryTerms = NULL;
    }

    if (vCurve->itsExpiryTerms)
       itsExpiryTerms = (ARM_Vector *) vCurve->itsExpiryTerms->Clone();

    if (itsVolatilities)
    {
       delete itsVolatilities;
       itsVolatilities = NULL;
    }

    if (vCurve->itsVolatilities)
       itsVolatilities = (ARM_Matrix *)
                     vCurve->itsVolatilities->Clone();

    if (itsFxVolSmileInterpolation)
    {
       delete itsFxVolSmileInterpolation;
       itsFxVolSmileInterpolation = NULL;
    }

    if (vCurve->itsFxVolSmileInterpolation)
       itsFxVolSmileInterpolation = (ARM_FXVolSmileInterpol *)
                                    vCurve->itsFxVolSmileInterpolation->Clone();

    memcpy(itsYearTermsX, vCurve->itsYearTermsX,
           sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);

    memcpy(itsYearTermsY, vCurve->itsYearTermsY,
           sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);

	// if (itsIndexName) 
    // {
	//  delete[] itsIndexName;
	// 
    //    itsIndexName = NULL;
    // }


	// if (vCurve->itsIndexName)
    // {
    //    itsIndexName = new char[strlen(vCurve->itsIndexName)+1];
	//
    //   strcpy(itsIndexName, vCurve->itsIndexName);
    // }

	itsIndexName=vCurve->itsIndexName; 

	if (vCurve->itsIndex) 
		SetIndex(vCurve->GetIndex()); 
}

 
 

/*----------------------------------------------------------------------------*/
/* Constructor(Copy)                                                          */
/*----------------------------------------------------------------------------*/

ARM_VolCurve::ARM_VolCurve(const ARM_VolCurve& volCurve) 
             : ARM_AbstractMarketClass(volCurve)
{
    Init();

    BitwiseCopy(&volCurve);
}



/*----------------------------------------------------------------------------*/
/* Assignment operator                                                        */
/*----------------------------------------------------------------------------*/

ARM_VolCurve& ARM_VolCurve::operator = (const ARM_VolCurve& volCurve)
{
    (*this).ARM_AbstractMarketClass::operator =(volCurve);

    BitwiseCopy(&volCurve);

    return (*this);
}



void ARM_VolCurve::SetCurrencyUnit(ARM_Currency* ccy) 
{
    if (itsCurrency && itsCurrency != ARM_DEFAULT_CURRENCY)
       delete itsCurrency;

    if (ccy)
       itsCurrency = (ARM_Currency *) ccy->Clone();
    else 
       itsCurrency = NULL;

    SetStrCurrency(itsCurrency->GetCcyName());
}



/******  THE FOLLOWING METHODS MUST BE OVERRIDEN BY SUBCLASSES  *******/



    
/*----------------------------------------------------------------------------*/
/* Returns the volatility for maturity yearTerm and strike.                   */
/*   This method MUST be overriden by subclasses.                             */
/*----------------------------------------------------------------------------*/

double ARM_VolCurve::ComputeVolatility(double yearTerm, double Strike)
{
    if ( yearTerm < 0.0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
              "Invalid maturity in <ComputeVolatility> method");
    }

    double vol = VolatilityFunction(yearTerm, Strike);

    return(vol);
}



double ARM_VolCurve::ComputeVolatility(double mat, double moneyness, double underlying, double strike)
{
    if ( mat < 0.0 )
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Invalid maturity in <ComputeVolatility> method");
    }

    double vol = VolatilityFunction(mat, moneyness, underlying);

    return(vol);
}



double ARM_VolCurve::ComputeVol(double mat, double underlying,
                                double moneyness)
{
    if ( mat < 0.0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Invalid maturity in <ComputeVol> method");
    }

    double vol;

    if ( GetName() == ARM_VOL_CUBE )
    {
       vol = VolatilityFunction(mat, moneyness, underlying);
    }
    else
    {
       vol = VolatilityFunction(mat, underlying);
    }

    return(vol);

}



double ARM_VolCurve::ComputeVolatility(ARM_Date& date, 
                                       double moneyness, double underlying)
{
   if ( date < itsLastKnownDate) 
   {
     throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Invalid maturity in <ComputeVolatility> method");
   }

   double mat = (date-itsLastKnownDate)/365.0;

   double vol = VolatilityFunction(mat, moneyness, underlying);

   return(vol);
}



double ARM_VolCurve::ComputeVol(ARM_Date& date,
								double underlying, double moneyness)
{
   if ( date < itsLastKnownDate) 
   {
     throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Invalid maturity in <ComputeVol> method");
   }

   double mat = (date-itsLastKnownDate)/365.0;

   double vol = ComputeVol(mat, underlying, moneyness);

   return(vol);
}


    
/*----------------------------------------------------------------------------*/
/* Returns the volatility for the date and strike.                            */
/*   This method MUST be overriden by subclasses.                             */
/*----------------------------------------------------------------------------*/

double ARM_VolCurve::ComputeVolatility(ARM_Date& date, double Strike)
{
   if ( date < itsLastKnownDate ) 
   {
     throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
              "Invalid maturity in <ComputeVolatility> method");
   }

   double yearTerm = (date-itsLastKnownDate)/365.0;

   double vol = VolatilityFunction(yearTerm, Strike);

   return(vol);
}


void ARM_VolCurve::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];
    char strDate[50];
 
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
 
       (void) unlink(fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

	fprintf(fOut, "\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
    fprintf(fOut, "\n     >>>>>>>>>>>> Generic Volatility matrix <<<<<<<<<<< \n");

    ARM_AbstractMarketClass::View(id, fOut);

    if ( itsInterpType == K_DIAG_INTERPOL )
    {
       fprintf(fOut, "\n\t Interpolation Mode : DIAG \n\n");
    }
    else
    {
       fprintf(fOut, "\n\t Interpolation Mode : LINEAR \n\n");
    }

     // Type d'option : K_IRG(CAP), K_SWOPT(SWAPTION)

    switch(itsOptionType)
    {
        case K_IRG:
        {
            fprintf(fOut, "\n OPTION type : IRG");
        };
        break;

        case K_SWOPT :
        {
            fprintf(fOut, "\n OPTION type : SWOPT");
        };
        break;

        default:
        {
             fprintf(fOut, "\n OPTION type : UNKNOWN");                       
        };
        break;
    };

    switch(GetVolType()) // K_ATMF_VOL(ATM), K_SMILE_VOL(SMILE), K_FX_VOL_SP_INTERP
    {
        case K_ATMF_VOL:
        {
            fprintf(fOut, "\n Volatility type : ATM");
        };
        break;

        case K_SMILE_VOL:
        {
            fprintf(fOut, "\n Volatility type : SMILE");
        };
        break;

        case K_FX_VOL_SP_INTERP :
        {
            fprintf(fOut, "\n Volatility type : classif FX SMILED Volatility");
        };
        break;

        default:
        {
            char buff[100];

            sprintf(buff, "\n Volatility type UNEXPECTED :----> %d", GetVolType());
        };
        break;
    }

    switch(GetStrikeType())
    {
        case K_STK_TYPE_PRICE:
        {
            fprintf(fOut, " (Price) \n");
        };
        break;

        default:
        {
            fprintf(fOut, " (Yield) \n");
        };
        break;
    }

    GetAsOfDate().JulianToStrDate(strDate);
    fprintf(fOut, "\n AsOfDate  : %s \n\n", strDate);

    if (!itsIndexName.empty()) 
       fprintf(fOut, "\n IndexName  : %s \n\n", itsIndexName.c_str());

	if (itsIndex) 
		fprintf(fOut, "\n Index is set: \n\n" );

	if ( itsAsOfDate != itsLastKnownDate )
    {
		itsLastKnownDate.JulianToStrDate(strDate);
		fprintf(fOut, "\n Last Known Date  : %s \n\n", strDate);
	}

    
	ARM_Vector* terms = GetExpiryTerms();

    if (terms)
    {
        ARM_Matrix* vols = GetVolatilities();
        int nbLg         = terms->GetSize();
        int nbCo         = vols->GetNumCols(); 

        int i, j;


		ARM_Vector* expDates = GetExpiryDates(); 
		if(expDates) 
		{
			fprintf(fOut,"          ");
			for (j = 0; j < nbCo; j++)
				fprintf(fOut,"%10.5lf", expDates->Elt(j));
			fprintf(fOut,"\n");
		}

		for (j = 0; j < nbCo; j++) 
			fprintf(fOut,"[%8s] ",itsYearTermsY[j]); 

		fprintf(fOut,"\n\n"); 

        for (i = 0; i < nbLg; i++)
        {
			fprintf(fOut," [%8s] ", itsYearTermsX[i] );

			fprintf(fOut," %10.5lf   ", terms->Elt(i));
            
            for (j = 0; j < nbCo; j++)
            {
               fprintf(fOut," %10.5lf", vols->Elt(i, j));
            }

            fprintf(fOut,"\n");
        }
    }
    else
    {
       // la vol est flat et on la calcule directement...
       fprintf(fOut," Vol Flat : %10.5lf\n", ComputeVolatility(1.0, 1.0));
    }


    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}



void ARM_VolCurve::ParallelShift(double value)
{
    int nLines = itsVolatilities->GetNumLines();
    int nCols  = itsVolatilities->GetNumCols();

    int i, j;


 
    for (i = 0; i < nLines; i++)
    {
        for (j = 0; j < nCols; j++)
        {
            itsVolatilities->Elt(i, j) += value;
        }
    }
}



void ARM_VolCurve::BumpVolatility(double value, int nthLine, int nthCol,
                                  int isCumul, int isAbsolute)
{
    int nLines = itsVolatilities->GetNumLines();
    int nCols  = itsVolatilities->GetNumCols();



    if (( nthLine > nLines ) || ( nthCol > nCols ))
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                     "Invalid nthLine or nthCol in <BumpVolatility> method");
    }

    int i, j;

    if ( nthLine <= 0 )
    {
       if ( nthCol <= 0 )
       {
          for (i = 0; i < nLines; i++)
          {
              for (j = 0; j < nCols; j++)
              {
				  if( strcmp(itsYearTermsY[j], "0.000000") )
				  {
					  if (isAbsolute == K_YES)
						itsVolatilities->Elt(i, j) += value;
					  else
						itsVolatilities->Elt(i, j) *= (1.+value/100.);
				  }
              }
          }
       }
       else
       {
          if ( isCumul == K_NO )
          {
              for (i = 0; i < nLines; i++)
              {
				  if (isAbsolute == K_YES)
					itsVolatilities->Elt(i, nthCol-1) += value;
				  else
					itsVolatilities->Elt(i, nthCol-1) *= (1.+value/100.);
              }
          }
          else
          {
              for (i = 0; i < nLines; i++)
              {
                  for (j = 0; j < nthCol; j++)
                  {
					  if (isAbsolute == K_YES)
						itsVolatilities->Elt(i, j) += value;
					  else
						itsVolatilities->Elt(i, j) *= (1.+value/100.);
                  }
              }
          }
       }
    }
    else
    {
       if ( nthCol <= 0 )
       {
          if (isCumul == K_NO)
          {
              for (j = 0; j < nCols; j++)
              {
				  if (isAbsolute == K_YES)
					itsVolatilities->Elt(nthLine-1, j) += value;
				  else
					itsVolatilities->Elt(nthLine-1, j) *= (1.+value/100.);
              }
          }
          else
          {
              for (i = 0; i < nthLine; i++)
              {
                  for (j = 0; j < nCols; j++)
                  {
					  if (isAbsolute == K_YES)
						itsVolatilities->Elt(i, j) += value;
					  else
						itsVolatilities->Elt(i, j) *= (1.+value/100.);
                  }
              }
          }
       }
       else
       {
          if ( isCumul == K_NO )
          {
			  if (isAbsolute == K_YES)
				itsVolatilities->Elt(nthLine-1, nthCol-1) += value;
			  else
				itsVolatilities->Elt(nthLine-1, nthCol-1) *= (1.+value/100.);
          }
          else if(isCumul == -1)  // NEW MODE: Cumulative on a matrix 
          {
              for (j = 0; j < nthCol-1; j++)
			  {
				  for(i=0; i<nLines; i++)
				  {
					  if (isAbsolute == K_YES)
						itsVolatilities->Elt(i, j) += value;
					  else
						itsVolatilities->Elt(i, j) *= (1.+value/100.);
				  }
			  }

              for (i = 0; i < nthLine; i++)
              {
				  if (isAbsolute == K_YES)
					itsVolatilities->Elt(i, nthCol-1) += value;
				  else
					itsVolatilities->Elt(i, nthCol-1) *= (1.+value/100.);
			  }
          } 
		  else
          {
              for (i = 0; i < nthLine; i++)
              {
                  for (j = 0; j < nthCol; j++)
                  {
                      if (j < nthCol - 1)
					  {
						  if (isAbsolute == K_YES)
							itsVolatilities->Elt(i, j) += value;
						  else
							itsVolatilities->Elt(i, j) *= (1.+value/100.);
					  }
                      else if (i <= nthLine - 1)
					  {
						  if (isAbsolute == K_YES)
							itsVolatilities->Elt(i, j) += value;
						  else
							itsVolatilities->Elt(i, j) *= (1.+value/100.);
					  }
				  }
              }
          }
       }
    }

    if ( itsVolType == K_FX_VOL_SP_INTERP )
    {
       ARM_FXVolSmileInterpol* fxVolInterpol = GetFxVolSmileInterpolation();

       if (fxVolInterpol)
          delete fxVolInterpol;

       fxVolInterpol = new ARM_FXVolSmileInterpol(this);

       SetFxVolSmileInterpolation(fxVolInterpol);
    }
}



void ARM_VolCurve::UpdateCol(ARM_Vector* pCol, double tenor)
{
    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
           "ARM_VolCurve::UpdateCol should never be called");
}



void ARM_VolCurve::UpdateLine(ARM_Vector* pLine, double yearterm)
{
    if (( pLine == NULL ) || ( GetVolatilities() == NULL ))
       return;

    if ( pLine->GetSize() != GetVolatilities()->GetNumCols() )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                       "Invalid nthCol in <UpdateLine> method");
    }

    long i = indexBeforeValue(GetExpiryTerms(), yearterm);

    long iToUpdate;

    if ( i == GetExpiryTerms()->GetSize()-1 )
       iToUpdate = GetExpiryTerms()->GetSize()-1;
    else if ( i == -1 )
       iToUpdate = 0;
    else if ( fabs(yearterm-GetExpiryTerms()->Elt(i) )
              < fabs(GetExpiryTerms()->Elt(i+1)-yearterm) 
            )
       iToUpdate = i;
    else
       iToUpdate = i+1;

    for (int j = 0; j < pLine->GetSize(); j++)
        GetVolatilities()->Elt(iToUpdate, j) = pLine->Elt(j);
}



double ARM_VolCurve::ComputeFxVol(ARM_Date& AsOf,
                                  ARM_Date& matuDate,
                                  double calcMatu,
                                  double fxSpot, // The FX spot, not the Fwd!
                                  double strike,
                                  ARM_ZeroCurve* discCrv, // JPY
                                  ARM_ZeroCurve* divCrv)  // USD
{
    double vol;

    vol = itsFxVolSmileInterpolation->ComputeFxVol(AsOf, matuDate,
                                                   calcMatu,
                                                   fxSpot,
                                                   strike,
                                                   discCrv,
                                                   divCrv);

    return(vol);
}


bool IS_SCALAR_EQUAL_ONE(double x) // at 1e-6
{
    long truncated_x = long(floor(x*100000.0));

    return( truncated_x == 100000L );
}


bool ARM_VolCurve::IsEqualToOne(void) // at 1e-6
{
    ARM_Matrix* values = GetVolatilities();

    if ( values == NULL ) 
    {
       return(true);
    }

    int nbLn  = values->GetNumLines();
    int nbCol = values->GetNumCols();

    for (int i = 0; i < nbLn; i++)
    {
        for (int j = 0; j < nbCol; j++)
        {
            if (!(IS_SCALAR_EQUAL_ONE(values->Elt(i, j))))
            {
               return(false);
            }
        }
    }

    return(true);
}


                         
                         
                         /*-------------------------------------*/
                         /*                                     */
                         /*            ARM_FxVolCurve           */
                         /*                                     */
                         /*-------------------------------------*/


ARM_FXVolCurve::ARM_FXVolCurve(ARM_Date& Asof,
                               const ARM_Vector& time2maturities, 
			                   const ARM_Vector& pivotVols, 
			                   const ARM_Vector& pivotTypes,
			                   const ARM_Vector& deltasCall,
			                   ARM_Matrix& volsCall,
			                   const ARM_Vector& deltasPut,
			                   ARM_Matrix& volsPut,
			                   const ARM_Vector& inFxFwds,
			                   ARM_Vector& interpolTypesVect,
			                   int whatIsInterpolated,
			                   int correctSplineWithLinear,
                               double spotFX,
                               ARM_ZeroCurve* domCurve,
                               ARM_ZeroCurve* ForCurve,
                               int inRRSTR,
                               int isATM) // ATM flag not used for now
{
    Init();

    try
    {
		itsSpotFX =	spotFX;
		
        if (domCurve)
           itsDomCcy = (ARM_Currency *) domCurve->GetCurrencyUnit()->Clone();
        
        if (ForCurve)
           itsForCcy =   (ARM_Currency *) ForCurve->GetCurrencyUnit()->Clone();

        if  (domCurve)
            itsDomCurve = (ARM_ZeroCurve *) domCurve->Clone();

        if (ForCurve)
           itsForCurve = (ARM_ZeroCurve *) ForCurve->Clone();


        SetExpiryTerms(static_cast<ARM_Vector*>(const_cast<ARM_Vector&>(time2maturities).Clone()));

        GenerateFXVols(Asof,
                       time2maturities, 
			           pivotVols, 
			           pivotTypes,
			           deltasCall,
			           volsCall,
			           deltasPut,
			           volsPut,
			           inFxFwds,
			           interpolTypesVect,
			           whatIsInterpolated,
			           correctSplineWithLinear,
                       spotFX,
                       domCurve,
                       ForCurve,
                       inRRSTR,
                       isATM);

        try
        {
            itsCalcFXATM = ComputeATMFxVol();
        }

        catch(Exception& exp)
        {
            throw exp;
        }

        catch(...)
        {
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                     "Error in : ARM_FXVolCurve::ARM_FXVolCurv and ComputeATMFxVol");
        }
    }

    catch(Exception& exp)
    {
        throw exp;
    }

    catch(...)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                     "Error in : ARM_FXVolCurve::ARM_FXVolCurve and GenerateFXVols");
    }
}



ARM_FXVolCurve::~ARM_FXVolCurve(void)
{   
    if (itsData)
	   delete[] itsData;

    itsData = NULL;

	if (itsDomCcy) 
    { 
        delete itsDomCcy; 
        itsDomCcy = NULL; 
    }
	
    if (itsForCcy) 
    { 
        delete itsForCcy; 
        itsForCcy = NULL; 
    }

    if (itsDomCurve)
    {
       delete itsDomCurve;
       itsDomCurve = NULL;
    }
	
    if (itsForCurve)
    {
       delete itsForCurve;
       itsForCurve = NULL;
    }

    nbTime2Maturities = 0;

    if (itsCalcFXATM)
    {
       delete itsCalcFXATM;
	   itsCalcFXATM = NULL;
	}
}



void ARM_FXVolCurve::GenerateFXVols(ARM_Date& Asof,
                                    const ARM_Vector& time2maturities, 
			                        const ARM_Vector& pivotVols, 
			                        const ARM_Vector& pivotTypes,
			                        const ARM_Vector& deltasCall,
			                        ARM_Matrix& volsCall,
			                        const ARM_Vector& deltasPut,
			                        ARM_Matrix& volsPut,
			                        const ARM_Vector& inFxFwds,
			                        ARM_Vector& interpolTypesVect,
			                        int whatIsInterpolated,
			                        int correctSplineWithLinear,
                                    double spotFX,
                                    ARM_ZeroCurve* domCurve,
                                    ARM_ZeroCurve* ForCurve,
                                    int inRRSTR,
                                    int isATM) // ATM flag not used for now)
{
       if (itsData)
       {
   	      delete[] itsData;
	                            
          itsData = NULL;
		                               
          nbTime2Maturities = 0;
       }

       itsRRSTRInputFlag = inRRSTR;

       SetAsOfDate(Asof);

       // Manage case of FX Vols == 1e-3 (ATM)
	   // WARNING : here we test if the FXVol is smiled
	   // in order to normalize the underlying price in the case of ARM_FXVOLAT(whatever is the model: Multi3F or 2IRFX)
	   int volIsSmiled = 0;
	   // If at least one vol value is > 10E-3, then the FXVol is smiled.
       for (int k = 0; k < volsPut.GetNumLines(); k++)
       {
           for (int h = 0; h < volsPut.GetNumCols(); h++)
           {
               if (( volsPut.Elt(k, h) < 0.0 ) || ( volsPut.Elt(k, h) > 1e-3 ))
               {
                  volIsSmiled = 1;
               }
			   else
			   {
                  volsPut.Elt(k, h) = 0.0;
			   }

               if (( volsCall.Elt(k, h) < 0.0 ) || ( volsCall.Elt(k, h) > 1e-3 ))
               {
                  volIsSmiled = 1;
               }
			   else
			   {
                  volsCall.Elt(k, h) = 0.0;
			   }
           }
       }

	   itsSmileFlag = volIsSmiled;

       if (inRRSTR)
       {
          itsRR  = volsPut;
          itsSTR = volsCall;

          ARM_Vector volsPivots = pivotVols;

          FromRRSTRtoVol(&volsPut, &volsCall, &volsPivots);
       }

       // Assign
       itsOptionsMatus  = time2maturities;

       itsPivotVols     = pivotVols;

       itsPivotTypes    = pivotTypes;
       itsInterpolTypes = interpolTypesVect;
       itsDeltasCall    = deltasCall;
       itsDeltasPut     = deltasPut;

 
       itsVolsPut       = volsPut;
       itsVolsCall      = volsCall;


       int sizeControlOK = (( time2maturities.GetSize() == pivotVols.GetSize() )
                            &&
                            ( time2maturities.GetSize() == pivotTypes.GetSize() )
                            &&
                            ( time2maturities.GetSize() == interpolTypesVect.GetSize() )
                           );

       if (!(sizeControlOK))
       {
          throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                           "Inconsistent vector sizes");
       }

       sizeControlOK = (( deltasCall.GetSize() == volsCall.GetNumCols() )
                         &&
                        ( deltasPut.GetSize()  == volsPut.GetNumCols() )
                       );

       if (!(sizeControlOK))
       {
          throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                           "Insconstent Delta vols matrixes sizes");
       }

       ARM_Vector fxFwds = inFxFwds;

       if ( inFxFwds.GetSize() == 0 )
       {
          // Calculate the Fwds

          int sz =  time2maturities.GetSize();

          fxFwds.Resize(sz);

          for (int i = 0; i < sz; i++)
          {
              ARM_Date maturity = Asof;
      
              double curMatu = time2maturities.Elt(i);

              int nbDays = int(curMatu*365.0);

              maturity.AddDays(nbDays); 

              fxFwds.Elt(i) = CalcFwdFXSpot(Asof,
										    spotFX,
										    maturity,
										    domCurve,
										    ForCurve);
          }
       }

       itsFxFwds = fxFwds;

       what_is_interpolated	= whatIsInterpolated;
       forceLinearWhenSplineFails = correctSplineWithLinear;
       nbTime2Maturities = time2maturities.GetNumLines();

       itsData = new DATA[nbTime2Maturities];

       // should test that all these matrixes have the same number of lines and that all the vectors are columns 
       // except the deltas and that the maturitties and deltas are ordered and that pivotTypes is 0 or 1
       double pivotStrike, pivotVol, pivotDeltaCall, pivotDeltaPut, fxFwd, maturity;
   
       int pivotType;
   
       for (int i = 0; i < nbTime2Maturities; i++)
       {
           int interpolType = int(interpolTypesVect.Elt(i));

	       pivotVol  = pivotVols.Elt(i);
	       fxFwd     = fxFwds.Elt(i);
	       maturity  = time2maturities.Elt(i);
	       pivotType = (int)pivotTypes.Elt(i);

	       itsData[i].time2maturity = maturity;
	       itsData[i].pivotVol      = pivotVol;
	       itsData[i].pivotType     = pivotType;
	       itsData[i].fxFwd         = fxFwd;

	       if ( pivotType == PIVOT_IS_ZDS )
           {
		      ARM_ComputeImpliedZDS_FWP(pivotVol, itsData[i].fxFwd, 
                                    itsData[i].time2maturity, 
                                    pivotStrike, pivotDeltaCall, pivotDeltaPut);
	       }
           else if ( pivotType == PIVOT_IS_ATMF )
	       {
		      ARM_ComputeImpliedATMF_FWP(pivotVol, itsData[i].fxFwd, itsData[i].time2maturity,
                                     pivotStrike, pivotDeltaCall, pivotDeltaPut);
	       }
	       else if ( pivotType == PIVOT_IS_ATMF_ANDNOTWP )
	       {
		      ARM_ComputeImpliedATMF(pivotVol, itsData[i].fxFwd, itsData[i].time2maturity, 
                                  pivotStrike, pivotDeltaCall, pivotDeltaPut);
	       }
           else if ( pivotType == PIVOT_IS_ZDS_NOPREM )
           {
               ARM_ComputeImpliedZDSWithoutPremium(pivotVol, itsData[i].fxFwd, 
                                               itsData[i].time2maturity, 
                                               pivotStrike, pivotDeltaCall, pivotDeltaPut);
           }
           else if ( pivotType == PIVOT_IS_ZDS_SPOT_WPREM )
           {
               ARM_ComputeImpliedZDSWithPremiumSpot(pivotVol, itsData[i].fxFwd, spotFX, 
                                                itsData[i].time2maturity,
                                                domCurve,
                                                ForCurve,
                                                pivotStrike, 
                                                pivotDeltaCall, 
                                                pivotDeltaPut);
           }
           else if ( pivotType == PIVOT_IS_ZDS_SPOT_NOPREM )
           {
               ARM_ComputeImpliedZDSWithoutPremiumSpot(pivotVol, itsData[i].fxFwd, spotFX, 
                                                   itsData[i].time2maturity,
                                                   domCurve,
                                                   ForCurve,
                                                   pivotStrike, 
                                                   pivotDeltaCall, 
                                                   pivotDeltaPut);
           }
           else
           {
              throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                              "Smiled FX vol.: Unexpected volatility Pivot type");
           }

           ARM_Vector* callLine = volsCall.GetRow(i);
           ARM_Vector* putLine  = volsPut.GetRow(i);

	       set_a_DATA(itsData[i], deltasCall, 
                      ARM_Vector(callLine), 
                      deltasPut, 
                      ARM_Vector(putLine), 
                      pivotStrike, pivotVol, pivotDeltaCall, pivotDeltaPut, interpolType);

           delete callLine;
           delete putLine;
       }
}



void ARM_FXVolCurve::GenerateFXVols(ARM_Matrix& volsCall,
									ARM_Matrix& volsPut,
									double spotFX,
									int inRRSTR)
{
		if (itsData)
		{
   			delete[] itsData;     
			itsData = NULL;                           
			nbTime2Maturities = 0;
		}

		itsRRSTRInputFlag = inRRSTR;

		// Manage case of FX Vols == 1e-3 (ATM)
		// WARNING : here we test if the FXVol is smiled
		// in order to normalize the underlying price in the case of ARM_FXVOLAT(whatever is the model: Multi3F or 2IRFX)
		int volIsSmiled = 0;
		// If at least one vol value is > 10E-3, then the FXVol is smiled.
		for (int k = 0; k < volsPut.GetNumLines(); k++)
		{
           for (int h = 0; h < volsPut.GetNumCols(); h++)
           {
               if (( volsPut.Elt(k, h) < 0.0 ) || ( volsPut.Elt(k, h) > 1e-3 ))
               {
                  volIsSmiled = 1;
               }
			   else
			   {
                  volsPut.Elt(k, h) = 0.0;
			   }

               if (( volsCall.Elt(k, h) < 0.0 ) || ( volsCall.Elt(k, h) > 1e-3 ))
               {
                  volIsSmiled = 1;
               }
			   else
			   {
                  volsCall.Elt(k, h) = 0.0;
			   }
           }
       }
	   itsSmileFlag = volIsSmiled;

       if (inRRSTR)
       {
          itsRR  = volsPut;
          itsSTR = volsCall;

          ARM_Vector volsPivots = itsPivotVols;

          FromRRSTRtoVol(&volsPut, &volsCall, &volsPivots);
       }
 
       itsVolsPut       = volsPut;
       itsVolsCall      = volsCall;


       int sizeControlOK = (( itsDeltasCall.GetSize() == itsVolsCall.GetNumCols() )
                           &&
                            ( itsDeltasPut.GetSize()  == itsVolsPut.GetNumCols() )
                            );

       if (!(sizeControlOK))
       {
          throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                           "Insconstent Delta vols matrixes sizes");
       }

	   nbTime2Maturities = itsOptionsMatus.GetNumLines();
       itsData = new DATA[nbTime2Maturities];

       // should test that all these matrixes have the same number of lines and that all the vectors are columns 
       // except the deltas and that the maturitties and deltas are ordered and that pivotTypes is 0 or 1
       double pivotStrike, pivotVol, pivotDeltaCall, pivotDeltaPut, fxFwd, maturity;
   
       int pivotType;
   
       for (int i = 0; i < nbTime2Maturities; i++)
       {
           int interpolType = int(itsInterpolTypes.Elt(i));

	       pivotVol  = itsPivotVols.Elt(i);
	       fxFwd     = itsFxFwds.Elt(i);
	       maturity  = itsOptionsMatus.Elt(i);
	       pivotType = (int)itsPivotTypes.Elt(i);

	       itsData[i].time2maturity = maturity;
	       itsData[i].pivotVol      = pivotVol;
	       itsData[i].pivotType     = pivotType;
	       itsData[i].fxFwd         = fxFwd;

	       if ( pivotType == PIVOT_IS_ZDS )
           {
		      ARM_ComputeImpliedZDS_FWP(pivotVol, itsData[i].fxFwd, 
                                        itsData[i].time2maturity, 
                                        pivotStrike, pivotDeltaCall, pivotDeltaPut);
	       }    
           else if ( pivotType == PIVOT_IS_ATMF )
	       {
		      ARM_ComputeImpliedATMF_FWP(pivotVol, itsData[i].fxFwd, itsData[i].time2maturity,
                                         pivotStrike, pivotDeltaCall, pivotDeltaPut);
	       }
	       else if ( pivotType == PIVOT_IS_ATMF_ANDNOTWP )
	       {
		       ARM_ComputeImpliedATMF(pivotVol, itsData[i].fxFwd, itsData[i].time2maturity, 
                                      pivotStrike, pivotDeltaCall, pivotDeltaPut);
	       }
           else if ( pivotType == PIVOT_IS_ZDS_NOPREM )
           {
               ARM_ComputeImpliedZDSWithoutPremium(pivotVol, itsData[i].fxFwd, 
                                                   itsData[i].time2maturity, 
                                                   pivotStrike, pivotDeltaCall, pivotDeltaPut);
           }
           else if ( pivotType == PIVOT_IS_ZDS_SPOT_WPREM )
           {
               ARM_ComputeImpliedZDSWithPremiumSpot(pivotVol, itsData[i].fxFwd, spotFX, 
                                                    itsData[i].time2maturity,
                                                    itsDomCurve,
                                                    itsForCurve,
                                                    pivotStrike, 
                                                    pivotDeltaCall, 
                                                    pivotDeltaPut);
           }
           else if ( pivotType == PIVOT_IS_ZDS_SPOT_NOPREM )
           {
               ARM_ComputeImpliedZDSWithoutPremiumSpot(pivotVol, itsData[i].fxFwd, spotFX, 
                                                       itsData[i].time2maturity,
                                                       itsDomCurve,
                                                       itsForCurve,
                                                       pivotStrike, 
                                                       pivotDeltaCall, 
                                                       pivotDeltaPut);
           }
           else
           {
              throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                              "Smiled FX vol.: Unexpected volatility Pivot type");
           }
	       
           ARM_Vector* callLine = volsCall.GetRow(i);
           ARM_Vector* putLine  = volsPut.GetRow(i);

	       set_a_DATA(itsData[i], itsDeltasCall, ARM_Vector(callLine), 
                      itsDeltasPut, ARM_Vector(putLine), pivotStrike, 
					  pivotVol, pivotDeltaCall, pivotDeltaPut, interpolType);

           delete callLine;
           delete putLine;
       }

	   delete itsCalcFXATM;
	   itsCalcFXATM = NULL;

	   itsCalcFXATM = ComputeATMFxVol();
}



void ARM_FXVolCurve::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];
    char strDate[50];
 
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
 
       (void) unlink(fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

    fprintf(fOut, "\n =====> Smiled FX Volatilty matrix Volatility matrix  \n");

 
    GetAsOfDate().JulianToStrDate(strDate);
    fprintf(fOut, "\n AsOfDate  : %s \n\n", strDate);


    
	int sz = itsOptionsMatus.GetSize();

    fprintf(fOut, "Maturities\tPivotsVols\tpivotTypes\tInterpTypes\tFxForwards\n\n");

    int i;
    for (i = 0; i < sz; i++)
    {
        char strInterpol[20];

        if ( int(itsInterpolTypes[i]) == K_SPLINE )
        {
           strcpy(strInterpol, "SPLINE");
        }
        else
        {
           strcpy(strInterpol, "LINEAR");
        }

        fprintf(fOut, "%10.5lf\t%10.4lf\t%10d\t%11s\t%14.8lf\n", 
                                                        itsOptionsMatus[i],
                                                        itsPivotVols[i],
                                                        int(itsPivotTypes[i]),
                                                        strInterpol,
                                                        itsFxFwds[i]);
    }

    fprintf(fOut, "\n\n ----> Delta PUT Inputs: \n\n");
    
    int nbLin = itsVolsPut.GetNumLines();
    int nbCol = itsVolsPut.GetNumCols();

    int j;

    // Delta Put

    fprintf(fOut,"Delta PUT\n");
    for (i = 0; i < nbCol; i++)
    {
        if ( i == 0 )
        {
           fprintf(fOut,"\t\t%10.4lf\t", itsDeltasPut[i]);    
        }
        else
        {
           fprintf(fOut,"%10.4lf\t", itsDeltasPut[i]);
        }
    }


    fprintf(fOut,"\n\n");

    for (i = 0; i < nbLin; i++)
    {
        for (j = 0; j < nbCol; j++)
        {
            if ( j == 0 )
            {
               fprintf(fOut,"%10.5lf\t", itsOptionsMatus[i]);
            
               fprintf(fOut,"\t%10.4lf", itsVolsPut.Elt(i, j));
            }
            else
            {
               fprintf(fOut,"\t%10.4lf", itsVolsPut.Elt(i, j));
            }   
        }

        fprintf(fOut, "\n");
    }

    fprintf(fOut, "\n\n ----> Delta CALL Inputs: \n\n");
    
    nbLin = itsVolsCall.GetNumLines();
    nbCol = itsVolsCall.GetNumCols();

    // Delta Call

    fprintf(fOut,"Delta CALL\n");
    for (i = 0; i < nbCol; i++)
    {
        if ( i == 0 )
        {
           fprintf(fOut,"\t\t%10.4lf\t", itsDeltasCall[i]);
        }
        else
        {
           fprintf(fOut,"%10.4lf\t", itsDeltasCall[i]);
        }
    }

    fprintf(fOut,"\n\n");

    for (i = 0; i < nbLin; i++)
    {
        for (j = 0; j < nbCol; j++)
        {
            if ( j == 0 )
            {
               fprintf(fOut,"%10.5lf\t", itsOptionsMatus[i]);
            
               fprintf(fOut,"\t%10.4lf", itsVolsCall.Elt(i, j));
            }
            else
            {
               fprintf(fOut,"\t%10.4lf", itsVolsCall.Elt(i, j));
            }   
        }

        fprintf(fOut, "\n");
    }

    fprintf(fOut, "\n\n----> ATM FX VOL\n");

    ARM_VolLInterpol* ATMFXVol = ComputeATMFxVol();

    ATMFXVol->View(id, fOut);
    delete ATMFXVol;

    if (itsRRSTRInputFlag)
    {
       fprintf(fOut, "\n\n ===============> INPUTS are in terms of : RR and STR\n\n");

       fprintf(fOut, "\n\n ====> RR \n\n");
 
       for (i = 0; i < itsDeltasPut.GetSize(); i++)
       {
           if ( i == 0 )
           {
              fprintf(fOut," %8.4lf\t", itsDeltasPut[i]);    
           }
           else
           {
              fprintf(fOut,"%8.4lf\t", itsDeltasPut[i]);
           }
       }

       fprintf(fOut, "\n\n");

       for (i = 0; i < itsRR.GetNumLines(); i++)
       {
           for (j = 0; j < itsRR.GetNumCols(); j++)
           {
               fprintf(fOut, "%8.4lf\t", itsRR.Elt(i, j));
           }

           fprintf(fOut, "\n");
       }

       fprintf(fOut, "\n\n ====> STR \n\n");
 
       for (i = 0; i < itsDeltasCall.GetSize(); i++)
       {
           if ( i == 0 )
           {
              fprintf(fOut," %8.4lf\t", itsDeltasCall[i]);    
           }
           else
           {
              fprintf(fOut,"%8.4lf\t", itsDeltasCall[i]);
           }
       }

       fprintf(fOut, "\n\n");

       for (i = 0; i < itsSTR.GetNumLines(); i++)
       {
           for (j = 0; j < itsSTR.GetNumCols(); j++)
           {        
               fprintf(fOut, "%8.4lf\t", itsSTR.Elt(i, j));
           }

           fprintf(fOut, "\n");
       }
    }
    else
    {
       fprintf(fOut, "\n\n ===============> INPUTS are in terms of : DELTA PUT, DELTA CALL\n\n");
    }

    fprintf(fOut, "\n\n*------------> TOTEM(depending on strikes) \n\n");

    ARM_Matrix TOTEMVols = get_info();

    fprintf(fOut, "\n              -10       -25         50        25       10\n\n");

    for (i = 0; i < TOTEMVols.GetNumLines(); i++)
    {
        for (j = 0; j < TOTEMVols.GetNumCols(); j++)
        {
            fprintf(fOut, "%8.4lf  ", TOTEMVols.Elt(i, j));
        }

        fprintf(fOut, "\n");
    }

    // Volatilities

    fprintf(fOut, "\n\n*------------> TOTEM(depending on volatilities) \n\n");

    int infoVol = 1;

    TOTEMVols = get_info(infoVol);

    fprintf(fOut, "\n              -10       -25         50        25       10\n\n");

    for (i = 0; i < TOTEMVols.GetNumLines(); i++)
    {
        for (j = 0; j < TOTEMVols.GetNumCols(); j++)
        {
            fprintf(fOut, "%8.4lf  ", TOTEMVols.Elt(i, j));
        }

        fprintf(fOut, "\n");
    }
    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}



ARM_VolLInterpol* ARM_FXVolCurve::ComputeATMFxVol(double lag)
{
    ARM_VolLInterpol* ATMVolFX = NULL;


    if (itsCalcFXATM)
    {
       ATMVolFX = (ARM_VolLInterpol *) itsCalcFXATM->Clone();

       return(ATMVolFX);
    }
    
    ARM_Vector ATMVols(nbTime2Maturities);

    
    double ATMFwdStrike = 1.0; // ATM
	ARM_Vector* ATMVolsYearTerms = (ARM_Vector *) itsOptionsMatus.Clone();
    for (int i = 0; i < nbTime2Maturities; i++)
    {
		(*ATMVolsYearTerms)[i] = (*ATMVolsYearTerms)[i]-lag;
        ATMVols[i] = computeVol(ATMFwdStrike, itsOptionsMatus[i]);
    }

    
    ARM_Vector* undMatu          = new ARM_Vector(1);

    ARM_Matrix* theVols          = new ARM_Matrix(ATMVols);
        
    ATMVolFX = new ARM_VolLInterpol(GetAsOfDate(), 
                                    ATMVolsYearTerms,
                                    undMatu,
                                    theVols);

    return(ATMVolFX);
}



double ARM_FXVolCurve::VolatilityFunction(double m1, double m2)
{
    ARM_VolLInterpol* ATMVolFX = NULL;
    double ATMVol;


    if (itsCalcFXATM)
    {
       ATMVolFX = itsCalcFXATM;

       ATMVol   = ATMVolFX->ComputeVolatility(m1, m2);
    }
    else
    {
       ATMVolFX = ComputeATMFxVol();
       ATMVol   = ATMVolFX->ComputeVolatility(m1, m2);

       if (ATMVolFX)
          delete ATMVolFX;
    }

    return(ATMVol);
}



double ARM_FXVolCurve::VolatilityFunction(double m1, double K, double m2)
{
    double volinterp;


    if (( K == 0.0 ) && (!( m2 == 0.0 ))) // Here we are AT THE MONEY
    {
       K = m2;
    }
    else
    {
       m2 = 0.0; // m2 : is not relevant
    }

    volinterp = ARM_FXVolCurve::VolatilityFunction(m1, K);

    return(volinterp);
}



void ARM_FXVolCurve::BitwiseCopy(const ARM_Object* srcVolCurve)
{
    ARM_FXVolCurve* vCurve = (ARM_FXVolCurve *) srcVolCurve;

    
    if (itsForCcy)
    {
       delete itsForCcy;

       itsForCcy = NULL;
    }

    if (vCurve->itsForCcy)
	   itsForCcy = (ARM_Currency *) vCurve->itsForCcy->Clone();
	
    if (itsDomCcy)
    {
       delete itsDomCcy;

       itsDomCcy = NULL;
    }

    if (vCurve->itsDomCcy)
       itsDomCcy		 = (ARM_Currency *) vCurve->itsDomCcy->Clone();
	
    if (itsDomCurve)
    {
       delete itsDomCurve;
       itsDomCurve = NULL;
    }

    if (vCurve->itsDomCurve)
       itsDomCurve = (ARM_ZeroCurve *) vCurve->itsDomCurve->Clone();

    if (itsForCurve)
    {
       delete itsForCurve;
       itsForCurve = NULL;
    }

    if (vCurve->itsForCurve)
       itsForCurve = (ARM_ZeroCurve *) vCurve->itsForCurve->Clone();

    itsOptionsMatus  = vCurve->itsOptionsMatus;
    itsPivotVols     = vCurve->itsPivotVols;
    itsPivotTypes    = vCurve->itsPivotTypes;
    itsInterpolTypes = vCurve->itsInterpolTypes;
    itsDeltasCall    = vCurve->itsDeltasCall;
    itsVolsCall      = vCurve->itsVolsCall;
    itsDeltasPut     = vCurve->itsDeltasPut;
    itsVolsPut       = vCurve->itsVolsPut;
    itsFxFwds        = vCurve->itsFxFwds;
	itsSpotFX	     = vCurve->itsSpotFX;
    itsRR            = vCurve->itsRR;
    itsSTR           = vCurve->itsSTR;

    nbTime2Maturities          = vCurve->nbTime2Maturities;
    what_is_interpolated       = vCurve->what_is_interpolated;
    forceLinearWhenSplineFails = vCurve->forceLinearWhenSplineFails;

    itsRRSTRInputFlag = vCurve->itsRRSTRInputFlag;

    if (itsData)
    {
       delete [] itsData;

       itsData = NULL;
    }

    if (vCurve->itsData)
    {
        itsData = new DATA[vCurve->nbTime2Maturities];
    
        for (int i = 0; i < vCurve->nbTime2Maturities; i++)
        {
            itsData[i].time2maturity = vCurve->itsData[i].time2maturity;
            itsData[i].pivotVol      = vCurve->itsData[i].pivotVol;
	        itsData[i].pivotType     = vCurve->itsData[i].pivotType;
            itsData[i].pivotStrike   = vCurve->itsData[i].pivotStrike;
	        itsData[i].fxFwd         = vCurve->itsData[i].fxFwd;
            itsData[i].volCall       = vCurve->itsData[i].volCall;
            itsData[i].deltaCall     = vCurve->itsData[i].deltaCall; 
            itsData[i].volPut        = vCurve->itsData[i].volPut;		
	        itsData[i].deltaPut      = vCurve->itsData[i].deltaPut; 
            itsData[i].strikeCall    = vCurve->itsData[i].strikeCall;
            itsData[i].strikePut     = vCurve->itsData[i].strikePut;
            itsData[i].interpol_type = vCurve->itsData[i].interpol_type;
	        itsData[i].strikeAll     = vCurve->itsData[i].strikeAll;
	        itsData[i].volAll        = vCurve->itsData[i].volAll;	
        }
    }    

    itsSmileFlag = vCurve->itsSmileFlag;

    if (itsCalcFXATM)
    {
       delete itsCalcFXATM;

       itsCalcFXATM = NULL;
    }

    if (vCurve->itsCalcFXATM)
    {
       itsCalcFXATM = (ARM_VolLInterpol *) vCurve->itsCalcFXATM->Clone();
    }
}



void ARM_FXVolCurve::get_relevant_maturities(double maturity, int& low, int& high) 
{
	int i = 0;

    int found = 0;

    while (( i < nbTime2Maturities )
           &&
           (!(found))
          )
    {
        if ( maturity < itsData[i].time2maturity )
        {
           found = 1;
        }
        else
        {
           i++;
        }
    }
	
    if ( i == 0 )
    {
       low=high=0;
       return;
    }
	
    if ( i == nbTime2Maturities )
	{
	   low = high=(nbTime2Maturities-1);
	}
	else
	{
	   low  = i-1;
	   high = i;
	}
}



void ARM_FXVolCurve::set_a_DATA(DATA& aData, 
						        const ARM_Vector& deltasCall, 
						        const ARM_Vector& volsCall, 
						        const ARM_Vector& deltasPut, 
						        const ARM_Vector& volsPut,
						        double pivotStrike,
						        double pivotVol,
						        double pivotDeltaCall,
						        double pivotDeltaPut,
						        int interpolType) 
{
	int nbC, nbP;

    nbC = deltasCall.GetNumLines();
	nbP = deltasPut.GetNumLines();
	
    aData.volCall   = ARM_Vector(nbC+1);
	aData.deltaCall = ARM_Vector(nbC+1);
	aData.volPut    = ARM_Vector(nbP+1);
	aData.deltaPut  = ARM_Vector(nbP+1);
	aData.interpol_type = interpolType;
	
	// THE CALLS, need mirroring
	for (int i = 0; i < nbC; i++)
    {
		aData.volCall.Elt(i)   = volsCall.Elt(nbC-1-i);
		aData.deltaCall.Elt(i) = deltasCall.Elt(nbC-1-i);
	}

	aData.volCall.Elt(nbC)   = pivotVol;
	aData.deltaCall.Elt(nbC) = pivotDeltaCall;
	
	// THE PUTS, need mirroring
	aData.volPut.Elt(0)   = pivotVol;
	aData.deltaPut.Elt(0) = pivotDeltaPut;
	
    for (int i = 0; i < nbP; i++)
	{
		aData.volPut.Elt(i+1)   = volsPut.Elt(nbP-1-i);
        aData.deltaPut.Elt(i+1) = deltasPut.Elt(nbP-1-i);
	}

	aData.pivotStrike = pivotStrike;
	aData.pivotVol = pivotVol;

	set_strikes(aData);


	if ( what_is_interpolated == FXINTERP_DELTA_SMOOTH ) //pseudo smoothing of the connection of left/put sida and right/call side
	{
		ARM_Vector extendedVolCall(nbC+2);
		ARM_Vector extendedVolPut(nbP+2);
		ARM_Vector extendedDeltaCall(nbC+2);
		ARM_Vector extendedDeltaPut(nbP+2);

		for (int i = 0; i <= nbP; i++)
		{
			extendedVolPut.Elt(i+1)   = aData.volPut.Elt(i);
			extendedDeltaPut.Elt(i+1) = aData.deltaPut.Elt(i);
		}

		extendedDeltaPut.Elt(0) = 0.99*pivotDeltaPut;
		
        for (int i = 0; i <= nbC; i++)
		{
			extendedVolCall.Elt(i)   = aData.volCall.Elt(i);
			extendedDeltaCall.Elt(i) = aData.deltaCall.Elt(i);
		}

		extendedDeltaCall.Elt(nbC+1) = 1.01*pivotDeltaCall;
		aData.volCall=extendedVolCall;
		aData.deltaCall=extendedDeltaCall;
		aData.volPut=extendedVolPut;
		aData.deltaPut=extendedDeltaPut;

		double extrapolC, extrapolP, deltaForceC, deltaForceP;
		aData.volCall.Elt(nbC+1) = extrapolP = aData.volPut.Elt(1);//start point
		aData.volPut.Elt(0)      = extrapolC = aData.volCall.Elt(nbC);//start point
		
		int nbIter = 0, nbIterMax = 1000;
		double tolerance=0.05;  //no need for a lot of precision, target is just to smooth the junction
		
		do
		{
			aData.volCall.Elt(nbC+1) = fabs(extrapolP+aData.volCall.Elt(nbC+1))/2;
			aData.volPut.Elt(0) = fabs(extrapolC+aData.volPut.Elt(0))/2; //is fabs valid???
			
            set_strikes(aData);

			if (( aData.pivotType == PIVOT_IS_ATMF_ANDNOTWP )  // not with premium
			    ||
                ( aData.pivotType == PIVOT_IS_ZDS_NOPREM )
                ||
                ( aData.pivotType == PIVOT_IS_ZDS_SPOT_NOPREM )
               )
            {
				deltaForceC = aData.deltaPut.Elt(0)+1;
				deltaForceP = aData.deltaCall.Elt(nbC+1)-1;
			}
			else // with premium
			{
				deltaForceC = aData.deltaPut.Elt(0)+aData.strikePut.Elt(0)/aData.fxFwd;
				deltaForceP = aData.deltaCall.Elt(nbC+1)-aData.strikeCall.Elt(nbC+1)/aData.fxFwd;
			}

			extrapolC = ComputeSigmaFromDeltaSIMPLEX(deltaForceC, &aData.deltaCall, &aData.volCall, aData.interpol_type);
			extrapolP = ComputeSigmaFromDeltaSIMPLEX(deltaForceP, &aData.deltaPut, &aData.volPut, aData.interpol_type);
			
            nbIter++;

		} 
        while ( nbIter < nbIterMax && (fabs(aData.volCall.Elt(nbC+1)-extrapolP) > tolerance || fabs(aData.volPut.Elt(0)-extrapolC) > tolerance));
   
        if ( nbIter == nbIterMax )
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "ARM_FXVolCurve::set_a_DATA, smoothing does not work");
        }

	}//end pseudo-smoothing

	if (true) // what_is_interpolated==FXINTERPOL_STRIKE)
	{
		aData.volAll = ARM_Vector(nbC+nbP+1);
		aData.strikeAll = ARM_Vector(nbC+nbP+1);

		for (int i = 0; i < nbP; i++) 
		{
			aData.volAll.Elt(i)    = aData.volPut.Elt(nbP-i);
			aData.strikeAll.Elt(i) = aData.strikePut.Elt(nbP-i);
		}

 		aData.volAll.Elt(nbP)    = aData.pivotVol;
		aData.strikeAll.Elt(nbP) = aData.pivotStrike;
		
        for (int i = 1; i <= nbC; i++)
		{
			aData.volAll.Elt(nbP+i)=aData.volCall.Elt(nbC-i);
			aData.strikeAll.Elt(nbP+i)=aData.strikeCall.Elt(nbC-i);
		}
	}
}



double ARM_FXVolCurve::computeVol(double moneyness, double maturity)
{
    try
    {
	    int low, high;

	    get_relevant_maturities(maturity, low, high);
	    
        double fwd_low, fwd_high, strike_low, strike_high; 
        double vol_low, vol_high, maturity_low, maturity_high, vol;
	    
	    fwd_low    = itsData[low].fxFwd;
	    strike_low = moneyness*fwd_low;
	    vol_low    = computeVol_givenMaturity(strike_low, low);

	    maturity_low = itsData[low].time2maturity;

        if ( low == high )
	    {
	       vol = vol_low;
	    }
	    else
	    {
	       fwd_high = itsData[high].fxFwd;
	       
           strike_high = moneyness*fwd_high;
	       
           vol_high = computeVol_givenMaturity(strike_high,high);
	       
           maturity_high = itsData[high].time2maturity;
	       
           vol = ((maturity_high-maturity)*vol_low+(maturity-maturity_low)*vol_high)
                 /(maturity_high-maturity_low);
	    }

	    return(vol);
    }

    catch(...)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in : ARM_FXVolCurve::computeVol");
    }
}



double ARM_FXVolCurve::computeVol_givenMaturity(double K, int iMaturity)
{
 	double res;

    // out of the range strikes are handled a la mano and extrapolated flat
    // (on the very small 
    // strikes any vol is good numerically (put delta fwd with premim null)
    // so that the algo does not give the 
    // expected solution

    int nbC, nbP;

    nbC = itsData[iMaturity].strikeCall.GetNumLines();
    nbP = itsData[iMaturity].strikePut.GetNumLines();
	
    if ( K < itsData[iMaturity].strikePut.Elt(nbP-1) )
    {
        res = itsData[iMaturity].volPut.Elt(nbP-1);
        
        return(res);
    }

    if ( K > itsData[iMaturity].strikeCall.Elt(0) )
    {
        res = itsData[iMaturity].volCall.Elt(0);
        
        return(res);
    }

    // in the range case
	if ( what_is_interpolated == FXINTERP_STRIKE )
	{
	   if ( itsData[iMaturity].interpol_type == K_LINEAR )
       {
          res = linInterpol2(&itsData[iMaturity].strikeAll, K, 
                             &itsData[iMaturity].volAll);
       }
       else
       {
          res = SplineInterpolateFunc(&itsData[iMaturity].strikeAll, &itsData[iMaturity].volAll,
                                      K,
                                      NULL, // SecondDerivCalc,
                                      1);   // keep2Der 
       }
	}
	else
	{
       if ( itsData[iMaturity].pivotType == PIVOT_IS_ZDS_SPOT_WPREM )
       {
          res = ComputeVolWithDeltaSpotAsInputWithPremium(K, 
                                                          itsData[iMaturity].fxFwd, 
                                                          itsSpotFX,
                                                          itsData[iMaturity].time2maturity,
                                                          &itsData[iMaturity].deltaCall, 
                                                          &itsData[iMaturity].volCall,
                                                          &itsData[iMaturity].deltaPut, 
                                                          &itsData[iMaturity].volPut,
                                                          itsData[iMaturity].pivotStrike, 
                                                          itsData[iMaturity].pivotVol,// guessed vol, 
                                                          itsData[iMaturity].interpol_type, 
                                                          forceLinearWhenSplineFails,
                                                          itsDomCurve,
                                                          itsForCurve);
          
       }
       else if ( itsData[iMaturity].pivotType == PIVOT_IS_ZDS_SPOT_NOPREM )
       {
           res = ComputeVolWithDeltaSpotAsInputWithoutPremium(K,
                                                              itsData[iMaturity].fxFwd,
                                                              itsSpotFX,
                                                              itsData[iMaturity].time2maturity,
                                                              &itsData[iMaturity].deltaCall,
                                                              &itsData[iMaturity].volCall,
                                                              &itsData[iMaturity].deltaPut, 
                                                              &itsData[iMaturity].volPut,
                                                              itsData[iMaturity].pivotStrike,
                                                              itsData[iMaturity].interpol_type,
                                                              itsDomCurve,
                                                              itsForCurve);
       }
       else if (( itsData[iMaturity].pivotType != PIVOT_IS_ATMF_ANDNOTWP )
                &&
                (itsData[iMaturity].pivotType != PIVOT_IS_ZDS_NOPREM )
               )
       {
		  res = computeVol_withDeltaFWPasInput_aux(K,
				                                   itsData[iMaturity].fxFwd,
					                               itsData[iMaturity].time2maturity, 
					                               &itsData[iMaturity].deltaCall,
					                               &itsData[iMaturity].volCall,
					                               &itsData[iMaturity].deltaPut,
					                               &itsData[iMaturity].volPut,
					                               itsData[iMaturity].pivotStrike,
										           itsData[iMaturity].pivotVol,//guessed vol
										           itsData[iMaturity].interpol_type,
                                                   forceLinearWhenSplineFails);
       }
	   else
       {
		  res = computeVol_withDeltaAsInput_aux(K,
					                            itsData[iMaturity].fxFwd,
					                            itsData[iMaturity].time2maturity, 
					                            &itsData[iMaturity].deltaCall,
					                            &itsData[iMaturity].volCall,
					                            &itsData[iMaturity].deltaPut,
					                            &itsData[iMaturity].volPut,
					                            itsData[iMaturity].pivotStrike,
											    itsData[iMaturity].interpol_type);
       }
	}

	return(res);
}



void ARM_FXVolCurve::set_strikes(DATA& aData)

{
    // First call strikes

    int nbC, nbP, k;

    nbC = aData.deltaCall.GetNumLines();
    aData.strikeCall = ARM_Vector(nbC);


    if ( aData.pivotType == PIVOT_IS_ZDS_SPOT_WPREM )
    {
       for (k = 0; k < nbC; k++)
       {
           ARM_ComputeImpliedStrikeFromSpotDeltaWithPremium(aData.volCall.Elt(k), 
                                                            aData.fxFwd,
                                                            itsSpotFX,
                                                            aData.time2maturity, 
                                                            aData.deltaCall.Elt(k), 
                                                            K_CALL,
                                                            itsDomCurve,
                                                            itsForCurve,
                                                            aData.strikeCall.Elt(k));
       }
    }
    else if ( aData.pivotType == PIVOT_IS_ZDS_SPOT_NOPREM )
    {
       for (k = 0; k < nbC; k++)
       {
           ARM_ComputeImpliedStrikeFromSpotDeltaWithoutPremium(aData.volCall.Elt(k), 
                                                               aData.fxFwd,
                                                               itsSpotFX,
                                                               aData.time2maturity, 
                                                               aData.deltaCall.Elt(k), 
                                                               K_CALL,
                                                               itsDomCurve,
                                                               itsForCurve,
                                                               aData.strikeCall.Elt(k));
       }
    }
    else if (( aData.pivotType == PIVOT_IS_ATMF_ANDNOTWP )
             ||
             ( aData.pivotType == PIVOT_IS_ZDS_NOPREM )
            )
    {
       for (k = 0; k < nbC; k++)
       {
           ARM_ComputeImpliedStrike(aData.volCall.Elt(k), 
                                    aData.fxFwd,
                                    aData.time2maturity,
                                    aData.deltaCall.Elt(k),
                                    K_CALL,
                                    aData.strikeCall.Elt(k));
       }
    }
    else
    {
       for (k = 0; k < nbC; k++)
       {
           ARM_ComputeImpliedStrike_FWP(aData.volCall.Elt(k), 
                                     aData.fxFwd,aData.time2maturity,
                                     aData.deltaCall.Elt(k),
                                     K_CALL,
                                     aData.strikeCall.Elt(k));
       }
    }

    // idem for puts
    nbP = aData.deltaPut.GetNumLines();
    aData.strikePut = ARM_Vector(nbP);

    if ( aData.pivotType == PIVOT_IS_ZDS_SPOT_WPREM )
    {
       for (k = 0; k < nbP; k++)
       {
           ARM_ComputeImpliedStrikeFromSpotDeltaWithPremium(aData.volPut.Elt(k), 
                                                            aData.fxFwd,
                                                            itsSpotFX,
                                                            aData.time2maturity, 
                                                            aData.deltaPut.Elt(k), 
                                                            K_PUT,
                                                            itsDomCurve,
                                                            itsForCurve,
                                                            aData.strikePut.Elt(k));
       }
    }
    else if ( aData.pivotType == PIVOT_IS_ZDS_SPOT_NOPREM )
    {
       for (k = 0; k < nbP; k++)
       {
           ARM_ComputeImpliedStrikeFromSpotDeltaWithoutPremium(aData.volPut.Elt(k), 
                                                               aData.fxFwd,
                                                               itsSpotFX,
                                                               aData.time2maturity, 
                                                               aData.deltaPut.Elt(k), 
                                                               K_PUT,
                                                               itsDomCurve,
                                                               itsForCurve,
                                                               aData.strikePut.Elt(k));
       }
    }
    else if (( aData.pivotType == PIVOT_IS_ATMF_ANDNOTWP )
             ||
             ( aData.pivotType == PIVOT_IS_ZDS_NOPREM )
            )
    {
       for (k = 0; k <nbP; k++)
       {
          ARM_ComputeImpliedStrike(aData.volPut.Elt(k), aData.fxFwd, 
                                   aData.time2maturity,
                                   aData.deltaPut.Elt(k),
                                   K_PUT,
                                   aData.strikePut.Elt(k));
       }
    }
    else
    {
       for (k = 0; k < nbP; k++)
       {
          ARM_ComputeImpliedStrike_FWP(aData.volPut.Elt(k),
                                       aData.fxFwd,aData.time2maturity,
                                       aData.deltaPut.Elt(k),
                                       K_PUT,
                                       aData.strikePut.Elt(k));
       }
    }
}



ARM_Matrix ARM_FXVolCurve::get_info(int infoVol)
{
	int nbK, nbM;
	nbM = nbTime2Maturities;
	nbK = itsData[0].volAll.GetNumLines(); 

    ARM_Matrix res(nbM,nbK+1);
	
    for (int i = 0; i < nbM; i++)
	{
		for (int j = 0; j < nbK; j++)
        {
            if (infoVol)
            {
               res.Elt(i, j+1) = itsData[i].volAll.Elt(j);
            }
            else
            {
			   res.Elt(i, j+1) = itsData[i].strikeAll.Elt(j);
            }
        }

		res.Elt(i,0) = itsData[i].time2maturity;
	}

	return(res);
}



void ARM_FXVolCurve::BumpFxVol(	ARM_Vector			newVol,
								ARM_ZeroCurve*		domCurve,
								ARM_ZeroCurve*		ForCurve){

	if ( newVol.GetSize() != itsPivotVols.GetSize() )	
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "the Input Vector is not of the appropriate size");
			
	for (int i = 0; i < itsPivotVols.GetSize(); i++)
		itsPivotVols[i] = newVol[i];

	BumpFxVol(0.0,0,0,itsSpotFX,domCurve,ForCurve);
         
}



void ARM_FXVolCurve::BumpFxVol(double shiftValue,
                               int firstLine,
                               int lastLine,
                               double spotFX,
                               ARM_ZeroCurve* domCurve,
                               ARM_ZeroCurve* ForCurve)
{
     try
     {
         if (itsCalcFXATM)
         {
            delete itsCalcFXATM;

            itsCalcFXATM = NULL;
         }

         for (int i = firstLine; i <= lastLine; i++)
         {
             itsPivotVols[i] += shiftValue;
         }

         ARM_Vector time2maturities   = itsOptionsMatus;
         ARM_Vector pivotVols         = itsPivotVols;
         ARM_Vector pivotTypes        = itsPivotTypes;
         ARM_Vector deltasCall        = itsDeltasCall;
         ARM_Matrix volsCall          = itsVolsCall;
         ARM_Vector deltasPut         = itsDeltasPut;
         ARM_Matrix volsPut           = itsVolsPut;
         ARM_Vector inFxFwds          = itsFxFwds;
         ARM_Vector interpolTypesVect = itsInterpolTypes;

         ARM_Matrix RR                = itsRR;
         ARM_Matrix STR               = itsSTR;

         int RRSTRInputFlag = itsRRSTRInputFlag;

         int inNbTime2Maturities      =  nbTime2Maturities;

         int inWhatIsInterpolated     = what_is_interpolated;
         int correctSplineWithLinear  = forceLinearWhenSplineFails;
     
         GenerateFXVols(GetAsOfDate(),
                        time2maturities, 
					    pivotVols, 
					    pivotTypes,
					    deltasCall,
                        ( RRSTRInputFlag ? STR : volsCall ) ,
					    deltasPut,
                        ( RRSTRInputFlag ? RR : volsPut ), 
					    inFxFwds,
					    interpolTypesVect,
					    inWhatIsInterpolated,
					    correctSplineWithLinear,
                        spotFX,
                        domCurve,
                        ForCurve,
                        RRSTRInputFlag);
     }

	 catch(Exception& e)
	 {
		 throw e;
	 }

     catch(...)
     {
         throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Error in : ARM_FXVolCurve::BumpFxVol");
     }

}



// Here, UP TO NOW!, we must not regenerate the ATM Vol
// contrary to : BumpFxVol where the ATM is bumped so must be 
// generated in ComputeATMFxVol() (after being deleted in BumpFxVol)

void ARM_FXVolCurve::BumpRRorSTR(const vector<int>& matusX,
                                 const vector<int>& deltasY,
                                 const vector<double>& shiftValues,
                                 double spotFX,
                                 ARM_ZeroCurve* domCurve,
                                 ARM_ZeroCurve* ForCurve)
{
     try
     {
         int sz = shiftValues.size();
         
         int RR_size = itsRR.GetNumCols();

         for (int i = 0; i < sz; i++)
         {
             if ( deltasY[i] == 0 ) // Bumping ATM Vol is not allowed in this mode!
             {
                throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Error in : ARM_FXVolCurve::BumpRRorSTR: Bumping ATM Vol is not allowed in this mode!");
             }
           
             if ( deltasY[i] <= RR_size ) // Bump RR 
             {
                int shiftIndex = deltasY[i]-1;

                itsRR.Elt(matusX[i], shiftIndex) += shiftValues[i];
             }
             else // Bump STR
             {
                int shiftIndex = deltasY[i]-RR_size-1; 

                itsSTR.Elt(matusX[i], shiftIndex) += shiftValues[i];
             }
         }

         ARM_Vector time2maturities   = itsOptionsMatus;
         ARM_Vector pivotVols         = itsPivotVols;
         ARM_Vector pivotTypes        = itsPivotTypes;
         ARM_Vector deltasCall        = itsDeltasCall;
         ARM_Matrix volsCall          = itsSTR;
         ARM_Vector deltasPut         = itsDeltasPut;
         ARM_Matrix volsPut           = itsRR;
         ARM_Vector inFxFwds          = itsFxFwds;
         ARM_Vector interpolTypesVect = itsInterpolTypes;

         ARM_Matrix RR                = itsRR;
         ARM_Matrix STR               = itsSTR;

         int RRSTRInputFlag           = itsRRSTRInputFlag;

         int inNbTime2Maturities      =  nbTime2Maturities;

         int inWhatIsInterpolated     = what_is_interpolated;
         int correctSplineWithLinear  = forceLinearWhenSplineFails;

     
         GenerateFXVols(GetAsOfDate(),
                        time2maturities, 
					    pivotVols, 
					    pivotTypes,
					    deltasCall,
                        ( RRSTRInputFlag ? STR : volsCall ) ,
					    deltasPut,
                        ( RRSTRInputFlag ? RR : volsPut ), 
					    inFxFwds,
					    interpolTypesVect,
					    inWhatIsInterpolated,
					    correctSplineWithLinear,
                        spotFX,
                        domCurve,
                        ForCurve,
                        RRSTRInputFlag);
     }

	 catch(Exception& e)
	 {
		 throw e;
	 }

     catch(...)
     {
         throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Error in : ARM_FXVolCurve::BumpRRorSTR");
     }
}



void ARM_FXVolCurve::FXBumpRRorSTR(	double shiftValue, 
									int nbRow, 
									int nbCol, 
									double spotFX,
									int isRR,
									int isCumul,
									int isAbsolute)
{
	try
	{
		if (itsCalcFXATM)
		{
			delete itsCalcFXATM;

			itsCalcFXATM = NULL;
		}

		int nLines = itsRR.GetNumLines();
		int nCols  = itsRR.GetNumCols();

		if (( nbRow > nLines ) || ( nbCol > nCols ))
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                     "Invalid nthLine or nthCol in <FXBumpRRorSTR> method");
		}


		int first;
		int last;

		ARM_Matrix RRorSTR = isRR == K_YES ? itsRR : itsSTR ;

		last = RRorSTR.GetNumLines() - 1;
		first = (nbRow == 0) ? -1 : last - nbRow;
			
		if ( isCumul == K_YES )
		{
			for (int i = last; i > first; i--)
			{ 
				if ( nbCol != 0 )
				{
					if (isAbsolute == K_YES)
						RRorSTR.Elt(i,nbCol-1) += shiftValue;
					else
						RRorSTR.Elt(i,nbCol-1) *= (1+shiftValue/100);
				} 
				else
				{
					if ( isAbsolute == K_YES )
					{
						RRorSTR.Elt(i,0) += shiftValue;
						RRorSTR.Elt(i,1) += shiftValue;
					}
					else
					{
						RRorSTR.Elt(i,0) *= (1+shiftValue/100);
						RRorSTR.Elt(i,1) *= (1+shiftValue/100);
					}
				}
			}
		}
		else
		{
			if ( nbCol != 0 )
			{
				if (isAbsolute == K_YES)
					RRorSTR.Elt(last - nbRow + 1,nbCol-1) += shiftValue;
				else
					RRorSTR.Elt(last - nbRow + 1,nbCol-1) *= (1+shiftValue/100);

			} 
			else
			{
				if ( isAbsolute == K_YES )
				{
					RRorSTR.Elt(last - nbRow + 1,0) += shiftValue;
					RRorSTR.Elt(last - nbRow + 1,1) += shiftValue;
				}
				else
				{
					RRorSTR.Elt(last - nbRow + 1,0) *= (1+shiftValue/100);
					RRorSTR.Elt(last - nbRow + 1,1) *= (1+shiftValue/100);
				}
			}
		}

		ARM_Matrix RR  = (isRR == K_YES) ? RRorSTR : itsRR;
		ARM_Matrix STR = (isRR == K_YES) ? itsSTR  : RRorSTR;

		GenerateFXVols(STR, RR, spotFX);
     }

	 catch(Exception& e)
	 {
		 throw e;
	 }

     catch(...)
     {
         throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Error in : ARM_FXVolCurve::BumpFxVol");
     }
}



void ARM_FXVolCurve::FXBumpRRorSTR(	ARM_Matrix mR, 	int isRR)
{
	ARM_Matrix* tmp= NULL;

	try
    {
		if ( isRR == K_YES) 
        {
			if ( mR.GetNumLines() != itsRR.GetNumLines() ) 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Error in : The input volatilty do not have the rigth dimension");

			if ( mR.GetNumCols() != itsRR.GetNumCols() ) 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Error in : The input volatilty do not have the rigth dimension");


			tmp = (ARM_Matrix *) (&itsSTR)->Clone();
			GenerateFXVols(*tmp, mR, itsSpotFX);
		}
		else
        {
			if ( mR.GetNumLines() != itsSTR.GetNumLines()) 
			   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
               "Error in : The input volatilty do not have the rigth dimension");

			if ( mR.GetNumCols() != itsSTR.GetNumCols()) 
			   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
               "Error in : The input volatilty do not have the rigth dimension");

			tmp = (ARM_Matrix *) (&itsRR)->Clone();
			
            GenerateFXVols(mR, *tmp, itsSpotFX);
		}

		if (tmp) 
        { 
           delete tmp; 
           tmp = NULL;
        }
	}

	 catch(Exception& e)
	 {
		 throw e;
	 }

     catch(...)
     {
         throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Error in : The input volatilty do not have the right dimension");
     }
}

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/