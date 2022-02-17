/*
 * $Log: volint.cpp,v $
 * Revision 1.28  2004/03/29 10:20:55  mab
 * Correction in vol conversion: LogNor <---> Nor
 *
 * Revision 1.27  2004/03/18 15:26:40  mab
 * ApproximatedForwardRate : displaced in module : utile: fromto.h(.cpp)
 *
 * Revision 1.26  2004/02/20 15:19:57  jpriaudel
 * SetCurrencyUnit added in ConvertToNormalVol
 *
 * Revision 1.25  2003/11/20 15:20:29  ebenhamou
 * change view for lastknown date
 *
 * Revision 1.24  2003/11/14 11:29:19  mab
 * Take in account Splined FX VOL
 *
 * Revision 1.23  2003/10/22 16:39:35  jpriaudel
 * modif pour Olivier
 *
 * Revision 1.22  2003/10/22 16:19:45  jpriaudel
 * correction dans convertToBsVol
 *
 * Revision 1.21  2003/09/25 10:36:05  mab
 *  A flag added to : ConvertToBSVol()
 *
 * Revision 1.20  2003/09/23 09:33:50  mab
 * Added: ConvertToBSVol(..)
 *
 * Revision 1.19  2003/09/22 12:39:35  mab
 * added: const in constructors parameters
 * corrections in : ApproximatedForwardRate()
 * ImpliedLogNorIRVol()
 *
 * Revision 1.15  2003/08/18 09:49:51  mab
 * improvement in the principal method :
 * virtual double VolatilityFunction(double m1, double K, double m2)
 * virtual double VolatilityFunction(double m1, double m2) :
 * No default parameter but 2 methods!
 *
 * Revision 1.14  2003/02/11 15:01:27  mab
 * void UpdateCol(ARM_Vector* pCol, double tenor);
 * ARM_Vector* GetCol(double tenor);
 *
 * Revision 1.13  2002/08/01 07:48:14  mab
 * Added : fprintf(fOut, "\n\t Currency : %s \n\n\n",
 *                      GetCurrency()->GetCcyName());
 * in The View() Method
 *
 * Revision 1.12  2001/04/27 09:26:22  smysona
 * makeup
 *
 * Revision 1.11  2001/04/20 13:48:08  mab
 * Amelioration du View
 *
 * Revision 1.10  2001/04/20 13:43:53  mab
 * message en Anglais dans le View
 *
 * Revision 1.9  2001/04/20 07:52:17  mab
 * amelioration du View
 *
 * Revision 1.8  2001/04/20 07:49:33  mab
 * Amelioration du view
 *
 * Revision 1.7  2001/04/18 18:10:41  abizid
 * Modif View
 *
 * Revision 1.6  2001/04/18 16:33:59  abizid
 * Ajout du View
 *
 * Revision 1.5  2001/04/03 11:58:53  nicolasm
 * unused variables.
 *
 * Revision 1.4  1999/10/01 13:39:38  nicolasm
 * Suppression SetVolType dans le constructeur
 *
 * Revision 1.3  1999/09/07 08:15:48  nicolasm
 * Ajout champ volType dans le constructeur principal
 *
 * Revision 1.2  1999/04/14 13:28:53  mab
 * Rajout Log pour RCS
 *
 */


/*----------------------------------------------------------------------------*
 
    volint.cpp
 
    This file implements the ARM_VolLInterpol class, a class for 
         computing a ARM_VolCurve linear interpolation.

*----------------------------------------------------------------------------*/



#ifdef unix
#include <sys/types.h>
#include <unistd.h>
#endif

#include "firsttoinc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util.h"
#include "volint.h"
#include "volfxspinterp.h"
#include "expt.h"


#include "volflat.h"

//#include "swap.h"


//#include <ICMKernel/util/icm_macro.h>

/*----------------------------------------------------------------------------*/

void ARM_VolLInterpol::Init(void)
{
    SetName(ARM_VOL_LIN_INTERPOL);

    itsStrikes   = NULL;
}


ARM_VolLInterpol::ARM_VolLInterpol(void)
{
    Init();

    itsStrikes   = NULL;
}


ARM_VolLInterpol::ARM_VolLInterpol(const ARM_Date& asOf,
                                   ARM_Vector* yearTerms,
                                   ARM_Vector* strikes,
                                   ARM_Matrix* volatilities,
                                   int KType, int volType,
                                   ARM_Currency* ccy):ARM_VolCurve(asOf, 
                                                                   KType,
                                                                   volType)
{
    // check variables
    Init();

    SetCurrencyUnit(ccy);

    if ( yearTerms->GetSize() != volatilities->GetNumLines() )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
     "The input year terms size and volatility number of lines must be equal");
    }

    if ( strikes->GetSize() != volatilities->GetNumCols() )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
      "The input strikes size and volatility number of columns must be equal");
    }

    // set variables
    SetVolatilities(volatilities);

    SetExpiryTerms(yearTerms);

    SetStrikes(strikes);

    if ( volType == K_FX_VOL_SP_INTERP )
    {
       ARM_FXVolSmileInterpol* fxVolInterpol = NULL;

       fxVolInterpol = new ARM_FXVolSmileInterpol(this);

       SetFxVolSmileInterpolation(fxVolInterpol);
    }
}



ARM_VolLInterpol::ARM_VolLInterpol(const ARM_Date& asOf,
                                   ARM_Vector* yearTerms,
                                   ARM_Vector* volatilities,
                                   int KType,
                                   int volType,
                                   ARM_Currency* ccy) 
                                   : ARM_VolCurve(asOf, KType, volType)
{
    Init();

    // check variables

    SetCurrencyUnit(ccy);

    if ( yearTerms->GetSize() != volatilities->GetSize() )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
     "The input year terms size and volatility number of lines must be equal");
    }


    // Set variables

    itsStrikes = new ARM_Vector(1, 1.0);

    ARM_Matrix* vols = new ARM_Matrix(*volatilities);

    SetVolatilities(vols);

    SetExpiryTerms(yearTerms);

    if ( volType == K_FX_VOL_SP_INTERP )
    {
       ARM_FXVolSmileInterpol* fxVolInterpol = NULL;

       fxVolInterpol = new ARM_FXVolSmileInterpol(this);

       SetFxVolSmileInterpolation(fxVolInterpol);
    }
}



ARM_VolLInterpol::ARM_VolLInterpol(ARM_VolFlat* volFlat, 
                                   ARM_Vector* expiries,
                                   ARM_Vector* undTenors)
{
    int expiriesSize = 0;
    int tenorsSize   = 0;


    if ( expiries == NULL )
    {
       expiriesSize = 2;

       ARM_Vector* tmpExp = new ARM_Vector(expiriesSize);

       tmpExp->Elt(0) = 0.0;
       tmpExp->Elt(1) = 1000.0;
   
       SetExpiryTerms(tmpExp);
    }
    else
    {
       expiriesSize = expiries->GetSize();

       SetExpiryTerms((ARM_Vector *) expiries->Clone());
    }

    if ( undTenors == NULL )
    {
       tenorsSize = 2;

       ARM_Vector* tmpStrikes = new ARM_Vector(tenorsSize);

       tmpStrikes->Elt(0) = 0.0;
       tmpStrikes->Elt(1) = 1000.0;

       itsStrikes = tmpStrikes;
    }
    else
    {
       tenorsSize = undTenors->GetSize();

       itsStrikes = (ARM_Vector *) undTenors->Clone();
    }

    ARM_Matrix* vols = new ARM_Matrix(expiriesSize, tenorsSize, volFlat->GetVolatility());

    SetVolatilities(vols);

    SetAsOfDate(volFlat->GetAsOfDate());

    SetCurrencyUnit(volFlat->GetCurrency());
}



ARM_VolLInterpol::ARM_VolLInterpol(const ARM_VolLInterpol& volLInterpol)
                                  :ARM_VolCurve(volLInterpol)
{
    Init();

    BitwiseCopy(&volLInterpol);
}



ARM_VolLInterpol& ARM_VolLInterpol::operator = 
                                    (const ARM_VolLInterpol& volLInterpol)
{
    (*this).ARM_VolCurve::operator = (volLInterpol);

    BitwiseCopy(&volLInterpol);

    return(*this);
}

    
ARM_VolLInterpol::~ARM_VolLInterpol(void)
{
    if (itsStrikes)
       delete itsStrikes;
}


void ARM_VolLInterpol::BitwiseCopy(const ARM_Object* srcVollint)
{
    ARM_VolLInterpol* vollint = (ARM_VolLInterpol *) srcVollint;

    if (itsStrikes)
    {
       delete itsStrikes;
       itsStrikes = NULL;
    }

    if (vollint->itsStrikes)
       itsStrikes = (ARM_Vector *) vollint->itsStrikes->Clone();
}


void ARM_VolLInterpol::Copy(const ARM_Object* vollint)
{
    ARM_VolCurve::Copy(vollint);

    BitwiseCopy(vollint);
}


ARM_Object* ARM_VolLInterpol::Clone(void)
{
    ARM_VolLInterpol* theClone = new ARM_VolLInterpol();

    theClone->Copy(this);

    return(theClone);
}


void ARM_VolLInterpol::Set(int szyt, double* yt, int szstk, double* stk, 
							 double* vol, char* date, int strikeType)
{
    ARM_Vector* terms = new ARM_Vector(szyt, yt);

    SetExpiryTerms(terms);

    ARM_Vector* strikes = new ARM_Vector(szstk, stk);

    SetStrikes(strikes);
     
    ARM_Matrix* vols = new ARM_Matrix(szyt, szstk, vol);

    SetVolatilities(vols);
     
    SetAsOfDate((ARM_Date) date);

    SetStrikeType(strikeType);
}


ARM_Vector* ARM_VolLInterpol::GetStrikes(void)
{
    return(itsStrikes);
}

ARM_Vector* ARM_VolLInterpol::GetStrikes() const
{
    return itsStrikes;
}

void ARM_VolLInterpol::SetStrikes(ARM_Vector* strikes)
{
    if ( itsStrikes == strikes )
       return;

    if (itsStrikes)
    {
       delete itsStrikes;
       itsStrikes = NULL;
    }

    if (strikes)
       itsStrikes = strikes;
}

inline double 
linInterpol2Matrix(const ARM_Vector&strike,const ARM_Vector&terms,const ARM_Matrix& vols,double m,double K)
{
	if (strike.empty()) throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,"linInterpol2Matrix: no strikes"); 
	if (terms.empty()) throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,"linInterpol2Matrix: no terms"); 

    double yval, lambda;
    int    index=0;

    
    // Get index of point right before xval
    // or the appropriate index if xval does not belong to (xmin, xmax)

	if(strike.size()==1)
		// return vc.Elt(0); 
		return linInterpol2Col(terms, m, vols,index); 

    if ( K < strike[0]-K_DOUBLE_TOL ) 
    {
       //index = 0;
	   return linInterpol2Col(terms, m, vols,0);
    }
    else if ( K > strike[strike.size()-1]+K_DOUBLE_TOL) 
    {
//       index = strike.size()-2;
	   return linInterpol2Col(terms, m, vols,strike.size()-1);
    }
    else 
    {
       index = indexBeforeValue(strike.GetElt(),strike.size(),K);
    }

    //  do linear interpolation

    if ( fabs (K - strike[index]) < K_DOUBLE_TOL) 
    {
       return // vc(index) ((yElt)[index]);
		   linInterpol2Col(terms, m, vols,index); 
    }

	
    lambda = (strike[index+1]-K)/(strike[index+1]-strike[index]);

    // yval = lambda*(yElt)[index]+(1.0-lambda)*(yElt)[index+1];
	yval = lambda*  linInterpol2Col(terms, m, vols,index) 
		+ (1.0-lambda)* linInterpol2Col(terms, m, vols,index+1) ; 

    return(yval);

}

double ARM_VolLInterpol::VolatilityFunction(double m, double K, double m2)
{
    double volinterp;

    int volType = GetVolType();

    ARM_Matrix* vols = GetVolatilities();

    ARM_Vector* terms = GetExpiryTerms();

    int nCols = vols->GetNumCols();

    if (( K == 0.0 ) && (!( m2 == 0.0 ))) // Here we are AT THE MONEY
    {
       K = m2;
    }
    else
    {
       m2 = 0.0; // m2 : is not relevant
    }

    if ( GetInterpType() == K_DIAG_INTERPOL )
    {
       volinterp = triangularInterpol(terms, 
                                      itsStrikes,
                                      vols, 
                                      m, 
                                      K);

       return(volinterp);
    }


	ARM_Vector vc (nCols, 0.0);

	/* First interpolate for each strike  */
	/** JLA : optimized. 
	for (int i = 0; i < nCols; i++)
	{
		ARM_Vector vectcol (*vols, i);

	   vc.Elt(i) = linInterpol2(terms, m, &vectcol); 
	}

	if (itsStrikes->GetSize() > 1 )
	{
	   volinterp = linInterpol2(itsStrikes, K , &vc);
	}
	else
	{
	   volinterp = vc.Elt(0);
	}
	*/ 
	volinterp = linInterpol2Matrix(*itsStrikes,*terms,*vols,m,K); 
	

    return(volinterp);
}


double ARM_VolLInterpol::VolatilityFunction(double m1, double m2)
{
    // We are in the case of a Matrix
    // the Third param is not relevant
    return(VolatilityFunction(m1, m2, 0.0));
}


void ARM_VolLInterpol::View(char* id, FILE* ficOut)
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
    fprintf(fOut, "\n Volatility matrix \n\n");

    ARM_AbstractMarketClass::View(id, fOut);


    if ( GetInterpType() == K_DIAG_INTERPOL )
    {
       fprintf(fOut, "\n\t Interpolation Mode : DIAG \n\n");
    }
    else
	{
		if ( GetInterpType() == K_SPLINE )
		{
			fprintf(fOut, "\n\t Interpolation Mode : SPLINE \n\n");
		}
		else
		{
			fprintf(fOut, "\n\t Interpolation Mode : LINEAR \n\n");
		}
	}
    switch(GetVolType())
    {
        case K_ATMF_VOL:
        {
            fprintf(fOut, " Volatility type : ATM");
        };
        break;

        case K_FX_VOL :
        {
            fprintf(fOut, " Volatility type : FX");
        };
        break;

        default:
        {
            fprintf(fOut, " Volatility type : SMILE");
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

    GetAsOfDate().JulianToStrDateDay(strDate);
    fprintf(fOut, " AsOfDate\t : %s \n", strDate);
	if (!GetIndexName().empty()) fprintf(fOut, "\n IndexName  : %s \n", GetIndexName().c_str());
	if (hasIndex()) fprintf(fOut, "\n Index is defined. \n\n" );
	if( GetAsOfDate() != GetLastKnownDate() )
    {
		GetLastKnownDate().JulianToStrDateDay(strDate);
		fprintf(fOut, " Last Known Date : %s \n", strDate);
	}

    fprintf(fOut, "\n\n" );

    ARM_Vector* terms = GetExpiryTerms();
    ARM_Matrix* vols  = GetVolatilities();
    int nbLg  = terms->GetSize();
    int nbCol = itsStrikes->GetSize();

    int i, j;

    fprintf(fOut,"                     ");
    for (j = 0; j < nbCol; j++)
    {
        fprintf(fOut," %10.5lf", itsStrikes->Elt(j));
    }
    fprintf(fOut,"\n");
	fprintf(fOut,"                     ");
    for (j = 0; j < nbCol; j++)
		fprintf(fOut," [%8s]", itsYearTermsY[j]) ;
    fprintf(fOut,"\n");
	ARM_Vector* expDates = GetExpiryDates(); 
	if(expDates) 
	{
		fprintf(fOut,"                     ");
		for (j = 0; j < nbCol; j++)
			fprintf(fOut,"%10.5lf", expDates->Elt(j));
		fprintf(fOut,"\n");
	}
    
	fprintf(fOut,"\n");

    for (i = 0; i < nbLg; i++)
    {
		fprintf(fOut," [%8s] ", itsYearTermsX[i] );
        fprintf(fOut," %10.5lf  ", terms->Elt(i));

        for (j = 0; j < nbCol; j++)
        {
            fprintf(fOut," %10.5lf",vols->Elt(i, j));
        }

        fprintf(fOut,"\n");
    }


    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}



ARM_Vector* ARM_VolLInterpol::GetCol(double tenor)
{
    int sizeCol = GetVolatilities()->GetNumLines();

    ARM_Vector* tenorCol = new ARM_Vector(sizeCol);


    long i = indexBeforeValue(GetStrikes(), tenor);

    long iToUpdate;

    if ( i == GetStrikes()->GetSize()-1 )
       iToUpdate = GetStrikes()->GetSize()-1;
    else if (i == -1)
       iToUpdate = 0;
    else if ( fabs(tenor-GetStrikes()->Elt(i))
              < fabs(GetStrikes()->Elt(i+1)-tenor)
            )
        iToUpdate = i;
    else
        iToUpdate = i+1;

    for (int j = 0; j < sizeCol; j++)
        tenorCol->Elt(j) = GetVolatilities()->Elt(j,iToUpdate);

    return(tenorCol);
}



void ARM_VolLInterpol::UpdateCol(ARM_Vector* pCol, double tenor)
{
    if (( pCol == NULL ) || ( GetVolatilities() == NULL ))
       return;

    if ( pCol->GetSize() != GetVolatilities()->GetNumLines() )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                     "Invalid nthLine in <UpdateLine> method");
    }

    long i = indexBeforeValue(GetStrikes(), tenor);

    long iToUpdate;

    if ( i == GetStrikes()->GetSize()-1 )
       iToUpdate = GetStrikes()->GetSize()-1;
    else if (i == -1)
       iToUpdate = 0;
    else if ( fabs(tenor-GetStrikes()->Elt(i)) 
              < fabs(GetStrikes()->Elt(i+1)-tenor)
            )
       iToUpdate = i;
    else
       iToUpdate = i+1;

    for (int j = 0; j < pCol->GetSize(); j++)
        GetVolatilities()->Elt(j,iToUpdate) = pCol->Elt(j);
}



double ARM_VolLInterpol::ImpliedLogNorIRVol(double Mat,
                                            double Tenor,
                                            ARM_ZeroCurve* YieldCurve,
                                            int IsSwaptionVol,
                                            int InPct)
/*
 * This method returns the equivalent Black-Scholes volatility,
 * converted out of a normal interest rate volatility matrix.
 * The method interpolates linearly in the absolute volatility matrix, then converts
 * the interpolated volatility into lognormal BS volatility.
 */
{
    ARM_VolLInterpol* This = this;

    double Sigma = This->ComputeVolatility(Mat, Tenor);

    double CutoffMat = 1.0/(double) YieldCurve->GetCurrencyUnit()->GetFixedPayFreq();

    int IsSwapFwd = (IsSwaptionVol) ? (Tenor >= CutoffMat) : (Tenor > 1.0);

    double F = ApproximatedForwardRate(Mat, Tenor, YieldCurve, IsSwapFwd);

    F /= 100.0; // In ARM the fwd are in base 100
    
    if ( F > 0 )
    {
       double z = ((Sigma/F*sqrt(Mat*0.5/PI))+1.0)/2.0;

       z = CDFInvNormal(z)*2.0/sqrt(Mat);

       return((InPct) ? z * 100.0 : z);
    }

    return(0.0);
}



ARM_VolLInterpol* ARM_VolLInterpol::ConvertToNormalVol(ARM_ZeroCurve* YieldCurve,
                                                       int IsSwaptionVol,
                                                       int InPct)
/*
 * This method returns a vol surface containing a normal (absolute) volatility surface.
 * The normal surface is converted out of the surface contained in this instance
 * and of the yield curve passed in argument (to obtain the forward rates)
 * This instance is supposed to contain a lognormal ATM Interest Rate Volatility surface.
 * If the InPct argument is nonzero, then the lognormal surface is assumed in percent,
 * i.e., 1 means 1%. Otherwise, the surface is assumed unitary, i.e., 1 means 100%.
 */
{
    int Rows, Cols;
    ARM_Matrix* BSVolSurface;

    ARM_Date asOfDate = YieldCurve->GetAsOfDate();

    int spotDays  = this->GetCurrency()->GetSpotDays(); 

    char payCal[50];

    this->GetCurrency()->CalcFloatPayCal(payCal);

    char* ccyName = this->GetCurrency()->GetCcyName(); 


    // If we have abscissae, ordinates, the surface that goes with those, and
    // if we passed a yield curve:

    if ((this->GetExpiryTerms()) 
        && (Rows = this->GetExpiryTerms()->GetSize())
        && (this->GetStrikes())
        && (Cols = this->GetStrikes()->GetSize()) 
        && (BSVolSurface = this->GetVolatilities())
        && (YieldCurve)
       )
    {
        // Local copy of this->Absicssae and this->Ordinates:
        ARM_Vector* Abscissae = (ARM_Vector *) this->GetExpiryTerms()->Clone();
        ARM_Vector* Ordinates = (ARM_Vector *) this->GetStrikes()->Clone();

        // The matrix for the normal vol:
        ARM_Matrix* NormalSurface = new ARM_Matrix(Rows, Cols);

        unsigned int i, j, IsSwapFwd;
        double T, Tenor, F, SigBS, SigHW;

        double CutoffMat = 1.0
                       /(double) YieldCurve->GetCurrencyUnit()->GetFixedPayFreq();

        for (i = 0; ( i < Rows ); ++i)
        {
            for (j = 0; ( j < Cols ); ++j)
            {
                // Maturity: In fact a swap Start Date en years fraction
                T = Abscissae->Elt(i);

                // Convert T to a swaption Expiry

                T = T*365.0;

                ARM_Date tmpDate = asOfDate;

                int nbDays = ROUND(T);

                ARM_Date newDate = tmpDate.AddDays(nbDays);

                // old code newDate.PreviousBusinessDay(spotDays, ccyName);
                // In order to be more coherent with : ApproximatedForwardRate

                ARM_Date spotDate = asOfDate;

                spotDate = spotDate.NextBusinessDay(spotDays, payCal);

                T = (newDate-spotDate)/365.0;
                
                Tenor = Ordinates->Elt(j);

                // Type of the forward rate (MM or Par Swap):
                IsSwapFwd = (IsSwaptionVol) ? (Tenor >= CutoffMat) : (Tenor > 1.0);

                // BS Sigma:
                if (InPct)
                   SigBS = BSVolSurface->Elt(i, j)/100.0;
                else
                   SigBS = BSVolSurface->Elt(i, j);

                // Fwd Rate:
                F = ApproximatedForwardRate(T, Tenor, YieldCurve, IsSwapFwd);

                F /= 100.0;

                // Absolute Sigma:
                SigHW = F*sqrt(2.0*PI/T)*(2.0*cdfNormal(SigBS*sqrt(T)/2.0)-1.0);

                // Matrix Update:
                NormalSurface->Elt(i, j) = SigHW;
            }
        }

        // Creation of the returned instance:
        ARM_VolLInterpol* RV = new ARM_VolLInterpol(GetAsOfDate(),
                                                    Abscissae, Ordinates,
                                                    NormalSurface);

		memcpy(RV->itsYearTermsX, itsYearTermsX,
			sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);

		memcpy(RV->itsYearTermsY, itsYearTermsY,
			sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);

        RV->SetCurrencyUnit(GetCurrency());

		// Update Mkt data characteristics
		string	vType("VOL SWOPT ABS");
		string	vIndex = GetExternalIndex();
		string	vCurrency = GetStrCurrency();
		string	vCrvName = GetExternalCrvId();

		RV->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvName, vType);

        return(RV);
    }

    return((ARM_VolLInterpol *) NULL);
}



ARM_VolLInterpol* ARM_VolLInterpol::ConvertToBSVol(ARM_ZeroCurve* YieldCurve,
                                                   int IsSwaptionVol,
                                                   int InPct,
                                                   int Post100)
/*
 * by Olivier Croissant
 * inspired by ConvertToHWVol
 */

{
    int Rows, Cols;
    ARM_Matrix* HWVolSurface;

    ARM_Date asOfDate = YieldCurve->GetAsOfDate();

    int spotDays  = this->GetCurrency()->GetSpotDays();

    char* ccyName = this->GetCurrency()->GetCcyName();

    char payCal[50];

    this->GetCurrency()->CalcFloatPayCal(payCal);

    // Pasted code gets hereunder:
    // If we have abscissae, ordinates, the surface that goes with those, and
    // if we passed a yield curve:

    if ((this->GetExpiryTerms())
        && (Rows = this->GetExpiryTerms()->GetSize())
        && (this->GetStrikes())
        && (Cols = this->GetStrikes()->GetSize())
        && (HWVolSurface = this->GetVolatilities())
        && (YieldCurve)
       )
    {
        // Local copy of this->Absicssae and this->Ordinates:

        ARM_Vector* Abscissae = (ARM_Vector *) this->GetExpiryTerms()->Clone();

        ARM_Vector* Ordinates = (ARM_Vector *) this->GetStrikes()->Clone();


        // The matrix for the normal vol:

        ARM_Matrix* BSSurface = new ARM_Matrix(Rows, Cols);

        unsigned int i, j, IsSwapFwd;

        double T, F, SigBS, SigHW;

        double CutoffMat = 1.0/
                    (double) YieldCurve->GetCurrencyUnit()->GetFixedPayFreq();

        for (i = 0; ( i < Rows ); ++i)
        {
            for (j = 0; ( j < Cols ); ++j)
            {
                // Maturity: In fact a swap Start Date en years fraction

                T = Abscissae->Elt(i);

                // Convert T to a swaption Expiry
                T = T*365.0;

                ARM_Date tmpDate = asOfDate;

                int nbDays = ROUND(T);

                ARM_Date newDate = tmpDate.AddDays(nbDays);

                // old code: newDate.PreviousBusinessDay(spotDays, ccyName);
                // T = (newDate-asOfDate)/365.0;

                ARM_Date spotDate = asOfDate;

                spotDate = spotDate.NextBusinessDay(spotDays, payCal);

                T = (newDate-spotDate)/365.0;
                
                double Tenor = Ordinates->Elt(j);

                // Type of the forward rate (MM or Par Swap):

                IsSwapFwd = (IsSwaptionVol) ? ( Tenor >= CutoffMat) : (T > 1.0);


                // BS Sigma:

                if (InPct)
                   SigHW = HWVolSurface->Elt(i, j)/100.0;
                else
                   SigHW = HWVolSurface->Elt(i, j);

                // Fwd Rate:

                F = ApproximatedForwardRate(T, Tenor, 
                                            YieldCurve, IsSwapFwd);
                F /= 100.0;

                // BS Sigma:

                SigBS = (2.0/sqrt(T))
                        *CDFInvNormal((1.0+(SigHW/F)*sqrt(T/(2.0*PI)))/2.0);

                // Matrix Update:
                if (Post100)
                   SigBS *= 100.0;

                BSSurface->Elt(i, j) = SigBS;
            }
        }

        // Creation of the returned instance:

        ARM_VolLInterpol* RV = new ARM_VolLInterpol(GetAsOfDate(),
                                                    Abscissae, Ordinates,
                                                    BSSurface);

        return(RV);
    }

    return((ARM_VolLInterpol *) NULL);
}



double ARM_VolLInterpol::CalcNumericalObjectSignature(void)
{
    double signature = 0.0;

    ARM_Matrix* vols = GetVolatilities();

    int nbLines = vols->GetNumLines();
    int nbCols  = vols->GetNumCols();
    
    ARM_Vector* expiries = GetExpiryTerms(); 
    ARM_Vector* matus    = GetStrikes();

    for (int i = 0; i < nbLines; i++)
    {
        for (int j = 0; j < nbCols; j++)
        {
            signature += (expiries->Elt(i)*matus->Elt(j)*vols->Elt(i, j));

        }
    }

    signature += GetAsOfDate().GetJulian();

    return(signature);
}

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
