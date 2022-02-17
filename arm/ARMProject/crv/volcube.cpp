#include <stdio.h>
#include <stdlib.h>

#ifdef unix
#    include <unistd.h>
#    include <fcntl.h>
#    include <sys/stat.h>
#endif


#include "volcube.h"
//#include "frmutils.h"
#include "volfxspinterp.h"

#define cleanit(__ptr) \
    { if(__ptr)        \
     { delete __ptr;   \
       __ptr = NULL;   \
     }                 \
    }                


ARM_VolCube::ARM_VolCube(void)
{
    Init();
}




/*!
 * Constructor method
 */
ARM_VolCube::ARM_VolCube(ARM_VolCurve** inVols, int size, 
                         ARM_Vector* underlyings)
{
    int i;


    Init();

    // assume that all the vols have consistent parameters

    if ( size > 1 )
    {
        for (i = 1; i < size; i++)
        {
            if (inVols[0]->GetAsOfDate() != inVols[i]->GetAsOfDate())
            {
        		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
					    "Volatility matrices have different asOfDay");
            }

            if (inVols[0]->GetStrikeType() != inVols[i]->GetStrikeType())
            {
        		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
					    "Volatility matrices have different strikeType");
            }

            if (inVols[0]->GetOptionType() != inVols[i]->GetOptionType())
            {
        		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
					    "Volatility matrices have different optionType");
            }

            if (strcmp(inVols[0]->GetCurrency()->GetCcyName(), 
                       inVols[i]->GetCurrency()->GetCcyName()))
            {
        		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
					    "Volatility matrices have different currencies");
            }
        }
    }

    SetAsOfDate(inVols[0]->GetAsOfDate());    // Date of the  volcurve

    SetStrikeType(inVols[0]->GetStrikeType());// Type de Strike (Price or Delta)

    SetOptionType(inVols[0]->GetOptionType());// Type d'option : CAP, SWAPTION

    SetCurrencyUnit(inVols[0]->GetCurrency());

	if (inVols[0]->hasIndex()) SetIndex(inVols[0]->GetIndex()); 

    itsVols = new vector<ARM_VolCurve*>;

    itsVols->resize(size);

    for (i = size-1; i >= 0; i--)
    {
        itsVols->at(i) = (ARM_VolCurve*) inVols[i]->Clone();
    }

    itsUnderlyings = (ARM_Vector*) underlyings->Clone();


	// Update Mkt data characteristics
	ARM_VolCurve*	vVolCurve = itsVols->at(0);
	if(vVolCurve->GetName() == ARM_VOL_FLAT)
		vVolCurve = itsVols->at(1);
	if(vVolCurve)
	{
		string	vType = vVolCurve->GetStrType() + " SMILE";
		string	vIndex = vVolCurve->GetExternalIndex();
		string	vCurrency = vVolCurve->GetStrCurrency();
		string	vCrvName = vVolCurve->GetExternalCrvId();

		SetMktExternalCharacteristics(vIndex, vCurrency, vCrvName, vType);
	}
}


ARM_VolCube::ARM_VolCube(ARM_VolCurve* atmVol, ARM_VolCurve** inVols, 
                         int size, 
                         ARM_Vector* underlyings,
                         int volType, int checkCcy)
{
    int i;


    Init();

    SetVolType(volType);

	itsATMVol = (ARM_VolCurve*) atmVol->Clone();

	itsATMref = true;

    // assume that all the vols have consistent parameters

    if ( size > 1 )
    {
        for(i = 1; i < size; i++)
        {
            if ( inVols[0]->GetAsOfDate() != inVols[i]->GetAsOfDate() )
            {
        	   throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
					    "Volatility matrices have different asOfDay");
            }

            /*if ( inVols[0]->GetStrikeType() != inVols[i]->GetStrikeType() )
            {
        	   throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
					    "Volatility matrices have different strikeType");
            }

            if ( inVols[0]->GetOptionType() != inVols[i]->GetOptionType() )
            {
        	   throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
					    "Volatility matrices have different optionType");
            }*/

            if(checkCcy && (strcmp(inVols[0]->GetCurrency()->GetCcyName(), 
							inVols[i]->GetCurrency()->GetCcyName())))
            {
        		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
					    "Volatility matrices have different currencies");
            }
        }
    }

    SetAsOfDate(inVols[0]->GetAsOfDate());    // Date of the  volcurve

    SetStrikeType(inVols[0]->GetStrikeType());// Type de Strike (Price or Delta)

    SetOptionType(inVols[0]->GetOptionType());// Type d'option : CAP, SWAPTION

    SetCurrencyUnit(inVols[0]->GetCurrency());

	if (inVols[0]->hasIndex()) SetIndex(inVols[0]->GetIndex()); 

    itsVols = new vector<ARM_VolCurve*>;

    itsVols->resize(size);

    for (i = size-1; i >= 0; i--)
    {
        itsVols->at(i) = (ARM_VolCurve*) inVols[i]->Clone();
    }

    itsUnderlyings = (ARM_Vector*) underlyings->Clone();

    if ( volType == K_FX_VOL_SP_INTERP )
    {
       ARM_FXVolSmileInterpol* fxVolInterpol = NULL;

       fxVolInterpol = new ARM_FXVolSmileInterpol(this);

       SetFxVolSmileInterpolation(fxVolInterpol);
    }

	// Update Mkt data characteristics
	ARM_VolCurve*	vVolCurve = atmVol;
	if(vVolCurve->GetName() == ARM_VOL_FLAT)
		vVolCurve = inVols[0];
	if(vVolCurve)
	{
		string	vType = vVolCurve->GetStrType() + " SMILE";
		string	vIndex = vVolCurve->GetExternalIndex();
		string	vCurrency = vVolCurve->GetStrCurrency();
		string	vCrvName = vVolCurve->GetExternalCrvId();

		SetMktExternalCharacteristics(vIndex, vCurrency, vCrvName, vType);
	}
}



ARM_VolCube::ARM_VolCube(vector < ARM_VolCurve* > * inVols, ARM_Vector* underlyings, 
                         const ARM_Date& lastKnownDate )
{
    Init();

    int size=inVols->size();

	if ( size < 1 )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
                         "a volcube need at least one volcurve!");
	
	if ( size != underlyings->size() )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
    "the vector of volcurve should have same size as the vector of corresponding underlying!");

	vector < ARM_VolCurve* >::iterator it = inVols->begin();

    ARM_VolCurve* Vol0=*it;

	ARM_Date theLastKnownDate;

	if ( lastKnownDate == ARM_Date() )
	   theLastKnownDate = Vol0->GetLastKnownDate();
	else
	   theLastKnownDate  = lastKnownDate;
    
	vector < ARM_VolCurve* >::iterator it_end=inVols->end();

    /// check that all the vols have consistent parameters
    if ( size > 1 )
    {
		++it;

	    for (;it!=it_end;++it)
        {
            if (Vol0->GetAsOfDate() != (**it).GetAsOfDate())
            {
                throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
                                 "Volatility matrices have different asOfDay");
            }

            if ( Vol0->GetStrikeType() != (**it).GetStrikeType() )
            {
                throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
                               "Volatility matrices have different strikeType");
            }

            if (Vol0->GetOptionType() != (**it).GetOptionType())
            {
                throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
                               "Volatility matrices have different optionType");
            }

            if (strcmp(Vol0->GetCurrency()->GetCcyName(),
                       (**it).GetCurrency()->GetCcyName()))
            {
                throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
                               "Volatility matrices have different currencies");
            }
        }
    }

    itsVols = new vector<ARM_VolCurve*>(inVols->size());

    vector<ARM_VolCurve*>::iterator it_new= (*itsVols).begin();

    for (it = inVols->begin(); it != it_end; ++it, ++it_new)
    {
        *it_new = (ARM_VolCurve*) (**it).Clone();
    }


    SetAsOfDate(Vol0->GetAsOfDate());    // Date of the  volcurve

    SetStrikeType(Vol0->GetStrikeType());// Type de Strike (Price or Delta)

    SetOptionType(Vol0->GetOptionType());// Type d'option : CAP, SWAPTION

    SetCurrencyUnit(Vol0->GetCurrency());

	if (Vol0->hasIndex()) SetIndex(Vol0->GetIndex()); 

	SetLastKnownDate( theLastKnownDate ); // set the last known date

    itsUnderlyings = (ARM_Vector*) underlyings->Clone();


	// Update Mkt data characteristics
	ARM_VolCurve*	vVolCurve = Vol0;
	if(vVolCurve->GetName() == ARM_VOL_FLAT)
		vVolCurve = itsVols->at(1);
	if(vVolCurve)
	{
		string	vType = vVolCurve->GetStrType() + " SMILE";
		string	vIndex = vVolCurve->GetExternalIndex();
		string	vCurrency = vVolCurve->GetStrCurrency();
		string	vCrvName = vVolCurve->GetExternalCrvId();

		SetMktExternalCharacteristics(vIndex, vCurrency, vCrvName, vType);
	}
}


ARM_VolCube::~ARM_VolCube(void)
{
    cleanit(itsUnderlyings);

	cleanit(itsATMVol);

    if (itsVols)
    {
       vector<ARM_VolCurve*>::iterator it;
        
       for (it = itsVols->begin(); it < itsVols->end(); it++)
       {
           cleanit(*it);
       }
        
        cleanit(itsVols);
        
    }

    cleanit(itsStrikeLevels);
}



void ARM_VolCube::Init(void)
{
    SetName(ARM_VOL_CUBE);

    itsUnderlyings = NULL;
    itsVols      = NULL;
	itsATMVol    = NULL;
	itsATMref    = false;
	itsStickyDelta	= true;
	itsStrikeLevels	= NULL;

    SetVolType(K_CUBE_VOL);     
}


// If used with stickyStrike feature, moneyness will be
// the strike

double ARM_VolCube::VolatilityFunction(double expiry, double moneyness, 
                                       double underlying)
{
    double Vol = ComputeSmileOnly(expiry, moneyness, underlying);

	if (itsATMref)
	{
		// We have interpolated the smile
		// now we add the vol

		Vol += itsATMVol->ComputeVolatility(expiry, underlying);
	}

    return(Vol);
}

double ARM_VolCube::ComputeSmileOnly(double expiry, double moneyness, 
                                       double underlying)
{
    int size = itsUnderlyings->GetSize();

    double Vol, VolUp;


    if (!itsStickyDelta)
       return VolatilityFunctionByStrike(expiry, moneyness,  underlying);

    if ( underlying <= itsUnderlyings->Elt(0) )
    {
       Vol = itsVols->at(0)->ComputeVolatility(expiry, moneyness);
    }
    else if ( underlying >= itsUnderlyings->Elt(size-1) )
    {
       Vol = itsVols->at(size-1)->ComputeVolatility(expiry, moneyness);
    }
    else
    {
       int j = locateIndex(itsUnderlyings, underlying);


       Vol   = itsVols->at(j)->ComputeVolatility(expiry, moneyness);
       VolUp = itsVols->at(j+1)->ComputeVolatility(expiry, moneyness);

       Vol = linInterpol(underlying, itsUnderlyings->Elt(j), Vol,
                         itsUnderlyings->Elt(j+1), VolUp);
    }

    return(Vol);
}


double ARM_VolCube::VolatilityFunctionByStrike(double expiry, double strike, 
                                               double underlying)
{
    int size = itsUnderlyings->GetSize();

    double Vol, VolUp;

    if ( underlying <= itsUnderlyings->Elt(0) )
    {
        // We interpolate smile on the smallest underlying
        // We need to interpolate the strike level for the given expiry

        int k = locateIndex(((*itsVols)[0])->GetExpiryTerms(), expiry);

        strike -= linInterpol(expiry, ((*itsVols)[0])->GetExpiryTerms()->Elt(k),
                              itsStrikeLevels->Elt(0,k),
	                          itsUnderlyings->Elt(0,k+1), 
                              ((*itsVols)[0])->GetExpiryTerms()->Elt(k+1));
        
		Vol = itsVols->at(0)->ComputeVolatility(expiry, strike);
    }
    else if ( underlying >= itsUnderlyings->Elt(size-1) )
    {
        // We interpolate smile on the longest underlying
        // We need to interpolate the strike level for the given expiry

        int k = locateIndex(((*itsVols)[size-1])->GetExpiryTerms(), expiry);

        strike -= linInterpol(expiry, ((*itsVols)[size-1])->GetExpiryTerms()->Elt(k),
                              itsStrikeLevels->Elt(size-1,k),
	                          itsUnderlyings->Elt(size-1,k+1), 
                              ((*itsVols)[size-1])->GetExpiryTerms()->Elt(k+1));

		Vol = itsVols->at(size-1)->ComputeVolatility(expiry, strike);
    }
    else
    {
       int j = locateIndex(itsUnderlyings, underlying);

       int k = locateIndex(((*itsVols)[size-1])->GetExpiryTerms(), expiry);

       double levelDown = linInterpol(expiry, 
                             ((*itsVols)[j])->GetExpiryTerms()->Elt(k),
                             itsStrikeLevels->Elt(j,k),
	                         itsUnderlyings->Elt(j,k+1), 
                             ((*itsVols)[j])->GetExpiryTerms()->Elt(k+1));

       double levelUp   = linInterpol(expiry, 
                             ((*itsVols)[j+1])->GetExpiryTerms()->Elt(k),
                             itsStrikeLevels->Elt(j+1,k),
	                         itsUnderlyings->Elt(j+1,k), 
                             ((*itsVols)[j+1])->GetExpiryTerms()->Elt(k));

	   Vol   = itsVols->at(j)->ComputeVolatility(expiry, strike-levelDown);
	   VolUp = itsVols->at(j+1)->ComputeVolatility(expiry, strike-levelUp);

	   Vol = linInterpol(underlying, itsUnderlyings->Elt(j), Vol,
	                     itsUnderlyings->Elt(j+1), VolUp);
    }

	if (itsATMref)
	{
		// We have interpolated the smile
		// now we add the vol

		Vol += itsATMVol->ComputeVolatility(expiry, underlying);
	}

    return(Vol);
}


double	ARM_VolCube::ComputeCorrelByExpiry(double aExpiry, double aTenor1, double aTenor2)
{
	if(GetVolType() == K_CUBE_CORREL_DIAG)
	{
		return	ComputeVolatility(aTenor1, aTenor2, aExpiry);
	}
	else
	{
		return	ComputeVolatility(aExpiry, aTenor1, aTenor2);
	}
}
/*
double	ARM_VolCube::ComputeHyperCorrel(double aOptMat, double aTenor1,
										  double aTenor2, double aStrike)
{
	return	( ComputeVolatility(aOptMat, aTenor1, aTenor2) );
}
*/
void ARM_VolCube::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char  fOutName[200];

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
    fprintf(fOut, "\n\n >>>>>>>>>>>>>>>>>>>>>> VOLATILITY CUBE <<<<<<<<<<<<<<<<<<<<<<<\n\n");

    ARM_AbstractMarketClass::View(id, fOut);
	if (hasIndex()) fprintf(fOut, "\n Index Defined."); 

    fprintf(fOut, "\n\n\n ======================> THE COMPONENTS <========================\n\n\n");

    if (itsATMref)
    {
       if (itsATMVol)
       {
          fprintf(fOut, "\n\n ====================> ATM Vol.\n\n");

          itsATMVol->View(id, fOut);

          fprintf(fOut, "\n\n <==================== ATM Vol.\n\n");
       } 
    }

    if (itsVols)
    {
       int sz = itsUnderlyings->GetSize();

       for (int i = 0; i < sz; i++)
       {
           fprintf(fOut, "\n\n\n------------------> Tenor : %lf\n", 
                         (*itsUnderlyings)[i]);

           (*itsVols)[i]->View(id, fOut);

        }
    }

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}



void ARM_VolCube::ParallelShift(double value)
{
    if (itsATMVol)
       itsATMVol->ParallelShift(value);
}



void ARM_VolCube::ParallelShiftAll(double value)
{
    if (itsATMVol)
       itsATMVol->ParallelShift(value);

    if (itsVols)
    {
       int sz = itsUnderlyings->GetSize();

       for (int i = 0; i < sz; i++)
       {
           (*itsVols)[i]->ParallelShift(value);
       }
    }    
}



double ARM_VolCube::VolatilityFunction(double expiry, double underlying)
{
    return(VolatilityFunction(expiry, 0.0, underlying));
}


ARM_VolCube::ARM_VolCube(const ARM_VolCube& volCube) : ARM_VolLInterpol(volCube)
{
    Init();

    BitwiseCopy(&volCube);
}


void ARM_VolCube::BitwiseCopy(const ARM_Object* srcObject)
{
    const ARM_VolCube* srcVol = (const ARM_VolCube *) (srcObject);

    cleanit(itsVols);

    if (srcVol->itsVols)
    {
       itsVols = new vector<ARM_VolCurve*>;

       itsVols->resize(srcVol->itsVols->size());

       vector<ARM_VolCurve*>::iterator itA, itB;

       for (itA = itsVols->begin(), itB = srcVol->itsVols->begin();
            itB < srcVol->itsVols->end();
            itA++, itB++)
       {
            *itA = (ARM_VolCurve*) (*itB)->Clone();
       }
    }

    cleanit(itsATMVol);

    if (srcVol->itsATMVol)
       itsATMVol = (ARM_VolCurve*) srcVol->itsATMVol->Clone();

    itsATMref = srcVol->itsATMref;    

    cleanit(itsUnderlyings);

    if (srcVol->itsUnderlyings)
    {
       itsUnderlyings = (ARM_Vector*) srcVol->itsUnderlyings->Clone();
    }
}



void ARM_VolCube::Copy(const ARM_Object* VolIn)
{
    ARM_VolLInterpol::Copy(VolIn);

    BitwiseCopy(VolIn);
}



ARM_Object* ARM_VolCube::Clone(void)
{
    ARM_VolCube* theClone = new ARM_VolCube();

    theClone->Copy(this);

    return(theClone);
}



ARM_VolCurve* ARM_VolCube::GetNthTenorVolMatrix(int i)
{
    ARM_VolCurve* nThVolCurve = NULL;


    nThVolCurve =  itsVols->at(i);

    return(nThVolCurve);
}


ARM_VolCurve*	ARM_VolCube::GetVolCurve(string& aTenor)
{
	if( aTenor == "ATM")
	{
		return	GetATMVol();
	}

	double	vMatu = FromStrMatuToDouble(aTenor.c_str());

	int	vNthTenor = -1;
	ARM_Vector*	tenors = GetUnderLyings();
	if( tenors )
	{
		for( int i=0; i<tenors->GetSize(); i++ )
		{
			if( tenors->Elt(i) == vMatu )
			{
				vNthTenor = i;
				break;
			}
		}
	}

	if( vNthTenor == -1 )
	{
		char	vErrMsg[100];
		string	sErrMsg("ARM_VolCube::GetVolCurve(tenor), unable to get smile for tenor ");
		sErrMsg += aTenor;
		sprintf(vErrMsg, sErrMsg.c_str());

        throw	Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrMsg);
	}

	return	GetNthTenorVolMatrix(vNthTenor);
}


double ARM_VolCube::CalcNumericalObjectSignature(void)
{
    double signature = 0.0;

    if (itsATMVol)
       signature += itsATMVol->CalcNumericalObjectSignature();

    if (itsVols)
    {
       int sz = itsUnderlyings->GetSize();

       for (int i = 0; i < sz; i++)
       {
           signature += (*itsVols)[i]->CalcNumericalObjectSignature(); 
       }
    }

    return(signature);
}



void ARM_VolCube::BumpVolatility(double value,
                                 int nthLine,
                                 int nthCol,
                                 int isCumul,
								 int isAbsolute)
{
    if (itsATMVol)
       itsATMVol->BumpVolatility(value,nthLine,nthCol,isCumul,isAbsolute);
    else
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                "you need an ATM Vol to bump a cube");
}


void ARM_VolCube::BumpSmile(double value, 
							double tenor,
                            int nthLine, 
							int nthCol,
                            int isCumul,
							int isAbsolute)
{
	if(tenor)
	{
		int j = locateIndex(itsUnderlyings, tenor);

		if(j == -1)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
							"Cannot find VolCurve for specified tenor in smile in ARM_VolCube::BumpSmile()");

		itsVols->at(j)->BumpVolatility(value, nthLine, nthCol, isCumul, isAbsolute);
	}
	else
	{
		for( vector<ARM_VolCurve*>::iterator vIter = itsVols->begin(); vIter != itsVols->end(); vIter++ )
		{
			((ARM_VolCurve*)(*vIter))->BumpVolatility(value, nthLine, nthCol, isCumul, isAbsolute);
		}
	}
}


/*---------------------------------------------------------------*/
/*---- End Of File ----*/

