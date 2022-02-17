/*---------------------------------------------------------------------------------------*/
/*                                                                                       */
/*  ARM_HyperCube class methods                                                          */
/*                                                                                       */
/*---------------------------------------------------------------------------------------*/

#include "hypercube.h"
#include "fromto.h"





ARM_HyperCube::ARM_HyperCube(void)
{
	Init();
}


ARM_HyperCube::ARM_HyperCube(vector<ARM_VolCurve*>& aVolCurveList, vector<string>& aTenorList)
{
	Init();

    
	vector<string>::iterator vTenorIter = aTenorList.begin();
	vector<string>::iterator vTenorEnd  = aTenorList.end();

	vector<ARM_VolCurve*>::iterator	vVolCurveIter = aVolCurveList.begin();

	if ( aVolCurveList.size() != aTenorList.size() )
	{
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
			"The VolCurve list and the tenor list don't have the same size in hyper cube");
	}
 
	ARM_VolCurve* vVolCurve  = NULL;
	ARM_VolCurve* vVolCurve0 = *vVolCurveIter;
	
    if ( vVolCurve0 == NULL )
	{
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
					   "First vol cube is NULL in hyper cube");
	}

    // Set the currency

    SetCurrencyUnit(aVolCurveList[0]->GetCurrency());

	double	 vYearTerm;
	ARM_Date vAsOfDate = vVolCurve0->GetAsOfDate();

	for (; vTenorIter != vTenorEnd; vTenorIter++, vVolCurveIter++)
	{
		vVolCurve = *vVolCurveIter;

		if ( vVolCurve == NULL )
		{
		   throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
							"One of the vol cubes is NULL in hyper cube");
		}

		if ( vVolCurve->GetAsOfDate() != vAsOfDate )
		{
		   throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
							"Some vol cubes don't have the same AsOfDate in hyper cube");
		}

        ARM_VolCurve* theClone = (ARM_VolCurve *) vVolCurve->Clone();
	
		itsVolCurves.insert( pair<string, ARM_VolCurve*> (*vTenorIter, theClone) );

		// On convertit en nb de mois
		vYearTerm = 12. * StringMatuToYearTerm((char*)(*vTenorIter).c_str());

		itsTenorList.push_back((int)vYearTerm);
	}

	sort(itsTenorList.begin(), itsTenorList.end());
	vector<int>::iterator	vFirstDuplicate = unique(itsTenorList.begin(), itsTenorList.end());
	itsTenorList.erase(vFirstDuplicate, itsTenorList.end());

    SetAsOfDate(vAsOfDate);

	SetMktExternalCharacteristics(vVolCurve0->GetExternalIndex(), 
								  vVolCurve0->GetStrCurrency(),
								  vVolCurve0->GetExternalCrvId(),
								  string("Hyper Cube"));	// TMP : type à voir
}



ARM_HyperCube::~ARM_HyperCube(void)
{
	std::map<string, ARM_VolCurve*>::iterator vIter = itsVolCurves.begin();

	for (; vIter != itsVolCurves.end(); vIter++)
	{
		if (itsVolCurves[vIter->first])
        {
           delete itsVolCurves[vIter->first];

	       itsVolCurves[vIter->first] = (ARM_VolCurve *) NULL;
        }
	}
}



void ARM_HyperCube::Init(void)
{
	SetName(ARM_HYPER_CUBE);

	SetIntersurfaceInterpol(false);
}



ARM_Object* ARM_HyperCube::Clone(void)
{
    ARM_HyperCube*	theClone = new ARM_HyperCube();


    theClone->Copy(this);

    return	theClone;
}



void ARM_HyperCube::Copy(const ARM_Object* aCubeIn)
{
    ARM_VolCube::Copy(aCubeIn);

    BitwiseCopy(aCubeIn);
}



void ARM_HyperCube::BitwiseCopy(const ARM_Object* srcObject)
{
    ARM_HyperCube* srcCube = (ARM_HyperCube *) (srcObject);


	std::map<string, ARM_VolCurve*>::iterator vIter = srcCube->itsVolCurves.begin();

	for (; vIter != srcCube->itsVolCurves.end(); vIter++)
	{
		itsVolCurves[vIter->first] = (ARM_VolCurve*)(vIter->second)->Clone();
	}

	vector<int>::iterator	vIter2 = srcCube->itsTenorList.begin();

	for(; vIter2 != srcCube->itsTenorList.end(); vIter2++)
	{
		itsTenorList.push_back(*vIter2);
	}

	itsIntersurfaceInterpol = srcCube->itsIntersurfaceInterpol;
}



void ARM_HyperCube::View(char* id, FILE* ficOut)
{
    FILE*	fOut;
    char	fOutName[200];


    if( ficOut == NULL )
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
    fprintf(fOut, "\n\n >>>>>>>>>>>>>>>>>>>>>> HYPER CUBE <<<<<<<<<<<<<<<<<<<<<<<\n\n");

    ARM_AbstractMarketClass::View(id, fOut);

    fprintf(fOut, "\n\n ======================> HyperCube COMPONENTS <========================\n\n");

	ARM_VolCurve*	vVolCurve = NULL;

	string C_CurveTenor;
	map<string, ARM_VolCurve*>::iterator	vVolCurveIter = itsVolCurves.begin();
	map<string, ARM_VolCurve*>::iterator	vVolCurveEnd = itsVolCurves.end();

	for (; vVolCurveIter != vVolCurveEnd; vVolCurveIter++)
	{
		C_CurveTenor = vVolCurveIter->first;
		
        fprintf(fOut, "\n\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
		fprintf(fOut, "\n Cube Tenor: %s \n", C_CurveTenor.c_str());

		vVolCurve = vVolCurveIter->second;
		
        if (vVolCurve)
		   vVolCurve->View(id, fOut);
	}

	fprintf(fOut, "\n\n ======================< End of HyperCube COMPONENTS >========================\n\n");

	fprintf(fOut, "\n\n <<<<<<<<<<<<<<<<<<<<<< HYPER CUBE >>>>>>>>>>>>>>>>>>>>>>>\n\n");

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}



double	ARM_HyperCube::ComputeVolatility(double aOptMat, double aTenor1, double aTenor2, double aStrike)
{
	double	vTenor1 = aTenor1;
	double	vTenor2 = aTenor2;


	if ( aTenor2 < aTenor1 )
	{
	   vTenor1 = aTenor2;
	   vTenor2 = aTenor1;
	}

	string	vTenor1Str(ConvertYearTermToStringMatu(vTenor1));
	ARM_VolCurve*	vVolCurve = NULL;
	map<string, ARM_VolCurve*>::iterator	vIter = itsVolCurves.find(vTenor1Str);

	if ( vIter != itsVolCurves.end() )
	{
		vVolCurve = vIter->second;

		if ( vVolCurve->GetName() != ARM_VOL_CUBE )
		   aStrike = 0;

		return	( vVolCurve->ComputeVolatility(aOptMat, aStrike, vTenor2) );
	}

	if ( !itsIntersurfaceInterpol )
	{
		char vErrorMsg[50];

		sprintf(vErrorMsg, "Couldn't find Correl cube for Tenor = %s ", vTenor1Str.c_str());

		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
                        "ARM_HyperCube::ComputeVolatility(...): Tenor not found");
	}

	// Encadrement de Tenor1 par des tenors de la liste
	int	vLowerBound, vUpperBound;
	GetBounds(itsTenorList, (int)(12*vTenor1), vLowerBound, vUpperBound);

	int	vTenor1_d = itsTenorList.at(vLowerBound);
	int	vTenor1_u = itsTenorList.at(vUpperBound);

	string	sTenor1_d(ConvertYearTermToStringMatu((double)vTenor1_d / 12.));
	string	sTenor1_u(ConvertYearTermToStringMatu((double)vTenor1_u / 12.));
	
	ARM_VolCurve*	vVolCurve_d = GetVolCurve(sTenor1_d);
	ARM_VolCurve*	vVolCurve_u = GetVolCurve(sTenor1_u);

	double	vCorrel_d = vVolCurve_d->ComputeVolatility(aOptMat, aStrike, vTenor2);
	double	vCorrel_u = vVolCurve_u->ComputeVolatility(aOptMat, aStrike, vTenor2);

	double	vCorrel = linInterpol(12.0*vTenor1, vTenor1_d, vCorrel_d, vTenor1_u, vCorrel_u);

	return	vCorrel;
}



double ARM_HyperCube::ComputeCorrelByExpiry(double aExpiry, double aTenor1, double aTenor2)
{
	double	vCorrel = 0;

	if ( GetVolType() == K_CUBE_CORREL_DIAG )
	   vCorrel = ComputeVolatility(aTenor1, aTenor2, aExpiry);
	else
	   vCorrel = ComputeVolatility(aExpiry, aTenor1, aTenor2);

	if ( vCorrel > 100 )
	   vCorrel = 100;
	else if ( vCorrel < -100 )
	   vCorrel = -100;

	return	vCorrel;
}



double ARM_HyperCube::ComputeHyperCorrel(double aOptMat, double aTenor1,
										 double aTenor2, double aStrike)
{
	return	( ComputeVolatility(aOptMat, aTenor1, aTenor2, aStrike) );
}



vector<int>	ARM_HyperCube::GetTenorList()
{
	return itsTenorList;
}



ARM_VolCurve* ARM_HyperCube::GetVolCurve(string& aTenor)
{
	ARM_VolCurve* vVolCurve = NULL;

	map<string, ARM_VolCurve*>::iterator vIter = itsVolCurves.find(aTenor);
	
    if ( vIter != itsVolCurves.end() )
	{
		vVolCurve = vIter->second;
	}

	return	vVolCurve;
}



ARM_VolCube* ARM_HyperCube::CreateCorrelCubeByExpiry(vector<double>& aTenor1List, 
													 vector<double>& aTenor2List, 
													 vector<double>& aExpiryList)
{
	ARM_Vector*	vTenorList;

	string vExpiryStr;
	string vTenorStr;
	double vExpiry;

    ARM_Currency* currency = GetCurrency();

    ARM_VolCurve* vVolCurve = NULL;
	ARM_VolCube*  vVolCube  = NULL;

	map<string, ARM_VolCurve*>::iterator vVolCurvesIter;


	// Par défaut, la liste des expiries en sortie est égale à celle de l'hyperCube initial
	if (aExpiryList.empty())
	{
		for (vVolCurvesIter = itsVolCurves.begin(); vVolCurvesIter != itsVolCurves.end(); vVolCurvesIter++ )
		{
			vExpiryStr = vVolCurvesIter->first;
			
            vExpiry = StringMatuToYearTerm((char*)vExpiryStr.c_str());
			
            aExpiryList.push_back(vExpiry);
		}
	}

	sort(aExpiryList.begin(), aExpiryList.end());


	// Par défaut, la liste des tenors en sortie est égale à l'union des tenors des volCurve de l'hyperCube initial
	if (aTenor1List.empty())
	{
		for (vVolCurvesIter = itsVolCurves.begin(); vVolCurvesIter != itsVolCurves.end(); vVolCurvesIter++ )
		{
			vVolCurve = vVolCurvesIter->second;
			
            if ( vVolCurve->GetName() == ARM_VOL_CUBE )
			{
				vVolCube = (ARM_VolCube*)vVolCurve;
				vVolCurve = vVolCube->GetATMVol();
			}

			vTenorList = ((ARM_VolLInterpol *) vVolCurve)->GetStrikes();
			
            for (int i = 0; i < vTenorList->size(); i++)
			{
				aTenor1List.push_back(vTenorList->Elt(i));
			}
		}
	}

	sort(aTenor1List.begin(), aTenor1List.end());
	vector<double>::iterator vFirstDuplicate = unique(aTenor1List.begin(), aTenor1List.end());
	aTenor1List.erase(vFirstDuplicate, aTenor1List.end());

	aTenor2List = aTenor1List;

	
	vector<double>::iterator	vExpiryIter;
	vector<double>::iterator	vTenor1Iter;
	vector<double>::iterator	vTenor2Iter;
	
	vector<ARM_VolCurve*>		vCorrelCurves;
	vector<string>	vExpiryList;

	ARM_Vector*	vTenor1Terms = NULL;
	ARM_Vector*	vTenor2Terms = NULL;
	ARM_Matrix*	vVolatilities = NULL;

	int	vTenor1ListSize = aTenor1List.size();
	int	vTenor2ListSize = aTenor2List.size();

	double	vTenor1, vTenor2, vCorrel;
	int		i, j;
	string	vIndexName;

	for (vExpiryIter = aExpiryList.begin() ; vExpiryIter != aExpiryList.end(); vExpiryIter++ )
	{
		vExpiry = *vExpiryIter;

		vTenor1Terms  = new ARM_Vector(aTenor1List);
		vTenor2Terms  = new ARM_Vector(aTenor2List);
		vVolatilities = new ARM_Matrix(vTenor1ListSize, vTenor2ListSize);

		for( vTenor1Iter = aTenor1List.begin(), i=0; vTenor1Iter != aTenor1List.end(); vTenor1Iter++, i++ )
		{
			vTenor1 = *vTenor1Iter;

			for( vTenor2Iter = aTenor2List.begin(), j=0; vTenor2Iter != aTenor2List.end(); vTenor2Iter++, j++ )
			{
				vTenor2 = *vTenor2Iter;
				
				if ( vTenor1 == vTenor2 )
				   vCorrel = 100.;
				else
				   vCorrel = ComputeVolatility(vExpiry, vTenor1, vTenor2); 

				vVolatilities->Elt(i, j) = vCorrel;
			}
		}

		vVolCurve = new ARM_VolLInterpol(itsAsOfDate, vTenor1Terms, vTenor2Terms, vVolatilities);
 
        vVolCurve->SetCurrencyUnit(currency);

		vVolCurve->SetMktExternalCharacteristics(GetExternalIndex(), GetStrCurrency(), GetExternalCrvId(), GetStrType());
		vVolCurve->SetInterpType(K_DIAG_INTERPOL);

		vExpiryStr = ConvertYearTermToStringMatu(vExpiry);
		vIndexName = "Expiry = " + vExpiryStr;
		vVolCurve->SetIndexName((char*)vIndexName.c_str());

		vCorrelCurves.push_back(vVolCurve);
		vExpiryList.push_back(vExpiryStr);
	}


 	ARM_Vector	vUnderlyings(aExpiryList);

	vVolCube = new ARM_VolCube(&vCorrelCurves, &vUnderlyings, GetLastKnownDate());

    vVolCube->SetCurrencyUnit(currency);

	vVolCube->SetVolType(K_CUBE_CORREL_DIAG);

	vector<ARM_VolCurve*>::iterator	vCorrelCurvesIter;
	for (vCorrelCurvesIter = vCorrelCurves.begin(); vCorrelCurvesIter != vCorrelCurves.end(); vCorrelCurvesIter++)
	{
		delete *vCorrelCurvesIter;
	}

	return	vVolCube;
}


void ARM_HyperCube::BumpVolatility(double value, int nthLine, int nthCol, int isCumul, int isAbsolute)
{
	ARM_VolCurve* vCurve;
	ARM_Matrix* vMatrix;
	double valCorrel;
	int nbLine, nbCol, i,j;
	map<string, ARM_VolCurve*>::iterator iter = itsVolCurves.begin();

	for(; iter != itsVolCurves.end(); iter++)
	{
		vCurve = iter->second;
		vCurve->BumpVolatility(value, nthLine, nthCol, isCumul, isAbsolute);
		
		if(vCurve->GetName() == ARM_VOL_CUBE )
			vMatrix = ((ARM_VolCube*)vCurve)->GetATMVol()->GetVolatilities();	
		else
			vMatrix = vCurve->GetVolatilities();

		nbLine = vMatrix->GetNumLines();
		nbCol = vMatrix->GetNumCols();
		for(i = 0; i < nbLine; i++)
		{
			for(j=0;j<nbCol;j++)
			{
				valCorrel = vMatrix->Elt(i,j);
				if(valCorrel > 100.)
					vMatrix->Elt(i,j) = 100;
				if(valCorrel < -100.)
					vMatrix->Elt(i,j) = -100;
			}
		}
	}
}


void ARM_HyperCube::BumpSmile(double value, double aCubeTenor, double aSmileTenor, 
							  int nthLine, int nthCol, int isCumul, int isAbsolute)
{
	double			vCubeTenorDefaultValue = 0.;
	ARM_VolCurve*	vCurve;
	ARM_VolCube*	vVolCube;
	string			sCubeTenor(ConvertYearTermToStringMatu(aCubeTenor));
	map<string, ARM_VolCurve*>::iterator iter;

	if(aCubeTenor != vCubeTenorDefaultValue)
	{
		iter = itsVolCurves.find(sCubeTenor);
		if( iter != itsVolCurves.end() )
		{
			vCurve = iter->second;
			if(vCurve->GetName() == ARM_VOL_CUBE)
			{
				vVolCube = (ARM_VolCube*)vCurve;
				vVolCube->BumpSmile(value, aSmileTenor, nthLine, nthCol, isCumul, isAbsolute);
			}
		}
		else
		{
			char vErrorMsg[50];
			sprintf(vErrorMsg, "Couldn't find VolCube for Tenor = %s ", sCubeTenor.c_str());
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
							"ARM_HyperCube::BumpSmile(...): Tenor not found");
		}
	}
	else
	{
		for(iter = itsVolCurves.begin(); iter != itsVolCurves.end(); iter++)
		{
			vCurve = iter->second;
			if(vCurve->GetName() == ARM_VOL_CUBE)
			{
				vVolCube = (ARM_VolCube*)vCurve;
				vVolCube->BumpSmile(value, aSmileTenor, nthLine, nthCol, isCumul, isAbsolute);
			}
		}
	}
}


void	ARM_HyperCube::GetBounds(const vector<int>& aVector, int aTenor, int& aLowerBound, int& aUpperBound) const
{
	int	vSize = aVector.size();

	if(aTenor <= aVector.at(0))
	{
		aLowerBound = 0;
		aUpperBound = 0;
	}
	else if(aTenor >= aVector.at(vSize-1))
	{
		aLowerBound = vSize-1;
		aUpperBound = vSize-1;
	}
	else
	{
		int	i=0;
		while(aVector.at(i) < aTenor)
		{
			i++;
		}
		aUpperBound = i;
		if(aVector.at(i) == aTenor)
		{
			aLowerBound = i;
		}
		else
		{
			aLowerBound = i-1;
		}
	}
}


/*---------------------------------------------------------------------------------------*/
/*---- End Of File ----*/