#include "IndexIndexCorrelCube.h"
#include "fromto.h"



ARM_IndexIndexCorrelCube::ARM_IndexIndexCorrelCube(void)
{
	Init();
}


ARM_IndexIndexCorrelCube::ARM_IndexIndexCorrelCube(vector<ARM_VolCurve*>& aCorrelList, 
												   vector<string>& aTenor1List, vector<string>& aTenor2List)
{
	Init();

	if ( aTenor1List.size() != aCorrelList.size() )
	{
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
						"The Correl list and the tenor1 list don't have the same size in index-index correl cube");
	}

	if ( aTenor2List.size() != aTenor1List.size() )
	{
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
						"The 2 tenor lists don't have the same size in index-index correl cube");
	}

	vector<ARM_VolCurve*>::iterator	vVolCurveIter = aCorrelList.begin();

	ARM_VolCurve* vVolCurve  = NULL;
	ARM_VolCurve* vVolCurve0 = *vVolCurveIter;
	
    if ( vVolCurve0 == NULL )
	{
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
					    "First correl curve is NULL in index-index correl cube");
	}

	int		vTenor1, vTenor2;
	double	vYearTerm1, vYearTerm2;
	string	vCorrelStr;
	vector<string>::iterator vTenor1Iter = aTenor1List.begin();
	vector<string>::iterator vTenor1End  = aTenor1List.end();
	vector<string>::iterator vTenor2Iter = aTenor2List.begin();
	vector<string>::iterator vTenor2End  = aTenor2List.end();

	ARM_Date vAsOfDate = vVolCurve0->GetAsOfDate();

	for (; vTenor1Iter != vTenor1End; vTenor1Iter++, vTenor2Iter++, vVolCurveIter++)
	{
		vVolCurve = *vVolCurveIter;

		if ( vVolCurve == NULL )
		{
		   throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
							"One of the correl curves is NULL in index-index correl cube");
		}

		if ( vVolCurve->GetAsOfDate() != vAsOfDate )
		{
		   throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,
							"Some correl curves don't have the same AsOfDate in index-index correl cube");
		}

		// On convertit en nb de mois
		vYearTerm1 = 12. * StringMatuToYearTerm((char*)(*vTenor1Iter).c_str());
		vYearTerm2 = 12. * StringMatuToYearTerm((char*)(*vTenor2Iter).c_str());
				
		if( vYearTerm1 > vYearTerm2 )
		{
			char msg[255];
			sprintf( msg, "Tenor1 (%s) should be smaller than Tenor2 (%s)", 
					(char*)(*vTenor1Iter).c_str(), (char*)(*vTenor2Iter).c_str() );
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
		}
		else
		{
			vTenor1 = (int)vYearTerm2;
			vTenor2 = (int)vYearTerm1;
			vCorrelStr = (*vTenor1Iter) + "-" + (*vTenor2Iter);
		}

        ARM_VolCurve* theClone = (ARM_VolCurve *) vVolCurve->Clone();
		theClone->SetInterpType(K_DIAG_INTERPOL);
	
		itsVolCurves.insert( pair<string, ARM_VolCurve*> (vCorrelStr, theClone) );

		itsTenorList.push_back(vTenor1);
		itsTenorList.push_back(vTenor2);
	}

	sort(itsTenorList.begin(), itsTenorList.end());
	vector<int>::iterator	vFirstDuplicate = unique(itsTenorList.begin(), itsTenorList.end());
	itsTenorList.erase(vFirstDuplicate, itsTenorList.end());
		

    SetAsOfDate(vAsOfDate);

	SetMktExternalCharacteristics(vVolCurve0->GetExternalIndex(), 
								  vVolCurve0->GetStrCurrency(),
								  vVolCurve0->GetExternalCrvId(),
								  string("Index-Index Correl Cube"));	// TMP : type à voir
}


ARM_IndexIndexCorrelCube::~ARM_IndexIndexCorrelCube(void)
{
}


void ARM_IndexIndexCorrelCube::Init(void)
{
	SetName(ARM_INDEX_INDEX_CORREL_CUBE);

	SetIntersurfaceInterpol(true);
}


ARM_Object* ARM_IndexIndexCorrelCube::Clone(void)
{
    ARM_IndexIndexCorrelCube*	theClone = new ARM_IndexIndexCorrelCube();


    theClone->Copy(this);

    return	theClone;
}


void ARM_IndexIndexCorrelCube::Copy(const ARM_Object* aCubeIn)
{
    ARM_HyperCube::Copy(aCubeIn);

    //BitwiseCopy(aCubeIn);
}


void ARM_IndexIndexCorrelCube::BitwiseCopy(const ARM_Object* srcObject)
{
}


void ARM_IndexIndexCorrelCube::View(char* id, FILE* ficOut)
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


	fprintf(fOut, "\n  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
    fprintf(fOut, "\n\n >>>>>>>>>>>>>>>>>>>>>> Index-Index Correl Cube <<<<<<<<<<<<<<<<<<<<<<<\n\n");

    ARM_AbstractMarketClass::View(id, fOut);

    fprintf(fOut, "\n\n ======================> IndexIndexCorrelCube COMPONENTS <========================\n\n");

	ARM_VolCurve*	vVolCurve = NULL;
	string	vCurveKey;

	map<string, ARM_VolCurve*>::iterator	vVolCurveIter = itsVolCurves.begin();
	map<string, ARM_VolCurve*>::iterator	vVolCurveEnd = itsVolCurves.end();

	for(; vVolCurveIter != vVolCurveEnd; vVolCurveIter++)
	{
		vCurveKey = vVolCurveIter->first;
		fprintf(fOut, "\n\n xxxxxxxxxxxx\n Correl %s", vCurveKey.c_str());

		vVolCurve = vVolCurveIter->second;
		if(vVolCurve)
			vVolCurve->View(id, fOut);
	}

	fprintf(fOut, "\n\n ======================< End of IndexIndexCorrelCube COMPONENTS >========================\n\n");

	fprintf(fOut, "\n\n <<<<<<<<<<<<<<<<<<<<<< Index-Index Correl Cube >>>>>>>>>>>>>>>>>>>>>>>\n\n");
	fprintf(fOut, "\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}


double	ARM_IndexIndexCorrelCube::ComputeIndexIndexCorrel(string aTenor1, double aExpiry1,
														  string aTenor2, double aExpiry2) const
{
	// On convertit en nb de mois
	double	vYearTerm1 = 12. * StringMatuToYearTerm((char*)aTenor1.c_str());
	double	vYearTerm2 = 12. * StringMatuToYearTerm((char*)aTenor2.c_str());
	
	int	vTenor1 = (int)vYearTerm1;
	int	vTenor2 = (int)vYearTerm2;

	string	vTenor1Str = aTenor1;
	string	vTenor2Str = aTenor2;

	double	vExpiry1 = aExpiry1;
	double	vExpiry2 = aExpiry2;

	OrderTenorsAndExpiries(vTenor1, vTenor2, vTenor1Str, vTenor2Str, vExpiry1, vExpiry2);
	string	vCorrelStr = vTenor2Str + "-" + vTenor1Str;

	map<string, ARM_VolCurve*>::const_iterator	vIter = itsVolCurves.find(vCorrelStr);
	if ( vIter != itsVolCurves.end() )
	{
		ARM_VolCurve*	vVolCurve = vIter->second;

		return	( vVolCurve->ComputeVolatility(vExpiry2, vExpiry1) );
	}

	if ( !itsIntersurfaceInterpol )
	{
		char vErrorMsg[50];

		sprintf(vErrorMsg, "Couldn't find Correl cube for tenor couple = %s ", vCorrelStr.c_str());

		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
                        "ARM_IndexIndexCorrelCube::ComputeIndexIndexCorrel(...): Tenor couple not found");
	}

	// "Tenor2-Tenor1" n'est pas dans la liste des correl
	// On n'intervertit éventuellement les tenors et les expiries qu'en bout de chaîne
	vExpiry1 = aExpiry1;
	vExpiry2 = aExpiry2;
	vTenor1 = (int)vYearTerm1;
	vTenor2 = (int)vYearTerm2;

	// Encadrement de Tenor1 et Tenor2 par des tenors de la liste
	int	vLowerBound1, vUpperBound1, vLowerBound2, vUpperBound2;
	GetBounds(itsTenorList, vTenor1, vLowerBound1, vUpperBound1);
	GetBounds(itsTenorList, vTenor2, vLowerBound2, vUpperBound2);

	int	vTenor1_d = itsTenorList.at(vLowerBound1);
	int	vTenor1_u = itsTenorList.at(vUpperBound1);
	int	vTenor2_d = itsTenorList.at(vLowerBound2);
	int	vTenor2_u = itsTenorList.at(vUpperBound2);

	string	sTenor1_d(ConvertYearTermToStringMatu((double)vTenor1_d / 12.));
	string	sTenor1_u(ConvertYearTermToStringMatu((double)vTenor1_u / 12.));
	string	sTenor2_d(ConvertYearTermToStringMatu((double)vTenor2_d / 12.));
	string	sTenor2_u(ConvertYearTermToStringMatu((double)vTenor2_u / 12.));

	double	vCorrel_dd, vCorrel_du, vCorrel_ud, vCorrel_uu;
	
	ARM_VolCurve*	vVolCurve_dd = NULL;
	ARM_VolCurve*	vVolCurve_du = NULL;
	ARM_VolCurve*	vVolCurve_ud = NULL;
	ARM_VolCurve*	vVolCurve_uu = NULL;
	
	if(vLowerBound1 == vUpperBound1)
	{
		if(vLowerBound2 == vUpperBound2)
		{
			// Uniquement en cas d'extrapolation si on est au-dela de la derniere maturite sur les 2 axes	
			OrderTenorsAndExpiries(vTenor1_d, vTenor2_d, sTenor1_d, sTenor2_d, vExpiry1, vExpiry2);
			vCorrelStr = sTenor2_d + "-" + sTenor1_d;
			vVolCurve_dd = GetCorrelCurve(vCorrelStr);
			vCorrel_dd = vCorrel_du = vCorrel_ud = vCorrel_uu = vVolCurve_dd->ComputeVolatility(vExpiry2, vExpiry1);
		}
		else
		{
			vCorrel_dd = vCorrel_ud = ComputeCorrel(vTenor1_d, vTenor2_d, sTenor1_d, sTenor2_d, vExpiry1, vExpiry2);
			vCorrel_du = vCorrel_uu = ComputeCorrel(vTenor1_d, vTenor2_u, sTenor1_d, sTenor2_u, vExpiry1, vExpiry2);
		}
	}
	else
	{
		if(vLowerBound2 == vUpperBound2)
		{
			vCorrel_dd = vCorrel_du = ComputeCorrel(vTenor1_d, vTenor2_d, sTenor1_d, sTenor2_d, vExpiry1, vExpiry2);
			vCorrel_ud = vCorrel_uu = ComputeCorrel(vTenor1_u, vTenor2_d, sTenor1_u, sTenor2_d, vExpiry1, vExpiry2);
		}
		else
		{
			vCorrel_dd = ComputeCorrel(vTenor1_d, vTenor2_d, sTenor1_d, sTenor2_d, vExpiry1, vExpiry2);
			vCorrel_ud = ComputeCorrel(vTenor1_u, vTenor2_d, sTenor1_u, sTenor2_d, vExpiry1, vExpiry2);
			vCorrel_du = ComputeCorrel(vTenor1_d, vTenor2_u, sTenor1_d, sTenor2_u, vExpiry1, vExpiry2);
			vCorrel_uu = ComputeCorrel(vTenor1_u, vTenor2_u, sTenor1_u, sTenor2_u, vExpiry1, vExpiry2);
		}
	}

	ARM_Vector	vTenors1(2);
	vTenors1.Elt(0) = vTenor1_d;
	vTenors1.Elt(1) = vTenor1_u;

	ARM_Vector	vTenors2(2);
	vTenors2.Elt(0) = vTenor2_d;
	vTenors2.Elt(1) = vTenor2_u;

	ARM_Matrix	vCorrelValues(2, 2);
	vCorrelValues.Elt(0, 0) = vCorrel_dd;
	vCorrelValues.Elt(0, 1) = vCorrel_du;
	vCorrelValues.Elt(1, 0) = vCorrel_ud;
	vCorrelValues.Elt(1, 1) = vCorrel_uu;

	double	vCorrel = triangularInterpol(&vTenors2, &vTenors1, &vCorrelValues, vTenor2, vTenor1);

	return	vCorrel;
}


double	ARM_IndexIndexCorrelCube::ComputeIndexIndexCorrel(double aTenor1, double aExpiry1,
														  double aTenor2, double aExpiry2) const
{
	string	vTenor1Str(ConvertYearTermToStringMatu(aTenor1));
	string	vTenor2Str(ConvertYearTermToStringMatu(aTenor2));

	return	( ComputeIndexIndexCorrel(vTenor1Str, aExpiry1, vTenor2Str, aExpiry2) );
}


ARM_VolCurve*	ARM_IndexIndexCorrelCube::GetCorrelCurve(string& aCorrelStr) const
{
	map<string, ARM_VolCurve*>::const_iterator	vIter = itsVolCurves.find(aCorrelStr);

	if( vIter == itsVolCurves.end() )
	{
		char vErrorMsg[50];

		sprintf(vErrorMsg, "Couldn't find %s Correl cube", aCorrelStr.c_str());

		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, vErrorMsg, 
						"ARM_IndexIndexCorrelCube::GetCorrelCurve(...): Correl not found");
	}

	return	(vIter->second);
}

bool ARM_IndexIndexCorrelCube::isCorrelCurve(string& aCorrelStr) const
{
	map<string, ARM_VolCurve*>::const_iterator	vIter = itsVolCurves.find(aCorrelStr);

	return ( vIter != itsVolCurves.end() );
}


void ARM_IndexIndexCorrelCube::OrderTenorsAndExpiries(int& aTenor1, int& aTenor2, 
													  string& aTenor1Str, string& aTenor2Str,
													  double& aExpiry1, double& aExpiry2) const
{
	if(aTenor1 < aTenor2)
	{
		int	vTenorTmp = aTenor2;
		aTenor2 = aTenor1;
		aTenor1 = vTenorTmp;

		double	vExpiryTmp = aExpiry2;
		aExpiry2 = aExpiry1;
		aExpiry1 = vExpiryTmp;

		string	vTenorStrTmp = aTenor2Str;
		aTenor2Str = aTenor1Str;
		aTenor1Str = vTenorStrTmp;
	}
}


double	ARM_IndexIndexCorrelCube::ComputeCorrel(int aTenor1, int aTenor2, 
												string aTenor1Str, string aTenor2Str,
												double aExpiry1, double aExpiry2) const
{
	OrderTenorsAndExpiries(aTenor1, aTenor2, aTenor1Str, aTenor2Str, aExpiry1, aExpiry2);

	string	vCorrelKey = aTenor2Str + "-" + aTenor1Str;

	ARM_VolCurve*	vCorrelCurve = GetCorrelCurve(vCorrelKey);

	double	vCorrel = vCorrelCurve->ComputeVolatility(aExpiry2, aExpiry1);

	return	vCorrel;
}

ARM_VolCurve* ARM_IndexIndexCorrelCube::GetCorrelDiag(string& aTenor1, vector<string>& aTenor2) const
{
	int tenor1 = 12.*StringMatuToYearTerm((char*)aTenor1.c_str() );
	int tenor2;
	int i,j;
	int aLowerBound;
	int aUpperBound;
	vector<double> listTenor;
	int sizeList = itsTenorList.size();
	int IdTenor1;

	GetBounds(itsTenorList,tenor1,aLowerBound,aUpperBound);

	char vErrorMsg[50];

		
	if( (aLowerBound != aUpperBound) || (aUpperBound == 0 && tenor1 < itsTenorList[0]) 
		|| (aLowerBound == sizeList-1 && tenor1 > itsTenorList[sizeList-1]) ) 
	{
		sprintf(vErrorMsg, "Couldn't find Tenor1 %s ", aTenor1.c_str());
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,vErrorMsg,
					    "ARM_IndexIndexCorrelCube::GetCorrelDiag(...)");
	}

	IdTenor1 = aLowerBound;

	string strTenor2;
	string strCorrel;
	
	if(!aTenor2.empty())
	{
		vector<string>::iterator iterStr = aTenor2.begin();
		for( ;iterStr != aTenor2.end();iterStr++)
		{
			strTenor2 = (*iterStr);
			tenor2 = 12.*StringMatuToYearTerm( (char*)(iterStr->c_str()) );

			GetBounds(itsTenorList,tenor2,aLowerBound,aUpperBound);

			if( (aLowerBound != aUpperBound) || (aUpperBound == 0 && tenor2 < itsTenorList[0]) 
				|| (aLowerBound == sizeList-1 && tenor2 > itsTenorList[sizeList-1]) ) 
			{
				sprintf(vErrorMsg, "Couldn't find Tenor2 %s ", iterStr->c_str());
				throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,vErrorMsg,
								"ARM_IndexIndexCorrelCube::GetCorrelDiag(...)");
			}
			if(IdTenor1 > aLowerBound)
			{
				sprintf(vErrorMsg, "Tenor1: %s should be lower than Tenor2:%s", aTenor1.c_str(),iterStr->c_str());
				throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,vErrorMsg,
								"ARM_IndexIndexCorrelCube::GetCorrelDiag(...)");
			}
			strCorrel = aTenor1+"-"+strTenor2;
			if(!isCorrelCurve(strCorrel))
			{
				sprintf(vErrorMsg, "there is no correl Curve between %s and %s", aTenor1.c_str(),strTenor2.c_str());
				throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,vErrorMsg,
								"ARM_IndexIndexCorrelCube::GetCorrelDiag(...)");
			}
			listTenor.push_back(tenor2 / 12.);
		}
	}
	else
	{
		for(i = IdTenor1;i<sizeList; i++)
		{
			strTenor2 = ConvertYearTermToStringMatu(itsTenorList[i] / 12.);
			strCorrel = aTenor1+"-"+strTenor2;
			if(isCorrelCurve(strCorrel))
				listTenor.push_back(itsTenorList[i] / 12.);
		}
	}

	
	vector<double>::iterator iterInt = listTenor.begin();
	j=0;

	strTenor2 = ConvertYearTermToStringMatu((double)(*iterInt));
	strCorrel = aTenor1+"-"+strTenor2;

	ARM_VolLInterpol* CorrelCurve = (ARM_VolLInterpol*)GetCorrelCurve(strCorrel);
	ARM_Vector* ExpiryTermes1 = CorrelCurve->GetExpiryTerms();
	ARM_Vector* ExpiryDates1 = CorrelCurve->GetStrikes();
	ARM_Matrix* Volatilities = CorrelCurve->GetVolatilities();
	int sizeLine = ExpiryTermes1->GetNumLines();
	int sizeCol = ExpiryDates1->GetNumLines();

	if( sizeLine != sizeCol )
	{
		sprintf(vErrorMsg, "Curve Vol %s is not a symetric matrix", strCorrel.c_str());
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,vErrorMsg,
						"ARM_IndexIndexCorrelCube::GetCorrelDiag(...)");
	}

	ARM_Matrix* MatrixVol = new ARM_Matrix(sizeLine,listTenor.size(),0.);
	for( i = 0 ; i<sizeLine ; i++ )
	{
		if( fabs(ExpiryTermes1->Elt(i) - ExpiryDates1->Elt(i)) > K_NEW_DOUBLE_TOL)
		{
			delete MatrixVol;
			sprintf(vErrorMsg, "Curve Vol %s is not a symetric matrix", strCorrel.c_str());
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,vErrorMsg,
							"ARM_IndexIndexCorrelCube::GetCorrelDiag(...)");
		}

		MatrixVol->Elt(i,0) = Volatilities->Elt(i,i);
	}

	ARM_Vector* ExpiryTermes2;
	ARM_Vector* ExpiryDates2;
	iterInt++;
	j++;
	for( ; iterInt!=listTenor.end() ; iterInt++,j++ )
	{
		strTenor2 = ConvertYearTermToStringMatu((double)(*iterInt));
		strCorrel = aTenor1+"-"+strTenor2;
		CorrelCurve = (ARM_VolLInterpol*)GetCorrelCurve(strCorrel);
		ExpiryTermes2 = CorrelCurve->GetExpiryTerms();
		ExpiryDates2 = CorrelCurve->GetStrikes();
		Volatilities = CorrelCurve->GetVolatilities();
		
		if( ( sizeLine != ExpiryTermes2->GetNumLines() ) || ( sizeLine != ExpiryDates2->GetNumLines() ) )
		{
			delete MatrixVol;
			sprintf(vErrorMsg, "Curve Vol %s is not a symetric matrix", strCorrel.c_str());
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,vErrorMsg,
							"ARM_IndexIndexCorrelCube::GetCorrelDiag(...)");
		}

		for( i=0 ; i<sizeLine ; i++ )
		{
			if( fabs(ExpiryTermes2->Elt(i) - ExpiryDates2->Elt(i)) > K_NEW_DOUBLE_TOL)
			{
				delete MatrixVol;
				sprintf(vErrorMsg, "Curve Vol %s is not a symetric matrix", strCorrel.c_str());
				throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,vErrorMsg,
								"ARM_IndexIndexCorrelCube::GetCorrelDiag(...)");
			}
			if( fabs(ExpiryTermes1->Elt(i) - ExpiryTermes2->Elt(i)) > K_NEW_DOUBLE_TOL)
			{
				delete MatrixVol;
				sprintf(vErrorMsg, "not the same Expiry at %s with the first curve",strCorrel.c_str());
				throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,vErrorMsg,
								"ARM_IndexIndexCorrelCube::GetCorrelDiag(...)");
			}

			MatrixVol->Elt(i,j) = Volatilities->Elt(i,i);
		}
	}

	ARM_Vector*	vTenor2List = new ARM_Vector(listTenor);

	ARM_VolCurve* vVolcurve = new ARM_VolLInterpol(GetAsOfDate(),(ARM_Vector*)ExpiryTermes1->Clone(),
												   vTenor2List,MatrixVol);
	
	vVolcurve->SetCurrencyUnit((ARM_Currency*)GetCurrency()->Clone());
//	vVolcurve->SetVolatilities(MatrixVol);
//	vVolcurve->SetExpiryTerms((ARM_Vector*)ExpiryTermes1->Clone());
//	vVolcurve->SetExpiryDates(&ARM_Vector(listTenor));

	return vVolcurve;
}

void ARM_IndexIndexCorrelCube::BumpVolatility(double value, int nthLine, int nthCol, int isCumul, int isAbsolute)
{

	ARM_Matrix* vMatrix;
	vector<string> strVector;
	string strCorrel;
	int nbLine, i;
	char vErrorMsg[50];
	map<string, ARM_VolCurve*>::iterator iter = itsVolCurves.begin();

	ARM_HyperCube::BumpVolatility(value, nthLine, nthCol, isCumul, isAbsolute);

	for(; iter != itsVolCurves.end(); iter++)
	{
		strCorrel = iter->first;
		strVector = splitStringIntoPieces(strCorrel,"-");
		if( strVector[0] == strVector[1] )
		{
			if(iter->second->GetName() == ARM_VOL_CUBE )
				vMatrix = ((ARM_VolCube*)iter->second)->GetATMVol()->GetVolatilities();
			else
				vMatrix = iter->second->GetVolatilities();

			nbLine = vMatrix->GetNumLines();
			if( nbLine != vMatrix->GetNumCols() )
			{
				sprintf(vErrorMsg, "Curve Vol %s is not a symetric matrix", strCorrel.c_str());
					throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA ,vErrorMsg,
									"ARM_IndexIndexCorrelCube::BumpVolatility(...)");
			}
			for(i = 0; i < nbLine;i++)
				vMatrix->Elt(i,i) = 100;
		}
	}
}

