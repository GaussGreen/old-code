

#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ARMKernel\crv\correlmanager.h"
#include "ICMKernel\glob\icm_mktdatamng.h"
#include "ICMKernel\glob\icm_correlation.h"
#include <set>


void ICM_ModelMultiCurves::Init(void)
{
   SetName(ICM_MODELMULTICURVES);

	// itsNbDefCurves = 0;
	// itsDefaultCurves = NULL;
	itsRecoveryRates.Resize(0); 
	itsMktDataMng  = new ICM_MktDataMng();
 	itsCorrelation = NULL;
	itsCorrelManager = NULL;
       
	itsIndependantPart = 0.;
	itsFullCorrelPart = 0.;
}
// virtual 
void ICM_ModelMultiCurves::SetCorrelation(ICM_Correlation* correl) 
{
	if (GetFlgClone())
	{
		if (itsCorrelation)
			delete itsCorrelation; itsCorrelation=NULL; 
		if (correl) 
			itsCorrelation = (ICM_Correlation*) correl->Clone();
	}
	else
	{
		itsCorrelation =NULL; 
		if (correl) 
			itsCorrelation = correl;
	}
}
// void ICM_ModelMultiCurves::Copy(const ARM_Object* src)
// { ICM_DefaultCurveModel::Copy(src);
//   BitwiseCopy(src); }

ARM_Object* ICM_ModelMultiCurves::Clone(void)
{
//   ICM_ModelMultiCurves* theClone = new ICM_ModelMultiCurves();
//   theClone->Copy(this);
//   return(theClone);
	return new ICM_ModelMultiCurves(*this); 
}


// set pointer without copy. Destruct if pointer is null;
void ICM_ModelMultiCurves::SetCpnInfCurve(ARM::ARM_InfCurv* infcurve) 
{ 
	string lKeyNewInfla(INFCURV);
	if (infcurve) {	
		if (GetFlgClone())
		{
			itsMktDataMng->adopt((ARM_Object*) infcurve->Clone(), lKeyNewInfla);
		}
		else
			itsMktDataMng->associate((ARM_Object*) infcurve, lKeyNewInfla);
	}else {
		itsMktDataMng->erase(lKeyNewInfla);	
	}
}

ARM::ARM_InfCurv* ICM_ModelMultiCurves::GetCpnInfCurve(void) const {
	//return itsCpnInfCurve;
	return itsMktDataMng->GetUnicInflaCurve();
}


void ICM_ModelMultiCurves::SetCpnIRCurve(ARM_ZeroCurve* ircurve) 
{	
	string lKeyNewZC(IRCURVE);
	if (ircurve) {
		if (GetFlgClone())
		{
			itsMktDataMng->adopt((ARM_Object*) ircurve->Clone(), lKeyNewZC);
		}
		else
			itsMktDataMng->associate((ARM_Object*) ircurve, lKeyNewZC);
	}else {
		itsMktDataMng->erase(lKeyNewZC);
	}
}

 ARM_ZeroCurve* ICM_ModelMultiCurves::GetCpnIRCurve(void)  const {
	return itsMktDataMng->GetUnicZeroCurve();
}

// ------------------------------------------------------------------------------------------------------
//	Returns Union of tenors for a given vector of issuers labels ex:"1Y,3Y,5Y..."
// ------------------------------------------------------------------------------------------------------

 /** 
char** 
ICM_ModelMultiCurves::GetUnionTenorsForCollateral(const std::vector<std::string>& labels,
												  int& sizeout, 
												  bool ParallelFlag)
{
	int size =0;
	int maxsize = 0;
	int i=0,j=0,k=0;
	bool status = false;
	char* tmp = NULL;

	std::set<ORDER_CHAR> Tenors;

	int sizel = labels.size(); 
	for (k=0;k<sizel;k++) 
	{
	status = false;

	for (i=0;i<GetNbDefCurves();i++)
	{
		if (strcmp(GetDefaultCurves(i)->GetLabel().c_str(),labels[k].c_str())==NULL) 	{status = true;}

		if (status)
		{
		for (j=0;j<unconst(*GetDefaultCurves(i)).GetRates()->GetSize();j++)
		{
			if (GetDefaultCurves(i)->GetTerm(j) != string("X"))
			{
			tmp = new char[ARM_NB_MAX_CHAR_TERMS];
			memset(tmp,'\0',sizeof(char)*ARM_NB_MAX_CHAR_TERMS);
			char tmp2[ARM_NB_MAX_CHAR_TERMS];
			memset(tmp2,'\0',sizeof(tmp2));
			memcpy(tmp2,GetDefaultCurves(i)->GetTerm(j).c_str(),sizeof(char)*ARM_NB_MAX_CHAR_TERMS);
			strcpy(tmp,tmp2);
			ORDER_CHAR y(tmp,GetDefaultCurves(i)->GetYearTerms()->Elt(j+1));
			
			if (Tenors.find(y) == Tenors.end()) Tenors.insert(y);
			else delete[] tmp;
			}
		}
		break;
		}
	}

	}

	size = Tenors.size();
	if (ParallelFlag)
		size++;
	sizeout = size;
		
	char** output = new char*[size];
	i=0;

	for(set<ORDER_CHAR>::iterator iterateur = Tenors.begin(); iterateur != Tenors.end(); ++ iterateur)
	{
		output[i] = iterateur->id;
		i++;
	}

	if (ParallelFlag) output[i]	=	"NONE";

	Tenors.clear();

	return (output);
}

**/ 

// ------------------------------------------------------------------------------------------------------
//	Returns Union of Year Fractions for a given vector of issuers labels 
// ------------------------------------------------------------------------------------------------------
void ICM_ModelMultiCurves::GetUnionYFForCollateral(const std::vector<std::string>& labels,ARM_Vector& Output) const
{

	int i=0,j=0,k=0;
	bool status = false;

	set<ORDER_DOUBLE_WITH_ORDER> Tenors;

	int size = labels.size(); 
	for (k=0;k<size;k++) 
	{
		for (i=0;i<GetNbDefCurves();i++)
		{
			status = false;
			if (strcmp(GetDefaultCurves(i)->GetLabel().c_str(),labels[k].c_str())==NULL) {status = true;}

			if (status)
			{
				for (j=1;j<unconst(*GetDefaultCurves(i)).GetRates()->GetSize();j++)
				{ORDER_DOUBLE_WITH_ORDER y(GetDefaultCurves(i)->GetYearTerms()->Elt(j));
				if (Tenors.find(y) == Tenors.end()) Tenors.insert(y);}
				break;
			}
		}
	}

	size = Tenors.size();
	Output.Resize(size);
	i=0;

	for(set<ORDER_DOUBLE_WITH_ORDER>::iterator iterateur = Tenors.begin(); iterateur != Tenors.end(); ++ iterateur)
	{
		Output[i] = iterateur->m_value;
		i++;
	}

	Tenors.clear();
}
//	Take ownership (or not) of the container
//	Might not take ownership of the curves
void   ICM_ModelMultiCurves::SetDefaultCurves(const std::vector<const ICM_DefaultCurve*>&curves) 
{ 
	if (GetFlgClone())
	{
		for (int i=0;i<itsDefaultCurves.size();i++)
			if (itsDefaultCurves[i]) delete  (itsDefaultCurves[i]); 
	}
	itsDefaultCurves = curves;
}
// void   ICM_ModelMultiCurves::SetDefaultCurves(ICM_DefaultCurve** curves) 
// { 
// 	SetDefaultCurves((const ICM_DefaultCurve** )curves); 
// }

// LJ -> CC
void ICM_ModelMultiCurves::SetDefaultCurve(int Num, const ICM_DefaultCurve* data)
{
	if ((Num >= itsDefaultCurves.size()) || (Num < 0))
		ICMTHROW(ERR_INVALID_DATA,"Invalid Default Curve Id: " << Num << ".")
	else{
		if (GetFlgClone())
		{
			if(itsDefaultCurves[Num]) delete itsDefaultCurves[Num];
		}
		itsDefaultCurves[Num] = data; 
	}
}
/** JLA : moved to copy ctor
void ICM_ModelMultiCurves::BitwiseCopy(const ARM_Object* src)
{
	int i = 0;
	ICM_ModelMultiCurves* srcMod = (ICM_ModelMultiCurves*) src;
	itsNbDefCurves = srcMod->itsNbDefCurves;
	if (GetFlgClone())
	{
		if (itsDefaultCurves)
		{
			for (i=0;i<itsNbDefCurves;i++)
			{
			 if (itsDefaultCurves[i])
				 delete  (itsDefaultCurves[i]);
			}
			delete[] itsDefaultCurves;
			itsDefaultCurves = NULL;
		}
		itsDefaultCurves = new  const ICM_DefaultCurve*[itsNbDefCurves];
		for (i=0;i<itsNbDefCurves;i++)
			itsDefaultCurves[i] = (ICM_DefaultCurve*) (((srcMod->itsDefaultCurves)[i])->Clone());
	}
	else
		itsDefaultCurves = srcMod->itsDefaultCurves;

	SetCorrelation((ICM_Correlation*) srcMod->GetCorrelation());
	itsRecoveryRates=srcMod->itsRecoveryRates; 
	
	if (itsMktDataMng) {
		delete itsMktDataMng; itsMktDataMng = NULL;
	}
	itsMktDataMng = (ICM_MktDataMng*) srcMod->GetMktDataMng()->Clone();


	SetCpnIRCurve((ARM_ZeroCurve*) srcMod->GetCpnIRCurve());
	SetCpnInfCurve((ARM::ARM_InfCurv*) srcMod->GetCpnInfCurve());

	if (GetFlgClone())
	{if (srcMod->itsCorrelManager)
	{itsCorrelManager = (ARM_CorrelManager*) srcMod->itsCorrelManager->Clone();}}
	else
		itsCorrelManager = srcMod->itsCorrelManager; // ??

}
**/ 

/** JLA TEMP 
void ICM_ModelMultiCurves::Set(int NbDefCurves, 
							   ICM_DefaultCurve** DefaultCurves,
							   ARM_ZeroCurve* DiscountCurve,
							   const ARM_Vector& RecoveryRates,
							   // double** CorrelationMatrix,
							   const ICM_QMatrix<double>& corrMatrix,
							   ARM::ARM_InfCurv* infcurve,
							   ARM_ZeroCurve* CpnIRCurve,
							   ARM_CorrelManager* CorrelManager)
{
	// AsOfDate of model is AsOfDate of discountCurve
	ARM_Date AsOf = DiscountCurve->GetAsOfDate();
	SetStartDate(AsOf);

	int i =0;
	itsNbDefCurves = NbDefCurves;

	//Check Symetric Matrix
	// if (CorrelationMatrix != NULL)
	if (!corrMatrix.IsEmpty() )
		for (i = 0; i< NbDefCurves; i++)
		{
			for (int j = 0; j< NbDefCurves; j++)
			{
				if (corrMatrix(i,j) != corrMatrix(j,i) ) 
				{
					throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
						"ERROR : The Correlation Matrix is not symetric !");
				}

				if (fabs(corrMatrix(i,j) )>1.)
				{
					throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
						"ERROR : a value in correlation matrix is >1 in absolute value");
				}

			}
		}

	if (GetFlgClone())
	{
	if (itsDefaultCurves)
	{
		for (i=0;i<itsNbDefCurves;i++)
		{
		 if (itsDefaultCurves[i])
			 delete  (itsDefaultCurves[i]);
		}
		delete[] itsDefaultCurves;
		itsDefaultCurves = NULL;
	}

	itsDefaultCurves = new ICM_DefaultCurve*[itsNbDefCurves];

	for (i=0;i<itsNbDefCurves;i++)
		itsDefaultCurves[i] = (ICM_DefaultCurve*) (((DefaultCurves)[i])->Clone());
	}
	else
		itsDefaultCurves = DefaultCurves;

	// if (itsRecoveryRates)
	// 	delete[] itsRecoveryRates;
	// itsRecoveryRates = NULL;

	// if (RecoveryRates)
	// {	
	// 	itsRecoveryRates = new double[itsNbDefCurves];
	// 	memcpy(itsRecoveryRates, RecoveryRates, itsNbDefCurves*sizeof(double));
	// }
	itsRecoveryRates=RecoveryRates; 

	std::vector<std::string> labels ( NbDefCurves ); 
	for (i=0;i<itsNbDefCurves;i++)
		labels[i] = DefaultCurves[i]->GetLabel() ; 
	

	if (itsMktDataMng)
		delete itsMktDataMng;
	itsMktDataMng = new ICM_MktDataMng();

	if (!corrMatrix.IsEmpty() ){
		ICM_CorrMatrix* lpCorrelation = new ICM_CorrMatrix(AsOf,"DBLCORRMATRIX",labels,corrMatrix);
		SetCorrelation((ICM_Correlation*) lpCorrelation);
	}

	
	SetCpnIRCurve((ARM_ZeroCurve*) CpnIRCurve);
	SetCpnInfCurve((ARM::ARM_InfCurv*) infcurve);

	if (GetFlgClone())
	{
	if (itsCorrelManager)
		delete itsCorrelManager;
	itsCorrelManager = (ARM_CorrelManager*) CorrelManager->Clone();
	}
	else
		itsCorrelManager = CorrelManager;
}
**/ 
void ICM_ModelMultiCurves::Set( /** int NbDefCurves, 
							   const ICM_DefaultCurve** DefaultCurves, **/ 
							   const std::vector<const ICM_DefaultCurve*>& defCurves,
							   ARM_ZeroCurve* DiscountCurve,
							   const ARM_Vector& RecoveryRates,
							   ICM_Correlation* Correlation,
							   ARM::ARM_InfCurv* infcurve,
							   ARM_ZeroCurve* CpnIRCurve,
							   ARM_CorrelManager* CorrelManager)
{
	// AsOfDate of model is AsOfDate of discountCurve
	ARM_Date AsOf = DiscountCurve->GetAsOfDate();
	SetStartDate(AsOf);

	int i =0;
	// itsNbDefCurves = NbDefCurves;

	SetCorrelation((ICM_Correlation*) Correlation);

	
	if (GetFlgClone())
	{
		// if (itsDefaultCurves)
		// {
			for (i=0;i<itsDefaultCurves.size();i++)
			{
			 if (itsDefaultCurves[i])
				 delete  (itsDefaultCurves[i]);
			}
			// delete[] itsDefaultCurves;
			// itsDefaultCurves = NULL;
		// }
		// itsDefaultCurves = new const ICM_DefaultCurve*[itsNbDefCurves];
		itsDefaultCurves.resize(defCurves.size()); 
		for (i=0;i<itsDefaultCurves.size();i++)
			itsDefaultCurves[i] = (ICM_DefaultCurve*)  defCurves[i]->Clone();
	}
	else
		itsDefaultCurves = defCurves;
	
	
	if (itsMktDataMng)
		delete itsMktDataMng;
	itsMktDataMng = new ICM_MktDataMng();

	
	SetCpnIRCurve( CpnIRCurve);
	SetCpnInfCurve( infcurve);

	/**if (itsRecoveryRates)
		delete[] itsRecoveryRates;
	itsRecoveryRates = NULL;

	if (RecoveryRates)
	{
		itsRecoveryRates = new double[itsNbDefCurves];
		memcpy(itsRecoveryRates, RecoveryRates, itsNbDefCurves*sizeof(double));
	}**/
	itsRecoveryRates=RecoveryRates; 
	

	if (GetFlgClone())
	{
	if (itsCorrelManager)
		delete itsCorrelManager;
	if (CorrelManager) itsCorrelManager = (ARM_CorrelManager*) CorrelManager->Clone();
	}
	else
		itsCorrelManager = CorrelManager;

}

void ICM_ModelMultiCurves::Set(/** int NbDefCurves,
				const ICM_DefaultCurve** DefaultCurves, **/ 
				const std::vector<const ICM_DefaultCurve*>& defCurves,
				ARM_ZeroCurve* DiscountCurve,
				const ARM_Vector& RecoveryRates,
				ICM_Correlation*   lpCorrelation,
				ICM_MktDataMng * pMktDataMng,
				ARM_CorrelManager* CorrelManager)
{
	// AsOfDate of model is AsOfDate of discountCurve
	ARM_Date AsOf = DiscountCurve->GetAsOfDate();
	SetStartDate(AsOf);

	int i =0;
	// itsNbDefCurves = NbDefCurves;


	if (GetFlgClone())
	{
		// if (itsDefaultCurves)
		// {
			for (i=0;i<itsDefaultCurves.size();i++)
			{
			 if (itsDefaultCurves[i])
				 delete  (itsDefaultCurves[i]);
			}
		// 	delete[] itsDefaultCurves;
		// 	itsDefaultCurves = NULL;
		// }
		// itsDefaultCurves = new const ICM_DefaultCurve*[itsNbDefCurves];
		itsDefaultCurves.resize(defCurves.size()); 
		for (i=0;i<itsDefaultCurves.size();i++)
			itsDefaultCurves[i] = (ICM_DefaultCurve*) defCurves[i]->Clone();
	}
	else
		itsDefaultCurves = defCurves;
	

	SetCorrelation(lpCorrelation);

	if (itsMktDataMng)
		delete itsMktDataMng;
	itsMktDataMng = NULL;
	if (pMktDataMng){
		if (GetFlgClone())
			itsMktDataMng = dynamic_cast<ICM_MktDataMng*>(pMktDataMng->Clone());
		else 
			itsMktDataMng = new ICM_MktDataMng();		
		SetCpnIRCurve((ARM_ZeroCurve*) pMktDataMng->GetUnicZeroCurve());
		SetCpnInfCurve((ARM::ARM_InfCurv*) pMktDataMng->GetUnicInflaCurve());
	}else {
		itsMktDataMng = new ICM_MktDataMng();
	}

	/** if (itsRecoveryRates)
		delete[] itsRecoveryRates;
	itsRecoveryRates = NULL;

	if (RecoveryRates)
	{
		itsRecoveryRates = new double[itsNbDefCurves];
		memcpy(itsRecoveryRates, RecoveryRates, itsNbDefCurves*sizeof(double));
	} **/ 
	itsRecoveryRates=RecoveryRates; 
	//useless
	if (GetFlgClone())
	{
	if (itsCorrelManager)
		delete itsCorrelManager;
	if (CorrelManager) itsCorrelManager = (ARM_CorrelManager*) CorrelManager->Clone();
	}
	else
		itsCorrelManager = CorrelManager;
}
/*----------------------------------------------------------------------------*
    Constructor    (copy).
*----------------------------------------------------------------------------*/

ICM_ModelMultiCurves::ICM_ModelMultiCurves(const ICM_ModelMultiCurves& ref) : ICM_DefaultCurveModel(ref)
{
    Init();
    // BitwiseCopy(&ycModelCredit);
	int i = 0;
	// itsNbDefCurves = ref.itsNbDefCurves;
	// itsDefaultCurves = new  const ICM_DefaultCurve*[itsNbDefCurves];
	itsDefaultCurves.resize(ref.itsDefaultCurves.size()) ; 
	for (i=0;i<itsDefaultCurves.size();i++)
		itsDefaultCurves[i] = (ICM_DefaultCurve*) ref.itsDefaultCurves[i]->Clone();
	
	SetCorrelation((ICM_Correlation*) ref.GetCorrelation());
	itsRecoveryRates=ref.itsRecoveryRates; 
	
	itsMktDataMng = (ICM_MktDataMng*) ref.GetMktDataMng()->Clone();

	SetCpnIRCurve( ref.GetCpnIRCurve());
	SetCpnInfCurve(ref.GetCpnInfCurve());

	if (ref.itsCorrelManager)
		itsCorrelManager = (ARM_CorrelManager*) ref.itsCorrelManager->Clone();  
}


/*----------------------------------------------------------------------------*
    Destructor
*-----------------------------------------------------------------------------*/

ICM_ModelMultiCurves::~ICM_ModelMultiCurves(void)
{
	int i =0;	
	if (GetFlgClone())
	{
		// if (itsDefaultCurves)
		// {
			for (i=0;i<itsDefaultCurves.size();i++)	{ if (itsDefaultCurves[i]) delete  (itsDefaultCurves[i]);}
			// delete[] itsDefaultCurves;
		// 	itsDefaultCurves = NULL;
		// }
	} 
	// else 
	// {
	// 	//	we do not delete the underlying objs, but we delete the container. 
	// 	if (itsDefaultCurves) delete [] itsDefaultCurves; itsDefaultCurves=0; 
	// }


	if (GetFlgClone())
	{
		
		if (itsMktDataMng) {
			itsMktDataMng->clear();
			delete itsMktDataMng ;
			itsMktDataMng =NULL;
		}
		itsMktDataMng = NULL;
		if (itsCorrelManager)
			delete itsCorrelManager;
		itsCorrelManager = NULL;
	}
	if (itsMktDataMng) { delete itsMktDataMng; itsMktDataMng =NULL;} 

	if (GetFlgClone())
	{
	if (itsCorrelation)
		delete itsCorrelation;
	itsCorrelation = NULL;
	}
}


void ICM_ModelMultiCurves::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];

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

	char* Date = GetZeroCurve()->GetAsOfDate().GetStrDate();

    fprintf(fOut, "\t\t\t ----------------- Model Multi Curves ----------------- \n");
	fprintf(fOut,"\nInterest Rate Curve\n");
	fprintf(fOut,"Date Of The Discount Curve : %s\n", Date);

	if (GetZeroCurve()) GetZeroCurve()->View(id, fOut);

	

	if (itsCorrelation ) 
	{
		fprintf(fOut,"\n\t\t\t ----------------- Correlation Matrix ----------------- \n\n");
		itsCorrelation->View("",fOut);
	} else 
	{
		fprintf(fOut,"\n\t\t\t ----------------- No Correlation Defined ----------------- \n\n");
	}

	delete[] Date;

	if (itsMktDataMng) itsMktDataMng->View(id,fOut);

	fprintf(fOut,"\n");
	fprintf(fOut,"IndependantPart : %f\n",itsIndependantPart);
	fprintf(fOut,"FullCorrelPart : %f\n",itsFullCorrelPart);

	fprintf(fOut,"\tLabel\t\t\tRecovery\n");
	fprintf(fOut,"	-----------------------------------------\n");
	fprintf(fOut,"\n	VolCurve -----------------------------------------\n");
    if (GetVolCurve()) GetVolCurve()->View(id, fOut);
	fprintf(fOut,"\nn\n");
    
	for (int i =0; i<itsDefaultCurves.size(); i++)
	{
		fprintf(fOut,"\t%s\t\t\t",itsDefaultCurves[i]->GetLabel().c_str());
		if (!itsRecoveryRates.empty()) fprintf(fOut,"\t\t%.5lf",itsRecoveryRates[i]);
		fprintf(fOut,"\n");
	}

	fprintf(fOut,"\n \t\t\t ----------------- Default Curves ----------------- \n\n");
	fprintf(fOut,"Number of default curves : %ld\n\n", (long)itsDefaultCurves.size());
	for ( i =0; i<itsDefaultCurves.size(); i++)
	{		
		if (itsDefaultCurves[i]){
			fprintf(fOut,"\t%s\t\t\t",itsDefaultCurves[i]->GetLabel().c_str());
			itsDefaultCurves[i]->View(id, fOut);
			fprintf(fOut,"\n");
		}
	}

	if (itsCorrelManager)
	{
		fprintf(fOut, "\n\t\t\t ----------------- Correlation Manager ----------------- \n");
		itsCorrelManager->View(id, fOut);
	}		

	fprintf(fOut,"\n");

    if ( ficOut == NULL )
    {fclose(fOut);}

}

void ICM_ModelMultiCurves::SetMktDataMng(ICM_MktDataMng* pMktDataMng) { 
			if (itsMktDataMng) delete itsMktDataMng;
			if (pMktDataMng){
				if (GetFlgClone())
					itsMktDataMng = dynamic_cast<ICM_MktDataMng*>(pMktDataMng->Clone());
				else {
					itsMktDataMng = pMktDataMng;
					SetCpnIRCurve((ARM_ZeroCurve*) pMktDataMng->GetUnicZeroCurve());
					SetCpnInfCurve((ARM::ARM_InfCurv*) pMktDataMng->GetUnicInflaCurve());
				}
			}else {
				itsMktDataMng = NULL;
			}
}
 
 

// -------------------------------------------------------------------------
// Bump a package of Default Curves according a given bump profile
// -------------------------------------------------------------------------
ICM_ModelMultiCurves* 
ICM_ModelMultiCurves::GenerateShiftModel(qSENSITIVITY_TYPE typesensi,
										 const std::string& plot, 
										 const std::string& label,
										 int& outNoCurve,
										 double epsvalue)
{
	if (GetFlgClone() == false)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Function Not available for uncloned parameters");

    double sensitivity =0.;
	int NoCurve = 0;
	int i=0,h=0;
	// const ICM_DefaultCurve** portfolio = itsDefaultCurves;
	const std::vector<const ICM_DefaultCurve*>& portfolio = itsDefaultCurves; 
	// int nbcurves = GetNbDefCurves();
	int nbcurves = portfolio.size(); // GetNbDefCurves();

	std::vector<const ICM_DefaultCurve*>  ModifVectorDefCurves  (nbcurves) ; // = itsDefaultCurves; 
	// ICM_DefaultCurve** ModifVectorDefCurves = new ICM_DefaultCurve*[nbcurves];	
	for (i=0;i<nbcurves;i++) ModifVectorDefCurves[i] = NULL;

	double result = 0.0;
	bool Parallelshift = false,CurvesImpact = false; 

	vector<string> Term;
	char TermZC[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];
	memset(TermZC,'\0',sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);
	// ARM_Vector* epsilon = NULL;
	ICM_MktDataMng* pMNG = dynamic_cast<ICM_MktDataMng*>(this->GetMktDataMng());
	ICM_DefaultCurve* pDefCurveIndex = NULL;
	
				

	double EPSL = 0.01;		//Bump sur les taux fixé à 1bp	
	double EPSL_DEF = 0.01; //Bump sur le spread fixé à 1 bp.
	double EPSL_REC = 0.1;  //Bump fixe sur le recovery de 0.1
	double EPSL_CORR = 0.1;  //Bump fixe sur la correl de 0.1
	double EPSL_BETA = 0.1;  //Bump fixe sur la correl de 0.1
	double EPSL_SPREL = 0.1 ;//Bump relatif des spreads de 0.1
	double EPSL_DTR = CREDIT_DEFAULT_VALUE;//Bump relatif des spreads de 0.1

	if (epsvalue != CREDIT_DEFAULT_VALUE)
		EPSL_DTR = EPSL_BETA = EPSL_CORR = EPSL_REC = EPSL_DEF = EPSL = EPSL_SPREL = epsvalue;

	const ICM_DefaultCurve* InitDefCurve = NULL;
	ICM_DefaultCurve* ModifDefCurve = NULL;
	ARM_ZeroCurve* InitShortCurve = NULL;
	ARM_ZeroCurve* ModifShortCurve = NULL;
	ARM_ReferenceValue* Recovery = NULL;
	ARM_Vector* DiscreteValues = NULL,DiscreteDates = NULL;
	ICM_Correlation* Correlation = NULL;

	ARM_Vector epsilon (1,epsvalue);

	// if (!strcmp(plot,"NONE"))  //Detection d'un parallel Shift
	if (plot=="NONE")
	{
		Term.push_back("ALL");
		Parallelshift = true;
	}
	else
	{
		Term.push_back(plot);
		strcpy(TermZC[0],plot.c_str());
		switch (typesensi)
		{
			case ICMRECOVERY_TYPE :
			{
				epsilon.InitElt(0,EPSL_REC);
			}
			break;
			case ICMCORRELATION_TYPE :
			{
				epsilon.InitElt(0,EPSL_CORR);
			}
			break;
			case ICMIRCURVE_TYPE :
			case ICM_GREEK_RHO_TYPE :
			{
				epsilon.InitElt(0,EPSL);
			}
			break;
			default:
			case ICM_INDX_SPREAD_RESCALING :
			case ICMSPREAD_TYPE :
			{
				epsilon.InitElt(0,EPSL_DEF);
			}
			break;
			case ICMBETA_TYPE :
			{
				epsilon.InitElt(0,EPSL_BETA);
			}
			break;
			case ICMSPRELSHIFT_TYPE :
			{
				epsilon.InitElt(0,EPSL_SPREL);
			}
			break;
			case ICM_DTR_TYPE :
			{
				epsilon.InitElt(0,EPSL_DTR);
			}
			break;
		}
	}

	if	((typesensi == ICM_DTR_TYPE) || 
		(typesensi == ICMRECOVERY_TYPE) || 
		(typesensi == ICMIRCURVE_TYPE) ||
		(typesensi == ICMSPRELSHIFT_TYPE) || 
		(typesensi == ICM_INDX_SPREAD_RESCALING) || 
		(typesensi == ICMSPREAD_TYPE) ||
		(typesensi == ICM_GREEK_RHO_TYPE))  CurvesImpact = true;

	if (  label!="NONE"  && CurvesImpact )  //Detection d'un label particulier
	{
		// case of DefCurve in MarketDataMng
		if (pMNG) {
				pDefCurveIndex = dynamic_cast<ICM_DefaultCurve*>(pMNG->GetDefaultCurve(label,this->GetAsOfDate()));
				if (pDefCurveIndex) {	NoCurve = -1; }
		} 
		if (NoCurve == 0){		
			NoCurve = -1;
			for (h=0; h<nbcurves; h++)
			{
				// if (strcmp(label,portfolio[h]->GetLabel().c_str()) == NULL)
				if ( label==portfolio[h]->GetLabel() )	{
					NoCurve = h;
					break;
				}
			}
			if (NoCurve == -1)
				ICMTHROW(ERR_INVALID_ARGUMENT,"The curve label is unknown !");
		}
		for (h=0; h<nbcurves; h++)
			ModifVectorDefCurves[h] = (ICM_DefaultCurve*) portfolio[h]->Clone();

	} else {
		NoCurve = nbcurves;
		/*for (h=0; h<nbcurves; h++)
			ModifVectorDefCurves[h] = (ICM_DefaultCurve*) portfolio[h]->Clone();*/
	}
	outNoCurve = NoCurve;

// new model
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) Clone();
	ICM_MktDataMng* pnewMNG = dynamic_cast<ICM_MktDataMng*>(model->GetMktDataMng());
		// ---------------------------------------------------------------------------
		// Bump Recovery / IrCurve / Spread 
		// ---------------------------------------------------------------------------
		switch (typesensi)
		{
			case ICM_DTR_TYPE :
			{

				if (( NoCurve>=0) && (NoCurve != nbcurves))
				{
					if (ModifVectorDefCurves[NoCurve])
						delete ModifVectorDefCurves[NoCurve];
					ICM_DefaultCurve* item = (ICM_DefaultCurve*) portfolio[NoCurve]->GenerateShiftCurve(Term,epsilon,ICM_DTR_TYPE);
					item->SetLabel(portfolio[NoCurve]->GetLabel());

					ModifVectorDefCurves[NoCurve] =item; // ->SetLabel(portfolio[NoCurve]->GetLabel());
				}
				else
				{
					for (i = 0; i <nbcurves; i++)
					{
					if (ModifVectorDefCurves[i])
						delete ModifVectorDefCurves[i];
					ICM_DefaultCurve* item = (ICM_DefaultCurve*) portfolio[i]->GenerateShiftCurve(Term,epsilon,ICM_DTR_TYPE);
					item->SetLabel(portfolio[i]->GetLabel());
					// ModifVectorDefCurves[i] = (ICM_DefaultCurve*) portfolio[i]->GenerateShiftCurve(Term,epsilon,ICM_DTR_TYPE);

					ModifVectorDefCurves[i] =item; // ->SetLabel(portfolio[i]->GetLabel());
					}
				}

			}
			break;
			case ICMRECOVERY_TYPE :
			{

				if (( NoCurve>=0) && (NoCurve != nbcurves))
				{
					if (ModifVectorDefCurves[NoCurve])
						delete ModifVectorDefCurves[NoCurve];

					// if (Parallelshift)
					// 	ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) portfolio[NoCurve]->GenerateShiftCurve(EPSL_REC,ICMRECOVERY_TYPE);
					// else
					ICM_DefaultCurve*item = (ICM_DefaultCurve*) portfolio[NoCurve]->GenerateShiftCurve(Term,epsilon,ICMRECOVERY_TYPE);
					item ->SetLabel(portfolio[NoCurve]->GetLabel());

					// ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) portfolio[NoCurve]->GenerateShiftCurve(Term,epsilon,ICMRECOVERY_TYPE);

					ModifVectorDefCurves[NoCurve]=item; // ->SetLabel(portfolio[NoCurve]->GetLabel());
					if (!itsRecoveryRates.empty()) model->itsRecoveryRates[NoCurve] = MIN(model->itsRecoveryRates[NoCurve] + EPSL_REC,0.999);
				}
				else
				{
					for (i = 0; i <nbcurves; i++)
					{
					if (ModifVectorDefCurves[i])
						delete ModifVectorDefCurves[i];

					// if (Parallelshift)
					// 	ModifVectorDefCurves[i] = (ICM_DefaultCurve*) portfolio[i]->GenerateShiftCurve(EPSL_REC,ICMRECOVERY_TYPE);
					// else
					ICM_DefaultCurve* item = (ICM_DefaultCurve*) portfolio[i]->GenerateShiftCurve(Term,epsilon,ICMRECOVERY_TYPE);
					// ModifVectorDefCurves[i] = (ICM_DefaultCurve*) portfolio[i]->GenerateShiftCurve(Term,epsilon,ICMRECOVERY_TYPE);
					item ->SetLabel(portfolio[i]->GetLabel());
	
					ModifVectorDefCurves[i]=item; // ->SetLabel(portfolio[i]->GetLabel());
					if (!itsRecoveryRates.empty()) model->itsRecoveryRates[i] = MIN(model->itsRecoveryRates[i] + EPSL_REC,0.999);
					}
				}

			}
			break;
			case ICMIRCURVE_TYPE :
			case ICM_GREEK_RHO_TYPE :
			{

				if (( NoCurve>=0) && (NoCurve != nbcurves))
				{
					if (ModifVectorDefCurves[NoCurve])
						delete ModifVectorDefCurves[NoCurve];

					ICM_DefaultCurve* item = (ICM_DefaultCurve*) portfolio[NoCurve]->Clone(); 

					// ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) portfolio[NoCurve]->Clone();
					InitShortCurve = (ARM_ZeroCurve*) portfolio[NoCurve]->GetZeroCurve()->Clone();

					if (Parallelshift)
					{	
						ModifShortCurve = (ARM_ZeroCurve*)(InitShortCurve->Clone());
						ModifShortCurve->ParallelShift(EPSL);
					}
					else
						ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->GenerateShiftCurve(TermZC,&epsilon);


					// Impact sur la DefProbCurve  ********************************
					item/** ModifVectorDefCurves[NoCurve]**/ ->SetZeroCurve((ARM_ZeroCurve*)ModifShortCurve->Clone());
					item/**ModifVectorDefCurves[NoCurve]**/->Calibrate();
					item/** ModifVectorDefCurves[NoCurve]**/->SetLabel(portfolio[NoCurve]->GetLabel());
					ModifVectorDefCurves[NoCurve]=item; 

					if (InitShortCurve)
						delete InitShortCurve;
					InitShortCurve = NULL;

					if (ModifShortCurve)
						delete ModifShortCurve;
					ModifShortCurve = NULL;
				}
				else
				{
					for (i = 0; i <nbcurves; i++)
					{

					if (ModifVectorDefCurves[i])
						delete ModifVectorDefCurves[i];

					ICM_DefaultCurve* item =(ICM_DefaultCurve*) portfolio[i]->Clone();
					// ModifVectorDefCurves[i] = (ICM_DefaultCurve*) portfolio[i]->Clone();
					InitShortCurve = (ARM_ZeroCurve*) item->GetZeroCurve()->Clone();

					if (Parallelshift)
					{
						ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->Clone();
						ModifShortCurve->ParallelShift(EPSL);
					}
					else
						ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->GenerateShiftCurve(TermZC,&epsilon);


					// Impact sur la DefProbCurve  ********************************
					item->SetZeroCurve((ARM_ZeroCurve*) ModifShortCurve->Clone());
					item->Calibrate();
					item->SetLabel(portfolio[i]->GetLabel());
					ModifVectorDefCurves[i]=item; 

					if (InitShortCurve)
						delete InitShortCurve;
					InitShortCurve = NULL;

					if (ModifShortCurve)
						delete ModifShortCurve;
					ModifShortCurve = NULL;

					}
				}
			}
			break;
			case ICM_INFLATION_CPN_CURVE :
			{
				if (model->GetCpnInfCurve()==NULL) break; //pas de courbe inflation

				vector<double> vepsilon;
				vector<string> mktterms;
				vepsilon.resize(epsilon.GetSize());
				mktterms.resize(epsilon.GetSize());

				for (int il_ = 0; il_<epsilon.GetSize(); il_++)
				{vepsilon[il_]=epsilon.Elt(il_);
				 mktterms[il_]= (string)TermZC[il_];}

				ARM::ARM_InfCurv* ModifInfCurve = (ARM::ARM_InfCurv*) GetCpnInfCurve()->GenerateShiftCurve(mktterms,vepsilon);
				model->SetCpnInfCurve(ModifInfCurve);
				if (ModifInfCurve)
					delete ModifInfCurve;
				ModifInfCurve = NULL;
			}
			break;
			case ICM_IRCURVE_WITH_CPN:
			{
				//InitShortCurve = (ARM_ZeroCurve*) GetZeroCurve()->Clone();

				if (Parallelshift)
				{
					ModifShortCurve = (ARM_ZeroCurve*) GetZeroCurve()->Clone();
					ModifShortCurve->ParallelShift(EPSL);
				}
				else
					ModifShortCurve = (ARM_ZeroCurve*) GetZeroCurve()->GenerateShiftCurve(TermZC,&epsilon);

				//if (InitShortCurve) delete InitShortCurve;
				if (ModifShortCurve) model->SetZeroCurve(ModifShortCurve);
				

				if (model->GetCpnIRCurve()==NULL) break; //pas de coupon variable

			
				ARM_ZeroCurve* CpnIRCurve = (ARM_ZeroCurve*) model->GetCpnIRCurve()->Clone();
				ARM_ZeroCurve* ModifCpnIRCurve = NULL;
				
				/*if ( string(CpnIRCurve->GetCurrencyUnit() ->GetCcyName()) != string(ModifShortCurve->GetCurrencyUnit() ->GetCcyName())) 
				{
					if (CpnIRCurve) delete CpnIRCurve;
					CpnIRCurve = NULL;
					break;
				}*/

				if (Parallelshift)
				{
					ModifCpnIRCurve = (ARM_ZeroCurve*) CpnIRCurve->Clone();
					ModifCpnIRCurve->ParallelShift(EPSL);
				}
				else
					ModifCpnIRCurve = (ARM_ZeroCurve*) CpnIRCurve->GenerateShiftCurve(TermZC,&epsilon);

				if (CpnIRCurve) delete CpnIRCurve;
				CpnIRCurve = NULL;
				model->SetCpnIRCurve(ModifCpnIRCurve);
				if (ModifCpnIRCurve)
					delete ModifCpnIRCurve;
				ModifCpnIRCurve = NULL;
				
			}
			break;
			case ICM_INTEREST_CPN_CURVE :
			{
				if (model->GetCpnIRCurve()==NULL) break; //pas de courbe inflation
				ARM_ZeroCurve* CpnIRCurve = (ARM_ZeroCurve*) model->GetCpnIRCurve()->Clone();
				ARM_ZeroCurve* ModifCpnIRCurve = NULL;

				if (Parallelshift)
				{
					ModifCpnIRCurve = (ARM_ZeroCurve*) CpnIRCurve->Clone();
					ModifCpnIRCurve->ParallelShift(EPSL);
				}
				else
					ModifCpnIRCurve = (ARM_ZeroCurve*) CpnIRCurve->GenerateShiftCurve(TermZC,&epsilon);

				if (CpnIRCurve) delete CpnIRCurve;
				CpnIRCurve = NULL;
				model->SetCpnIRCurve(ModifCpnIRCurve);
				if (ModifCpnIRCurve)
					delete ModifCpnIRCurve;
				ModifCpnIRCurve = NULL;
			}
			break;
			case ICMIRCURVE_WITHOUT_DEFCURVE_TYPE :
			{
				InitShortCurve = (ARM_ZeroCurve*) GetZeroCurve()->Clone();

				if (Parallelshift)
				{
					ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->Clone();
					ModifShortCurve->ParallelShift(EPSL);
				}
				else
					ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->GenerateShiftCurve(TermZC,&epsilon);

				if (InitShortCurve) delete InitShortCurve;
				if (ModifShortCurve) model->SetZeroCurve(ModifShortCurve);
			}
			break;
			case ICMSPRELSHIFT_TYPE :
			{

				if (( NoCurve>=0) && (NoCurve != nbcurves) )
				{
					if (ModifVectorDefCurves[NoCurve])
						delete ModifVectorDefCurves[NoCurve];

					InitDefCurve = (ICM_DefaultCurve*)portfolio[NoCurve]->Clone();

					// if (Parallelshift)
					// 	ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(EPSL_DEF,ICMSPRELSHIFT_TYPE);
					// else
					ICM_DefaultCurve* item = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,epsilon,ICMSPRELSHIFT_TYPE);

					// ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,epsilon,ICMSPRELSHIFT_TYPE);

					item->SetLabel(portfolio[NoCurve]->GetLabel());
					ModifVectorDefCurves[NoCurve]=item; 

					if (InitDefCurve)
						delete InitDefCurve;
					InitDefCurve = NULL;

				}
				else
				{
					for (i = 0; i <nbcurves; i++)
					{
					if (ModifVectorDefCurves[i])
						delete ModifVectorDefCurves[i];

					InitDefCurve = (ICM_DefaultCurve*)portfolio[i]->Clone();

					// if (Parallelshift)
					// 	ModifVectorDefCurves[i] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(EPSL_DEF,ICMSPRELSHIFT_TYPE);
					// else
					ICM_DefaultCurve* item  = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,epsilon,ICMSPRELSHIFT_TYPE);
					item->SetLabel(portfolio[i]->GetLabel());	
					ModifVectorDefCurves[i]=item; 

					if (InitDefCurve)
						delete InitDefCurve;
					InitDefCurve = NULL;
					}
				}
				
			}
			break;
			case ICM_INDX_SPREAD_RESCALING :
			{
				if (ModifVectorDefCurves[NoCurve])
					delete ModifVectorDefCurves[NoCurve];

				InitDefCurve = (ICM_DefaultCurve*)portfolio[NoCurve]->Clone();

				// if (Parallelshift)
				// 	ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(EPSL_DEF,ICMSPREAD_TYPE);
				// else
				ICM_DefaultCurve* item = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,epsilon,ICMSPREAD_TYPE);
				// ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,epsilon,ICMSPREAD_TYPE);
				item->SetLabel(portfolio[NoCurve]->GetLabel());
				ModifVectorDefCurves[NoCurve]=item; 

				if (InitDefCurve)
					delete InitDefCurve;
				InitDefCurve = NULL;
			}
			break;
			case ICMSPREAD_TYPE :
			{
				// case of DefaultCurve In MarketDataMng for CorridorLeg
				if (NoCurve == -1) {
					// case of INDEX spread in CorridorLeg
					//pDefCurveIndex = dynamic_cast<ICM_DefaultCurve*>(pMNG->GetDefaultCurve(label,this->GetAsOfDate()));
					if (pDefCurveIndex){
						ICM_DefaultCurve* newDefaultCurve = NULL;
						// if (Parallelshift)
						// 	newDefaultCurve = dynamic_cast<ICM_DefaultCurve*>(pDefCurveIndex->GenerateShiftCurve(EPSL_DEF,ICMSPREAD_TYPE));
						// else
						newDefaultCurve = dynamic_cast<ICM_DefaultCurve*>(pDefCurveIndex->GenerateShiftCurve(Term,epsilon,ICMSPREAD_TYPE));
						pnewMNG->adopt(newDefaultCurve);
					}
				}
				else if (NoCurve == nbcurves) { // case all curves in the model ( not in the marketDataMng
					for (i = 0; i <nbcurves; i++)	{
						if (ModifVectorDefCurves[i])
							delete ModifVectorDefCurves[i];

						InitDefCurve = portfolio[i];

						// if (Parallelshift)
						// 	ModifVectorDefCurves[i] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(EPSL_DEF,ICMSPREAD_TYPE);
						// else
						ICM_DefaultCurve* item = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,epsilon,ICMSPREAD_TYPE);
						// ModifVectorDefCurves[i] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,epsilon,ICMSPREAD_TYPE);
						item->SetLabel(portfolio[i]->GetLabel());	
						ModifVectorDefCurves[i]=item; 
						InitDefCurve = NULL;
					}
				}
				else if (NoCurve>=0 ) 	{ // one Curve
					if (ModifVectorDefCurves[NoCurve])
						delete ModifVectorDefCurves[NoCurve];

					InitDefCurve = (ICM_DefaultCurve*)portfolio[NoCurve]->Clone();

					// if (Parallelshift)
					// 	ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(EPSL_DEF,ICMSPREAD_TYPE);
					// else

					 ICM_DefaultCurve* item = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,epsilon,ICMSPREAD_TYPE);
					 // ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,epsilon,ICMSPREAD_TYPE);
					item->SetLabel(portfolio[NoCurve]->GetLabel());
					ModifVectorDefCurves[NoCurve]=item; 

					if (InitDefCurve)
						delete InitDefCurve;
					InitDefCurve = NULL;
				}	
			}
			break;
			default :
			;
		}

	if (CurvesImpact) model->SetDefaultCurves(ModifVectorDefCurves);

	// ---------------------------------------------------------------------------
	// Bump Beta / Correlation
	// ---------------------------------------------------------------------------

	switch (typesensi)
	{
		case ICMBETA_TYPE :
		{
			Correlation = GetCorrelation()->GenerateShiftBetas(label,typesensi,EPSL_BETA);
		}
		break;
		case ICM_INDX_SPREAD_RESCALING :
		{
			Correlation = (ICM_Correlation*)GetCorrelation()->Clone();
			Correlation->ResetBetas();
		}
		break;
		case ICMCORRELATION_TYPE :
		{
			Correlation = GetCorrelation()->GenerateShiftBetas(label,typesensi,EPSL_CORR);
		}
		break;
		case ICMCORREL_STRIKE_DOWN_TYPE :
		{
			Correlation = GetCorrelation()->GenerateShiftBetas(label,typesensi,-EPSL_CORR*10.);
		}
		break;
		case ICMCORREL_STRIKE_UP_TYPE :
		{
			Correlation = GetCorrelation()->GenerateShiftBetas(label,typesensi,EPSL_CORR*10.);
		}
		break;
		default :
		;
	}

	if (Correlation)
	{
		model->SetCorrelation(Correlation);
		delete Correlation;
	}

	// if (epsilon) delete epsilon;
	if (InitDefCurve) delete InitDefCurve;
	if (ModifDefCurve) delete ModifDefCurve;
#ifdef	_DEBUG
	FILE * pFile = NULL;
	pFile = fopen("C:\\temp\\testModefShifted.txt","w");
	fprintf(pFile, "\n\n MODEL SHIFTED \n\n");
	model->View("",pFile);
	if (pFile) fclose(pFile);
#endif	

	return (model);     
}


// -----------------------------------------------------------------
// Returns the N° of curve in a Defcurve Package for a given Label
// -----------------------------------------------------------------
void ICM_ModelMultiCurves::CptNoCurve(const std::string& label,
									  int& outNoCurve) const 
{
	int NoCurve = -1;
	int nbcurves = GetNbDefCurves();

	for (int h=0; h<nbcurves; h++)
	{
		if (label==GetDefaultCurves(h)->GetLabel())
		{
		NoCurve = h;
		break;
		}
	}

	if (NoCurve == -1)
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_ModelMultiCurves::CptNoCurve: not found "<<label); 

	// }
	outNoCurve = NoCurve;
}

//
double ICM_ModelMultiCurves::GetRecoveryRate(const std::string& label) const
{
	for (int i = 0; i<itsDefaultCurves.size(); i++) 
	{
		if (itsDefaultCurves[i]->GetLabel()==label)
		{
			if (!itsRecoveryRates.empty()) 
				return itsRecoveryRates[i];
			else
				return (itsDefaultCurves[i])->GetRecovery();
		}
	}
	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_ModelMultiCurves::GetRecoveryRate: not found "<<label); 
}
//

void 
ICM_ModelMultiCurves::SetCorrelManager(ARM_CorrelManager* CorrelManager) 
{ if (GetFlgClone())
	{if (itsCorrelManager)
		delete itsCorrelManager;
	itsCorrelManager = (ARM_CorrelManager*) CorrelManager->Clone();}
	else
		itsCorrelManager = CorrelManager;
}
