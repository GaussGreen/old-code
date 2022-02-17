#include "ARMKernel\glob\firsttoinc.h"
//#include "ARMKernel\util\merge.h"
#include "ICMKernel\inst\ICM_CDO2.h"
//#include "ICMKernel\util\icm_utils.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\util\icm_matrix.h"

// ----------------------------
//	Copy of members data
// ----------------------------
void ICM_Cdo2::BitwiseCopy(const ARM_Object* src)
{
	int i = 0;
    ICM_Cdo2* cdo = (ICM_Cdo2 *) src;

	if (cdo->itsPortfolio)
	{
	if (itsPortfolio)
		delete itsPortfolio;
	itsPortfolio = (ICM_Portfolio*) cdo->itsPortfolio->Clone();
	}

	itsIsCrossSubordinate = cdo->itsIsCrossSubordinate;

}

// -------------
//	Copy Method 
// -------------
void ICM_Cdo2::Copy(const ARM_Object* srccds)
{
     ICM_Mez::Copy(srccds);
 
     BitwiseCopy(srccds);
}

// --------------
//	Clone Method
// --------------
ARM_Object* ICM_Cdo2::Clone(void)
{
     ICM_Cdo2* theClone = new ICM_Cdo2();

     theClone->Copy(this);
 
     return(theClone);
}

// ---------------------
//	Init of members data
// ---------------------

void ICM_Cdo2::Init()
{
	SetName(ICM_CDO2);

	itsPortfolio = NULL;
	itsIsCrossSubordinate = false;
}

// ----------------------------------------------
//	Constructor of CDS (Reference Date is char*)
// ----------------------------------------------

ICM_Cdo2::ICM_Cdo2(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const double& SubAmount,
				const double& MezzAmount,
				const int& Frequency,
				const int& DayCountBasis,
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string & Ccy, 
				const int& stubrule,
				ICM_Portfolio* ptf,
				const double& CreditLag,
				const int& FreqDefLeg,
				const double& Binary,
				const std::string &  payCalName,
				const bool& CrossSub,
				const bool& IncludeMaturity)
{

	Init();

	Set(EffectiveDate,
		ScheduleTerminateDate,
		FirstPeriodReferenceDate,
		FstCpnEffDate,
		FixedRate,
		intRule,
		adjStartDate,
		SubAmount,
		MezzAmount,
		Frequency,
		DayCountBasis,
		AccruedOnDefault,
		Ccy, 
		stubrule,
		ptf,
		CreditLag,
		FreqDefLeg,
		Binary,
		payCalName,
		CrossSub,
		IncludeMaturity);

}


// ----------------------------
//	Set Method of members data (Reference Date is char*)
// ----------------------------

void ICM_Cdo2::Set(const ARM_Date& EffectiveDate,
				const ARM_Date& ScheduleTerminateDate,
				const ARM_Date* FirstPeriodReferenceDate,
				const ARM_Date* FstCpnEffDate,
				const double& FixedRate,
				int intRule,
				int adjStartDate,
				const double& SubAmount,
				const double& MezzAmount,
				const int& Frequency,
				const int& DayCountBasis,
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				const std::string & Ccy, 
				const int& stubrule,
				ICM_Portfolio* ptf,
				const double& CreditLag,
				const int& FreqDefLeg,
				const double& Binary,
				const std::string &  payCalName,
				const bool& CrossSub,
				const bool& IncludeMaturity)
{
	// int nbissuers =0;	
	
	std::vector<std::string> issuers ;
	ARM_Vector notionals ;


	ptf->GetIssuersDatas(issuers,notionals);
	std::vector<double> notios(notionals.size()); 
	for(int i=0;i<notionals.size();i++) notios[i]=notionals.Elt(i);
	if (itsPortfolio)
		delete itsPortfolio;
	itsPortfolio = (ICM_Portfolio*) ptf->Clone();

	itsIsCrossSubordinate = CrossSub;
	
	ICM_Mez::Set(EffectiveDate,
				ScheduleTerminateDate,
				FirstPeriodReferenceDate,
				FstCpnEffDate,
				FixedRate,
				intRule,
				adjStartDate,
				SubAmount,
				fabs(MezzAmount),
				issuers,
				notios,
				//issuers.size(),
				Frequency,
				DayCountBasis,
				MezzAmount,
				AccruedOnDefault,
				Ccy, 
				MezzAmount,
				stubrule,
				CreditLag,
				FreqDefLeg,
				Binary,
				payCalName,
				qRunning_Leg,
				qStandart_Recovery_Leg,
				IncludeMaturity);

// 	delete[] notionals;
// 	FreePointerTabChar(issuers,nbissuers);	
}


void ICM_Cdo2::View(char* id, FILE* ficOut)
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

	fprintf(fOut, "\n");
	fprintf(fOut, "-----------------------------------------------------------------\n");
	fprintf(fOut, "-----					Cdo Square Viewer				   -----\n");
	fprintf(fOut, "-----------------------------------------------------------------\n\n");
	
	if (itsIsCrossSubordinate) fprintf(fOut, "Cross Subination : YES\n\n");

	ICM_Mez::View(id, fOut);

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}

//	-------------------------------------------------------------------------------------------
ICM_QMatrix<int>* 
ICM_Cdo2::GenAppMatrix(std::vector<std::string>& _UnionIssuers)
{
	ARM_Vector notionals ;

	GetPortfolio()->GetIssuersDatas(_UnionIssuers,notionals) ;
	// if (notionals) delete notionals;
	int sizeptf = GetPortfolio()->GetNbSec();
	int nbnames = _UnionIssuers.size(); 
	ICM_QMatrix<int>* MatAppartenance = new ICM_QMatrix<int>(nbnames,sizeptf,0);
	ICM_Mez* cdo = NULL;
#ifdef _DEBUG
	FILE* pFile;
	pFile = fopen("c:\\temp\\GenAppMatrix.txt", "w");
#endif 
	for (int i=0; i<nbnames;i++) 
		for (int j=0; j<sizeptf;j++)
		{
			cdo = (ICM_Mez*) GetPortfolio()->GetSecurities()[j];
			if (cdo->GetCollateral()->ContainIssuer(_UnionIssuers[i])) {
				(*MatAppartenance)(i,j) = 1;
#ifdef _DEBUG
				fprintf(pFile, "_UnionIssuers[i]= %s appartient à la tranche %d\n", _UnionIssuers[i].c_str(), j);
#endif 
			}
		}
#ifdef _DEBUG
	if (pFile) fclose(pFile);
#endif
	return (MatAppartenance);
}

ICM_QMatrix<int>*  ICM_Cdo2::AppMatrix(const vector<string>& vUnionIssuers){
	int sizeptf = GetPortfolio()->GetNbSec();
	ICM_QMatrix<int>* MatAppartenance = new ICM_QMatrix<int>(vUnionIssuers.size(),sizeptf,0);
	ICM_Mez* cdo = NULL;

	for (int i=0; i<vUnionIssuers.size();i++)
		for (int j=0; j<sizeptf;j++)
		{
			cdo = (ICM_Mez*) GetPortfolio()->GetSecurities()[j];
			if (cdo->GetCollateral()->ContainIssuer(vUnionIssuers[i])) 
			{	(*MatAppartenance)(i,j) = 1;
			}
		}
	return (MatAppartenance);
}
	

//	----------------------------------------------------------------------------------
void 
ICM_Cdo2::CptIssuersData(void)
{
	std::vector<std::string> issuers; 
	ARM_Vector notionals; 
	GetPortfolio()->GetIssuersDatas(issuers,notionals);
	if (GetCollateral()==NULL) 
	{
		SetCollateral(ICM_Collateral());
	} 
	GetCollateral()->SetIssuersInfos(issuers,notionals);
}
//	----------------------------------------------------------------------------------
ICM_Cdo2::ICM_Cdo2(ICM_Cds* cds,
		double SubAmount,
		ICM_Portfolio* pf,
		const double& Binary ) : ICM_Mez(cds,SubAmount,ICM_Collateral(),Binary)
{ Init();
  SetPortfolio((ICM_Portfolio*)pf->Clone());
}
