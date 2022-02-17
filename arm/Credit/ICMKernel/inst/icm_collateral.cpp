#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel\util\icm_utils.h"
#include <set>
#include <algorithm>
#include <deque>

// ---------------------
//	Init of members data
// ---------------------

//	-------------------------------------------------------------------------------------------------------------
ICM_Collateral::ICM_Collateral(const ICM_Collateral&ref) : ARM_Object(ref)
	
{
	itsNbIssuersShort=ref.itsNbIssuersShort; 
	itsIssuersLabels=ref.itsIssuersLabels; 
	itsIssuersNotionals=ref.itsIssuersNotionals;
	itsRecovCoef=ref.itsRecovCoef; 
	itsFlgFullHomogeneous=ref.itsFlgFullHomogeneous; 
	itsCollateralInDefault=0; 
	if (ref.itsCollateralInDefault) 
		itsCollateralInDefault=dynamic_cast<ICM_Collateral*>(ref.itsCollateralInDefault->Clone()) ;

}
//	-------------------------------------------------------------------------------------------------------------
ARM_Object* ICM_Collateral::Clone()
{
	return new ICM_Collateral(*this); 
}
//	-------------------------------------------------------------------------------------------------------------
void ICM_Collateral::Init()
{
	itsNbIssuersShort = 0;
	itsFlgFullHomogeneous = false;
	itsRecovCoef = CREDIT_DEFAULT_VALUE;
	if (itsCollateralInDefault) delete itsCollateralInDefault; itsCollateralInDefault=NULL; 
	itsIssuersLabels.resize(0); 
	itsIssuersNotionals.resize(0); 
}

//	-------------------------------------------------------------------------------------------------------------
ICM_Collateral::ICM_Collateral() 
{	
	SetName(ICM_COLLATERAL);
	itsCollateralInDefault=NULL; 
	Init();	
}	


//	-------------------------------------------------------------------------------------------------------------
// ICM_Collateral::ICM_Collateral(const std::vector<std::string>& Labels,
// 			   const ARM_Vector&  Notionals,
// 			   const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/&	IsInDefault,
// 			   const double& RecovCoef )
// {
// 	Init();
// 	Set(Labels,Notionals,IsInDefault,RecovCoef);
// }
 
//	-------------------------------------------------------------------------------------------------------------
ICM_Collateral::ICM_Collateral(const std::vector<std::string>&Labels,
			   const ARM_Vector&  Notionals )
{
	SetName(ICM_COLLATERAL);
	itsCollateralInDefault=NULL; 
	Init();
	Set(Labels,Notionals);
}

//	-------------------------------------------------------------------------------------------------------------

ICM_Collateral::~ICM_Collateral() 
{
	if (itsCollateralInDefault) delete itsCollateralInDefault;
	itsCollateralInDefault=0; 
}
//	-------------------------------------------------------------------------------------------------------------
// void 
// ICM_Collateral::Set(const std::vector<std::string>& Labels,
// 				   const ARM_Vector&  Notionals,
// 				   const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/&	IsInDefault,
// 				   const double& RecovCoef)
// {
// 	// do not init ?
// 	itsIssuersLabels=Labels; 
// 	itsIssuersNotionals.resize(itsIssuersLabels.size());
// 	for(int i=0;i<itsIssuersNotionals.size();i++) 
// 		itsIssuersNotionals[i]=Notionals[i]; 
// 	itsRecovCoef = RecovCoef;
// 	CptInternals(); 
// }	

//	-------------------------------------------------------------------------------------------------------------
void 
ICM_Collateral::Set(const std::vector<std::string>& Labels,
						 const ARM_Vector& Notionals)
{
	// do not init ? 
	itsIssuersLabels=Labels; 
	itsIssuersNotionals.resize(itsIssuersLabels.size());
	for (int i=0; i<itsIssuersNotionals.size(); i++) itsIssuersNotionals[i] = Notionals[i];
	CptInternals(); 
}
//	-------------------------------------------------------------------------------------------------------------
void
ICM_Collateral::Set(const std::vector<std::string>& labels,
					const ICM_Vector<ARM_ReferenceValue>& notionals)
{
	// do not init ?
	itsIssuersLabels=labels; 
	itsIssuersNotionals=notionals; 
	CptInternals(); 
}

//	-------------------------------------------------------------------------------------------------------------
void ICM_Collateral::View(char* id, FILE* ficOut)
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

	fprintf(fOut, "\t\t\t ----------------------------------------------------- \n");
    fprintf(fOut, "\t\t\t ----------------- Collateral Viewer ----------------- \n");
	fprintf(fOut, "\t\t\t ----------------------------------------------------- \n");
	//Global Recovery coefficient
	if (GetRecovCoefFlg()) {fprintf(fOut, "Recovery Coefficient :"); fprintf(fOut, "\titsRecovCoef");}

	
	fprintf(fOut, "\tIssuer\tAmount");
	fprintf(fOut, "\n\n");

	// Discret case
	if ( itsIssuersNotionals[0].IsConstant())
	{
		for (int i=0;i<itsIssuersLabels.size();i++)
		{
			fprintf(fOut,"%14s\t%14.0lf\t ",itsIssuersLabels[i].c_str(),itsIssuersNotionals[i].CptReferenceValue(0));
		
			fprintf(fOut,"\n");
		}
	}
	else
	{
		ICM_Vector<ARM_Date> notionalDates;
		GetNotionalDates(notionalDates);
		fprintf(fOut,"\t\t") ;

		for ( int j=0;j<notionalDates.size();j++)
		{
			char d[20] ;
			notionalDates[j].JulianToStrDate(d);
			fprintf(fOut,"%s\t",d);
		}
		
		fprintf(fOut,"\n");

		for (int i=0;i<itsIssuersLabels.size();i++)
		{
			fprintf(fOut,"%14s\t ",itsIssuersLabels[i].c_str()) ;

			for (int j=0;j<notionalDates.size();j++)
				fprintf(fOut,"%14.0lf\t ",itsIssuersNotionals[i].CptReferenceValue(notionalDates[j]));

			fprintf(fOut,"\n");
		}


	}

	if (itsCollateralInDefault){
	fprintf(fOut, "\n\t\t\t ----------------- Default Collateral ----------------- \n");
	itsCollateralInDefault->View(id,ficOut);
	}

	//End
	fprintf(fOut, "\n");

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}

}

//	-------------------------------------------------------------------------------------------------------------
void ICM_Collateral::ExcludeIssuer(const std::string&label)
{
	std::vector<std::string>::iterator it = std::find(itsIssuersLabels.begin(),itsIssuersLabels.end(),label); 
	if (it==itsIssuersLabels.end()) return; 
	unsigned int position = it - itsIssuersLabels.begin(); 
	itsIssuersLabels.erase(it); 
	itsIssuersNotionals.erase(position); 
	

	/** 
	std::vector<std::string> IssuersLabels(GetNbIssuers()-1); 
	ARM_Vector IssuersNotionals(GetNbIssuers()-1); 
	bool status = true;
	int i=0,k=0;
	vector<string> CcyNotional;		//Vector of Notional Ccy
	vector<string> CcyDeliver;		//Vector of Deliver Ccy
	vector<bool>  IsInDefault;		//Vector of Defaulted Issuers
	vector<ARM_Date> DefaultDate;	//Vector of Default Date
	vector<double> ObsRecovery;		//Vector of Obs Recovery

	if (CcyNotional.size()>0) CcyNotional.resize(GetNbIssuers()-1);
	if (CcyDeliver.size()>0) CcyDeliver.resize(GetNbIssuers()-1);
	if (IsInDefault.size()>0) IsInDefault.resize(GetNbIssuers()-1);
	if (DefaultDate.size()>0) DefaultDate.resize(GetNbIssuers()-1);
	if (ObsRecovery.size()>0) ObsRecovery.resize(GetNbIssuers()-1);

	for (i = 0; i<GetNbIssuers(); i++) 
	{
		if (GetIssuersLabels(i)!=label) // strcmp(GetIssuersLabels(i).c_str(),label))
		{	
			IssuersLabels[k]=GetIssuersLabels(i);
			IssuersNotionals[k] = itsIssuersNotionals[i];

			if (IsInDefault.size()>0) IsInDefault[k] = itsIsInDefault[i];
			k++;
		}
		else
			status= false;
	}


	SetIssuersInfos(IssuersLabels,IssuersNotionals);

	itsIsInDefault=IsInDefault;

	CheckIssuersInDefault();
	**/ 
}
//	-------------------------------------------------------------------------------------------------------------

// void ICM_Collateral::OrderFollowRank(const vector<double>& absLossRates,bool& AlreadyOrder)
// {
// 	AlreadyOrder = false;
// 
// 	if (itsOrderFollowLossUnit || itsFlgFullHomogeneous) {AlreadyOrder=true;return;}
// 
// 	int i=0,k=0,l=0;
// 	vector<int> InitialRank;InitialRank.resize(absLossRates.size());
// 	vector<int> NormalRank;NormalRank.resize(absLossRates.size());
// 	vector<int> FinalRank;FinalRank.resize(absLossRates.size());
// 	
// 	ARM_Matrix Matrix(absLossRates.size(),2);
// 
// 	for (i=0;i<absLossRates.size();i++) 
// 	{	Matrix.Elt(i,0) = absLossRates[i];
// 		Matrix.Elt(i,1) = (double)i;	}
// 
// 	Matrix.SortLines(0);	
// 
// 	for (i=0;i<absLossRates.size();i++)
// 	{InitialRank[i] = (int)Matrix.Elt(i,1);
// 	 NormalRank[i] = i;}
// 
// 	if (NormalRank==InitialRank){AlreadyOrder=true;itsOrderFollowLossUnit=true;return;}
// 
// 	for (k=0;k<absLossRates.size();k++)
// 	for (i=0;i<absLossRates.size();i++)
// 		{ if (InitialRank[i]==k) 
// 			{FinalRank[k] = i;break;}}
// 
// 	vector<string>	m_itsIssuersLabels;	m_itsIssuersLabels.resize(itsNbIssuers);
// 	vector<double>  m_itsIssuersNotionals;m_itsIssuersNotionals.resize(itsNbIssuers);
// 	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/    m_itsIsInDefault;m_itsIsInDefault.resize(itsIsInDefault.size());		
// 	vector<int>		m_itsShortIndices;m_itsShortIndices.resize(0);
// 	vector<int>		m_itsLongIndices;m_itsLongIndices.resize(0);
// 
// 	for (k=0;k<absLossRates.size();k++)
// 	{
// 		m_itsIssuersLabels[FinalRank[k]] =  itsIssuersLabels[k];
// 		m_itsIssuersNotionals[FinalRank[k]]=itsIssuersNotionals[k];
// 
// 		if ((!itsIsInDefault.empty()) && (itsIsInDefault.size() == itsNbIssuers))
// 			m_itsIsInDefault[FinalRank[k]]=itsIsInDefault[k];
// 
// 	}
// 
// 	//If Long Short : Indices Vector
// 	for (k=0; k<itsNbIssuers; k++) 
// 	{
// 		if (lt(m_itsIssuersNotionals[k],0.))
// 			m_itsShortIndices.push_back(k);
// 		else 
// 			m_itsLongIndices.push_back(k);
// 	}
// 
// 	itsIsInDefault = m_itsIsInDefault;		
// 
// 	for (k=0;k<absLossRates.size();k++)
// 	{	
// 		itsIssuersLabels[k]=m_itsIssuersLabels[k] ;
// 		itsIssuersNotionals[k] = m_itsIssuersNotionals[k];
// 	}
// 
// 	itsOrderFollowLossUnit=true;
// }
//	-------------------------------------------------------------------------------------------------------------

double ICM_Collateral::SumNotionals(const ARM_Date&date) const 
{
	double result =0.;
	for (int i=0; i<itsIssuersNotionals.size(); i++) 
	{
		// this special line is for "ALL NOTIONALS <0 case" . 
		if ((IsLongShort())&&(itsIssuersNotionals[i].CptReferenceValue(date)<0)) continue;
		result += itsIssuersNotionals[i].CptReferenceValue(date);
	}
	return (result);
}

//	-------------------------------------------------------------------------------------------------------------

/** void ICM_Collateral::BitwiseCopy(const ARM_Object* src)
{
	ICM_Collateral* info = (ICM_Collateral*) src;
//	itsNbIssuers = info->itsNbIssuers;
	itsNbIssuersShort = info->itsNbIssuersShort;
	itsIssuersLabels = info->itsIssuersLabels; 

	itsIssuersNotionals = info->itsIssuersNotionals ;
	// itsIssuersNotionals.resize(info->itsNbIssuers);
	
	// for (int i=0;i<info->itsNbIssuers;i++) 
	// { 
	// 	itsIssuersNotionals[i] = info->itsIssuersNotionals[i];
	// }

	itsIsInDefault = info->itsIsInDefault;

	itsRecovCoef = info->itsRecovCoef;

	itsFlgFullHomogeneous = info->itsFlgFullHomogeneous;
}
//	-------------------------------------------------------------------------------------------------------------

void ICM_Collateral::Copy(const ARM_Object* src)
{
 ARM_Object::Copy(src);

 BitwiseCopy(src);
}

ARM_Object* ICM_Collateral::Clone(void)
{
 ICM_Collateral* theClone = new ICM_Collateral();

 theClone->Copy(this);

 return(theClone);
} **/ 
//	-------------------------------------------------------------------------------------------------------------
void 
ICM_Collateral::CptFullHomog(void) 
{
	itsFlgFullHomogeneous=false;

	for (int i=1;i<itsIssuersLabels.size();i++)
	{ 
		// if (strcmp(itsIssuersLabels[i-1].c_str(),itsIssuersLabels[i].c_str())) {return;}
		if (itsIssuersLabels[i-1] != itsIssuersLabels[i] ) return;
	}
	itsFlgFullHomogeneous=true;
}

//	-------------------------------------------------------------------------------------------------------------

 bool ICM_Collateral::ContainIssuer(const std::string& label) const 
{
	int i=0;

	for (i = 0; i<GetNbIssuers(); i++) 
		if (GetIssuersLabels(i)==label) return (true);
	return (false);
}

//	-------------------------------------------------------------------------------------------------------------

 void ICM_Collateral::SetIssuersInfos(const std::vector<std::string>&issuersLabels, const ARM_Vector& IssuersNotionals) 
{
	int i =0;
	// itsNbIssuers = issuersLabels.size(); 
	itsIssuersLabels=issuersLabels; 
 
	itsIssuersNotionals.resize(itsIssuersLabels.size());
	
	//Long Short
	itsNbIssuersShort = 0;
	
	for (i=0; i<itsIssuersLabels.size(); i++) 
	{
		itsIssuersNotionals[i] = IssuersNotionals[i];
		if (lt(itsIssuersNotionals[i].CptReferenceValue(0),0.))itsNbIssuersShort++; 
	}


}
//	-------------------------------------------------------------------------------------------------------------

/** void ICM_Collateral::CheckIssuersInDefault()
{
	int i=0;
	string s;
	
	if (itsIsInDefault.size() != itsIssuersLabels.size())
	{
		itsIsInDefault.resize(itsIssuersLabels.size());
	}

	for (i=0; i<itsIssuersLabels.size(); i++)
	{
		if (itsIssuersLabels[i].find(ISSUER_IN_DEFAULT)!=-1) {itsIsInDefault[i]=true;}
	}
}
**/ 
//	-------------------------------------------------------------------------------------------------------------

 void ICM_Collateral::GetIssuersInDefault(vector<string>& issuers,const bool& indef)
{
	issuers.clear();
	// if (itsIsInDefault.size()==0) return;
	for (int i=0; i<itsIssuersLabels.size(); i++) 
	{ 
		bool isIndefault = itsIssuersLabels[i].find(ISSUER_IN_DEFAULT)!=-1; 
		if (isIndefault && indef) 
			issuers.push_back(itsIssuersLabels[i]);
		else if ( isIndefault==false  &&  indef==false ) 
			issuers.push_back(itsIssuersLabels[i]);
	}
}	
//	-------------------------------------------------------------------------------------------------------------

 void ICM_Collateral::TransferIssuerInDefaultCollateral(const std::string& label)
{
	if (itsCollateralInDefault==NULL) itsCollateralInDefault = new ICM_Collateral();	

	for (int k = 0; k<itsIssuersLabels.size(); k++) 
	{	
		if (GetIssuersLabels(k)!=label)
		{
			itsCollateralInDefault->PushIssuer(label,
											   itsIssuersNotionals[k]
											 );

			ExcludeIssuer(label);
			break;
		}	
	}
}
//	-------------------------------------------------------------------------------------------------------------

 void ICM_Collateral::PushIssuer(const std::string& label,
									   const ARM_ReferenceValue& notional
									   )
{
	itsIssuersLabels.push_back(label); 
	itsIssuersNotionals.push_back(notional); 
	CptInternals(); 

	/** 
	vector<string> CcyNotional;		//Vector of Notional Ccy
	vector<string> CcyDeliver;		//Vector of Deliver Ccy
    vector<bool> IsInDefault;		//Vector of Defaulted Issuers
	vector<ARM_Date> DefaultDate;	//Vector of Default Date
	vector<double> ObsRecovery;		//Vector of Obs Recovery
	vector<int>	ShortIndices;		//Vector of Short Index
	vector<int>	LongIndices;		//Vector of Long Index

	if (itsIsInDefault.size()>0) IsInDefault.resize(itsNbIssuers+1);

	std::vector<std::string> IssuersLabels(itsNbIssuers+1); 
 	ARM_Vector IssuersNotionals(itsNbIssuers+1); 
 
	for (int i = 0; i<itsNbIssuers; i++) 
	{	 
		IssuersLabels[i]=itsIssuersLabels[i];	
		IssuersNotionals[i]=itsIssuersNotionals[i];
 
		if (IsInDefault.size()>0) IsInDefault[i] = itsIsInDefault[i];
	}

	if (!_IssuersLabel.empty())
	{	 
		IssuersLabels[itsNbIssuers]=_IssuersLabel ;
	}
		
	if (_notional) IssuersNotionals[itsNbIssuers] = *_notional;
 
	if ((IsInDefault.size()>0)&&(_isindefault)) IsInDefault[itsNbIssuers] = *_isindefault;

	SetIssuersInfos(IssuersLabels,IssuersNotionals);
 
	itsIsInDefault=IsInDefault; **/ 

}
//	-------------------------------------------------------------------------------------------------------------
const std::string& 
ICM_Collateral::GetIssuersLabels(unsigned int i) const
{
	if (i>=itsIssuersLabels.size()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Collateral::GetIssuersLabels: out of bond "<<i); 
	return itsIssuersLabels[i]; 
}
//	-------------------------------------------------------------------------------------------------------------
bool 
ICM_Collateral::isNotionalConstant() const
{
	for(int i=0;i<itsIssuersNotionals.size();i++) 
	{
		if (!itsIssuersNotionals[i].IsConstant()) return false ;
	}
	return true; 
}
//	-------------------------------------------------------------------------------------------------------------
void 
ICM_Collateral::CptInternals()
{
	itsNbIssuersShort=0; 
	for (int i=0; i<itsIssuersNotionals.size(); i++) 
		if (lt(itsIssuersNotionals[i].CptReferenceValue(0),0.)) itsNbIssuersShort++ ; 
	CptFullHomog();
}
//	-------------------------------------------------------------------------------------------------------------
void 
ICM_Collateral::GetConstantNotionals(ARM_Vector&notionals) const 
{
	notionals.Resize(itsIssuersNotionals.size()); 
	for(int i=0;i<itsIssuersNotionals.size();i++) 
	{
		if (!itsIssuersNotionals[i].IsConstant()) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Collateral::GetConstantNotionals: not constant"); 
		notionals[i]=itsIssuersNotionals[i].CptReferenceValue(0);  
	}
}
//	-------------------------------------------------------------------------------------------------------------
void 
ICM_Collateral::GetNotionalDates(ICM_Vector<ARM_Date>&ret) const
{
	std::set<double> alldates; 
	for(int i=0;i<itsIssuersNotionals.size();i++) 
	{
		if (!itsIssuersNotionals[i].IsConstant()) 
		{
			const ARM_Vector* dates = itsIssuersNotionals[i].GetDiscreteDates();
			if (dates) 
				alldates.insert(dates->begin(),dates->end()); 
		}
	}
	std::set<double>::const_iterator it = alldates.begin(); 
	ret.resize(alldates.size()); 
	i=0; 
	while (it!=alldates.end()) 
	{
		ret[i++]=*it; 
		++it; 
	}
}