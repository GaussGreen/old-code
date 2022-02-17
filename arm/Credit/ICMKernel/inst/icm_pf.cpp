
#include "ICMKernel\inst\icm_pf.h"
#include "ICMKernel\inst\icm_collateral.h"
#include "ICMKernel/util/icm_matrix.h"

void ICM_Portfolio::Init(void)
{
	SetName(ICM_PF);

	itsParams = NULL;
	itsSecurities = NULL;  
	itsNbSec = 0;
}
// virtual 
ICM_Portfolio::~ICM_Portfolio()
{

	if (itsNbSec)
	{
		delete[] itsSecurities;
		itsSecurities = NULL;
	}

	if (itsParams)
		delete itsParams;
	itsParams = NULL;

	itsNbSec = 0;
}	
void ICM_Portfolio::SetParams(ICM_Parameters* params)
{
	if (itsParams)
		delete itsParams;
	itsParams = (ICM_Parameters*) params->Clone();
}
void ICM_Portfolio::Set(ARM_Security** securities,int nbofsecurities,ICM_Parameters* params)
{

	int i =0;

	if (itsSecurities)
	{
		delete[] itsSecurities;
		itsSecurities = NULL;
	}

	itsSecurities = new ARM_Security*[nbofsecurities];
	itsNbSec = nbofsecurities;

	for (i=0;i<itsNbSec;i++)
		itsSecurities[i] = (securities)[i];

	if (params)
	{
	if (itsParams)
		delete itsParams;
	itsParams = (ICM_Parameters*) params->Clone();
	}

}

// ----------------------------
//	Copy of members data
// ----------------------------
void ICM_Portfolio::BitwiseCopy(const ARM_Object* src)
{
	int i = 0;
    ICM_Portfolio* pf = (ICM_Portfolio*) src;

	if (pf->itsSecurities)
	{
		if (itsSecurities)
			delete[] itsSecurities;
		itsSecurities = NULL;

		itsSecurities = new ARM_Security*[pf->itsNbSec];

		for (i = 0; i<pf->itsNbSec; i++)
				itsSecurities[i] = (pf->itsSecurities)[i];

	}

	if (pf->itsParams)
	{
		if (itsParams)
			delete itsParams;
		itsParams = (ICM_Parameters*) pf->itsParams->Clone();
	}		

	itsNbSec = pf->itsNbSec;

}

// -------------
//	Copy Method 
// -------------
void ICM_Portfolio::Copy(const ARM_Object* src)
{
    ARM_Security::Copy(src);

    BitwiseCopy(src);
}

// --------------
//	Clone Method
// --------------
ARM_Object* ICM_Portfolio::Clone(void)
{
ICM_Portfolio* theClone = new ICM_Portfolio();

theClone->Copy(this);
 
return(theClone);
}

void ICM_Portfolio::View(char* id, FILE* ficOut)
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

	for (int i=0; i<itsNbSec ; i++)
	{
	// Affichage de la matrice
	fprintf(fOut, "\t\t\t ----------------- Security N°%ld ----------------- \n",i);

	itsSecurities[i]->View(id, fOut);
	}

	if ( ficOut == NULL )
	{
	fclose(fOut);
	}
}

void ICM_Portfolio::GetIssuersDatas(std::vector<std::string>& names,ARM_Vector&notional)
{
	names.resize(0); 
	notional.Resize(0); 

	//	no failure if the portfolio is empty. 
	//
	if (itsNbSec==0) return; 

	//	first security..  
	//
	int nbissuers = ((ICM_Ftd*)(itsSecurities[0]))->GetCollateral()->GetNbIssuers();	
	names = ((ICM_Ftd*)(itsSecurities[0]))->GetCollateral()->GetIssuersLabels();
	// double*  tmpnotio =  ((ICM_Ftd*)(itsSecurities[0]))->GetIssuersNotionals();
	notional.Resize(nbissuers); 
	for (int i=0; i<nbissuers; i++) 
	{
		notional.Elt(i) = ((ICM_Ftd*)(itsSecurities[0]))->GetCollateral()->GetIssuersNotional(i)
			* ((ICM_Ftd*)(itsSecurities[0]))->GetTradedCoef(); 
	}
	if (itsNbSec==1) return; 

	//	more than one security.. 
	//	we add remaining issuers, 
	// 	and update notionals with the last found in the basket. 
	// 	
	std::vector<double>  tmp_notionals; 
	notional.populate(tmp_notionals); 
	for (i=1; i<itsNbSec; i++)
	{
		int nbissuers = ((ICM_Ftd*)(itsSecurities[i]))->GetCollateral()->GetNbIssuers();	
		const std::vector<std::string>& labels2 =((ICM_Ftd*)(itsSecurities[i]))->GetCollateral()->GetIssuersLabels();	
		// double * tmpnotio = ((ICM_Ftd*)(itsSecurities[i]))->GetIssuersNotionals();
		for(int j=0;j<nbissuers;j++) 
		{
			std::vector<std::string>::const_iterator it = 
				std::find(names.begin(),names.end(),labels2[j]); 
			if (it!=names.end()) 
			{
				// name already present.... we do just update the notional... 
				tmp_notionals[ it - names.begin() ] = ((ICM_Ftd*)(itsSecurities[i]))->GetCollateral()->GetIssuersNotional(j) * ((ICM_Ftd*)(itsSecurities[i]))->GetTradedCoef(); 
			}
			else
			{
				// adding name 
				names.push_back(labels2[j]); 
				tmp_notionals.push_back(((ICM_Ftd*)(itsSecurities[i]))->GetCollateral()->GetIssuersNotional(j)  * ((ICM_Ftd*)(itsSecurities[i]))->GetTradedCoef() ); 
			}
		}
	}
	notional = tmp_notionals ;
}


