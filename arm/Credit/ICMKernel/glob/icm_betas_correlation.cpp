#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\glob\icm_betas_correlation.h"


//	----------------------------------------------------------------------------------------------------
ICM_Beta_Correlation::ICM_Beta_Correlation(const ARM_Date& AsOf,
										   const std::string& name, 
					const ARM_Vector& betas, 
					const std::vector<std::string> &labels,
					const ARM_IRIndex*index1,
					const ARM_IRIndex*index2)
{
	Init();
	Set(AsOf,name,betas,labels,index1,index2);
}//	----------------------------------------------------------------------------------------------------
ICM_Beta_Correlation::ICM_Beta_Correlation(const ARM_Date& AsOf,
										   const double& beta_fixed,
										   const std::string& name,
										   const ARM_IRIndex*index1,
										   const ARM_IRIndex*index2)
{
	Init();
	std::vector<std::string> labels; 
	ICM_Correlation::Set(AsOf,labels,name,index1,index2 ); 

	if (itsBetas)
		delete itsBetas;
	itsBetas = new ARM_Vector(1,beta_fixed);

	itsFixedBeta = beta_fixed;
}
//	----------------------------------------------------------------------------------------------------
void ICM_Beta_Correlation::Init(void)
{
	ICM_Correlation::Init(); 
	SetName(ICM_BETA_CORRMATRIX);
	itsBetas = NULL;
	itsFixedBeta = CREDIT_DEFAULT_VALUE;
}
//	----------------------------------------------------------------------------------------------------
void ICM_Beta_Correlation::Set(const ARM_Date& AsOf,
		 const string& name,
		 const ARM_Vector& betas, 
		 const std::vector<std::string>& labels,
		 const ARM_IRIndex*index1,
		 const ARM_IRIndex*index2) 
{
	
	ICM_Correlation::Set(AsOf,labels,name,index1,index2 ) ;

	itsFixedBeta = CREDIT_DEFAULT_VALUE;

	if (!itsBetas) itsBetas=new ARM_Vector; 
	*itsBetas=betas; 
}
//	----------------------------------------------------------------------------------------------------
void ICM_Beta_Correlation::SetBetas(const ARM_Vector& betas) 
{
	if (!itsBetas) itsBetas=new ARM_Vector(); 
	*itsBetas=betas ; 
}

//	----------------------------------------------------------------------------------------------------
void ICM_Beta_Correlation::SetBetas(const ARM_Vector&betas,const std::vector<std::string>& labels) 
{
	int i=0,k=0;
	if (!itsBetas) itsBetas=new ARM_Vector(); 
	itsBetas->Resize(betas.size()); 
	for (k=0;k<betas.GetSize();k++)
	{
	i = GetLabelNo(labels[k]);	
	itsBetas->Elt(i)=betas.Elt(k);
	}
	itsFixedBeta = CREDIT_DEFAULT_VALUE;
}

//	----------------------------------------------------------------------------------------------------
ARM_Vector* ICM_Beta_Correlation::GetBetas(const std::vector<std::string>&labels,int size) const 
{
	int i=0,k=0;
	ARM_Vector * newbetas = new ARM_Vector (size,0.);

	for (k=0;k<size;k++)
	{
	i = GetLabelNo(labels[k].c_str());	

	if (itsFixedBeta == CREDIT_DEFAULT_VALUE)
		newbetas->Elt(k)=itsBetas->Elt(i);
	else
		newbetas->Elt(k)=itsFixedBeta;
	}

	return (newbetas);
}

// ----------------------------
//	Copy of members data
// ----------------------------
void ICM_Beta_Correlation::BitwiseCopy(const ARM_Object* src)
{
	ICM_Beta_Correlation* BetaCorrelation = (ICM_Beta_Correlation *) src;

	if (BetaCorrelation->itsBetas)
	{
		if (itsBetas)
			delete itsBetas;
		itsBetas = (ARM_Vector*)BetaCorrelation->itsBetas->Clone();
	}

	itsFixedBeta = BetaCorrelation->itsFixedBeta;
}
// -------------
//	Copy Method 
// -------------
void ICM_Beta_Correlation::Copy(const ARM_Object* src)
{
	ICM_Correlation::Copy(src);

	BitwiseCopy(src);
}

// --------------
//	Clone Method
// --------------
ARM_Object* ICM_Beta_Correlation::Clone(void)
{
 ICM_Beta_Correlation* theClone = new ICM_Beta_Correlation();

 theClone->Copy(this);

 return(theClone);
}


// --------------------------------------------------------------------
// View Mathod
// --------------------------------------------------------------------
void ICM_Beta_Correlation::View(char* id, FILE* ficOut)
{
	FILE* fOut;
	char  fOutName[200];
	int i = 0;

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

		fprintf(fOut, "\n ======> Betas Correlation Matrix :\n\n");
		
		if (!GetLabels().empty())
		{

		int size = itsBetas->GetSize();
		int k =0;

		for (i = 0; i<size; i++)
		{
			fprintf(fOut, "%s\t", GetLabel(i).c_str()); 
			fprintf(fOut, "\t %f ", itsBetas->Elt(i)); 
			fprintf(fOut, "\n");
		}

		}

		if (itsFixedBeta != CREDIT_DEFAULT_VALUE)
			fprintf(fOut, "Fixed Beta : %f\t", itsFixedBeta); 

		ICM_Correlation::View(id,ficOut);

		if ( ficOut == NULL )
		{
			fclose(fOut);
		}
}
//	----------------------------------------------------------------------------------------------------

//	virtual 
double ICM_Beta_Correlation::GetCorrelation(const std::string& issuer1,
											const std::string& issuer2,
							  double maturity ,
							  double strike ,
							  double actualYF )
{
	double correlation = GetBeta(issuer1) * GetBeta(issuer2);

	return (correlation);
}
//	----------------------------------------------------------------------------------------------------

// virtual 
double ICM_Beta_Correlation::GetBeta(const std::string&issuer,
					   double maturity ,
					   double strike ,
					   double actualYF )
{
	if (itsFixedBeta != CREDIT_DEFAULT_VALUE) return itsFixedBeta;

	bool find = false;

	for (int i=0; i<itsBetas->GetSize(); i++)
		// if (strcmp(issuer,GetLabel(i).c_str()) == NULL) { find = true; break; }
		if (issuer==GetLabel(i)) { find = true; break; }

	if (find) return itsBetas->Elt(i);

	ICMTHROW(ERR_INVALID_MODEL,"ICM_Beta_Correlation::GetBeta: not found "<<issuer); 
}
//	----------------------------------------------------------------------------------------------------

// virtual 
ICM_QMatrix<double> ICM_Beta_Correlation::ComputeCholeskyMatrix()
{
	ICM_QMatrix<double> matrix(1,1,0.);
	return matrix;
}

//	----------------------------------------------------------------------------------------------------

ICM_Beta_Correlation::~ICM_Beta_Correlation() 
{

if (itsBetas)
	delete itsBetas;
itsBetas = NULL;

};

// virtual 
ICM_Correlation* ICM_Beta_Correlation::GenerateShiftCorrel(const std::string&label,
											qSENSITIVITY_TYPE typesensi,
											double epsilon )
{ 
	ICM_Beta_Correlation* correl = (ICM_Beta_Correlation*) Clone(); 
	double defbeta = 0.999;

	if (itsFixedBeta != CREDIT_DEFAULT_VALUE)
	{
		correl->SetFixedBeta(MIN(sqrt(GetFixedBeta()*GetFixedBeta()+epsilon),defbeta));
		return correl;
	}
		
	int nolabel = GetLabelNo(label);

	ARM_Vector  Betas ( correl->GetBetas() );
	

	if (nolabel == CREDIT_DEFAULT_VALUE) 
	{
	for (int i=0; i<Betas.GetSize(); i++)
		{
			double value = MIN(sqrt(MAX(epsilon + Betas.Elt(i)*Betas.Elt(i),0.)),defbeta);

			Betas.Elt(i)=value;
		}	
	}
	else
	{
		Betas.Elt(nolabel) =
			MIN(sqrt(MAX(epsilon + Betas.Elt(nolabel)*Betas.Elt(nolabel),0.)),defbeta);
	}

	correl->SetBetas(Betas);

	return (correl);
}

// virtual 
ICM_Correlation* ICM_Beta_Correlation::GenerateShiftBetas(const std::string& label,
											qSENSITIVITY_TYPE typesensi,
											double epsilon )
{
	ICM_Beta_Correlation* correl = (ICM_Beta_Correlation*) Clone(); 

	if (itsFixedBeta != CREDIT_DEFAULT_VALUE)
	{
		correl->SetFixedBeta(GetFixedBeta()+epsilon);
		return correl;
	}
		
	if ((typesensi != ICMBETA_TYPE) && (typesensi != ICMCORRELATION_TYPE)) return (correl);
	
	int nolabel = GetLabelNo(label);
	double defbeta = 0.999;

	ARM_Vector  Betas (correl->GetBetas());

	if (nolabel == CREDIT_DEFAULT_VALUE) 
	{
	for (int i=0; i<Betas.GetSize(); i++)
		{
			double value = MIN(epsilon + Betas.Elt(i),defbeta);
			value = MAX(value,-defbeta);

			Betas.Elt(i)=value;
		}	
	}
	else
	{
		Betas.Elt(nolabel)=epsilon + Betas.Elt(nolabel);
	}

	correl->SetBetas(Betas);

	return (correl);
}