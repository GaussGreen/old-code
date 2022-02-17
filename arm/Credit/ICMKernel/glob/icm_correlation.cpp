#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\glob\icm_correlation.h"
#include "irindex.h"

//	---------------------------------------------------------------------
void 
ICM_Correlation::Init(void) 
{
	SetName(ICM_CORRELATION);
 	itsStructName = "CORR_NO_NAME";
	itsIndex1=0; 
	itsIndex2=0; 
}
//	---------------------------------------------------------------------
void 
ICM_Correlation::Set(const ARM_Date& AsOf,
		const std::vector<std::string>& labels,
		const std::string& name, 
		const ARM_IRIndex* index1,
		const ARM_IRIndex* index2)
{ 	
	itsLabels=labels; // SetLabels(label);
	SetStructName(name);
	itsAsOf = AsOf ;
	if (itsIndex1) delete itsIndex1; itsIndex1=0; 
	if (index1) itsIndex1 = dynamic_cast<ARM_IRIndex*>(unconst(*index1).Clone()); 
	if (itsIndex2) delete itsIndex2; itsIndex2=0; 
	if (index2) itsIndex2 = dynamic_cast<ARM_IRIndex*>(unconst(*index2).Clone()); 
}
 // ----------------------------
//	Copy of members data
// ----------------------------
void ICM_Correlation::BitwiseCopy(const ARM_Object* src)
{
	const ICM_Correlation* Correlation = (const ICM_Correlation *) src;

	itsLabels = Correlation->itsLabels ;
	
 	itsAsOf = Correlation->itsAsOf;
	itsStructName = Correlation->itsStructName;
	if (itsIndex1) delete itsIndex1; itsIndex1=0; 
	if (Correlation->itsIndex1) itsIndex1 = dynamic_cast<ARM_IRIndex*>(Correlation->itsIndex1->Clone()); 
	if (itsIndex2) delete itsIndex2; itsIndex2=0; 
	if (Correlation->itsIndex2) itsIndex2 = dynamic_cast<ARM_IRIndex*>(Correlation->itsIndex2->Clone()); 

}

// -------------
//	Copy Method 
// -------------
void ICM_Correlation::Copy(const ARM_Object* src)
{
	ARM_Object::Copy(src);

	BitwiseCopy(src);
}

// --------------
//	Clone Method
// --------------
ARM_Object* ICM_Correlation::Clone(void)
{
 ICM_Correlation* theClone = new ICM_Correlation();

 theClone->Copy(this);

 return(theClone);
}

//	---------------------------------------------------------------------
//virtual 
ARM_CLASS_NAME 
ICM_Correlation::GetRootName(void)
{
	return(ICM_CORRELATION);
}
//	---------------------------------------------------------------------
ICM_Correlation::ICM_Correlation(const ARM_Date& AsOf,
					const std::vector<std::string>& labels,
					const std::string& name, 
					const ARM_IRIndex* index1,
					const ARM_IRIndex* index2)
{
	Init();
	Set(AsOf,labels, name,index1,index2);
}
 //	---------------------------------------------------------------------
ICM_Correlation::~ICM_Correlation() 
{
 	if (itsIndex1) delete itsIndex1; itsIndex1=0; 
	if (itsIndex2) delete itsIndex2; itsIndex2=0; 
}


//	---------------------------------------------------------------------
int 
ICM_Correlation::GetLabelNo(const std::string&issuer) const 
{ 
	if (issuer=="NONE") { return CREDIT_DEFAULT_VALUE;}

	for (int i =0; i<itsLabels.size();i++)
		if (issuer==itsLabels[i]) { return i;}

	ICMTHROW(ERR_INVALID_MODEL,"ICM_Correlation::GetLabelNo: not found "<<issuer);
	
	return CREDIT_DEFAULT_VALUE;
}

//	---------------------------------------------------------------------
int		
ICM_Correlation::Get_IssuersRankFromLabel(const std::string& Label) const
{	
	int Index=0; 
	while (Label!=GetLabel(Index) && Index<itsLabels.size())
	{
		Index++;
	}
	
	return Index; 
}
//	---------------------------------------------------------------------
//	not virtual
ARM_Vector
ICM_Correlation::GetBetaVector(const std::vector<std::string>&labels,int size,double maturity,double strike)
{
	int j = 0;
	ARM_Vector vector(size);
	for ( j=0; j<size; j++) 
	{
		double beta = GetBeta(labels[j],5.,strike);	// virtual call 
		vector.Elt(j)=beta;				
	}
	return (vector);
}	

//	---------------------------------------------------------------------
void 
ICM_Correlation::View(char* id, FILE* ficOut)
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

		fprintf(fOut, "\n ======> Correlation Object :\n\n");
		char * tempDate= itsAsOf.GetStrDate(); 
		fprintf(fOut, "AsOfDate:%s\t\n", tempDate); 
		delete [] tempDate ; 
		fprintf(fOut, "StructName:%s\t\n", (char*)itsStructName.c_str()); 
		
		if (!itsLabels.empty())
		{

		int size = itsLabels.size();
		int k =0;

		for (i = 0; i<size; i++)
		{
			fprintf(fOut, "%s\t", itsLabels[i].c_str()); 
			fprintf(fOut, "\n");
		}

		}

		if ( ficOut == NULL )
		{
			fclose(fOut);
		}
}

void ICM_Correlation::SetIndex1(const ARM_IRIndex* index){
	if (itsIndex1) delete itsIndex1; itsIndex1=0; 
	if (index) itsIndex1 = dynamic_cast<ARM_IRIndex*>(unconst(*index).Clone()); 
}
void ICM_Correlation::SetIndex2(const ARM_IRIndex* index){
	if (itsIndex2) delete itsIndex2; itsIndex2=0; 
	if (index) itsIndex2 = dynamic_cast<ARM_IRIndex*>(unconst(*index).Clone()); 
}
void ICM_Correlation::SetAsOf(const ARM_Date& aDate){
	itsAsOf = aDate ;
}

void ICM_Correlation::SetLabels(const std::vector<std::string>& labels)
{
	itsLabels = labels;
}