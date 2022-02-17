#include "ICMKernel\util\icm_pch.h"
#include "irindex.h"
#include "ICMKernel\crv\icm_fixing_curve.h"
#include "ICMKernel/util/icm_utils.h"

ICM_Fixing_Curve::~ICM_Fixing_Curve()
{
	if (itsIndex) delete itsIndex; 
	itsIndex=0;
}
void ICM_Fixing_Curve::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[40];
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

   fprintf(fOut, "\n\n--------------------------------------------------------------------------------\n");
   fprintf(fOut, "----------------  Fixing Curve -----------------------------\n");
   fprintf(fOut, "--------------------------------------------------------------------------------\n\n");
   fprintf(fOut, "Index Name : %s\n", itsIndexName.c_str());
   if (hasIndex()) fprintf(fOut, "Index Is defined."); 
   else fprintf(fOut, "Index is not defined."); 
   fprintf(fOut, "AsOfDate : %ld.%ld.%ld\n", itsAsOfDate.GetDay(), itsAsOfDate.GetMonth(), itsAsOfDate.GetYear());
   std::map<ARM_Date, double>::iterator it ;
   fprintf(fOut, "Date  \t Value \n");
   // if (itspMFixing) {
	   for ( it = itspMFixing.begin(); it != itspMFixing.end(); it++)
	   {
			fprintf(fOut, "%ld.%ld.%ld \t %lf\n", it->first.GetDay(), it->first.GetMonth(), it->first.GetYear(), it->second);
	   }
	// }	
	if ( ficOut == NULL ) {   fclose(fOut); }
}

ICM_Fixing_Curve::ICM_Fixing_Curve(const std::string& Index,
						const ARM_Date& AsOfDate,
						const std::map<ARM_Date, double>& aMap , const ARM_IRIndex*index)
{
	Init();
	Set(Index, AsOfDate, aMap,index) ;
}


ICM_Fixing_Curve::ICM_Fixing_Curve(const std::string& Index,
						const ARM_Date& AsOfDate,
						const ARM_Date& aDate, double aValue,const ARM_IRIndex*index )
{
	Init();
	Set(Index, AsOfDate , aDate, aValue,index) ;
}


void ICM_Fixing_Curve::Set (const std::string& Index,
						const ARM_Date& AsOfDate,
						const std::map<ARM_Date, double>& aMap ,const ARM_IRIndex*index)
{
	itsIndexName = Index;
	itsAsOfDate = AsOfDate;
	// if (itspMFixing) {itspMFixing->clear(); delete itspMFixing; }
	itspMFixing = aMap ; // new std::map<ARM_Date, double>;
	// *itspMFixing = aMap;
	if (itsIndex) delete itsIndex; itsIndex=0; 
	if (index) itsIndex=dynamic_cast<ARM_IRIndex*> ( unconst(*index).Clone() ); 
}

void ICM_Fixing_Curve::Set(const std::string& Index,
						const ARM_Date& AsOfDate,
						const ARM_Date& aDate, double aValue, const ARM_IRIndex*index)
{
	itsIndexName = Index;
	itsAsOfDate = AsOfDate;
	ARM_Date date = aDate;
	// if (itspMFixing) {itspMFixing->clear(); delete itspMFixing; }
	itspMFixing.clear(); 
	// itspMFixing = new std::map<ARM_Date, double>;
	itspMFixing[date] = aValue;
	if (itsIndex) delete itsIndex; itsIndex=0; 
	if (index) itsIndex=dynamic_cast<ARM_IRIndex*> ( unconst(*index).Clone() ); 

}



void ICM_Fixing_Curve::Init(void)
{
	SetName(ICM_FIXING_CURVE);
	itspMFixing.clear(); 
	itsIndex= NULL; 
}

void ICM_Fixing_Curve::BitwiseCopy(const ARM_Object* src) {
	ICM_Fixing_Curve* apICM_Fixing_Curve = (ICM_Fixing_Curve *)src;
	itsIndexName = apICM_Fixing_Curve->itsIndexName;
	itsAsOfDate = apICM_Fixing_Curve->itsAsOfDate;
	// if (itspMFixing) { itspMFixing->clear(); delete itspMFixing;}
	itspMFixing = apICM_Fixing_Curve->itspMFixing ; // std::map<ARM_Date, double>(*(apICM_Fixing_Curve->itspMFixing));
	if (apICM_Fixing_Curve->itsIndex) 
		itsIndex=dynamic_cast<ARM_IRIndex*> ( unconst(apICM_Fixing_Curve->itsIndex)->Clone() ); 
}


void ICM_Fixing_Curve::AddDateValue(const ARM_Date& aDate, double value) 
{
	itspMFixing[aDate]=value ; 
}

double ICM_Fixing_Curve::getValue(const ARM_Date&valueDate) const
{
	std::map<ARM_Date,double>::const_iterator it = itspMFixing.find(valueDate); 
	if (it==itspMFixing.end()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't find fixing at date "<<valueDate); 
	return it->second; 
}

bool ICM_Fixing_Curve::hasValue(const ARM_Date&valueDate) const
{
	std::map<ARM_Date,double>::const_iterator it = itspMFixing.find(valueDate); 
	return it!=itspMFixing.end(); 
}