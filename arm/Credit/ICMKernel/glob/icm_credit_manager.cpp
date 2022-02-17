
#include "ICMKernel\glob\icm_credit_manager.h"

// GLOBAL Variables
ICM_Credit_Manager	The_Credit_Manager;

//------------------------------------------------
// Initialization of datas members
//------------------------------------------------

void ICM_Credit_Manager::Init(void)
{
	itsValDate	=	"01/01/2001";
	itsMaxCDSDate = "01/01/2001";
	itsZCValuesForCalibration = NULL;
	itsElapsed_Time	=	0.0;
	
}


void ICM_Credit_Manager::BitwiseCopy(const ARM_Object* src)
{
	ICM_Credit_Manager* dc = (ICM_Credit_Manager*) src;

	itsValDate = dc->itsValDate;
	itsMaxCDSDate = dc->itsMaxCDSDate;

	if (dc->itsZCValuesForCalibration)
		itsZCValuesForCalibration = (ARM_Vector*) dc->itsZCValuesForCalibration->Clone();

	itsElapsed_Time = dc->itsElapsed_Time;
}

void ICM_Credit_Manager::Copy(const ARM_Object* src)
{
   ARM_Object::Copy(src);

   BitwiseCopy(src);
}


void ICM_Credit_Manager::Reset()
{
	if (itsZCValuesForCalibration)
		delete itsZCValuesForCalibration;
	itsZCValuesForCalibration	=	NULL;
}


// some tests to be added?
double	ICM_Credit_Manager::DiscountPrice(int i)
{
	return (*itsZCValuesForCalibration)[i];
}



//--------------------------------------------
//	Display Method: TEMPORARY
//--------------------------------------------
void ICM_Credit_Manager::Display()
{
    FILE* fOut;
    char fOutName[100]	=	"C:\\Credit\\Excel\\display_Credit_Manager.txt";
 

    if ((fOut = fopen(fOutName, "w")) == NULL)
		// pb.
		return ;

	fprintf(fOut, "\n\n------------------------------------------------------------------\n");
	fprintf(fOut, "---------------- Credit Manager -----------------------------\n");
	fprintf(fOut, "------------------------------------------------------------------\n\n\n");

	fprintf(fOut, "Val Date\t:\t%ld.%ld.%ld\n\n", itsValDate.GetDay(), itsValDate.GetMonth(), itsValDate.GetYear());
	fprintf(fOut, "Max CDS Date\t:\t%ld.%ld.%ld\n\n", itsMaxCDSDate.GetDay(), itsMaxCDSDate.GetMonth(), itsMaxCDSDate.GetYear());
	fprintf(fOut, "Elapsed Time (in s)\t:\t%.lf\n\n", itsElapsed_Time);

	fprintf(fOut, "\n\nIndex\t\tZC Value\n");

	int NbPoints;
	NbPoints	=	itsZCValuesForCalibration->GetSize();

	for (int i = 0; i < NbPoints ; i++)
	{
		fprintf(fOut, "\n");	
		fprintf(fOut, "\t%lu\t\t%f", i, (*itsZCValuesForCalibration)[i]);
	}		 

	fclose(fOut);
 
}