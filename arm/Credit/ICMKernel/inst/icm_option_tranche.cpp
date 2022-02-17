#include "ARMKernel\glob\firsttoinc.h"
#include "icm_option_tranche.h"
#include "ICMKernel/inst/icm_mez.h"

void ICM_Option_Tranche::Init(void)
{
	ICM_Option::Init();
	itsUnderlying = NULL;
	itsExerciceType = K_BERMUDAN ;
	itsRehauss = 0.0;
	itsDiffCDO = 0 ;
	itsObservationFreq=0.0 ;
	itsInitSpread = 0.0 ;
	itsIsCMSpread = 0 ;
	itsCMSpreadMatu = 5.0 ;


}

/*void ICM_Option_Tranche::Set( const ARM_Date& underMaturity,
						const ARM_Date& ExpiryDate,
						double strike,
						int optionType,
						ICM_Security* Underlying,
						double rehauss,
						double Freq,
						int ExecType,
						int DiffCDO)
{
	ICM_Option::Set(underMaturity,ExpiryDate,strike,optionType,q_KO,1.0);
	itsUnderlying = Underlying ;
	itsObservationFreq = Freq ;
	itsExerciceType = ExecType ;
	itsRehauss = rehauss;
	itsDiffCDO = DiffCDO;
	ICM_Security::SetStartDateNA(Underlying->GetStartDateNA()); // default value


}*/
void ICM_Option_Tranche::Set( const ARM_Date& TrigerStartDate,
						const ARM_Date& ExpiryDate,
						double strike,
						double InitSpread,
						int optionType,
						ICM_Security* Underlying,
						double rehauss,
						double Freq,
						int ExecType,
						int DiffCDO,
						int IsCMSpread,
						double CMSpreadMatu)
{

	ICM_Option::Set(Underlying ->GetEndDateNA(),ExpiryDate,strike,optionType,q_KO,1.0);
	itsInitSpread = InitSpread;
	itsUnderlying = Underlying ;
	itsObservationFreq = Freq ;
	itsExerciceType = ExecType ;
	itsRehauss = rehauss;
	itsDiffCDO = DiffCDO;
	itsIsCMSpread = IsCMSpread;
	itsCMSpreadMatu = CMSpreadMatu;
	ICM_Security::SetStartDateNA(TrigerStartDate); // default value

}
void ICM_Option_Tranche::View(char* id , FILE* ficOut )
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
	
	fprintf(fOut, "\t\t\t ----------------- Restrikable Tranche ----------------- \n");
	

	ICM_Option::View(id,fOut);
	
	fprintf(fOut,"Trigger Frequency : %f\n",itsObservationFreq);
	fprintf(fOut,"Rehaussement : %f\n",itsRehauss);
	fprintf(fOut,"Init Spread : %f\n",itsInitSpread);
	fprintf(fOut,"Is DiffCDO : %i\n",itsDiffCDO);
	fprintf(fOut,"Is CM Spread : %i\n",itsIsCMSpread);
	fprintf(fOut," CM Spread Matu : %f\n",itsCMSpreadMatu);

	fprintf(fOut, "\t\t\t ----------------- Underlying  ----------------- \n");
	
	if ( itsUnderlying)
	{
		((ICM_Mez*)itsUnderlying) ->View("",fOut);
	}
	else
	{
		fprintf(fOut,"\n\t\t\t ----------------- No Underlying Defined ----------------- \n\n");
	}
	
	if ( ficOut == NULL )
    {
       fclose(fOut);
    }

}