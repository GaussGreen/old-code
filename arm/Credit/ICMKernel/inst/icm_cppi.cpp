#include "ICMKernel\inst\icm_cppi.h"



void ICM_Cppi::Init()
{
	SetName(ICM_CPPI);

	itsSecurity = NULL;
	itsCurrency = NULL;
	itsRangeFactor = NULL;

	// Caracteristique du produit
	itsNotional = 1.0;
	itsProtectedAmount = 1.0;
	itsManagementCost = 0.0;
	itsAdditionalLeverage = 0.0;
	itsDesactivateCushion = 0.0;

	// Output
	itsDefTime = 0.0;
	itsAvgLossValue = 0.0;
	itsNbLoss = 0.0;
	itsAvgFinalValue = 0.0;

	itsCorrelName = "TRAXX";
}


void ICM_Cppi::Set(	const ARM_Date& startDate,
					const ARM_Date& endDate,
					ARM_Security* sec,
					const double& notional,
					const double& ProtectedNotional,
					const double& AdditionalLeverage,
					ICM_RangeFactor* riskyFactor,
					const double& DesactivateCushion,
					const double& ManagementCost,
					ARM_Currency* ccy,
					string CorrelName)
{
	itsStartDate = startDate;
	itsEndDate   = endDate;
	if (itsSecurity)
		delete itsSecurity;
	itsSecurity  = dynamic_cast<ARM_Security*>( sec->Clone() );
	itsNotional  = notional;
	itsProtectedAmount = ProtectedNotional;
	itsManagementCost  = ManagementCost;
	itsAdditionalLeverage = AdditionalLeverage;
	itsDesactivateCushion = DesactivateCushion;
	if (itsCurrency)
		delete itsCurrency;
	itsCurrency = (ARM_Currency*)ccy->Clone();
	if (itsRangeFactor)
		delete itsRangeFactor;
	itsRangeFactor = (ICM_RangeFactor*) riskyFactor->Clone();

	itsCorrelName = CorrelName;
}

	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void ICM_Cppi::BitwiseCopy (const ARM_Object* src)
	{
	    ICM_Cppi* cppi= (ICM_Cppi*) src;

		itsStartDate	= cppi->itsStartDate;
		itsEndDate		= cppi->itsEndDate;
		itsSecurity		= cppi->itsSecurity;
		itsNotional		= cppi->itsNotional;
		itsProtectedAmount	= cppi->itsProtectedAmount;
		itsManagementCost	= cppi->itsManagementCost;
		itsAdditionalLeverage = cppi->itsAdditionalLeverage;
		itsDesactivateCushion = cppi->itsDesactivateCushion;
		
		itsDefTime = cppi->itsDefTime;
		itsAvgLossValue = cppi->itsAvgLossValue;
		itsNbLoss = cppi->itsNbLoss;
		itsAvgFinalValue = cppi->itsAvgFinalValue;
		itsCorrelName = cppi->itsCorrelName;
		
		if (itsCurrency /** && ( itsCurrency != ARM_DEFAULT_CURRENCY )**/  )
           delete itsCurrency;
 
        if (cppi->itsCurrency)
           itsCurrency = (ARM_Currency *) cppi->itsCurrency->Clone();
        else
           itsCurrency = NULL;

		itsRangeFactor	= cppi->itsRangeFactor;
	}

	// -------------
	//	Copy Method 
	// -------------
	void ICM_Cppi::Copy(const ARM_Object* src)
	{
		ARM_Object::Copy(src);
 
		BitwiseCopy(src);
	}

	// --------------
	//	Clone Method
	// --------------
	ARM_Object* ICM_Cppi::Clone(void)
	{
	ICM_Cppi* theClone = new ICM_Cppi();

	theClone->Copy(this);
 
    return(theClone);
	}

	// --------------
	//	View Method
	// --------------
	void ICM_Cppi::View(char* id, FILE* ficOut)
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

		char d1[20];
		char d2[20];

		itsStartDate.JulianToStrDate(d1);
		itsEndDate.JulianToStrDate(d2);

		fprintf(fOut, "\t\t\t ------------------ CPPI Informations ----------------- \n\n\n");


		fprintf(fOut, " Start Date : %14s \n",d1);
		fprintf(fOut, " End Date :   %14s \n",d2);
		fprintf(fOut, " Notional : %f \n",itsNotional);
		fprintf(fOut, " ProtectedAmount : %f \n",itsProtectedAmount);
		fprintf(fOut, " ManagementCost : %f \n",itsManagementCost);
		fprintf(fOut, " AdditionalLeverage : %f \n",itsAdditionalLeverage);
		fprintf(fOut, " DesactivateCushion : %f \n",itsDesactivateCushion);
		fprintf(fOut, " Correl Name: %s \n",itsCorrelName.c_str());
		fprintf(fOut, " DefTime : %f \n",itsDefTime);
		fprintf(fOut, " Average Loss Value: %f \n",itsAvgLossValue);
		fprintf(fOut, " Nb Loss: %f \n",itsNbLoss);
		fprintf(fOut, " Average Final Value: %f \n",itsAvgFinalValue);
		
		fprintf(fOut, "\n\n");
		fprintf(fOut, "\t\t\t ----------------- Security Informations ----------------- \n\n");
		itsSecurity->View(id,fOut);

		fprintf(fOut, "\n\n");
		fprintf(fOut, "\t\t\t ----------------- RangeFactor Informations ----------------- \n\n");
		itsRangeFactor->View(id,fOut);

		 if ( ficOut == NULL )
		{
			fclose(fOut);
		}
	}
