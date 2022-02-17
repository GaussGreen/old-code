#include "ARMKernel\glob\firsttoinc.h"
#include "icm_option.h"
#include "ICMKernel/crv/icm_defaultcurve.h"

void ICM_Option::Init(void)
{
	SetName(ICM_OPTION);
	itsStrike = 0.0;
	itsOptionType = 1;
	itsKoType = (qDEF_MAT)1;
	itsUnderlyingMaturityTerm = "";
	itsNotional = 1;
}

ICM_Option::ICM_Option( const string&  underMaturity,
					   const ARM_Date& underMaturityDate,
					   const ARM_Date&	ExpiryDate,
					   const string&	ccy,
					   const qCDS_ADJ cds_adj,
					   bool endAdj,
					   double strike,
					   int optionType,
				       qDEF_MAT KO,
					   double notional)
{
	Init();
	SetCcy(ccy);
	if ( underMaturity != "") {
		itsUnderlyingMaturityTerm = underMaturity;
		itsUnderlyingMaturityDate = AddPeriod( ExpiryDate, underMaturity, ccy, endAdj, cds_adj);
	}else {
		itsUnderlyingMaturityTerm = "";
		itsUnderlyingMaturityDate = underMaturityDate;
	}
	Set( itsUnderlyingMaturityDate,ExpiryDate,strike,
		optionType, KO,notional);
}

ICM_Option::ICM_Option( const ARM_Date&  underMaturity,
					   const ARM_Date&	ExpiryDate,
					   double strike,
					   int optionType,
				       qDEF_MAT KO )
{
	Init();
	itsUnderlyingMaturityTerm = "";
	Set( underMaturity,ExpiryDate,strike,
		optionType, KO,1);

}

void ICM_Option::BitwiseCopy(const ARM_Object* srcOption)
{
	ICM_Option* option = (ICM_Option *) srcOption;

	itsStrike        = option->itsStrike;
	itsOptionType    = option->itsOptionType;
	SetExpiry(option->GetExpiry());
	itsKoType		 = option->itsKoType ;
	itsUnderlyingMaturityTerm = option->itsUnderlyingMaturityTerm;
	itsUnderlyingMaturityDate = option->itsUnderlyingMaturityDate;
	itsNotional = option->itsNotional;
}

// -------------
//	Copy Method 
// -------------
void ICM_Option::Copy(const ARM_Object* srccds)
{
	ARM_Object::Copy(srccds);
 
     BitwiseCopy(srccds);
}

// --------------
//	Clone Method
// --------------
ARM_Object* ICM_Option::Clone(void)
{
     ICM_Option* theClone = new ICM_Option();
     theClone->Copy(this);
     return(theClone);
}

void ICM_Option::Set( const ARM_Date&  underMaturity,	 
							const ARM_Date& ExpiryDate,
							double Strike,
							int optiontype,
							qDEF_MAT KO,
							double N)
{
	SetExpiry(ExpiryDate) ; // expiry of the option.
	SetStrike(Strike) ;
	SetOptionType(optiontype);
	SetUnderMaturityDate(underMaturity);
	SetKoType(KO);
	SetNotional(N);
}


void ICM_Option::View(char* id, FILE* ficOut)
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

    fprintf(fOut, "\t\t\t ----------------- Option ----------------- \n");
	fprintf(fOut,"Strike : %f\n",itsStrike);
	fprintf(fOut,"OptionType : %i\n",itsOptionType);
	fprintf(fOut,"KO Type : %i\n",itsKoType);
	fprintf(fOut,"Notional of option : %d\n",itsNotional);
	char* tmp = GetExpiry().GetStrDate() ;
	fprintf(fOut,"ExpiryDate : %s\n",tmp);
	delete [] tmp; 
	tmp = itsUnderlyingMaturityDate.GetStrDate() ;
	if(itsUnderlyingMaturityTerm != "") {
		fprintf(fOut,"Underlying Maturity Term : %s\n",itsUnderlyingMaturityTerm.c_str());
		fprintf(fOut,"Underlying Maturity Date from calculation : %s\n",tmp);
	} else { 
		fprintf(fOut,"Underlying Maturity Term : Not given");
		fprintf(fOut,"Underlying Maturity in Date : %s\n",tmp);		
	}
	delete [] tmp; 

	if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}