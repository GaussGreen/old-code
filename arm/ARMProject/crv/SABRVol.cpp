#include "SABRVol.h"
#include "volcurv.h"

//--------------------------------------------------------------//
// Initialization of parameters									//
//--------------------------------------------------------------//
void ARM_SABRVol::Init()
{
	SetName(ARM_VOL_SABR);

	// Volatility curves

	itsSigmaOrAlpha		=	NULL;

	itsBeta				=	NULL;
	itsRho				=	NULL;
	itsNu				=	NULL;

	// Flags

	itsSigmaOrAlphaFlag	=	1;

	itsModelType		=	K_LD;

	// Weight

	itsWeight			=	0.5;
}


//--------------------------------------------------------------//
// Copy of parameters											//
//--------------------------------------------------------------//
void ARM_SABRVol::BitWiseCopy(const ARM_Object* srcSABRVol)
{
	ARM_SABRVol* tmpSABRVol = (ARM_SABRVol*) srcSABRVol;

	if ( itsSigmaOrAlpha )
	{
		delete itsSigmaOrAlpha;
		itsSigmaOrAlpha = NULL;
	}

	itsSigmaOrAlpha = new ARM_VolCurve(*(tmpSABRVol->itsSigmaOrAlpha));

	if ( itsBeta )
	{
		delete itsBeta;
		itsBeta = NULL;
	}

	itsBeta = new ARM_VolCurve(*(tmpSABRVol->itsBeta));

	if ( itsRho )
	{
		delete itsRho;
		itsRho = NULL;
	}

	itsRho = new ARM_VolCurve(*(tmpSABRVol->itsRho));

	if ( itsNu )
	{
		delete itsNu;
		itsNu = NULL;
	}

	itsNu = new ARM_VolCurve(*(tmpSABRVol->itsNu));

	itsSigmaOrAlphaFlag = tmpSABRVol->itsSigmaOrAlphaFlag;

	itsModelType		= tmpSABRVol->itsModelType;

	itsWeight			= tmpSABRVol->itsWeight;
}


//--------------------------------------------------------------//
// Global copy													//
//--------------------------------------------------------------//
void ARM_SABRVol::Copy(const ARM_Object* srcSABRVol)
{
	ARM_VolCurve::Copy(srcSABRVol);

	BitWiseCopy(srcSABRVol);
}


//--------------------------------------------------------------//
// Clone method													//
//--------------------------------------------------------------//
ARM_Object* ARM_SABRVol::Clone()
{
	ARM_SABRVol* theClone = new ARM_SABRVol();

	theClone->Copy(this);

	return theClone;
}


//--------------------------------------------------------------//
// Assignment operator											//
//--------------------------------------------------------------//
ARM_SABRVol& ARM_SABRVol::operator = (const ARM_SABRVol& SABRVol)
{
	(*this).ARM_VolCurve::operator = (SABRVol);

	BitWiseCopy(&SABRVol);

	return (*this);
}


//--------------------------------------------------------------//
// Default constructor											//
//--------------------------------------------------------------//
ARM_SABRVol::ARM_SABRVol()
{
	Init();
}


//--------------------------------------------------------------//
// Copy constructor												//
//--------------------------------------------------------------//
ARM_SABRVol::ARM_SABRVol(const ARM_SABRVol& SABRVol)
			:ARM_VolCurve(SABRVol)
{
	Init();

	BitWiseCopy(&SABRVol);
}


//--------------------------------------------------------------//
// Main constructor												//
//--------------------------------------------------------------//
ARM_SABRVol::ARM_SABRVol(ARM_VolCurve*	aSigmaOrAlpha,
						 ARM_VolCurve*	aBeta,
						 ARM_VolCurve*	aRho,
						 ARM_VolCurve*	aNu,
						 int			aSigmaOrAlphaFlag,
						 int			aModelType,
						 double			aWeight)
{
	// Initialization

	Init();

	
	// Put the flags to their right value

	itsSigmaOrAlphaFlag = aSigmaOrAlphaFlag;
	itsModelType		= aModelType;

	
	// Load volatility curves by using clone methods

	if(aSigmaOrAlpha)
	{
		itsSigmaOrAlpha = (ARM_VolCurve*) aSigmaOrAlpha->Clone();
	}
	
	if(aBeta)
	{
		itsBeta = (ARM_VolCurve*) aBeta->Clone();
	}

	if(aRho)
	{
		itsRho = (ARM_VolCurve*) aRho->Clone();
	}

	if(aNu)
	{
		itsNu = (ARM_VolCurve*) aNu->Clone();
	}

	// Load Option Type parameter

	SetOptionType( (int) aSigmaOrAlpha->GetOptionType() );


	// Load Weight

	itsWeight = aWeight;
}


//--------------------------------------------------------------//
// Destructor													//
//--------------------------------------------------------------//
ARM_SABRVol::~ARM_SABRVol()
{
	if ( itsSigmaOrAlpha )
	{
		delete itsSigmaOrAlpha;
	}

	itsSigmaOrAlpha	= NULL;

	if ( itsBeta )
	{
		delete itsBeta;
	}

	itsBeta = NULL;

	if ( itsRho )
	{
		delete itsRho;
	}

	itsRho = NULL;

	if ( itsNu )
	{
		delete itsNu;
	}

	itsNu = NULL;
}


//--------------------------------------------------------------//
// Viewer														//
//--------------------------------------------------------------//
void ARM_SABRVol::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];
	char* ModelType;


    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

    fprintf(fOut, "\n ====> SABR Volatility Structure :\n");

	if ( GetOptionType() == K_IRG )
	{
		fprintf(fOut, "\n ====> Volatility type: IRG \n");
	}
	else
	{
		fprintf(fOut, "\n ====> Volatility type: SWOPT \n");
	}


	// Smiled Model Type
	//-------------------------------------------------------
	
	switch ( GetModelType() )
	{
		case K_LD:
			ModelType = "LD";
			break;
		case K_SABR_GEO:
			ModelType = "SABR_G";
			break;
		case K_SABR_ARITH:
			ModelType = "SABR_A";
			break;
		case K_SABR_IMPLNVOL:
			ModelType = "SABR_IMPLNVOL";
			break;
		case K_SABR_WEIGHT:
			ModelType = "SABR_WEIGHT";
			break;
		default:
			ModelType = "Unknown type";
			break;
	}

	fprintf(fOut, "\n ====> Smiled model type: %s\n", ModelType);


	// Weight if ModelType == SABR_WEIGHT
	//-------------------------------------------------------
	if ( ModelType == "SABR_WEIGHT" )
	{
		fprintf(fOut, "\n ====> SABR Weight: %d\n", itsWeight);
	}


	// Volatility curves printing
	//-------------------------------------------------------
    
	fprintf(fOut, "\n\n\n ====> SIGMA :\n");
    itsSigmaOrAlpha->View(id, fOut);

    fprintf(fOut, "\n\n\n ====> RHO :\n");
    itsRho->View(id, fOut);

    fprintf(fOut, "\n\n\n ====> NU :\n");
    itsNu->View(id, fOut);

    fprintf(fOut, "\n\n\n ====> BETA :\n");
    if (itsBeta) 
       itsBeta->View(id, fOut);
    else
       fprintf(fOut, "\n ====> NO BETA <> 1");


	// Sigma Or Alpha
	//-------------------------------------------------------

    if (itsSigmaOrAlphaFlag && itsBeta)
    {
       fprintf(fOut, "\n ====> SigmaOrAlphaInput : SIGMA\n");
    }
    else
    {
       if (itsBeta)
       {
          fprintf(fOut, "\n ====> SigmaOrAlphaFlag : ALPHA considered\n");
       }
       else
       {
          fprintf(fOut, "\n ====> SigmaOrAlphaFlag : SIGMA considered\n");
       }
    }

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}


//--------------------------------------------------------------//
// Calculation of volatility within SABR model					//
//--------------------------------------------------------------//
double ARM_SABRVol::computeVol(double aMaturity,
							   double aVolTenor,
							   double aFwd,
							   double aStrike)
{
	return -1;
}