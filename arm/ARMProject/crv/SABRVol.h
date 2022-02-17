/*--------------------------------------------------------------------------------------*/
/*																						*/
/* FILE        : SABRVol.h																*/
/*																						*/
/* DESCRIPTION : Header for the SABRVol class, a class									*/
/*               for dealing with options implicit volatility							*/ 
/*																						*/
/* DATE        : Thu Feb 15 2007														*/
/*																						*/
/*--------------------------------------------------------------------------------------*/

#ifndef _SABR_VOL_H
#define _SABR_VOL_H


#include "volcurv.h"

#include <stdlib.h>



class ARM_SABRVol : public ARM_VolCurve
{
	private :

		ARM_VolCurve* itsSigmaOrAlpha;

		ARM_VolCurve* itsBeta;
		ARM_VolCurve* itsRho;
		ARM_VolCurve* itsNu;

		int           itsSigmaOrAlphaFlag;	// Default: Sigma: 1; Alpha: 0

		int			  itsModelType;			// Smiled Model Type: LD, SABR_G, SABR_A,
											//					  SABR_IMPLNVOL, SABR_WEIGHT

		double		  itsWeight;			// Used only with SABR_WEIGHT case

		void Init(void);

	public :

		// Constructors
        ARM_SABRVol(void);
   
		ARM_SABRVol(ARM_VolCurve* aAlpha, 
					ARM_VolCurve* aBeta, 
					ARM_VolCurve* aRho, 
					ARM_VolCurve* aNu,
					int           SigmaOrAlphaFlag = 1,		// Default: Sigma
					int           ModelType        = 0,		// Default: K_LD
					double        Weight		   = 0.5);	// Default: no weight

		ARM_SABRVol(const ARM_SABRVol& copy);

		// Destructor
	    virtual ~ARM_SABRVol(void);

		// Copy methods
	    void Copy		(const ARM_Object* srcSABRVol);
		void BitWiseCopy(const ARM_Object* srcSABRVol);

		virtual ARM_Object* Clone();

		// Assignment operator
		ARM_SABRVol& operator = (const ARM_SABRVol& SABRVol);

	public:

		// Accessors
		ARM_VolCurve* GetSigmaOrAlpha()	{	return itsSigmaOrAlpha;		}
		ARM_VolCurve* GetBeta()			{	return itsBeta;				}
		ARM_VolCurve* GetRho()			{	return itsRho;				}
		ARM_VolCurve* GetNu()			{	return itsNu;				}

		int GetSigmaOrAlphaFlag()		{	return itsSigmaOrAlphaFlag;	}
		int GetModelType()				{	return itsModelType;		}

		double GetWeight()				{	return itsWeight;			}

		double computeVol(double aMaturity, double aVolTenor, double aFwd, double aStrike);

		// Viewer
		void View(char* id = NULL, FILE* ficout = NULL);
};

#endif


/*---------------------------------------------------------------------------------------*/
/*---- End Of File ----*/