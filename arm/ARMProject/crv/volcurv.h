/*
 $Log: volcurv.h,v $
 Revision 1.22  2004/05/07 14:22:49  jpriaudel
 stockage des MarketDatas bruts

 Revision 1.21  2004/01/13 11:50:09  mab
 Added : ComputeVol(double mat, double underlying,
                                 double moneyness)

 Revision 1.20  2003/11/20 15:24:37  ebenhamou
 added last known date

 Revision 1.19  2003/11/14 11:24:30  mab
 Added itsFxVolSmileInterpolation
 destructor displaced to .cpp

 Revision 1.18  2003/10/07 15:09:36  ebenhamou
 added const accessor for const correctness

 Revision 1.17  2003/09/24 15:48:30  ebenhamou
 remove the ambiguous dynamic_cast

 Revision 1.16  2003/09/23 09:15:25  mab
 Added: ARM_VolCurve * ProductBy(double x)

 Revision 1.15  2003/09/19 15:15:25  jpriaudel
 cumulative bump added

 Revision 1.14  2003/09/18 17:10:29  ebenhamou
 change const

 Revision 1.13  2003/09/18 17:07:02  mab
 Added: virtual ARM_VolLInterpol* ConvertToNormalVol(..)

 Revision 1.12  2003/09/09 12:56:52  mab
 Added : double GetNthMaturity(int i)

 Revision 1.11  2003/08/18 09:46:49  mab
 improvement in the principal method :
 virtual double VolatilityFunction(double m1, double K, double m2)
 virtual double VolatilityFunction(double m1, double m2) :
 No default parameter but 2 methods!

 Revision 1.10  2003/04/03 08:58:09  mab
 Correction of View prototype

 Revision 1.9  2003/02/11 14:52:46  mab
 Added : virtual void ParallelShift(double value);
 (key word virtual added)
 virtual void UpdateCol(ARM_Vector* pCol, double tenor);
 void UpdateLine(ARM_Vector* pLine, double yearterm);
 virtual ARM_Vector* GetCol(double tenor)

 Revision 1.8  2002/10/11 08:27:04  mab
 Improvements

 Revision 1.7  2002/09/17 16:33:45  mab
 Added : SetCurrencyUnit(vCurve->itsCurrency);
 For the object ARM_VolCurve

 Revision 1.6  2001/07/23 17:08:26  smysona
  ComputeVolatility a trois parametres

 Revision 1.5  2001/04/03 11:55:41  nicolasm
 Print

 Revision 1.4  2001/03/21 15:05:56  nicolasm
 *** empty log message ***

 */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : volcurv.h                                                    */
/*                                                                            */
/* DESCRIPTION : Header for the ARM_VolCurve class, a class for dealing with  */
/*               vol curves. Always define subclasses of this class.          */
/*               Subclasses should override the VolatilityFunction            */
/*               and D1DiscountFunction methods.                              */
/*                                                                            */
/* DATE        : Tue Jan  1 1997                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#ifndef _VOL_CURV_H
#define _VOL_CURV_H



#include <stdio.h>


#include "AbstractMarketClass.h"
#include "dates.h"
#include "currency.h"
#include "bsxtic.h"




class ARM_Date;

class ARM_VolLInterpol;
class ARM_ZeroCurve; 
class ARM_FXVolSmileInterpol;
class ARM_Currency;
class ARM_VolCurve; 

// #include "irindex.h"

class ARM_IRIndex; 

class ARM_VolCurve : public ARM_AbstractMarketClass
{
    private:

        ARM_Date itsLastKnownDate; 

        // the last known date is by default the asOfDate
        // only in the case of the Inflation, it is different and 
        // represents the date of the last known index

        int itsStrikeType;       // Type de Strike (Price or Delta)

        int itsVolType;          // Type de Vol K_ATMF_VOL(ATM), K_SMILE_VOL(SMILE), K_FX_VOL_SP_INTERP

        int itsOptionType;       // Type d'option : K_IRG(CAP), K_SWOPT(SWAPTION)

		int itsInterpType;       // 0 : comme avant,
		                         // 1 : on cherche l'index dans la liste des ténors et on interpole sur la colonne

        // JLA usesless .
		// ARM_INDEX_TYPE itsIndexType;   // Index Type : IRIndex type  

        ARM_Currency* itsCurrency;

        ARM_Vector* itsExpiryDates;

        ARM_Vector* itsExpiryTerms;

        ARM_Matrix* itsVolatilities; // Corresponding volatility values

        
        ARM_FXVolSmileInterpol* itsFxVolSmileInterpolation;

		std::string	itsIndexName;		// Summit Index Name
		ARM_IRIndex* itsIndex;		// optional Index, aggregated if exists. 

        virtual double VolatilityFunction(double m1, double K, double m2)
        {
            return(0.0);
        }

 public:
        virtual double VolatilityFunction(double m1, double m2)
        {
            return(0.0); 
        }



        ARM_VolCurve(const ARM_Date& asOf, ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        ARM_VolCurve(const ARM_Date& asOf, int KType, int volType = K_ATMF_VOL,
                     ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        ARM_VolCurve(const ARM_VolCurve& volCurve);

        char itsYearTermsX[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];
        char itsYearTermsY[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];

        ARM_VolCurve(void)
        {
           //  set default values of other variables

           Init();
        }

       virtual ~ARM_VolCurve(void);

        ARM_VolCurve& operator = (const ARM_VolCurve& volCurve);
  
        void Init(void)
        {
           SetName(ARM_VOL_CURVE);

           itsStrikeType   = K_STK_TYPE_PRICE;

           itsVolType      = K_ATMF_VOL;

           itsOptionType   = K_SWOPT;

           //JLA useless itsIndexType    = LIBOR3M;

           itsInterpType   = K_LINEAR;

           itsCurrency     = new ARM_Currency(ARM_DEFAULT_COUNTRY);

           itsExpiryDates  = NULL;

           itsExpiryTerms  = NULL;

           itsVolatilities = NULL; 

           itsFxVolSmileInterpolation = NULL;
		   
		   // itsIndexName = NULL;

           memset(itsYearTermsX, '\0',
                  sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);

           memset(itsYearTermsY, '\0',
                  sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);

		   itsIndex=NULL; 
        }

        ARM_FXVolSmileInterpol* GetFxVolSmileInterpolation(void)
        {
            return(itsFxVolSmileInterpolation);
        }

        void SetFxVolSmileInterpolation(ARM_FXVolSmileInterpol* fxVolSmileInt)
        {
            itsFxVolSmileInterpolation = fxVolSmileInt;
        }

        int GetInterpType(void)
        {
            return(itsInterpType);
        }

virtual void SetInterpType(int interpType)
        {
			itsInterpType = interpType;
        }

        double GetNthMaturity(int i)
        {
             if ( i < 0 )
             {
                return(itsExpiryTerms->GetSize());
             }
             else if ( i < itsExpiryTerms->GetSize() )
             {
                return(itsExpiryTerms->Elt(i));
             }
             else
             {
                throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                 "ARM_VolCurve::GetNthMaturity(int i): Invalid Index");
             }
        }

        // Function that multiplies a volcurve by x

        ARM_VolCurve* ProductBy(double x)
        {
            ARM_VolCurve* newvolcurve= (ARM_VolCurve*)( Clone() );

            ARM_Matrix* M  = newvolcurve->GetVolatilities();

            ARM_Matrix* M1 = new ARM_Matrix(M);

            (*M1) *=x;

            newvolcurve->SetVolatilities(M1);

            return(newvolcurve);
        }

virtual void ParallelShift(double value);

virtual void BumpVolatility(double value, 
                            int nthLine = 0, int nthCol = 0,
                            int isCumul = K_NO,
							int isAbsolute = K_YES);
        
        // Conversion de la vol forward Black-Scholes en vol Normale:
        // Rq: utilise ImpliedLogNorIRVol
virtual ARM_VolLInterpol* ConvertToNormalVol(ARM_ZeroCurve* YieldCurve,
                                             int IsSwaptionVol = 1,
                                             int InPct = 1)
        {
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "ARM_VolCurve::ConvertToNormalVol: Should never be called");
        }

        double ComputeFxVol(ARM_Date& AsOf,
                            ARM_Date& matuDate,
                            double calcMatu,
                            double fxSpot, // The FX spot, not the Fwd!
                            double strike,
                            ARM_ZeroCurve* discCrv, // JPY
                            ARM_ZeroCurve* divCrv); // USD

        void BitwiseCopy(const ARM_Object* srcVolCurve);
 
        void Copy(const ARM_Object* vCurve)
        {
            ARM_AbstractMarketClass::Copy(vCurve);

            BitwiseCopy(vCurve);
        }

        virtual ARM_Object* Clone(void)
        {
            ARM_VolCurve* theClone = new ARM_VolCurve();


            theClone->Copy(this);

            return(theClone);
        }

        ARM_CLASS_NAME GetRootName(void)
        {
           return(ARM_VOL_CURVE);
        }

        void View(char* id = NULL, FILE* ficOut = NULL);


        double ComputeVolatility(double mat, double strike); // returns volatility 

        double ComputeVolatility(ARM_Date& MatDate , 
                                 double strike);  // returns volatility 

		// strike utilisé uniquement dans l'hypercube
        virtual double ComputeVolatility(double mat, double moneyness, double underlying, double strike = 0.0);

		virtual double computeVol(double moneyness, double maturity) { return 0; };

        double ComputeVolatility(ARM_Date& date, double moneyness, double underlying);

		virtual double ComputeCorrelByExpiry(double aExpiry, double aTenor1, double aTenor2)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							"ARM_VolCurve::ComputeCorrelByExpiry: Should only be called for a VolCurve");
        }

		virtual double ComputeHyperCorrel( double aTenor1, double aTenor2, double aExpiry, double strike = 0.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							"ARM_VolCurve::ComputeHyperCorrel: Should only be called for a HyperCube");
        }

        // new methods more suitable
        double ComputeVol(double mat, double underlying,
                          double moneyness = 0.0);

        double ComputeVol(ARM_Date& date, double underlying,
                          double moneyness = 0.0);


        void SetCurrencyUnit(ARM_Currency* ccy);

        // when we set the asOfDate
        // we also initialise the last known date 
        // to the same value
        void SetAsOfDate(ARM_Date& asOf) 
        {
            itsAsOfDate = asOf;
            itsLastKnownDate = asOf;
        }

        /// Beware that the SetLastKnownDate should always
        /// ALWAYS be set up after SetAsOfDate since the latter
        /// gives by default the same value to asOf and lastKnownDate!
        void SetLastKnownDate(const ARM_Date& lastKnownDate) 
        {
           itsLastKnownDate = lastKnownDate;
        }
      
        inline ARM_Date GetLastKnownDate() const
        {
            return(itsLastKnownDate);
        }
      
        int GetStrikeType(void) const
        {
            return(itsStrikeType);
        }

        void SetStrikeType(int strikeType)
        {
            itsStrikeType = strikeType;
        }

        int GetVolType(void) const
        {
            return(itsVolType);
        }

        void SetVolType(int volType)
        {
            itsVolType = volType;
        }

        int GetOptionType(void) const
        {
           return(itsOptionType);
        }

        void SetOptionType(int OptionType)
        {
            itsOptionType = OptionType;
        }

        ARM_Currency* GetCurrency(void)
        {
            return(itsCurrency);
        }

        /// const version for const correctness!
        ARM_Currency* GetCurrency(void) const
        {
            return itsCurrency;
        }

        void SetCurrency(ARM_Currency* ccy)
        {
            itsCurrency = ccy;

            SetStrCurrency(itsCurrency->GetCcyName());
        }

        inline virtual ARM_Matrix* GetVolatilities() const
        {
           return(itsVolatilities);
        }

        void SetVolatilities(ARM_Matrix* volatilities)
        {
            if ( itsVolatilities == volatilities )
               return;

            if (itsVolatilities)
            {
               delete itsVolatilities;
               itsVolatilities= NULL;
            }

            if (volatilities)
               itsVolatilities = volatilities;
        }


virtual void UpdateCol(ARM_Vector* pCol, double tenor);
        void UpdateLine(ARM_Vector* pLine, double yearterm);

        virtual ARM_Vector* GetCol(double tenor)
        {
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "ARM_VolCurve::GetCol should never be called");

            return(NULL);
        }

        virtual ARM_Vector* GetExpiryTerms(void) const
        {
           return(itsExpiryTerms);
        }

        void SetExpiryTerms(ARM_Vector* terms)
        {
           if ( itsExpiryTerms == terms )
              return;

           if (itsExpiryTerms)
           {
              delete itsExpiryTerms;
              itsExpiryTerms= NULL;
           }

           if (terms)
              itsExpiryTerms = terms;
        }

		virtual ARM_Vector* GetStrikes() const
        {
            return NULL;
        }
		ARM_Vector* GetExpiryDates(void) const
        {
           return(itsExpiryDates);
        }

        void SetExpiryDates(ARM_Vector* dates)
        {
           if ( itsExpiryDates == dates )
              return;

           if (itsExpiryDates)
           {
              delete itsExpiryDates;
              itsExpiryDates= NULL;
           }

           if (dates)
              itsExpiryDates = dates;
        }

		inline void SetIndexName(const std::string& name) 
		{
			itsIndexName=name; 
		}

		inline const std::string& GetIndexName(void) const 
        {
            return itsIndexName;
        }
		//	once set, index can't be cleared.
		void SetIndex(const ARM_IRIndex&ref) ;
		const ARM_IRIndex& GetIndex() const ;
		bool hasIndex() const { return itsIndex!=NULL; } 
		// 
		bool IsEqualToOne(void); // at 1e-6

		virtual ARM_VolLInterpol* ComputeATMFxVol(double lag=0.0) {return NULL;};
		virtual int GetSmileFlag(void) const  { return true;   }
		virtual void SetATMVolFX(ARM_VolLInterpol*) {};
};



                    /*------------------------------*/
                    /*                              */
                    /*         FX Vol Curve         */
                    /*                              */
                    /*------------------------------*/







class ARM_FXVolCurve : public ARM_VolCurve
{
    private:

        // Inputs


		ARM_Currency*  itsDomCcy;
        ARM_Currency*  itsForCcy;

        ARM_ZeroCurve* itsDomCurve;
        ARM_ZeroCurve* itsForCurve;

        ARM_Vector itsOptionsMatus;
        ARM_Vector itsPivotVols;
        ARM_Vector itsPivotTypes;
        ARM_Vector itsInterpolTypes;

        ARM_Vector itsDeltasCall;
        ARM_Matrix itsVolsCall;

        ARM_Vector itsDeltasPut;
        ARM_Matrix itsVolsPut;

        ARM_Vector itsFxFwds;

        ARM_Matrix itsRR;
        ARM_Matrix itsSTR;

		int nbTime2Maturities ;
		int what_is_interpolated; // FXINTERP_STRIKE or FXINTERP_DELTA or FXINTER_DELTA_SMOOTH
		
        int forceLinearWhenSplineFails;

        int itsRRSTRInputFlag;
        
		double itsSpotFX;

		typedef struct 
        {
			double time2maturity;
			double pivotVol;
			int pivotType;
			double pivotStrike;
			double fxFwd;
			ARM_Vector volCall;
			ARM_Vector deltaCall; // would allow to have different delta for diffent time2maturity
			ARM_Vector volPut;		
			ARM_Vector deltaPut;  // idem
            ARM_Vector strikeCall;
            ARM_Vector strikePut;

			int interpol_type;    // would allow to have different interpolation for diffent time2maturity K_LINEAR or K_SPLINE
			ARM_Vector strikeAll; // ordered strikes for the delta inputs
			ARM_Vector volAll;	  // vols

		} 
        DATA;
		
        DATA* itsData;

        // Smile flag setted in the constructor

        int itsSmileFlag;

        ARM_VolLInterpol* itsCalcFXATM;

    public:


        ARM_FXVolCurve(void)
        {
           //  set default values of other variables

           Init();
        }


        ARM_FXVolCurve(ARM_Date& Asof,
                       const ARM_Vector& time2maturities, 
					   const ARM_Vector& pivotVols, 
					   const ARM_Vector& pivotTypes,
					   const ARM_Vector& deltasCall,
					   ARM_Matrix& volsCall,
					   const ARM_Vector& deltasPut,
					   ARM_Matrix& volsPut,
					   const ARM_Vector& inFxFwds,
					   ARM_Vector& interpolTypesVect,
					   int whatIsInterpolated,
					   int correctSplineWithLinear,
                       double spotFX = 0.0,
                       ARM_ZeroCurve* domCurve = NULL,
                       ARM_ZeroCurve* ForCurve = NULL,
                       int inRRSTR = 0,
                       int isATM = 0);

       ~ARM_FXVolCurve(void);

	   	/// Accessors
        virtual int GetSmileFlag(void) const  { return itsSmileFlag;   }
		inline ARM_Vector GetOptionsMatus() const { return itsOptionsMatus; }
		
        inline ARM_Vector GetDeltasCall()	const { return itsDeltasCall;	}
		inline ARM_Vector GetDeltasPut ()	const { return itsDeltasPut;	}
		
        inline const ARM_Vector* GetPivotVols() const { return &itsPivotVols; };
	
        inline ARM_Matrix GetRR() const { return itsRR; }
		inline ARM_Matrix GetSTR() const { return itsSTR; }


        // Bump ATM VOL (pivots bump)
        void BumpFxVol(double shiftValue,
                       int firstLine,
                       int lastLine,
                       double spotFX = 0.0,
                       ARM_ZeroCurve* domCurve = NULL,
                       ARM_ZeroCurve* ForCurve = NULL);


		void BumpFxVol(ARM_Vector newVol,
                       ARM_ZeroCurve* domCurve = NULL,
                       ARM_ZeroCurve* ForCurve = NULL);



        // This method is just artificial for now in order to use the Excel interface
        // for testing
        void BumpVolatility(double value, 
                            int nthLine = 0, int nthCol = 0,
                            int isCumul = K_NO,
							int isAbsolute = K_YES)
        {
            BumpFxVol(value, nthLine, nthCol);
        }

        // Bump Risk reversal and strangles
        void BumpRRorSTR(const vector<int>& matusX,
                         const vector<int>& deltasY,
                         const vector<double>& shiftValues,
                         double spotFX = 0.0,
                         ARM_ZeroCurve* domCurve = NULL,
                         ARM_ZeroCurve* ForCurve = NULL);

		void FXBumpRRorSTR(	double shiftValue, 
							int nbRow, 
							int nbCol, 
							double spotFX,
							int isRR,
							int isCumul,
							int isAbsolute);
   
		void FXBumpRRorSTR(	ARM_Matrix mR, 	int isRR=K_YES);

        void GenerateFXVols(ARM_Date& Asof,
                            const ARM_Vector& time2maturities, 
					        const ARM_Vector& pivotVols, 
					        const ARM_Vector& pivotTypes,
					        const ARM_Vector& deltasCall,
					        ARM_Matrix& volsCall,
					        const ARM_Vector& deltasPut,
					        ARM_Matrix& volsPut,
					        const ARM_Vector& inFxFwds,
					        ARM_Vector& interpolTypesVect,
					        int whatIsInterpolated,
					        int correctSplineWithLinear,
                            double spotFX = 0.0,
                            ARM_ZeroCurve* domCurve = NULL,
                            ARM_ZeroCurve* ForCurve = NULL,
                            int inRRSTR = 0,
                            int isATM = 0);

		void GenerateFXVols(ARM_Matrix& volsCall,
							ARM_Matrix& volsPut,
					        double spotFX = 0.0,
                            int inRRSTR = 1);

        void View(char* id = NULL, FILE* ficOut = NULL);
        void BitwiseCopy(const ARM_Object* srcVolCurve);
 
        void Copy(const ARM_Object* vCurve)
        {
           ARM_VolCurve::Copy(vCurve);

           BitwiseCopy(vCurve);
        }

        ARM_Object* Clone(void)
        {
            ARM_FXVolCurve* theClone = new ARM_FXVolCurve();
            
            theClone->Copy(this);
            
            return(theClone);
        }
        
        void Init(void)
        {
            SetName(ARM_FX_VOLAT);

            itsDomCcy   = NULL;
            itsForCcy   = NULL;

            itsDomCurve = NULL;
            itsForCurve = NULL;

            itsRRSTRInputFlag = 0;
			itsData           = NULL;
			nbTime2Maturities = 0;
            itsSmileFlag      = 1;
            itsCalcFXATM      = NULL;
        }

        virtual ARM_VolLInterpol* ComputeATMFxVol(double lag=0.0);
		virtual void SetATMVolFX(ARM_VolLInterpol* vol) { itsCalcFXATM = vol;};

		void set_a_DATA(DATA& aData, 
						const ARM_Vector& deltasCall, 
						const ARM_Vector& volsCall, 
						const ARM_Vector& deltasPut, 
						const ARM_Vector& volsPut,
						double pivotStrike,
						double pivotVol,
						double pivotDeltaCall,
						double pivotDeltaPut,
						int interpolType);

		double computeVol_givenMaturity(double K, int iMaturity);

		virtual double computeVol(double moneyness, double maturity);

        virtual double VolatilityFunction(double m1, double m2);

        // here m2 is not relevant
        virtual double VolatilityFunction(double m1, double k, double m2);
		
        void get_relevant_maturities(double maturity, int& low, int& high);

        void set_strikes(DATA& aData);

		ARM_Matrix get_info(int infoVol = 0);

		ARM_Currency* GetDomCcy(void)	{	return itsDomCcy;	}
        ARM_Currency* GetForCcy(void)	{	return itsForCcy;	}

		
		void SetDomCcy(string domCcy) 
        {
            if (itsDomCcy)
            {
               delete itsDomCcy;

               itsDomCcy = NULL;
            }

            itsDomCcy = new ARM_Currency(domCcy.c_str()); 
        }
		
        void SetForCcy(string forCcy) 
        { 
            if (itsForCcy)
            {
               delete itsForCcy;

               itsForCcy = NULL;
            }

            itsForCcy = new ARM_Currency(forCcy.c_str()); 
        }
};





#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/