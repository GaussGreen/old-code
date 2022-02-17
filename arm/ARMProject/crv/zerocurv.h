/*
 * $Log: zerocurv.h,v $
 * Revision 1.25  2004/05/11 10:45:05  jpriaudel
 * virtual added in discountyield and fwdyield
 *
 * Revision 1.24  2004/05/07 14:09:53  jpriaudel
 * modif rene
 *
 * Revision 1.23  2004/04/22 13:27:27  jpriaudel
 * MINDOUBLE added
 *
 * Revision 1.22  2003/12/15 17:55:12  arm
 * Syntax correction
 *
 * Revision 1.21  2003/12/15 14:59:43  arm
 * SetMktData now in .h : because of strange compilation PB
 *
 * Revision 1.20  2003/11/13 13:31:00  jpriaudel
 * modif pour le credit
 *
 * Revision 1.19  2003/10/22 16:39:13  jpriaudel
 * modif pour Olivier
 *
 * Revision 1.18  2003/09/03 18:03:59  jpriaudel
 * ajout de SetMktData
 *
 * Revision 1.17  2003/07/01 19:56:15  jpriaudel
 * suppression of include armglob
 *
 * Revision 1.16  2003/02/11 14:57:00  mab
 * Added : double ForwardYieldAdj(ARM_Date& fwdDate,
 * >                            ARM_Date& zeroTerm,
 * >                            int compMeth = 0)
 *
 * Revision 1.15  2002/10/11 08:22:28  mab
 * Improvements
 *
 * Revision 1.14  2002/09/24 13:15:44  mab
 * Added : TOY management (See TOY section)
 *
 * Revision 1.13  2002/09/17 16:38:47  mab
 * in : CptFuturesZeroRates parameter int SizeS added
 *
 * Revision 1.12  2002/08/08 10:20:41  mab
 * Added const in declaration : ARM_ZeroCurve(const ARM_ZeroCurve& zeroCurve);
 * ARM_ZeroCurve& operator = (const ARM_ZeroCurve& zeroCurve);
 *
 * Revision 1.11  2002/05/30 13:39:22  mab
 * in char Terms[ARM_NB_TERMS][6] : 6 replaced by 12
 *
 * Revision 1.10  2001/04/23 09:23:20  smysona
 * inlining
 *
 * Revision 1.9  2001/04/03 12:00:01  nicolasm
 * Modif d'un set, inlining des accesseurs
 *
 * Revision 1.8  2001/03/12 19:21:27  smysona
 * Ajout d'un set
 *
 * Revision 1.7  2000/12/11 16:11:16  mab
 * Rajout de : extern int GOTO_END_OF_MONTH;
 *
 * Revision 1.6  2000/10/25 10:15:07  smysona
 * Rajout d'un test sur itsParameters dans le Print
 *
 * Revision 1.5  2000/05/02 10:03:17  mab
 * Rajout de comments. RCS et passage par reference de
 * ARM_Matrix ds les generations de courbes
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : zerocurv.h                                                   */
/*                                                                            */
/* DESCRIPTION : Header for the ARM_ZeroCurve class, a class for dealing with */
/*               zero curves. Always define subclasses of this class.         */
/*               Subclasses should override the DiscountFunction              */
/*               and D1DiscountFunction methods.                              */
/*                                                                            */
/* DATE        : Thu Aug  1 1996                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#ifndef _ZERO_CURV_H
#define _ZERO_CURV_H



#include <stdio.h>
#include <string>

#include "AbstractMarketClass.h"
#include "bkshift.h"
#include "dates.h"
#include "matlib.h"


extern int GOTO_END_OF_MONTH;

#define ANNUA            1
#define CONTINUOUS        0
#define MONEY            -1

#define MINDOUBLE        -10000.0
 
class ARM_Date;

class ARM_Portfolio;

class ARM_YCModel;
class ARM_Container;
class ARM_Currency;



double EvalEcart(ARM_Matrix* , ARM_Vector*  , void **);
            

class ARM_MarketData : public ARM_Object 
{
    public :

        ARM_CRV_TERMS itsMktTerms;

        ARM_Vector* itsMktValue;

        int itsMMVsFut;
        
        int itsSwapVsFut;
        
        int itsraw;
        
        int itsCont_Lin ;
        
        int itsConstructionMeth ;


        ARM_MarketData()
        {
            Init();
        }


        ARM_MarketData(ARM_CRV_TERMS& MktTerms, 
                       ARM_Vector* MktData, int MMVsFut, 
                       int SwapVsFut, int raw, 
                       int Cont_Lin, int ConstructionMeth)
        {
            Init();

            memcpy(itsMktTerms, MktTerms, 
                   sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);

            itsMktValue = (ARM_Vector*)MktData->Clone();
            itsMMVsFut = MMVsFut;
            itsSwapVsFut = SwapVsFut;
            itsraw = raw;
            itsCont_Lin = Cont_Lin;
            itsConstructionMeth = ConstructionMeth;
        }    

       ~ARM_MarketData()
        {
            if (itsMktValue)
                delete itsMktValue;

            itsMktValue = NULL;
        }
        
        void Init(void)
        {
            SetName(ARM_MARKET_DATA);

            memset(itsMktTerms, '\0',
                   sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);

            itsMktValue = NULL;
            itsMMVsFut = K_FUT;
            itsSwapVsFut = K_FUT;
            itsraw = K_PAR;
            itsCont_Lin = K_LINEAR;
            itsConstructionMeth = -1;
        }

        void BitwiseCopy(const ARM_Object* srcMktData)
        {
            ARM_MarketData* mkt = (ARM_MarketData *) srcMktData;

            memcpy(itsMktTerms, mkt->itsMktTerms, 
                   sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);

            if (itsMktValue)
            {
               delete itsMktValue;
            }

            if (mkt->itsMktValue)
               itsMktValue = (ARM_Vector *) mkt->itsMktValue->Clone();
            else
               itsMktValue = NULL;

            itsMMVsFut = mkt->itsMMVsFut;

            itsSwapVsFut = mkt->itsSwapVsFut;
        
            itsraw = mkt->itsraw;
        
            itsCont_Lin = mkt->itsCont_Lin;
        
            itsConstructionMeth = mkt->itsConstructionMeth;
        }
        
        void Copy(const ARM_Object* srcmktData)
        {
            ARM_Object::Copy(srcmktData);
 
            BitwiseCopy(srcmktData);
        }
 
        ARM_Object* Clone(void)
        {
            ARM_MarketData* theClone = new ARM_MarketData();
 
 
            theClone->Copy(this);
 
            return(theClone);
        }
};


class ARM_ZeroCurve : public ARM_AbstractMarketClass
{
    private:
       
		string	itsSmoothingMethod;

        ARM_Vector* itsDateTerms;
        
        ARM_Vector* itsYearTerms;

        ARM_Vector* itsBPShifts;
        
        ARM_Vector* itsDiscountFactors;

        ARM_Vector* itsZeroRates;

        // Added FOR TOY

        ARM_Vector* itsForwardDates;

        ARM_Vector* itsForwardRates;

        ARM_Vector* itsForwardDF;


        // END TOY

        ARM_Vector* itsParameters;

        ARM_Currency* itsCurrency;

        double itsBPShift;           // shift of the yield curve in bp
        
        double itsBucketStartPeriod; // start period (in yearTerm) of the bucket 

        double itsBucketEndPeriod;   // end period (in yearTerm) of the bucket 

        double itsDFSpotFwd;


        // Storage of markets data parameters
        
        ARM_MarketData* itsMktData;

		int itsFixDayCount;   // To handle generating different FixedIndex Daycount


        // The following methods must be overriden by subclasses
	public:
virtual double DiscountFunction(double);  // returns discount price
	
	private:
virtual double D1DiscountFunction(double);//1st derivative of discount price

virtual double D2DiscountFunction(double);//2nd derivative of discount price
 

        void CptMMZeroRates(ARM_Matrix* MMRates, int SizeM);

        void CptMMZeroRatesFWD(ARM_Matrix* MMRates, int SizeM, int Cont_Lin);

        void CptFuturesZeroRates(ARM_Matrix& MMRates, int SizeM,
                                 ARM_Matrix* Future, int SizeF,
                                 int SizeS,
                                 ARM_Matrix& Swap, 
                                 int MonVsFut, int Cont_Lin);

		void CptAllForwardsZeroRates(ARM_Matrix& MMRates, int SizeM,
                                 ARM_Matrix* Future, int SizeF,
                                 int SizeS,
                                 ARM_Matrix& Swap, 
                                 int MonVsFut, int Cont_Lin);

        void CptFuturesZeroRatesGen(ARM_Matrix& MMRates, int SizeM,
                                    ARM_Matrix& Future, int SizeF,
                                    ARM_Matrix& Swap,
                                    int MonVsFut, int Cont_Lin);

        void CptSwapZeroRates(ARM_Matrix& MMRates, int SizeM, 
                              ARM_Matrix* Swap, 
                              int& SizeS, int fxPayFreq,
                              ARM_Matrix& Future, int SizeF, int MMVsFut, 
                              int SwapVsFut, int raw,
                              int Cont_Lin,
							  int& firstSwapIndex,
							  int indexAUD=0,
							  int fxDayCount=KNOBASE);

		void CptSwapZeroRates(ARM_Matrix& MMRates, int SizeM, 
                              ARM_Matrix* Swap, 
                              int& SizeS, int fxPayFreq,
                              ARM_Matrix& Future, int SizeF, int MMVsFut, 
                              int SwapVsFut, int raw,
                              int Cont_Lin,
							  int indexAUD=0, 
							  int fxDayCount=KNOBASE);

        void CptBondZeroRates(ARM_Matrix& MMRates, 
                              ARM_Container* bonds, 
                              ARM_Matrix* bondPrices);

		double InterpDfSwapMilieu(double dateIn, 
								  const ARM_Matrix& MMRates, int SizeM, 
								  const ARM_Matrix& Future, int SizeF, 
								  const ARM_Matrix& Swap, int SizeS,
								  int Cont_Lin);

		double CalculSwapFlowNoFreq(double dateJul, double txSwap, double startDF,
							        const ARM_Matrix& MMRates, int SizeM, 
							        const ARM_Matrix& Future, int SizeF, 
							        const ARM_Matrix& Swap, int SizeS,
							        int fxPayFreq,
							        int Cont_Lin,
							        int fxDayCount,
							        int spotdays,
							        char* ccyName);

        /*---------------- TOY METHODS ------------------*/

        int CptDF(ARM_Matrix& MMRates,
                  int SizeM,
                  ARM_Matrix& FRARates,
                  int SizeFRA,
                  int Cont_Lin);

        double InterpolForward(double date,ARM_Matrix& FRARates,
                               int SizeFRA,
                               ARM_Matrix& MMRates,
                               int SizeM,
                               int Cont_Lin);

        double InterpolDF(double date,ARM_Matrix& FRARates,
                          int SizeFRA,
                          ARM_Matrix& MMRates, 
                          int SizeM,
                          int Cont_Lin);

        int CptTOYFRA(ARM_Matrix& FRARates, int SizeFRA,
                      ARM_Matrix* Futures, int SizeF,
                      int Cont_Lin);

        int CptTOY(ARM_Matrix& MMRates,
                   int SizeM,
                   ARM_Matrix* Futures,
                   int SizeF,
                   int Cont_Lin);

        void CptSwapZeroRatesTOY(ARM_Matrix& MMRates, 
                                  int SizeM,
                                  ARM_Matrix& FRARates,
                                  int SizeFRA,
                                  ARM_Matrix* Swap, 
                                  int& SizeS, 
                                  int fxPayFreq,
                                  ARM_Matrix& Future, 
                                  int SizeF, int MMVsFut,
                                  int SwapVsFut, int raw,
                                  int Cont_Lin);

        void CptFuturesZeroRatesTOY(ARM_Matrix& MMRates, int SizeM,
                                    ARM_Matrix& FRARates, int SizeFRA,
                                    ARM_Matrix* Future, int SizeF,
                                    ARM_Matrix Swap,
                                    int MonVsFut, int Cont_Lin);

         /*------------- END OF TOY METHODS --------------*/ 


		/*-------------   FORWARD METHODS  --------------*/


		void FillFutSwapDisc(ARM_Matrix* LocDiscData, int& LocDiscSize,
							 double FixRate, double FixDaycount, double FixedFreq,
							 double StartFwd, double ValSwap, double FwdDaycount, double FwdFreq,
							 ARM_Date& DtAsofFwd,ARM_Date& DtRollFwd, double DiscStartFwd,
							 ARM_Date& DtThFirstSw, ARM_Date& DtThEndSw, double DiscStartSw,
							 char* calName, int Cont_Lin, int FwdAdj);

		double CptSwapFwdValue(ARM_Matrix *LocDiscData, int& LocDiscSize,
							   double FixRate, double FixDaycount, double FixedFreq,
							   double StartFwd, double PteFwd, double FwdDaycount, double FwdFreq,
							   ARM_Date &DtAsofFwd,ARM_Date &DtRollFwd, double DiscStartFwd,
							   ARM_Date &DtThFirstSw, ARM_Date &DtThEndSw, double DiscStartSw,
							   char* calName, int Cont_Lin, int FwdAdj);

		 
		void CptSwapZeroRatesFWD(ARM_Matrix& MMRates, int SizeM,
                                 ARM_Matrix* Swap, 
                                 int& SizeS, int fxPayFreq,
                                 ARM_Matrix& Future, int SizeF, 
                                 int MMVsFut,
                                 int SwapVsFut,
                                 int raw,
                                 int Cont_Lin,
								 int& firstIndexSwap,
								 int indexAUD,
								 int FwdMatu=3,
								 int FwdAdj=K_UNADJUSTED,
								 int fxDayCount=KNOBASE);


		
		double LocDiscInterp(ARM_Matrix *TabDiscData, int TabDiscSize,
										ARM_Date &AsofDate, ARM_Date &CurDate,
										int Cont_Lin, int IntFlatLeft = 0, double startDisc = 1.);

		double LocSensi(ARM_Matrix *TabDiscData, int TabDiscSize,
							ARM_Date &AsofDate, ARM_Date &DtThStartDate, ARM_Date &DtThEndDate,
							char* FreqUnit, int FreqQty, int FixDayCount,
							int Cont_Lin, char* calName);

		double LocFwdInterp(ARM_Matrix *TabDiscData, int TabDiscSize,
									ARM_Date &AsofDate, ARM_Date &DtStartFwd, ARM_Date &DtEndFwd,
									int Cont_Lin, int FwdDayCount);


		/*-----	Bump Forwards Methods ( AR - 10/10/04) --*/

		/* moved to public declaration

		ARM_ZeroCurve* GenerateShiftCurveFwd(
				char Term[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS] , ARM_Vector* epsilon);

		*/

		void DecodePlot ( char* Terms, ARM_Date &asofDate, ARM_Date &spotDate, char* ccyName, char &matu, ARM_Date& matDate, ARM_Date& startDate);

		double FillDiscFValToSpr(ARM_Matrix* LocDiscData, int& LocDiscSize, double ValSwap,
									  double FixRate, double FixDaycount, double FixedFreq,
									  double FwdDaycount, double FwdFreq,
									  ARM_Date &DtAsofFwd,ARM_Date &DtRollFwd, double DiscStartFwd,
									  ARM_Date &DtThFirstSw, ARM_Date &DtThEndSw, double DiscStartSw,
									  char* calName, int Cont_Lin, int FwdAdj);

		double spread_MM(double fwd, ARM_Date &asofDate, ARM_Date &spotDate, ARM_Date &startDate, ARM_Date &endDate,
						double TxEnd, double DiscStart, double DiscSpot, int FwdDaycount);

		double spread_MFUT(double fwd, ARM_Date &asofDate, ARM_Date &prevDate, ARM_Date &startFut, ARM_Date &endFut,
						double TxFut, double DiscPrev, double DiscStart, int FwdDaycount);

		
		double DiscExtrapol_MFUT(ARM_Date &asofDate, ARM_Date &prevDate, ARM_Date &startFut, ARM_Date &endFut,
				 double TxFut, double DiscPrev, int FwdDaycount);

		int InsereDisc(ARM_Date &curDate, double curDisc, ARM_Matrix* Tabdisc, int curSize, int priority=K_YES);

		double CptSwapSprFwdValue(ARM_Matrix *LocDiscData, int& LocDiscSize,
									  double FixRate, double FixDaycount, double FixedFreq,
									  ARM_Matrix TabFwd, double SprFwd, double FwdDaycount, double FwdFreq,
									  ARM_Date &DtAsofFwd,ARM_Date &DtRollFwd, double DiscStartFwd,
									  ARM_Date &DtThFirstSw, ARM_Date &DtThEndSw, double DiscStartSw,
									  char* calName, int Cont_Lin, int FwdAdj);
        
		
		
		/*---------	End of Bump Forwards Methods --------*/	



		/*------------   END FORWARD METHODS  ------------*/


    protected:

        void GenerateYearTermsAndRates(void);
        void GenerateDateTerms(void);
        void GenerateDiscountFactors(int compMeth = 0);

        void GenerateFields(void)
        {
            GenerateYearTermsAndRates();
            GenerateDateTerms();
            GenerateDiscountFactors();
        }

    public:

        ARM_ZeroCurve(ARM_Date& asOf, ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);
        
        ARM_ZeroCurve(const ARM_ZeroCurve& zeroCurve);

        ARM_ZeroCurve(ARM_Date& asOf, 
                      ARM_CRV_TERMS& terms, 
                      ARM_Vector* mktData, 
                      int MonVsFut, 
                      int SwaVsFut, 
                      int raw, 
                      int Cont_Lin,
                      ARM_Currency* ccy = ARM_DEFAULT_CURRENCY,
					  int swapFrqId = K_DEF_FREQ,
					  int fixDayCount = KNOBASE);

        ARM_ZeroCurve(ARM_CRV_TERMS& terms,
                      ARM_Date& asOf,
                      ARM_Vector* mktData,
                      int MonVsFut,
                      int SwaVsFut,
                      int raw,
                      int Cont_Lin,
                      ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        ARM_ZeroCurve(ARM_Date& asOf, 
                      ARM_CRV_TERMS& terms, 
                      ARM_Vector* mktData, 
                      ARM_Container* bonds, 
                      ARM_Vector* yields,
                      int MonVsFut, 
                      ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        ARM_ZeroCurve(ARM_Date& asOf,
                      char* Terms[ARM_NB_TERMS],
                      ARM_Vector* mktData,
                      double mean_rates,
                      int raw,
                      int Cont_Lin,
                      ARM_Currency* ccy);

        ARM_ZeroCurve(void)
        {
            //  set default values of other variables
            Init();
			
			SetName(ARM_ZERO_CURVE);
        }

		string GetSmoothingMethod(void) const
		{
			return	itsSmoothingMethod;
		}

		void SetSmoothingMethod(string aSmoothingMethod)
		{
			itsSmoothingMethod = aSmoothingMethod;
		}

        ARM_MarketData* GetMktData(void)
        {
            return itsMktData;
        }

        void SetMktData(ARM_MarketData* data)
        {
            itsMktData = data;
        }

		virtual double CalcNumericalObjectSignature(void);
        
		virtual ARM_ZeroCurve* GenerateShiftCurve(ARM_CRV_TERMS& Term, 
                                                  ARM_Vector* epsilon) 
        { 
            return(NULL);
        }

		
        ARM_ZeroCurve* GenerateShiftCurveFwd(
                      char Term[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS], 
                      ARM_Vector* epsilon,
					  int FreqFwd = 3,
					  int AdjFwd = K_UNADJUSTED); // Fit to Forward Stripping convention


        /*---------------- TOY METHODS ------------------*/

        ARM_ZeroCurve(ARM_Date& asOf,
                      ARM_CRV_TERMS& Terms,
                      ARM_Vector *mktData,
                      int MMVsFut,
                      int SwapVsFut,
                      int raw,
                      int Cont_Lin,
                      int Frequency,
                      ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);


        void ZCFromMarketRatesTOY(ARM_CRV_TERMS& Terms,
                                  ARM_Vector* data,
                                  int MMVsFut = K_FUT,
                                  int SwapVsFut = K_FUT,
                                  int raw = K_PAR,
                                  int Cont_Lin = K_LINEAR,
                                  ARM_Currency* ccy = ARM_DEFAULT_CURRENCY,
                                  int Frequency = 0);

        /*------------- END OF TOY METHODS --------------*/


        ARM_ZeroCurve& operator = (const ARM_ZeroCurve& zeroCurve);

        
        /*
         *   Fit to MM, Futures and Swap rates
         *   MonVsFut = 0 -> Privilegier les Futures   
         *   SwaVsFut = 0 -> Privilegier les Futures   
         *   raw = 0 -> Par (interpolation lineaire des taux de swap 
         *   manquant cf Summit) 
         *   Cont_Lin = 0 -> Cont
         *   Cont_Lin = 1 -> Lin
         *   Cont_Lin = 2 -> Spl cub
         */
        
        void ZCFromMarketRates(ARM_CRV_TERMS& Terms,
                               ARM_Vector* data, 
                               int MonVsFut=K_FUT,
                               int SwaVsFut=K_FUT,  
                               int raw = K_PAR,  
                               int Cont_Lin= K_LINEAR,  
                               ARM_Currency* ccy = ARM_DEFAULT_CURRENCY,
							   int swapFrqId = K_DEF_FREQ,
							   int fixDayCount = KNOBASE);

        void ZCFromMarketRatesGen(ARM_CRV_TERMS& Terms,
                                  ARM_Vector* data,
                                  int MonVsFut=K_FUT,
                                  int SwaVsFut=K_FUT,
                                  int raw = K_PAR,
                                  int Cont_Lin= K_LINEAR,
                                  ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);
        
        void ZCFromMarketRates(ARM_CRV_TERMS& Terms, 
                               ARM_Vector* data, 
                               ARM_Container* bonds, 
                               ARM_Vector* yields,
                               int MMVsFut, 
                               ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);
 
        
        void ZCFromMarketRates(char* Terms[ARM_NB_TERMS], 
                               ARM_Vector* data,
                               double mean_rates = 0.0,
                               int raw = K_PAR,
                               int Cont_Lin = K_LINEAR,
                               ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        


        // Returns the standard error of the OLS
        double FitToMarketPrices(ARM_Container* Assets, ARM_Matrix* data,
                                 int *CST=NULL, double *OP=NULL);

        void ParallelShift(double value);

        void BitwiseCopy(const ARM_Object* srcZc)
        {
            ARM_ZeroCurve* zc = (ARM_ZeroCurve *) srcZc;


            if (itsDateTerms)
            {
               delete itsDateTerms;
            }
 
            if (zc->itsDateTerms)
               itsDateTerms = (ARM_Vector *) zc->itsDateTerms->Clone();
            else
               itsDateTerms = NULL;

            if (itsYearTerms)
            {
               delete itsYearTerms;
            }
 
            if (zc->itsYearTerms)
               itsYearTerms = (ARM_Vector *) zc->itsYearTerms->Clone();
            else
               itsYearTerms = NULL;

            if (itsBPShifts)
               delete itsBPShifts;

            if (zc->itsBPShifts)
               itsBPShifts = (ARM_Vector *) zc->itsBPShifts->Clone();
            else
               itsBPShifts = NULL;


            if (itsParameters)
            {
               delete itsParameters;
            }
 
            if (zc->itsParameters)
               itsParameters = (ARM_Vector *) zc->itsParameters->Clone();
            else
                itsParameters = NULL;

            if (itsZeroRates)
            {
               delete itsZeroRates;
            }
 
            if (zc->itsZeroRates)
               itsZeroRates = (ARM_Vector *) zc->itsZeroRates->Clone();
            else
               itsZeroRates = NULL;

            if (itsDiscountFactors)
            {
               delete itsDiscountFactors;
            }
 
            if (zc->itsDiscountFactors)
               itsDiscountFactors = (ARM_Vector *) 
                        zc->itsDiscountFactors->Clone();
            else
               itsDiscountFactors = NULL;

            if (zc->itsCurrency)
                this->SetCurrencyUnit(zc->itsCurrency);
            
            if (itsMktData)
            {
               delete itsMktData;
            }
 
            if (zc->itsMktData)
               itsMktData = (ARM_MarketData *) zc->itsMktData->Clone();
            else
               itsMktData = NULL;

            itsBPShift = zc->itsBPShift;

            itsBucketStartPeriod = zc->itsBucketStartPeriod;

            itsBucketEndPeriod   = zc->itsBucketEndPeriod; 

            itsDFSpotFwd   = zc->itsDFSpotFwd;

			itsSmoothingMethod = zc->itsSmoothingMethod;

			itsFixDayCount = zc->itsFixDayCount;
        }

        void Copy(const ARM_Object* srcZc)
        {
            ARM_AbstractMarketClass::Copy(srcZc);
 
            BitwiseCopy(srcZc);
        }
 
        ARM_Object* Clone(void)
        {
            ARM_ZeroCurve* theClone = new ARM_ZeroCurve();
 
 
            theClone->Copy(this);
 
            return(theClone);
        }

        ARM_CLASS_NAME GetRootName(void)
        {
            return(ARM_ZERO_CURVE);
        }

       ~ARM_ZeroCurve(void);


        void Init(void);

        void ViewForward(char* id = NULL, FILE* ficOut = NULL);
        void View(char* id = NULL, FILE* ficOut = NULL);

        virtual double GetModelFactor(int factorId)
        {
            return(0.0);
        }

        double GetParameter(int paramId)
        {
            if (paramId < GetNbParams())
                return(itsParameters->Elt(paramId));
            else
                return(0.0);
        }
 


    /*

    Most of Yield Curve models use the instantaneous rate (DiscountYield(0.0)
    and the instantaneous forward rate SpotForwardRate(0.0) which are equal 
    for maturity=0.0. However, the usual formulae used to compute these rate
    for maturity > 0.0 do not apply for maturity=0.0. 
    
    To avoid problems raised by this issue, the 0.0 maturity is replaced by 
    YC_INSTANTANEOUS_MAT set to 1.0/365.0 in "glob.h" file which can be changed 
    
    The zero price error is less than 3x10-4 (the corresponding zero price is
    equal to 0.999726 instead of 1.0)
    
      !!! NOTE : The ForwardYield methods are given for convenience :
      the basis used to compute the rate term is ACTUAL/365. 
      So take care and use the Yield Curve Model forward yield method
      which compute the rate term with the right basis of the security
        
    */

    virtual void SetParams(double* para)
    {
        double* dummy;
 
        dummy = para;
    }

    virtual int GetNbParams(void)
    {
        return(0);
    }

    // The following methods shift the yield curve for hedge calculations
    
    inline double GetAsOfDateJul(void)
    {
        return(itsAsOfDate.GetJulian());
    }

    void SetDateTerms(ARM_Vector* dateTerms)
    {
        if ( itsDateTerms == dateTerms )
           return;

        if (itsDateTerms)
        {
           delete itsDateTerms;
           itsDateTerms = NULL;
        }

        if (dateTerms)
           itsDateTerms = dateTerms;
    }

    void SetYearTerms(ARM_Vector* yearTerms)
    {
        if ( itsYearTerms == yearTerms )
           return;

        if (itsYearTerms)
        {
           delete itsYearTerms;
           itsYearTerms = NULL;
        }

        if (yearTerms)
           itsYearTerms = yearTerms;
    }

    void SetZeroRates(ARM_Vector* zeroRates)
    {
        if ( itsZeroRates == zeroRates )
           return;

        if (itsZeroRates)
        {
           delete itsZeroRates;
           itsZeroRates = NULL;
        }

        if (zeroRates)
           itsZeroRates = zeroRates;
    }

    void SetDiscountFactors(ARM_Vector* discountFactors)
    {
        if ( itsDiscountFactors == discountFactors )
           return;

        if (itsDiscountFactors)
        {
           delete itsDiscountFactors;
           itsDiscountFactors = NULL;
        }

        if (discountFactors)
           itsDiscountFactors = discountFactors;
    }

    void SetCurrencyUnit(ARM_Currency* ccy);

	void SetFixDayCount(int fixDayCount)
	{
		itsFixDayCount = fixDayCount;
	}

	int GetFixDayCount()
	{
		return(itsFixDayCount);
	}

    inline ARM_Currency* GetCurrencyUnit(void) const
    {
        return(itsCurrency);
    }

    void InitFields(void)
    {
        SetDateTerms(NULL);
        SetYearTerms(NULL);
        SetZeroRates(NULL);
    }

    void SetBPShift(double bpShift) 
    {
        itsBPShift=bpShift;
    }
    
    inline double GetBPShift(void) const
    {
        return(itsBPShift);
    }

    void SetBPShifts(ARM_Vector* bpShifts)
    {
        itsBPShifts = bpShifts;
    }
 
    inline ARM_Vector* GetBPShifts(void) const
    {
        return(itsBPShifts);
    }

    void SetBucketAtMaturity(double maturity) 
    {
        itsBucketStartPeriod = maturity-0.00005;

        itsBucketEndPeriod = maturity+0.00005;
    }

    void SetParallelShift(double shift) 
    {
        itsBucketStartPeriod = -0.00005;
        itsBucketEndPeriod = 1E6;
        itsBPShift = shift;
    }

    void GetRidOfShift(void) 
    {
        itsBucketStartPeriod = 0.0;
        itsBucketEndPeriod = 0.0;
        itsBPShift = 0.0;
    }

    void SetBucketStartPeriod(double startPeriod) 
    {
        itsBucketStartPeriod=startPeriod;
    }

    inline double GetBucketStartPeriod(void) const
    {
        return(itsBucketStartPeriod);
    }

    void SetBucketEndPeriod(double endPeriod) 
    {
        itsBucketEndPeriod=endPeriod;
    }

    inline double GetBucketEndPeriod(void) const
    {
        return(itsBucketEndPeriod);
    }

    void SetBucket(ARM_BucketShift* bucket)
    {
        if (bucket)
        {
           SetBucketStartPeriod(bucket->GetBucketStartPeriod());

           SetBucketEndPeriod(bucket->GetBucketEndPeriod());

           SetBPShift(bucket->GetBPShift());

           if (GetBPShifts() && GetYearTerms())
           {
               int size = GetBPShifts()->GetSize();

               for (int i = 0; i < size; i++)
               {
                  if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                         && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
                     GetBPShifts()->Elt(i) = 0.01 * GetBPShift(); // unit is 1 BP
               }

           }
        }
    }

    inline double GetDFSpotFwd(void)
    {
        return itsDFSpotFwd;
    }

    void SetDFSpotFwd(double dfSpotFwd)
    {
       itsDFSpotFwd = dfSpotFwd;
    }


    // in the following methods,the reference date is the AsOfDate of the curve

    double DiscountPrice(ARM_Date& maturity);

    double DiscountYield(ARM_Date& maturity, int compMeth = 0);

    double ForwardPrice(ARM_Date& fwdDate, ARM_Date& zeroTerm);

    double ForwardYield(ARM_Date& fwdDate,
                        ARM_Date& zeroTerm, 
                        int compMeth = 0,
                        int adjFwd = 0);

    double ForwardYieldWithDayCount(ARM_Date& fwdDate,
                                    ARM_Date& zeroTerm, 
                                    int compMeth,
                                    int dayCount);

    double ForwardYieldAdj(ARM_Date& fwdDate,
                           ARM_Date& zeroTerm,
                           int compMeth = 0)
    {
        double fwdDependingOnCcy;

        fwdDependingOnCcy = ForwardYield(fwdDate, zeroTerm,
                                         compMeth, 1);

        return(fwdDependingOnCcy);
    }

    virtual double SpotForwardRate(ARM_Date &maturity);

    double DiscountPriceBPV(ARM_Date& maturity);

    double ForwardPriceBPV(ARM_Date& fwdDate, ARM_Date& zeroTerm);


    // the following methods overload the previous functions 
    // by using yearterms instead of dates

    double DiscountPrice(double yearTerm);

    virtual double DiscountYield(double yearTerm, int compMeth = 0);

    double ForwardPrice(double fwdDate, double zeroTerm);

    virtual double ForwardYield(double fwdDate, double zeroTerm,
                        int compMeth = 0);

    double DiscountPriceBPV(double yearTerm);

    double ForwardPriceBPV(double fwdDate, double zeroTerm);

    virtual double SpotForwardRate(double yearTerm);

    double Interpol(double date1, double date2, 
                    double date_interpol, double Df1,
                    double Df2, int meth_inter, int base);

    double bondInterpol(double date1, double date2, 
                        double date_interpol, double Df1, double Df2, 
                        int meth_inter);

    virtual void SetParameters(ARM_Vector* params)
    {
      if ( params == itsParameters )
         return;

      if (itsParameters)
         delete itsParameters;

      if (params)
         itsParameters = params;
      else
         itsParameters = NULL;
    }

 
    ARM_Vector* GetParameters(void)
    {
        return(itsParameters);
    }
    
    ARM_Vector* GetDateTerms(void)
    {
        return(itsDateTerms);
    }
    
    inline ARM_Vector* GetYearTerms(void)
    {
        return(itsYearTerms);
    }
    
    inline ARM_Vector* GetZeroRates(void)
    {
        return(itsZeroRates);
    }

    inline ARM_Vector* GetDiscountFactors(void)
    {
        return(itsDiscountFactors);
    }



    // Methods for TAm/T4M curve

    void CptSwapT4MZeroRates(double mean_rates, ARM_Date startSwapDate, 
                             double* first_fixing_df,
                             ARM_Matrix* SwapT4M, int &SizeT4M, 
                             double pseudo_T4M_swap_value, 
                             double pseudo_T4M_swap_mat,
                             int raw, int Cont_Lin, 
                             ARM_Currency* ccy);


    void CptSwapTAMZeroRates(double mean_rates, ARM_Date startSwapDate, 
                             double first_fixing_df,
                             ARM_Matrix SwapT4M, int SizeT4M,
                             ARM_Matrix *SwapTAM, int &SizeTAM, 
                             int raw, int Cont_Lin, 
                             ARM_Currency* ccy);

    double EquivalentRate(double MktSwaprates, int nber_of_days); 

    double InvertEquivalentRate(double MktSwaprates, int nber_of_days); 

    double DecapitalizedTAMRate(double MktSwapRates, int nber_of_days, 
                                double baseACTUAL);

    double AddOneDayBis(double julian_date, 
                        ARM_Currency* ccy);

    double RefDate(double julian_date, ARM_Currency* ccy);
        
    double KeepFourFigures(double rate);

    double FirstFloatingT4MNPV(double last_ref_df,
                               double start_date,
                               double first_df,
                               double first_date,
                               double mean,
                               int fwd_bwd,
                               int Cont_Lin, 
                               ARM_Currency* ccy,
                               double* first_ref,
                               double* first_fixing,
                               double* last_df);

    double FirstFloatingTAMNPV(double last_ref_df,
                               double start_date,
                               double first_swap_fixing,
                               double lastT4M_ref_df,
                               double lastT4M_ref_date,
                               double mean,
                               int fwd_bwd,
                               int Cont_Lin, 
                               ARM_Currency* ccy,
                               double *last_fixing_df);

    double FloatingTAMNPVWithRompu(double last_ref_df, // pilier cherche
                                   double start_date,
                                   double first_swap_fixing,
                                   double first_ref_df,// premier pilier
                                                       // de la courbe TAM
                                   double first_ref_date, 
                                   ARM_Matrix SwapT4M,
                                   int SizeT4M,
                                   double mean,
                                   int nb_months,
                                   int fwd_bwd,
                                   int Cont_Lin,
                                   ARM_Currency* ccy);


    double FirstFixedNPV(double last_ref_df,
                         double last_ref_date,
                         double first_df,
                         double first_date,
                         double swap_start_date,
                         int fwd_bwd,
                         int T4M_TAM,
                         int Cont_Lin, 
                         ARM_Currency* ccy);


        double FixedTAMNPVWithRompu(double last_ref_df,  // pilier cherche
                                    double start_date,
                                    double first_ref_df, // premier pilier 
                                                         // de la courbe TAM
                                    double first_ref_date, 
                                    double swap_rate,
                                    ARM_Matrix SwapT4M,
                                    int SizeT4M,
                                    int nb_months,
                                    int Cont_Lin,
                                    ARM_Currency* ccy);

        double FloatingNPV(double last_ref_df,
                           double last_ref_date,
                           double first_ref_df,
                           double first_ref_date,
                           double first_fixing_df,
                           double first_fixing_date,
                           double last_NPV,
                           int nb_missing_flows,
                           int T4M_TAM,
                           int Cont_Lin, 
                           ARM_Currency* ccy);


        double FixedNPV(double last_ref_df,
                        double last_ref_date,
                        double first_ref_df,
                        double first_ref_date, 
                        double last_swap_reset_date,
                        double last_NPV,
                        int nb_missing_flows,
                        int T4M_TAM,
                        int Cont_Lin, 
                        ARM_Currency* ccy);


        virtual void Set(ARM_Vector* terms, ARM_Vector* rates,
                         ARM_Date& date, int meth = 0)
        {
            SetYearTerms(terms);

            SetZeroRates(rates);

            VectorFill(itsBPShifts, itsYearTerms->GetSize(),0.0);

            SetAsOfDate(date);
        }

        ARM_ZeroCurve* ProductBy(double x)
        {
            ARM_ZeroCurve* newCurve = (ARM_ZeroCurve*)(Clone());
            ARM_Vector* zeros = newCurve->GetZeroRates();

            (*zeros) *= x;

            newCurve->SetZeroRates(zeros);
            return newCurve;
        }

        void SetMktData(ARM_CRV_TERMS& terms,
                        ARM_Vector* mktData,
                        int MMVsFut,
                        int SwapVsFut,
                        int raw,
                        int Cont_Lin)
        {
           // Set variables Market Data
           if (itsMktData)
              delete itsMktData;
           itsMktData = NULL;

           itsMktData = new ARM_MarketData(terms, mktData,
                                    MMVsFut, SwapVsFut, raw, Cont_Lin, 0);
        }


/* special method to input in a standard way the real market data into the market curve*/

	void SetMktData(ARM_CRV_TERMS& terms, ARM_Vector* mktData)
		{
           // Set variables Market Data
		if (!itsMktData){
			delete itsMktData;
			itsMktData = NULL;
			}
		itsMktData = new ARM_MarketData(terms, mktData, 1, 0, 0, 1, 0);
        }

};

extern ARM_ZeroCurve* GenerateTwoCcyBSAdjusted(ARM_ZeroCurve* crv1, // Exp: EUR: The reference CCY
                                                ARM_ZeroCurve* adjBScrv1,
                                                ARM_ZeroCurve* crv2, // Exp: JPY
                                                ARM_ZeroCurve* adjBScrv2,
                                                int spreadMatuSize,
                                                ARM_CRV_TERMS& spreadsMatu,
												int inputAsSpread = 0,
                                                int spreadsOnly = 0);

extern double ComputeBasisSpread(char* matu, // "2D", "1M", ..., "1Y", ..., "10Y", ...
                                  ARM_ZeroCurve* crv1, 
                                  ARM_ZeroCurve* adjBScrv1,                       
                                  ARM_ZeroCurve* crv2, 
                                  ARM_ZeroCurve* adjBScrv2,
                                  int spotUSD = -1);

void ConstructDynamSwap(int NbMonthYear, int flagMonthOrYear, double dataIn,
						int flagAUDTri, //flagAUDTri==1 si matu==matu2 eg.1YY
						const ARM_Date& matDate,
						int lastMMTerm,
						const ARM_Matrix& Futures, int SizeF,
						ARM_Currency* ccy,						
						const ARM_Date& swapStartDate,
						ARM_Matrix* SwapRates,
						int fxPayFreq, int& SizeS, 
						int& compteurSwap,
						int& firstSwapIndex,
						int& lastAUDtriNb);	//lastAUDtriNb: last nb of years in the part of the AUD trimestre

void PaymentFlowVector(ARM_Matrix* FlowAll, int SizeS,
					   int fxPayFreq,
					   const ARM_Date& startDate,
					   char* ccyName,
					   int indexAUD,
					   ARM_Vector& FlowFreq, int& SizeFlowFreq,
					   ARM_Vector& FlowNoFreq, int& SizeNoFlowFreq,
					   int& SizeAUDtri);


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
