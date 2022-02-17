//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMInitialGuess.hpp
//
//   Description : Initial Guess for SRM
//
//   Author      : Louis A Moussu
//
//   Date        : 14 February 2006
//
//----------------------------------------------------------------------------

#ifndef EDR_SRMINITIALGUESS_HPP
#define EDR_SRMINITIALGUESS_HPP

#include "edginc/SRMRatesFactor.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/VolProcessedBSIR.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VanillaGrid.hpp"

DRLIB_BEGIN_NAMESPACE

class SRMInitialGuess: public CObject {

	public:
		
		
		static CClassConstSP const TYPE;
		SRMInitialGuess();
		// constructor for the calibration
		SRMInitialGuess(MarketDataSP				market,
						Calibrator::ObjFuncSP		objFunc,
						const Calibrator::InstanceIDArray& ids,
						const string&				methodType,
						const string&				volType);
		// constructor for the addin
		SRMInitialGuess(MarketDataSP			market,
						CAssetWrapper			asset,
						ExpiryArraySP			expiries,
						const DoubleArray&		strike,
						const string&			volType,
						const string&			methodType,
						const string&			volMap);		
		// save the guess in the objFunc for the calibration
		void saveGuess(Calibrator::ObjFuncSP objFunc);
		// get parameters
		ExpiryArraySP getExpiries();
		int			  getNbExpiries();
		string		  getVolMap();
		DoubleArraySP getAtmVol();
        void adjustAtmVol(DoubleArraySP guessAtmVal);
        void adjustSmileA1(DoubleArraySP guessSmileA1);
        void adjustSmileA2(DoubleArraySP guessSmileA2);
		DoubleArraySP getSmileA1();
		DoubleArraySP getSmileA2();
		DoubleArraySP getSmileA3();
        
        void validate(CDoubleMatrixSP result);

        int getCenterIndex(VanillaGrid::LeastSquareSimpleSP objFuncLS);

        void getAdjustedGuess(CDoubleMatrix sensitivities, CDoubleMatrix valGuess, DoubleArraySP& adjustedGuess);

        void tweakOneParameter(double shift, DoubleArraySP guessParam, 
                               DoubleArray& shiftGuess, DoubleArraySP& tweakGuess);

        void nbrOfExpiriesToImprove(CDoubleMatrix valGuess,double tol,int& count,DoubleArraySP& myRes);

        void getObjFuncArray(CDoubleMatrix valGuess,DoubleArraySP& objFuncVals, int& count1, int& count2,int& count3);
     
        int getNbStrikes();
        

	private:

		MarketDataSP market;
		CAssetWrapper asset;
		ExpiryArraySP expiries;
		DoubleArray strike;
		string volType;
		string methodType;
		string volMap;
		string volName; // $unregistered
		
		DoubleArraySP atmVol; // $unregistered
		DoubleArraySP smileA1; // $unregistered
		DoubleArraySP smileA2; // $unregistered
		DoubleArraySP smileA3; // $unregistered
        
        /*DoubleArraySP ASmile;
        DoubleArraySP BSmile;
        DoubleArraySP CSmile;*/
       

		//compute guess
		void computeGuess();

		// methods
		CDoubleMatrixSP smileImpVol();
		CDoubleMatrixSP smileLocVol();
		CDoubleMatrixSP smileDupireRates();
		CDoubleMatrixSP smileGen();
		DoubleArraySP	computeAtmVol();
		void validateStrike(VanillaGrid::LeastSquareSimpleSP objFuncLS);
		void validateGuess(CDoubleMatrixSP);
		void mapGuess(CDoubleMatrixSP);
        
		static void load(CClassSP& clazz);
		static IObject* defaultSRMInitialGuess();

	public:

		static const string SPOTVOL;
		static const string COMPVOL;
		static const string COMPMET;
		static const string IMPMET;
		static const string LOCMET;
		static const string DUPIREMET;
		static const string GENMET;
		static const string VOLPREF;

};

DRLIB_END_NAMESPACE
#endif