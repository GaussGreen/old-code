//----------------------------------------------------------------------------
//
//	 Group		 : QR
//
//	 Filename	 : SRMRatesHJMDiffuse.hpp
//
//	 Description : Markovian HJM / Ritchken-Sankarasubramanian path generation
//
//
//----------------------------------------------------------------------------

#ifndef	SRMRATESHJMDIFF_HPP
#define	SRMRATESHJMDIFF_HPP

#include <cstdio>
#include <cassert>
#include <set>

#include "edginc/DECLARE.hpp"
#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/Maths.hpp"

#include "edginc/QMCHelperBoundedDiffusion.hpp"

#include "edginc/SRMRatesHJMUtil.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/QMCFXBaseDiffuse.hpp"


DRLIB_BEGIN_NAMESPACE

class IQMCHelperTimeLogic;

//// base class	for	low	level IR path generator
class SRMRatesHJMDiffuse: public QMCRatesDiffuse, public DiffusionDriver
{

public:

	virtual	~SRMRatesHJMDiffuse();

	//// constructor
	SRMRatesHJMDiffuse();

	virtual	void setSRMRatesHJMDiffuse(
		int					   randomIndex,
		const DateTime&		   today,
		SRMRatesHJMUtilSP	   srmRatesHJMUtilSP,
		bool				   saveDomLnMoney,
		bool				   saveSigmaR,
		double				   NbSigmasMax,
		double				   NbSigmasMin,
		const vector<double>&  df, // for historic dates
		const double*		   irFxCorr,
		bool				   calibrateAgainstSwaptionVols);

	virtual	void setIrFxCorr(const double* irFxCorr) = 0;
	virtual	void setAlphaAndRho() =	0;


    size_t getDiscYCIdx(void) // /*TODO : fix the asserts*/ { assert(ratesHJMUtilSP.get()); return registerYCFlavor(ratesHJMUtilSP->getDiscYC());} // returns Idx of Discount YC
    {
        if (discYCIdx < 0)
        {
            assert(ratesHJMUtilSP.get());
            discYCIdx =  registerYCFlavor(ratesHJMUtilSP->getDiscYC());
        }
        return discYCIdx;

    }
    size_t getDiffYCIdx(void) // { assert(ratesHJMUtilSP.get()); return registerYCFlavor(ratesHJMUtilSP->getDiffYC());} // returns Idx of Diffusion YC
    {
        if (diffYCIdx < 0)
        {
            assert(ratesHJMUtilSP.get());
            diffYCIdx =  registerYCFlavor(ratesHJMUtilSP->getDiffYC());
        }
        return diffYCIdx;
    }

    virtual QMCRatesUtilSP    getRatesUtil(void) { return getRatesHJMUtil(); }
    SRMRatesHJMUtilSP   getRatesHJMUtil(void) { ASSERT(ratesHJMUtilSP.get()); return ratesHJMUtilSP;}

	/**	finalize the timeline, allocate	memory */
	// additional initialization	when allDates is know. no more new dates after this	point
	virtual	void finalize(DateTimeArrayConstSP allDates);

	/**	called after origSVol externally recalibrated via ICE */
	virtual	void recalibrate(QMCRatesUtilSP thisSRMRatesUtilSP);

protected:

    SRMRatesHJMUtilSP ratesHJMUtilSP;

	/**	populates fields required for calculating sigmaR0 during simulation.
	A subset of	irdiffuse::CalcSigmaR */
	void calcSigmaRParams();

	/**	Precompute X and Y variables at	specified dateIdx using	supplied
	parameters */
	virtual	void computeXAndY(
        int     dateIdx,
		const   vector<double>&	kT,
		const   vector<double>&	gT,
		double  rbarRatio) = 0;

    double          qLeft;
    double          qRight;
    vector<double>  QPivot;   // [Nall]
    bool            zeroQ; /* true if qLeft = qRight = 0   */
    double          pivotRatio; // invariably zero it seems
    int             diffYCIdx;
    int             discYCIdx;
};

//declare smartPtr versions
DECLARE(SRMRatesHJMDiffuse);

//// one factor	IR
class SRMRatesHJM1F: public SRMRatesHJMDiffuse{
	enum {nfact=1};
public:

	//// for capturing state of	diffusion
	struct DiffusedState{
		double PHI00, GAMMA0;
	};
	typedef	vector<DiffusedState>::const_iterator DiffusedStateIter;

	//// constructor
    SRMRatesHJM1F() : alpha(0), IrFxCorr0(0) {}

	virtual	void setIrFxCorr(const double* irFxCorr);
	virtual	void setAlphaAndRho();

	/**	generate path across all dates.	Essentially	irdiffuse::DiffuseIR_1F	plus
	a part of irdiffuse::CalcSigmaR	*/
    void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);
	/** finalize the timeline, allocate memory */
	void finalize(DateTimeArrayConstSP allDates); // addtional initialization when allDates	is know. no	more new dates after this point

	// performance critical	functions
	double getLnExpectedDiscFactor(size_t idx, FwdIdx i, FwdIdx	j);

	double getExpectedDiscFactor(size_t	idx, FwdIdx	i, FwdIdx j) {
		return exp(this->SRMRatesHJM1F::getLnExpectedDiscFactor(idx, i, j));
	}

	double getGFactor(FwdIdx i,	FwdIdx j, int factor = 0)
	{
		const int iEDF = getFwdIdx2EdfIdx(i);
		assert(iEDF	>= 0);
		return zeta[factor][iEDF]*(partialIntegral[factor][j] -	partialIntegral[factor][i]);
	}

protected:
	/**	Precompute X and Y variables at	specified dateIdx using	supplied
	parameters */
	virtual	void computeXAndY(
        int     dateIdx,
		const   vector<double>&	kT,
		const   vector<double>&	gT,
		double  rbarRatio);

private:
	//	  friend class Gen;
	double			alpha;
	double			IrFxCorr0;
	vector<double>	X00; //	[Nall] or [0]
	vector<double>	Y0;	 //	[Nall] or [0]
	vector<double>	kFactor; //	across sim dates [Nall]
	vector<double>	gFactor; //	across sim dates [Nall]
	vector<DiffusedState> expDF; //	for	computing expected df [Nedf]
	vector<double>	partialIntegral[nfact];	// [NunionEDF];	-\int_{\tau_i}^{T_max} exp(-beta s)	fwd(s) ds  // where	beta = beta[factor]
	vector<double>	zeta[nfact]; //	[Nedf].	It's simply	exp(-beta t) to	avoid recalc.


};

//// two factor	IR
class SRMRatesHJM2F: public SRMRatesHJMDiffuse{
	enum {nfact=2};
public:
	//// for capturing state of	diffusion
	struct DiffusedState{
		double PHI00;
		double PHI01;
		double PHI11;
		double GAMMA0;
		double GAMMA1;
	};
	typedef	vector<DiffusedState>::const_iterator DiffusedStateIter;

	//// constructor
    SRMRatesHJM2F() : alp0(0), alp1(0), rho01(0),	IrFxCorr0(0), IrFxCorr1(0) {}

	virtual	void setIrFxCorr(const double* irFxCorr);
	virtual	void setAlphaAndRho();

	/**	generate path across all dates.	Essentially	irdiffuse::DiffuseIR_1F	plus
	a part of irdiffuse::CalcSigmaR	*/
    void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);
	/** finalize the timeline,	allocate memory	*/
	void finalize(DateTimeArrayConstSP allDates); // addtional initialization when allDates	is know. no	more new dates after this point

	// performance critical	functions
	double getLnExpectedDiscFactor(size_t idx, FwdIdx i, FwdIdx	j);

	double getExpectedDiscFactor(size_t	idx, FwdIdx	i, FwdIdx j) {
		return exp(this->SRMRatesHJM2F::getLnExpectedDiscFactor(idx, i, j));
	}
	double getGFactor(FwdIdx i,	FwdIdx j, int factor = 0) {
		const int iEDF = getFwdIdx2EdfIdx(i);
		assert(iEDF	>= 0);
		return zeta[factor][iEDF]*(partialIntegral[factor][j] -	partialIntegral[factor][i]);
	}

protected:
	/**	Precompute X and Y variables at	specified dateIdx using	supplied
	parameters */
	virtual	void computeXAndY(
        int     dateIdx,
		const   vector<double>&	kT,
		const   vector<double>&	gT,
		double  rbarRatio);

private:
	struct Vars{
		double X00;
		double X11;
		double X01;
		double Y0;
		double Y1;
		double g0;
		double g1;
		double k0;
		double k1;
	};

	//	  friend class Gen;
	// same	names as in	SRM3
	double			alp0;
	double			alp1;
	double			rho01;
	double			IrFxCorr0;
	double			IrFxCorr1;
	vector<Vars>	varsPerTimePoint; // across	sim	dates
	vector<DiffusedState> expDF; //	for	computing expected df

	vector<double>	partialIntegral[nfact];	// [NunionEDF];	-\int_{\tau_i}^{T_max} exp(-beta s)	fwd(s) ds  // where	beta = beta[factor]
	vector<double>	zeta[nfact]; //	[Nedf].	It's simply	exp(-beta t) to	avoid recalc.


};


//// three factor IR
class SRMRatesHJM3F: public SRMRatesHJMDiffuse{
	enum {nfact=3};
public:
	//// for capturing state of	diffusion
	struct DiffusedState{
		double PHI00;
		double PHI01;
		double PHI02;
		double PHI11;
		double PHI12;
		double PHI22;
		double GAMMA0;
		double GAMMA1;
		double GAMMA2;
	};
	typedef	vector<DiffusedState>::const_iterator DiffusedStateIter;

	//// constructor
    SRMRatesHJM3F() : alp0(0), alp1(0), alp2(0), rho01(0), rho02(0), rho12(0),
		IrFxCorr0(0), IrFxCorr1(0),	IrFxCorr2(0) {}

		virtual	void setIrFxCorr(const double* irFxCorr);
		virtual	void setAlphaAndRho();

		/**	generate path across all dates.	Essentially	irdiffuse::DiffuseIR_1F	plus
		a part of irdiffuse::CalcSigmaR	*/
        void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);
		/** finalize the timeline,	allocate memory	*/
		void finalize(DateTimeArrayConstSP allDates); // addtional initialization when allDates	is know. no	more new dates after this point

		// performance critical	functions
		double getLnExpectedDiscFactor(size_t idx, FwdIdx i, FwdIdx	j);

		double getExpectedDiscFactor(size_t	idx, FwdIdx	i, FwdIdx j) {
			return exp(this->SRMRatesHJM3F::getLnExpectedDiscFactor(idx, i, j));
		}

		double getGFactor(FwdIdx i,	FwdIdx j, int factor = 0) {
			const int iEDF = getFwdIdx2EdfIdx(i);
			assert(iEDF	>= 0);
			double g = zeta[factor][iEDF]*(partialIntegral[factor][j] -	partialIntegral[factor][i]);
			return g;
		}



protected:
	/**	Precompute X and Y variables at	specified dateIdx using	supplied
	parameters */
	virtual	void computeXAndY(
        int     dateIdx,
		const   vector<double>&	kT,
		const   vector<double>&	gT,
		double  rbarRatio);

private:
	friend class SRMCreditHJM3F;		//FIXME: modify	later
	struct Vars{
		double X00;
		double X11;
		double X22;
		double X01;
		double X02;
		double X12;
		double Y0;
		double Y1;
		double Y2;
		double g0;
		double g1;
		double g2;
		double k0;
		double k1;
		double k2;
	};
	// friend class Gen;
	// same	names as in	SRM3
	double			alp0;
	double			alp1;
	double			alp2;
	double			rho01;
	double			rho02;
	double			rho12;
	double			IrFxCorr0;
	double			IrFxCorr1;
	double			IrFxCorr2;
	vector<Vars>	varsPerTimePoint; // across	sim	dates [Nall]
	vector<DiffusedState> expDF; //for computing expected df  [Nedf]
	// tabulated values for	calculating gFactor(t,T)
	// [NunionEDF]; -\int_{\tau_i}^{T_max} exp(-beta s)	fwd(s) ds
	// where beta = beta[factor]
	vector<double>	partialIntegral[nfact];

	vector<double>	zeta[nfact]; //[Nedf]. It's simply exp(-beta t) to avoid recalc.
};

DRLIB_END_NAMESPACE
#endif // SRMRATESHJMDIFF_HPP
