//	MlEqVolatilityStructure.h : Volatility structure & volatiltiy data
//
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQVOLATILITYSTRUCTURE_H_
#define _MLEQVOLATILITYSTRUCTURE_H_

#include <vector>
#include "cMatrix.h"
#include "MlEqInterpolator.h"
#include "utility.h"
#include "smart.h"
#include <math.h>
#include "solve.h"
#include "utility.h"
#include "MlEqBilinearInterpolator.h"
#include "MlEqDate.h"


/******************************************************************
**				
**				Volatility Structure Class
**				
******************************************************************/
class  MlEqVolatilityStructure : public RCObject
{					
public:	
	MlEqVolatilityStructure();

	MlEqVolatilityStructure(		
		const double					spot, 
		const CVector&					discountFactors, 
		const CVector&					forwards,
		const vector<DATE>&				terms, 
		MlEqDateHandle					hDate, 
		const vector < MlEqVolDataHandle> &		volBidData,
		const vector < MlEqVolDataHandle> &		volMidData,
		const vector < MlEqVolDataHandle> &		volAskData,
		const MlEqInterpolatorHandle&	interpAcrossTime,
		const int interpAcrossTimeFlag,
		const VolatilityDataTypeEnum& voltype,
		MlEqVolatilityStructureHandle pAsymptotic=NULL,
		double decayFactor = 1.0);

	void initialize( 
		const double					spot, 
		const CVector&					discountFactors, 
		const CVector&					forwards,
		const vector<DATE>&				terms, 
		MlEqDateHandle					hDate,
		const vector < MlEqVolDataHandle> &		volBidData,
		const vector < MlEqVolDataHandle> &		volMidData,
		const vector < MlEqVolDataHandle> &		volAskData,
		const MlEqInterpolatorHandle&	interpAcrossTime,
		const int interpAcrossTimeFlag,
		const VolatilityDataTypeEnum& voltype,
		MlEqVolatilityStructureHandle pAsymptotic= NULL,
		double decayFactor = 1.0);

	void initialize(
		const MlEqVolatilityStructure& spotVolStruct,
		const MlEqVolatilityStructure& asymptodicVolStruct,
		double decayFactor);
	
	MlEqVolatilityStructure(
		const MlEqVolatilityStructureHandle& spotVolStruct,
		const MlEqVolatilityStructureHandle& asymptodicVolStruct,
		double decayFactor);
	
	MlEqVolatilityStructure(const MlEqVolatilityStructure& rhs);
	virtual MlEqVolatilityStructureHandle copy() const {return new MlEqVolatilityStructure(*this);}
	virtual ~MlEqVolatilityStructure();
	virtual MlEqVolatilityStructure& operator=( const MlEqVolatilityStructure& rhs );		

public:
	DataSourceEnum								GetDataSource(void) const;
	void										PutDataSource(DataSourceEnum ds);
			
	long GetDateNumber(void) const;
	void PutDateNumber(long nDate);

//	void getVolDeriv(CVector& res,const MlEqStrike& strike,const DATE& date	) ;
	void getVolDeriv(CVector& res,const MlEqStrike& strike,const DATE& date, const CVector& mktForwards	) ;
	void getVolDeriv(CVector& res, const MlEqStrike& strike, int slice, double mktForward ) ;
	void getVolBumpDeriv(CVector& dvol, const MlEqStrike& strike, DATE date ) ;

	virtual double getVol(const MlEqStrike& strike, int islice, BidAskMidEnum bidOrAsk = Middle) ;
	virtual double getVol(const MlEqStrike& strike,  DATE date, BidAskMidEnum bidOrAsk = Middle) ;	
	virtual double getDensity(const MlEqStrike& strike, DATE date,double forward, BidAskMidEnum bidOrAsk = Middle, int cummulativeFlag = 0);

	virtual double StrikeInTermsOfStandartDeviations(const MlEqStrike& strike,const DATE& date,double forward, BidAskMidEnum bidOrAsk = Middle);

	double getSpot(){return m_spot;};
	const CVector& getDiscounts(){return m_discountFactors;};
	const CVector& getForwards(){return m_forwards;};
	MlEqConstDateHandle getDateToDouble(void) const;
	const vector<DATE>& getDates(){return m_Dates;};


	void putDateToDouble(MlEqDateHandle hDate);


	int getInterpAcrossTimeFlag(){return m_interpAcrossTimeFlag;};

	MlEqVolDataHandle   getAskVolData(int& numberStrikes,int islice,int istrike)const;
	MlEqVolDataHandle   getMidVolData(int& numberStrikes,int islice,int istrike)const;
	MlEqVolDataHandle   getBidVolData(int& numberStrikes,int islice,int istrike)const;

	virtual double getFutureVol(const MlEqStrike& strike,const DATE& startdate,const DATE& enddate,double spot,BidAskMidEnum bidOrAsk = Middle) ;
	virtual double getFutureVol(const MlEqStrike& strike,const DATE& startdate,const DATE& enddate,double spot,double fwdShort,double fwdLong,BidAskMidEnum bidOrAsk = Middle) ;
	virtual double getFutureVol(const MlEqAsset& asset,const MlEqStrike& strike,MlEqDate& startdate,long enddate,double volRef,const MlEqStrikeHandle& strikeInvariant,BidAskMidEnum bidOrAsk) ;
	virtual double getFutureVol(const MlEqAsset& asset,const MlEqStrike& strike,MlEqDate& startdate,long enddate,const MlEqStrikeHandle& strikeInvariant,BidAskMidEnum bidOrAsk) ;
	virtual double getFutureVol(const MlEqAsset& asset,const MlEqStrike& strike,MlEqDate& startdateHandle,long enddate,double volRef,BidAskMidEnum bidOrAsk) ;
	virtual double getFutureVol(const MlEqAsset& asset,const MlEqStrike& strike,MlEqDate& startdateHandle,long enddate,BidAskMidEnum bidOrAsk) ;

	double  getNaiveFutureVol(const MlEqStrike& shortStrike,const MlEqStrike& longStrike,const DATE& startdate,const DATE& enddate,BidAskMidEnum bidOrAsk = Middle) ;

	VolatilityDataTypeEnum				m_volType;
	
	vector <MlEqVolDataHandle >			getBidVolData(void) const {return m_BidVolData;}
	vector <MlEqVolDataHandle >			getMidVolData(void) const {return m_MidVolData;}
	vector <MlEqVolDataHandle >			getAskVolData(void) const {return m_AskVolData;}

	MlEqInterpolatorHandle				getInterpAcrossTime(void) const {return m_InterpAcrossTime;}

	void reinitializeStrikes(const MlEqAsset& asset);



public:
	MlEqVolatilityStructureHandle		m_asymtodicVolStruct;
	double								m_decayFactor;

	void setAsymtodicVolStruct(MlEqVolatilityStructureHandle volstruct);
	static int findATMFwdVol(int& result,vector< MlEqStrikeHandle >& Strikes, double fwd);
	static MlEqStrikeHandle setUpDefaultStrikeInvariant(MlEqStrikeHandle strikeInvariant,const MlEqAsset& asset,MlEqDateHandle startDateHandle,long enddate);

protected:	
	double getVol(const MlEqStrike& strike, int islice,	std::vector<MlEqVolDataHandle>& VolData);
	double getVol(const MlEqStrike& strike, DATE date, std::vector<MlEqVolDataHandle>& VolData);	
	void getYDataAcrossTime(const MlEqStrike& strike,CVector& yData,CVector& xData,int startPoint,int numPoints,vector <MlEqVolDataHandle >	& volData);
	double getVolFromYData(const MlEqStrike& strike,double xtime,double yData,CVector& xValues);
		
	// Member Data
	
	////////////////////////////////////////
	// Spot + Information to reconstruct Prices	
	
	double				m_spot;
	CVector				m_discountFactors;		// [idate]
	CVector				m_forwards;				// [idate]
	vector<DATE>		m_Dates;		
	MlEqDateHandle		m_hDate;				// contains the date of when these info were saved
	
	////////////////////////////////////////
	
	int									m_interpAcrossTimeFlag;		// [0] vol; [1] var; [2] parameters 
	vector <MlEqVolDataHandle >			m_BidVolData;				// [islice]	
	vector <MlEqVolDataHandle >			m_MidVolData;				// [islice]
	vector <MlEqVolDataHandle >			m_AskVolData;				// [islice]	
	MlEqInterpolatorHandle				m_InterpAcrossTime;			// this will interpolate across time
	DataSourceEnum						m_ds;
	int									m_asymptodicVolIsInitialized;

// Scenario Classes
public:
	enum scenario_type {
		stVolShift = 1,					// This scenario affects the get volatility
		stBetaShift = 2,				// This scenario affects the beta shift
		stForwardSkewShift = 3,			// This scenario affects MlEqVolatilityStructure::shiftFwdSkew
		stPostVolShift = 4				// These scenarios affect the get volatility function but are called after stVolShift scenarios
	};

protected:			
	// Base scenario class
	class Scenario : public RCObject
	{		
	public:		
		virtual double Apply(double f, MlEqConstStrikeHandle hStrike, DATE date1, DATE date2, double fReferenceVolatility) const = 0;
		virtual bool IsResetting(void) const {return false;}		// Set to true if the Apply method does not depend on the input value.
		virtual scenario_type GetScenarioType(void) const = 0;
	};
		
	// BetaShift type classes
	class Scenario_ShiftBeta : public Scenario
	{
	public:
		explicit Scenario_ShiftBeta(double bShift, estring szTenor);
		double Apply(double f, MlEqConstStrikeHandle, DATE dateStart, DATE dateEnd, double) const;
		scenario_type GetScenarioType(void) const;
	protected:
		MlEqInterpolatorHandle				m_hBetaShift;
	};

	// ForwardSkewShift type classes
	class Scenario_ShiftForwardSkew : public Scenario
	{
	public:
		explicit Scenario_ShiftForwardSkew(long nWindowStartDate, long nWindowEndDate, long nPeriodStartDays, long nPeriodEndDays, double fFwdSkewShiftAmount);
		double Apply(double f, MlEqConstStrikeHandle hStrikeInvariant, DATE dateStart, DATE dateEnd, double fReferenceVolatility) const;
		scenario_type GetScenarioType(void) const;

	protected:		
		long								m_nWindowStartDate;
		long								m_nWindowEndDate;
		long								m_nPeriodStartDays;
		long								m_nPeriodEndDays;
		double								m_fFwdSkewShiftAmount;
	};
	
	// VolShift type classes
	class Scenario_GetVolType : public Scenario
	{
	public:
		scenario_type GetScenarioType(void) const;
	};
	
	class Scenario_RelativeShift : public Scenario_GetVolType
	{
	public:
		explicit Scenario_RelativeShift(double Amount);
		double Apply(double f, MlEqConstStrikeHandle, DATE, DATE, double) const;

	protected:
		double								m_fAmount;
	};

	class Scenario_AbsoluteShift : public Scenario_GetVolType
	{
	public:
		explicit Scenario_AbsoluteShift(double Amount);
		double Apply(double f, MlEqConstStrikeHandle, DATE, DATE, double) const;

	protected:
		double								m_fAmount;
	};

	class Scenario_Set : public Scenario_GetVolType
	{
	public:
		explicit Scenario_Set(double Amount);
		double Apply(double f, MlEqConstStrikeHandle, DATE, DATE, double) const;
		bool IsResetting(void) const;

	protected:
		double								m_fAmount;
	};

	class Scenario_ShiftSkew : public Scenario_GetVolType
	{
	public:
		explicit Scenario_ShiftSkew(const std::vector<MlEqStrikeHandle>& Strikes, const CVector& Maturities, const CMatrix& Shift, MlEqConstAssetHandle hAsset);
		double Apply(double f, MlEqConstStrikeHandle hStrike, DATE date, DATE, double) const;

	protected:
		std::vector<MlEqStrikeHandle>		m_shiftStrikes;
		MlEqConstAssetHandle				m_hAsset;
		MlEq2DInterpolatorHandle			m_ShiftInterpolator;

	private:
		mutable bool						m_bRecursionBlock;
	};

	class Scenario_ShiftStrikeSkew : public Scenario_GetVolType
	{
	public:
		explicit Scenario_ShiftStrikeSkew(MlEqInterpolatorHandle hVolShiftInterpolator, MlEqInterpolatorHandle hVolSkewShiftInterpolator, MlEqConstAssetHandle hAsset, bool bIsRelative, bool bStrikeSkewInvariant);
		double Apply(double f, MlEqConstStrikeHandle hStrike, DATE date, DATE, double) const;

	protected:
		MlEqInterpolatorHandle				m_hShiftStrikeInvariantVolInterpolator;
		MlEqInterpolatorHandle				m_hShiftStrikeInvariantSkewInterpolator;
		MlEqConstAssetHandle				m_hAsset;		
		bool								m_bStrikeSkewInvariant;
		bool								m_bRelativeShift;

	private:
		mutable bool						m_bRecursionBlock;
	};

	// PostVolShift type classes
	class Scenario_FsaLiquidityShift : public Scenario		
	{
	public:
		explicit Scenario_FsaLiquidityShift(double fAverageVolume, double fCoreVolatility, double fNetGamma, double fPositionDirection, bool bIsStock);
		double Apply(double f, MlEqConstStrikeHandle, DATE date, DATE, double) const;
		scenario_type GetScenarioType(void) const;

	protected:
		double m_fAverageVolume;
		double m_fCoreVolatility;
		double m_fNetGamma;
		double m_fPositionDirection;
		bool m_bIsStock;				
	};

	// Scenarios storage class
	class Scenarios : private std::map<scenario_type, std::vector<RCPtr<Scenario> > >
	{
	public:		
		void push_back(const RCPtr<Scenario>& h);
		double Apply(double f, MlEqConstStrikeHandle hStrike, DATE date1, DATE date2, double fReferenceVolatility, scenario_type st) const;
		bool IsFirstResetting(scenario_type st) const;
		long size(scenario_type st) const;			
	} m_scenarios, m_scenarios_original;
				
public:		
	void FsaLiquidityShift(double fAverageVolume, double fCoreVolatility, double fNetGamma, double fPositionDirection, bool bIsStock);
	void parallelShiftVol(double VolShift, bool shiftIsRelative = false);
	void Reset();
	void SetVol(double vol);
	void shiftSkew(std::vector<MlEqStrikeHandle>& Strikes, CVector& Maturities, CMatrix& Shift, MlEqConstAssetHandle hAsset);
	void shiftStrikeInvariantSkew(MlEqInterpolatorHandle VolShiftInterpolator, MlEqInterpolatorHandle VolSkewShiftInterpolator, MlEqAssetHandle hAsset, bool shiftIsRelative = false);
	void shiftStrikeSkew(MlEqInterpolatorHandle VolShiftInterpolator, MlEqInterpolatorHandle VolSkewShiftInterpolator, MlEqAssetHandle hAsset);
	void Stick(void);
	
	void setUpFwdSkewShift(double shiftAmount, long WindowStartDate, long WindowEndDate, long PeriodStartDays, long PeriodEndDays);		
	virtual void setUpBetaShift(double bshift, const std::string& szTenor);
		
protected:
	double shiftFwdSkew(long startdate,long enddate,double volRef,const MlEqStrikeHandle& strikeInvariant);	
};		





/******************************************************************
**				
**		Stochastic Beta Volatility Structure Class
**				
******************************************************************/

	
class  MlEqStochBetaVolatilityStructure:  public MlEqVolatilityStructure
{	
public:
	MlEq2DInterpolatorHandle					m_hbiLevel;
	MlEq2DInterpolatorHandle					m_hbiBeta;
	MlEq2DInterpolatorHandle					m_hbiVolVol;
	MlEq2DInterpolatorHandle					m_hbiMeanReversion;	
	MlEq2DInterpolatorHandle					m_excessVolOfVol;
	MlEqInterpolatorHandle						m_hermiteElasticity;

	virtual double getFutureVol(const MlEqStrike& strike, const DATE& startdate, const DATE& enddate, double spot, BidAskMidEnum bidOrAsk = Middle);
	virtual double getFutureVol(const MlEqStrike& strike, const DATE& startdate, const DATE& enddate, double spot, double fwdLong, double fwdShort, BidAskMidEnum bidOrAsk = Middle);
	virtual double getFutureVol(const MlEqAsset& asset, const MlEqStrike& strike, MlEqDate& startdate, long enddate, double volRef, const MlEqStrikeHandle& strikeInvariant, BidAskMidEnum bidOrAsk);
	virtual double getFutureVol(const MlEqAsset& asset, const MlEqStrike& strike, MlEqDate& startdate, long enddate, const MlEqStrikeHandle& strikeInvariant, BidAskMidEnum bidOrAsk);
	virtual double getFutureVol(const MlEqAsset& asset, const MlEqStrike& strike, MlEqDate& startdateHandle, long enddate, double volRef, BidAskMidEnum bidOrAsk);
	virtual double getFutureVol(const MlEqAsset& asset, const MlEqStrike& strike, MlEqDate& startdateHandle, long enddate, BidAskMidEnum bidOrAsk);

	static double getFrequencyRoundUp(DATE dateStart, DATE dateEnd);
		
	double getLevel(DATE dateStart, DATE dateEnd)
	{
		double frequencyRU = getFrequencyRoundUp(dateStart, dateEnd);
		double f = m_hbiLevel->getValue(dateStart, frequencyRU);
		return f;
	}
	
	double getBeta(DATE dateStart, DATE dateEnd)
	{
		double frequencyRU = getFrequencyRoundUp(dateStart, dateEnd);
		double f = m_hbiBeta->getValue(dateStart, frequencyRU);		
		f = m_scenarios.Apply(f, NULL, dateStart, dateEnd, 0.0, stBetaShift);
		return f;
	}
	
	double getVolVol(DATE dateStart, DATE dateEnd)
	{
		double frequencyRU = getFrequencyRoundUp(dateStart, dateEnd);
		double f = m_hbiVolVol->getValue(dateStart, frequencyRU);
		return f;
	}
	
	double getMeanReversion(DATE dateStart, DATE dateEnd)
	{
		double frequencyRU = getFrequencyRoundUp(dateStart, dateEnd);
		double f = m_hbiMeanReversion->getValue(dateStart, frequencyRU);
		return f;
	}

	void setUpBetaShift(double bshift = 0.0, const std::string& szTenor = "");
};





















/******************************************************************
**
**	Generic VolData Class
**
******************************************************************/
	
class  MlEqVolData :  public RCObject
{		
	
public:
	
	MlEqVolData(){};
	
	MlEqVolData( const MlEqVolData& rhs );
	virtual ~MlEqVolData();
	
	MlEqVolData& operator=(const MlEqVolData& rhs);

	virtual MlEqVolDataHandle copy() const {return new MlEqVolData(*this);}
	
	MlEqVolData( vector<MlEqInterpolatorHandle>&	InterpAcrossStrike,
				 vector< vector< MlEqStrikeHandle > >&	strikes	);

	MlEqVolData( MlEqInterpolatorHandle&	InterpAcrossStrike,
				 vector< MlEqStrikeHandle > &	strikes	);

	MlEqVolData( vector< MlEqStrikeHandle  >&	strikes	,int numberSlices);


	void initialize(vector< MlEqStrikeHandle  >&	strikes	,int numberSlices);
	void initialize(const vector < vector< MlEqStrikeHandle  > >&	strikes	,int numberSlices);

	void initialize(	vector<MlEqInterpolatorHandle>&	InterpAcrossStrike,
						vector< vector< MlEqStrikeHandle > >&	strikes	);



public:


//	virtual void getVolDeriv(CVector& res,const MlEqStrike& strike,int whichSlice) ;
	virtual void getVolDeriv(CVector& res,const MlEqStrike& strike,int whichSlice, double mktForward) ;

	virtual int		getNumberOfSlices() const;
	virtual int		getNumberOfStikes(int idate)const;
	virtual double	getVol(const MlEqStrike& strike,int islice) ;

	virtual void	PutNumberOfSlices(int n)			{m_numberOfSlices = n;}

	double			InterpAcrossStrikeGetValue(double xval,int whichSlice) const;		

	const vector< vector< MlEqStrikeHandle > >&  getStrikes()const{return m_strikes;} ;


	virtual void reinitializeStrikes(const MlEqAsset& asset);

	
protected:
	virtual void getVolSpaceDerivs(CVector& res, MlEqStrike& strike, int slice) ;

	const MlEqStrikeHandle& getpStrike(int idate,int istrike);

	int	m_numberOfSlices;
	
	vector<MlEqInterpolatorHandle>	m_InterpAcrossStrike;	//[idate]
	vector< vector< MlEqStrikeHandle > >	m_strikes;				//[idate][istrike]
};	
	
		
/******************************************************************
**
**				JumpWing VolData Class
**
******************************************************************/

	
struct JumpWingCoeffs
{	
	double	ATMVolF;
	double	a;
	double	CallWing;
	double	PutWing;
	double	ASym;
	double  t;
	long	currentDate;

	double	skew;
	double	volMin;

	JumpWingCoeffs();

};	
	
	
class MlEqJumpWingVolData: public MlEqVolData
{
	
public:	


	MlEqJumpWingVolData(){};
	MlEqJumpWingVolData(const MlEqJumpWingVolData& rhs);
	MlEqJumpWingVolData( vector< vector< MlEqStrikeHandle > > &strikes, vector<JumpWingCoeffs> & jmp);
	MlEqJumpWingVolData( vector < MlEqStrikeHandle  > &strikes,JumpWingCoeffs& jmp);


	void initialize( vector < MlEqStrikeHandle  > &strikes,JumpWingCoeffs& jmp);
	void initialize( vector < MlEqStrikeHandle  > &strikes,vector < JumpWingCoeffs >& jmp);
	void initialize( vector< vector< MlEqStrikeHandle > > &strikes, vector<JumpWingCoeffs> & jmp);
	void initialize(const MlEqJumpWingVolData& rhs);

	MlEqVolDataHandle copy() const{ return new MlEqJumpWingVolData( *this ); }


//	int getNumberOfSlices()const { return m_numberOfSlices;};
	virtual double getVol(const MlEqStrike& strike,int islice) ;

	const 	vector<JumpWingCoeffs>&  getJumpWingParameters()const{return m_jw;} ;

	
protected:

	
	void setVolMin(JumpWingCoeffs& jmp);
	void setSkewFactor(JumpWingCoeffs& jmp);
	void setSkewStrenght_and_Asymmetry(JumpWingCoeffs& jmp);


	double getVol(const JumpWingCoeffs& jmp,double xstrike) ;

	// Jump Wing order is:	
	// pars			= ATMVolF, a, Call Wing slope, Put Wing slope, asymmetry param,
	//				   volO	   a         c                  p           b
	//					0	   1	     2					3			4		
	// vol^2 = f^2 = vol0^2 * ( 1 + a * log( (exp(c*z) + b*exp(-p*z)) / (1+b) ) )
	
	vector<JumpWingCoeffs> m_jw;//[islice]

};	


/******************************************************************
**
**				JumpWingFitVolData Fit Class
**
******************************************************************/


class MlEqJumpWingFitVolData: public MlEqJumpWingVolData, public LSGRGSolver						  
{		
			

protected:

		MlEqVolatilityStructureHandle	m_mktVols;	

		int m_totalNumVariables;		
		int m_currentidate;
		int m_fittingDimension;		

		vector < vector< MlEqStrikeHandle  > >	m_FittingStrikes;


public:
		
		int m_useVega;
		double m_fitSpread;
		double m_fitMidVols;

		MlEqJumpWingFitVolData(){};
		MlEqJumpWingFitVolData(MlEqVolatilityStructureHandle& mktVols,vector < vector< MlEqStrikeHandle  > >&	FittingStrikes	,vector < long > fittingDates,CMatrix& BoundsOnConstraints,int useVega,double fitSpread,double fitMidVols);
		MlEqJumpWingFitVolData(vector<JumpWingCoeffs>& JWGuess,MlEqVolatilityStructureHandle& mktVols,vector < vector< MlEqStrikeHandle  > >&	FittingStrikes	,vector < long > fittingDates,CMatrix& BoundsOnConstraints,int useVega,double fitSpread,double fitMidVols);
		MlEqJumpWingFitVolData(MlEqJumpWingVolDataHandle guess,MlEqVolatilityStructureHandle& mktVols,vector < long > fittingDates,CMatrix& BoundsOnConstraints,int useVega,double fitSpread,double fitMidVols);

		void initialize(vector<JumpWingCoeffs>& JWGuess,MlEqVolatilityStructureHandle& mktVols,const vector < vector< MlEqStrikeHandle  > >&	FittingStrikes	,vector < long > fittingDates,CMatrix& BoundsOnConstraints,int useVega,double fitSpread,double fitMidVols,bool initilaizeBaseClass=true);	
		void initialize(MlEqJumpWingVolDataHandle guess,MlEqVolatilityStructureHandle& mktVols,vector < long > fittingDates,CMatrix& BoundsOnConstraints,int useVega,double fitSpread,double fitMidVols);

		void  ObjectiveFcn(double* gVals,double* xVals);		
		
};		
			

/******************************************************************
**
**				Hull VolData Class
**
******************************************************************/
	

class HullCoeffs
{
  public:
		double	VoltMin;
		double	NsAtmVoltMin;
		double	RightInflectionPoint;
		double	CallWing;
		double	LeftInflectionPoint;
		double	PutWing;
};
		
double getHullVol(double normStrike, const HullCoeffs& sc);

		
class MLHullVolData: public MlEqVolData
{

public:

	MLHullVolData( vector< vector< MlEqStrikeHandle > > &strikes, vector<HullCoeffs> & hull);
	MLHullVolData():MlEqVolData(){};
	MLHullVolData(const MLHullVolData& rhs);
	MlEqVolDataHandle copy() const{ return new MLHullVolData( *this ); }
	MLHullVolData& operator=(const MLHullVolData& rhs);
	void initialize( vector< vector< MlEqStrikeHandle > > &strikes, vector<HullCoeffs> & hull);

	int getNumberOfSlices()const				{return m_numberOfSlices;};
	const HullCoeffs& getHullData(int whichSlice)const;

	virtual double getVol(const MlEqStrike& strike,int islice) ;


protected:

	// Hull order is:	(L = Left = put, R = Right = call)
	// pars			= volMin, nsMin, nsInflecR, volR-volMin, nsInflecL, volL-volMin, 
	//					0	   1	     2	        3	     4		5

	vector < HullCoeffs >  m_HullData; //idate]
	int m_numberOfSlices;

};


/******************************************************************
**
**				RAM VolData Class
**
******************************************************************/
	

class RamCoeffs
{
  public:

		double Vol90Strike;
		double Vol110Strike;

		
		double VolLowStrike;
		double VolHighStrike;
		double t;
		DATE currentDate;
		MlEqStrikeHandle m_strike;

  public:
		double	BaseVol;
		double	Slope;
		double	Curvature;
		MlEqStrikeHandle	UpCutoff;
		MlEqStrikeHandle	DownCutoff;

};
		
double getRamVol(double spotBasedStrike, const RamCoeffs& sc);
void   getRamVolDerivs(CVector& res, double spotBasedStrike, const RamCoeffs& sc);

		
class MLRamVolData: public MlEqVolData
{

public:

	MLRamVolData(vector< vector< MlEqStrikeHandle > >& strikes, vector<RamCoeffs>& rm);
	MLRamVolData():MlEqVolData(){}
	MLRamVolData(const MLRamVolData& rhs);
	MlEqVolDataHandle copy() const{ return new MLRamVolData( *this ); }
	MLRamVolData& operator=(const MLRamVolData& rhs);
	
	void initialize( vector< vector< MlEqStrikeHandle > > &strikes, vector<RamCoeffs> & rm);

	const RamCoeffs& getRamData(int whichSlice)const;

	const vector<RamCoeffs>&  getRamParameters()const {return m_RamData;}

	virtual double getVol(const MlEqStrike& strike,int islice) ;

	virtual void getVolSpaceDerivs(CVector& res, MlEqStrike& strike,int whichSlice) ;

protected:


	vector < RamCoeffs >  m_RamData; //[idate]


};


/******************************************************************
**
**				Ram Fit Class
**
******************************************************************/


class RamFit: public LSGRGSolver,
			  public MLRamVolData
{		
			
/*		vector<CVector> m_bidVols;//[idate]
		vector<CVector> m_askVols;//[idate]
		vector<CVector> m_midVols;//[idate]
*/

protected:

		MlEqVolatilityStructureHandle	m_mktVols;	

		int m_totalNumVariables;		
		int m_currentidate;
		int m_fittingDimension;		
public:
		
		int m_fitBoundaries;
		int m_useVega;
		double m_fitSpread;
		double m_fitMidVols;

		RamFit(){};

		RamFit(MlEqVolatilityStructureHandle& mktVols,vector < vector< MlEqStrikeHandle  > >&	FittingStrikes	,vector < long > fittingDates,CMatrix& BoundsOnConstraints,int fitBoundaries,int useVega,double fitSpread,double fitMidVols,vector < vector < MlEqStrikeHandle > >& upCutoff,vector < vector < MlEqStrikeHandle > > & downCutoff);

		void initialize(MlEqVolatilityStructureHandle& mktVols,vector < vector< MlEqStrikeHandle  > >&	FittingStrikes	,vector < long > fittingDates,CMatrix& BoundsOnConstraints,int fitBoundaries,int useVega,double fitSpread,double fitMidVols,vector < vector < MlEqStrikeHandle > > & upCutoff,vector < vector < MlEqStrikeHandle > > & downCutoff);
		void  ObjectiveFcn(double* gVals,double* xVals);
		
		
};		
			

/******************************************************************
**
**				HermiteCoeff VolData Class
**
******************************************************************/
	

class HermiteCoeff
{
  public:

	double mat;
	double fwd;
	MlEqDateHandle nToday;
	long	maturity;

	CVector hCoeff;

	double getHermiteCoeff(int n);
};
		
void getHermiteVol(CVector& vols, GVector< MlEqStrikeHandle >& strikes, const HermiteCoeff& sc);
void getHermitePrice(CVector& callPrice,GVector< MlEqStrikeHandle >& strikes, const HermiteCoeff& sc);
double getHermiteVol(const  MlEqStrike& strike, const HermiteCoeff& sc);

		
class MLHermiteVolData: public MlEqVolData
{

public:

	MLHermiteVolData( vector<HermiteCoeff>& rm);
	MLHermiteVolData():MlEqVolData(){}
//	MLHermiteVolData(const MLHermiteVolData& rhs);
	MlEqVolDataHandle copy() const{ return new MLHermiteVolData( *this ); }
	MLHermiteVolData& operator=(const MLHermiteVolData& rhs);
	
	void initialize( vector<HermiteCoeff> & rm);

	const HermiteCoeff& getHermiteData(int whichSlice)const;

	virtual double getVol(const MlEqStrike& strike,int islice) ;

protected:

	vector < HermiteCoeff >  m_HermiteCoeff; //[idate]

};


/******************************************************************
**
**				Svl VolData Class
**
******************************************************************/


class SvlCoeffs
{	
public:

	double  t;//time as double in order to re-construct vols
	bool	bNewValid;		// true if the new (svl) parameterisation is valid.
	bool	bOldValid;		// true if the old (svi) parameterisation is valid

	SvlCoeffs();
	SvlCoeffs(CVector& x);

// new parameterisation

	double ATMVolF;	 
	double Skew;			
	double CallWing;
	double PutWing;
	double VolMin;

// old parameterisation

	double	a;
	double	b;
	double	m;
	double	rho;
	double	sig;



	
};	



double getSvlVar(double xstrike,const SvlCoeffs& svl);// applies svl-formula whatever strike definition is
double getSvlVol(double xstrike,const SvlCoeffs& svl,double t);// applies svl-formula whatever strike definition is
void getSvlVolDeriv(CVector& res,double xstrike,const SvlCoeffs& svl);



/******************************************************************
**
**				Svl Vol Data
**
******************************************************************/

class MlSvlVolData: public MlEqVolData
{	

public:
	
	static void initNewSvl(SvlCoeffs& svl, bool bThrow = true);
	static void initOldSvl(SvlCoeffs& svl, bool bThrow = true);


public:
	
	MlSvlVolData(){m_numSvlParams=5;};
	MlSvlVolData(const MlSvlVolData& rhs);
	MlEqVolDataHandle copy() const{ return new MlSvlVolData( *this ); }
	MlSvlVolData( vector< vector< MlEqStrikeHandle > > &strikes, vector<SvlCoeffs> & svl);
	void initialize(const vector< vector< MlEqStrikeHandle > > &strikes, vector<SvlCoeffs> & svl);	
	void initialize( ){m_numSvlParams=5;};

	int getNumberOfSlices()const { return m_numberOfSlices;};	
	virtual double getVol(const MlEqStrike& strike,int islice) ;
//	void getVolDeriv(CVector& res,const MlEqStrike& strike,int whichSlice) ;
	void getVolDeriv(CVector& res,const MlEqStrike& strike,int whichSlice, double mktForward) ;

	int m_numSvlParams;

	const 	vector<SvlCoeffs>&  getSvlParameters()const{return m_svl;} ;

	
protected:

	vector<SvlCoeffs> m_svl;//[idate]
	
};	

	
/******************************************************************
**
**				Svl Fit Class
**
******************************************************************/


class SvlFit: public LSGRGSolver,
			  public MlSvlVolData
{		
protected:

		vector<CVector> m_bidVols;//[idate]
		vector<CVector> m_askVols;//[idate]
		vector<CVector> m_midVols;//[idate]
		
		
		CVector m_doubleTime;
		CMatrix m_ideals;//[iparam][3]
		
		double getIdeal(int ifitParam,const double time,const double parmconst,const double parmlong,const double parmlambda);
		void initialize(CVector& initialGuess,vector < vector< MlEqStrikeHandle  > >&	strikes	,CVector& doubleTimes,double fitIdeals,double fitMidPoints,double fitSpread,CMatrix& BoundsOnConstraints,CVector& initialGuess2);
		
		int m_totalNumVariables;
		
		
public:
		
		SvlFit(){};
		SvlFit(vector<CVector>& bidVols,vector<CVector>& midVols,vector<CVector>& askVols,vector< MlEqStrikeHandle  >&	strikes	,CVector& doubleTimes,CVector& initialGuess,double fitIdeals,double fitMidPoints,double fitSpread);
		void initialize(vector<CVector>& bidVols,vector<CVector>& midVols,vector<CVector>& askVols,vector < vector< MlEqStrikeHandle  > >&	strikes	,CVector& doubleTimes,double fitIdeals,double fitMidPoints,double fitSpread,CMatrix& initial_Guess,CMatrix& BoundsOnConstraints);
		void initialize(CMatrix& vols,vector< MlEqStrikeHandle  >&	strikes	,CVector& doubleTimes,CVector& initialGuess,double fitIdeals,double fitMidPoints,double fitSpread);
		
		
		const vector<CVector>&	GetBidVols(void) const {return m_bidVols;}
		const vector<CVector>&  GetAskVols(void) const {return m_askVols;}
		const vector<CVector>&	GetMidVols(void) const {return m_midVols;}
		const CVector&          GetDoubleTime(void) const {return m_doubleTime;}


		double m_fitIdeals;
		double m_fitMidVols;
		double m_fitSpread;
		
		
		
		void  ConstraintFcn(double* contraintVals, double* xVals);
		void  PartialDerivatives(CVector& xVal,CVector& partialDerivatives){};
		void  ObjectiveFcn(double* gVals,double* xVals);
		
		
};		
			
		
/******************************************************************
**
**				AR (Adil Reghai) Prop VolData Class
**
******************************************************************/


	
class MlARPropVolData: public MlEqVolData
{
	
public:	

	MlARPropVolData(){};
	MlARPropVolData(const MlARPropVolData& rhs);

	MlARPropVolData(
								double	AtmVolLamda,
								double	AtmVolS0,
								double	AtmVolSInfty,
								double	SkewVolLambda,
								double	SkewVolS0,
								double  SkewVolInfty,	
								double	curvature,
								int		timeCutoff,
								CVector	fwds,
								MlEqDateHandle hDate,
								vector<long> dates	
								);

	void initialize(
								double	AtmVolLamda,
								double	AtmVolS0,
								double	AtmVolSInfty,
								double	SkewVolLambda,
								double	SkewVolS0,
								double  SkewVolInfty,	
								double	curvature,
								int		timeCutoff,
								CVector	fwds,
								MlEqDateHandle hDate, 
								vector<long> dates
								);


	MlEqVolDataHandle copy() const{ return new MlARPropVolData( *this ); }

	int getNumberOfSlices()const { return m_numberOfSlices;};
	virtual double getVol(const MlEqStrike& strike,int islice) ;
	
protected:
	
	double  S(double T);
	double  getAtmVol(double T);
	double	m_AtmVolLamda;
	double	m_AtmVolS0;
	double	m_AtmVolSInfty;
	double	m_SkewVolLambda;
	double	m_SkewVolS0;
	double  m_SkewVolInfty;	
	double	m_curvature;
	double  m_timeCutoff;
	CVector m_fwds;//[idate]
	CVector m_t;//[idate]
	CVector m_VolAtm;

};	





class MlARPropFit: public LSGRGSolver,
				   public MlARPropVolData
{		
			

protected:

		MlEqVolatilityStructureHandle	m_mktVols;	
		int m_totalNumVariables;		
		CVector m_VolAtmMkt;//[idate]
		CVector m_VolSkewMkt;

		int m_fitSkew;

		vector < vector< MlEqStrikeHandle  > >	m_FittingStrikes;

public:
		

		MlARPropFit(){};

		MlARPropFit(MlEqVolatilityStructureHandle& mktVols,vector < vector< MlEqStrikeHandle  > >&	FittingStrikes	,vector < long > fittingDates,CVector& fittingFwds, MlEqDateHandle hDate,CMatrix& BoundsOnConstraints,int timeCutoff,double curvature);

		void initialize(MlEqVolatilityStructureHandle& mktVols,vector < vector< MlEqStrikeHandle  > >&	FittingStrikes	,vector < long > fittingDates,CVector& fittingFwds,MlEqDateHandle hDate,CMatrix& BoundsOnConstraints,int timeCutoff,double curvature);
		void  ObjectiveFcn(double* gVals,double* xVals);
		
		
};		

		
/******************************************************************
**			
**				Multiplier Data Class
**
******************************************************************/
	
	
			
class  MlEqVolMultiplierData: public MlEqVolData
{	
		
public:
		
	MlEqVolMultiplierData(){};	
	MlEqVolMultiplierData(const MlEqVolMultiplierData& Data);
	MlEqVolDataHandle copy() const{ return new MlEqVolMultiplierData( *this ); }
	MlEqVolMultiplierData& operator=(const MlEqVolMultiplierData& rhs);
	
	int getNumberOfSlices() const					{return m_numberOfSlices;}
	double getReferenceVol(int islice) const;  // returns the atm vol for islice
	virtual double getVol(const MlEqStrike& strike,int islice) ;

	void			PutReferenceVol(CVector& refVols);

	const CVector&	GetReferenceVol(void) const		{return m_referenceVol;}
		
protected:
	
	CVector			m_referenceVol;//[islice]	

};	
	

#endif 
	