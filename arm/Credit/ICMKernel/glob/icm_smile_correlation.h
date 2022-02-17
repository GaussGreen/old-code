
#ifndef _ICM_SMILE_CORRELATION_H_
#define _ICM_SMILE_CORRELATION_H_

#include "ICMKernel/util/icm_matrix.h"
#include "ICMKernel/glob/icm_correlation.h"

/*********************************************************************************/
/*! \class  ICM_Smile_Correlation icm__smile_correlation.h "icm__smile_correlation.h"
 *  \author D Pouponneau
 *	\version 1.0
 *	\date   November 2004
 *	\brief Creates a correlation by strike object 
/***********************************************************************************/

 

class ICM_Smile_Correlation;
class ICM_Parameters; 
class ICM_ModelMultiCurves ;
class ICM_DefaultCurve; 
class ICM_Credit_Index; 
class ICM_VolInterpol; 

class TSBaseCorrelCallib
{
	
public :
	double itsPrice;
	ICM_Pricer* itsPricer;
	double itsMaturity;
	double itsK1;
	double itsK2;
	ICM_Smile_Correlation* itsCorrel;
	
	ICM_ModelMultiCurves* itsMmc;
	ARM_CLASS_NAME itsCname;
	ICM_Parameters  itsParams;

	ARM_Date its_Global_StartDate;
	ARM_Date its_Global_Maturity;
	// TypeCollateral its_Global_Collateral;
	std::vector<std::string> its_Global_Collateral; 
	int	its_Global_NBnames;
	bool its_Global_IncMatu;
	bool its_Global_AdjStart;
	int	its_Global_CreditLag;
	
	double its_Global_FeePV;
	double its_Global_DefPV;
	double its_Global_NPV;

	TSBaseCorrelCallib()
	{
		itsPrice = 0.;
		itsPricer = NULL;
		itsMaturity = 0.;
		itsK1 = 0.;
		itsK2 = 0.;
		itsCorrel = NULL;
		itsMmc = NULL;
		itsCname = ICM_PRICER_HOMOGENEOUS_SMILE; 
		its_Global_FeePV=0.;
		its_Global_DefPV=0.;
		its_Global_NPV=0.;
	}

	double ValuationCDOIndex(const double& x);

	ICM_QMatrix<double>* GenerateTSBaseCorrelation(ARM_ZeroCurve* ircurve,
							   ICM_DefaultCurve* defcurveindex,
							   ICM_VolInterpol& basecorrel,
							   int creditlag);
};


class _Slice_Correl: public ARM_Object
{
private:
	ARM_VolLInterpol* itsVolCurve;		//Correlation By Strike & Maturities

	double	itsProportion;				//Proportion for each Index
	ICM_Credit_Index* itsIndex;			//Vector of Index

public:
	vector<double> itsMaturities;//Vector of maturities
	vector<double> itsSmileStrikeLow;//Vector of Low strikes following maturity
	vector<double> itsSmileStrikeHigh;//Vector of Hight strikes following maturity

public:
	_Slice_Correl() {Init();}
	~_Slice_Correl() ;

	void Init() ; 

	_Slice_Correl(const _Slice_Correl& input) ; 

	_Slice_Correl& operator= (const _Slice_Correl& ref) 
	{if (this!=&ref) 
		{	this->~_Slice_Correl(); 
			new(this)_Slice_Correl(ref); }
		return *this;}

	inline ARM_VolLInterpol * GetVolCurve() { return itsVolCurve;}
	void SetVolCurve(const ARM_VolCurve* volcurve);


	inline ICM_Credit_Index* GetIndex() { return itsIndex;}
	void SetIndex(const ICM_Credit_Index* index);

	inline void SetProportion(const double& prop) 
	{
		itsProportion = prop;
	}
	inline double GetProportion(void) 
	{
		return itsProportion;
	}

	vector<double>& GetSmileStrikeLow() 
	{
		return itsSmileStrikeLow;
	}
	void SetSmileStrikeLow(const vector<double>& SmileStrikeLow) 
	{
		itsSmileStrikeLow=SmileStrikeLow;
	}
	void SetSmileStrikeLow(const double& SmileStrikeLow) 
	{
		for (int i=0; i<itsSmileStrikeLow.size();i++)
		{itsSmileStrikeLow[i]=SmileStrikeLow;}
	}

	vector<double>& GetSmileStrikeHigh() 
	{
		return itsSmileStrikeHigh;
	}
	void SetSmileStrikeHigh(const vector<double>& SmileStrikeHigh) 
	{
		itsSmileStrikeHigh=SmileStrikeHigh;
	}
	void SetSmileStrikeHigh(const double& SmileStrikeHigh) 
	{
		for (int i=0; i<itsSmileStrikeHigh.size();i++)
		{itsSmileStrikeHigh[i]=SmileStrikeHigh;}
	}

	vector<double>& GetMaturities() 
	{
		return itsMaturities;
	}
	void SetMaturities(const vector<double>& Maturities) 
	{
		itsMaturities=Maturities;
	}

};

class ICM_Mez;
class ICM_Smile_Correlation : public ICM_Correlation
{

private:
	
	qCorrel_By_Strike itsForcedStrikeType;   //Forced Strike Type
	vector<_Slice_Correl>	itsSlices;		 //Correlation By Currency

	//Calibration of equivalent strikes only
	ARM_Security*	its_Callib_PtfBespoke;
	ARM_Security*	its_Callib_Index;

	ARM_Model*		its_Callib_model;
	double			its_Callib_K1_PtfBespoke;
	double			its_Callib_K2_PtfBespoke;
	double			its_Callib_SizeEq_Index;
	int				its_Callib_NoIndex;
	double			its_Callib_ActiveMaturity;
	
	double			its_Callib_ELBespoke;
	double			its_Callib_ELIndex;

	bool			its_Already_rescal;		//Rescaling has already be done
	bool			its_never_rescal;		//Rescaling is allowed 
	qTERM_STRUCTURE	its_TermStructureRescaling;//Rescaling using Term Structure
	bool			its_Rescaling_Full;		//NOT ONLY ZEROCOUPON (default false --> ZEROCOUPON)

	bool			its_Interpolation_At_Maturity;  // optimisation only at maturity
	bool			its_NormalizeByEL;				// renormalisation par l'EL 
	bool			its_IndexJoinOptim;			// optimisation jointe sur les indices
	qRescalType		itsRescalType;
	double			its_BaseCorrelShift;

	
private: // copy const & = are disabled
	ICM_Smile_Correlation& operator= (const ICM_Smile_Correlation& );
	ICM_Smile_Correlation (const ICM_Smile_Correlation&) ;
private: 
	//Pricer Type : set in ComputeStrikeEq
	ARM_CLASS_NAME	its_PricerType;

	inline void Init(void) ;
	/*
	{
		ICM_Correlation::Init(); 
		SetName(ICM_SMILE_CORRMATRIX);
		itsForcedStrikeType = qStrike_NONE;

		its_Callib_PtfBespoke=NULL;
		its_Callib_Index=NULL;
		its_Callib_model=NULL;
		its_Callib_K1_PtfBespoke=0.;
		its_Callib_K2_PtfBespoke=0.;
		its_Callib_NoIndex=0;
		its_Callib_SizeEq_Index=0.;
		its_Callib_ActiveMaturity =0.;
		itsSlices.clear();
		its_Already_rescal = false;
		its_never_rescal = false;
		its_TermStructureRescaling = qNoTermStructure;
		its_PricerType = ICM_PRICER_HOMOGENEOUS_SMILE;
		its_Rescaling_Full = false;
		its_Interpolation_At_Maturity=false;

		its_Callib_ELBespoke=1.;
		its_Callib_ELIndex=1.;
		its_NormalizeByEL=false;
		its_IndexJoinOptim=false;
		itsRescalType = qRescal_Std;
		its_BaseCorrelShift=0.;
	}
	*/ 
	void Set(const ARM_Date& AsOf,
		const std::string& name, 
		const ARM_VolCurve** VolCurves, 
		const std::vector<std::string>& labels,
		const ARM_Vector* Proportions,
		const ARM_Vector* SmileStrikeLow ,
		const ARM_Vector* SmileStrikeHigh )  ; 

 
	void Set(const ARM_Date& AsOf,
					const std::string& name,  
					const ARM_VolCurve** VolCurves, 
					const std::vector<std::string>& labels,
					const ARM_Vector* Proportions,
					const ICM_Credit_Index** IndexVector = NULL) ; 
	/*
	{
		ICM_Correlation::Set(AsOf,labels,name,(ARM_IRIndex*)0,(ARM_IRIndex*)0);

		// SetName(ICM_SMILE_CORRMATRIX);

		for (int i=0;i<labels.size();i++)
		{
			_Slice_Correl S;
			S.SetProportion((*Proportions)[i]);
			S.SetVolCurve((VolCurves)[i]);
			S.SetIndex((IndexVector)[i]);
			itsSlices.push_back(S);
		}
	}*/

public: 

	ICM_Smile_Correlation() {Init();}		

	ICM_Smile_Correlation(const ARM_Date& AsOf,
						const string& name, 
						const ARM_VolCurve** matrix, 
						// char** label,
						const std::vector<std::string>& labels,
						const ARM_Vector* Proportions  ,
						const ARM_Vector* SmileStrikeLow ,
						const ARM_Vector* SmileStrikeHight  ) 
					
	{
		Init();

		Set(AsOf,name,matrix,labels,Proportions,SmileStrikeLow,SmileStrikeHight);
	};

	ICM_Smile_Correlation(const ARM_Date& AsOf,
						  const string& name, 
						  const ARM_VolCurve** matrix, 
						  const std::vector<std::string>& labels,
						  const ARM_Vector* Proportions = NULL,
						  const ICM_Credit_Index** IndexVector = NULL)
	{
		Init();

		Set(AsOf,name,matrix,labels,Proportions,IndexVector);
	};

	void Set(const ARM_Date& Asof,
				const std::string& name,
				const ARM_VolCurve** VolCurves, 
				const std::vector<std::string>& labels,
				const ARM_Vector& Proportions,
				const ICM_QMatrix<double>& fullStrikeLow,
				const ICM_QMatrix<double>& fullStrikeUp ) ;

	inline void SetAlready_rescal(bool Already_Rescal) {its_Already_rescal =  Already_Rescal;}

	inline void SetVolCurves(const ARM_VolCurve** volcurves, int size)
	{
		if (itsSlices.size() != size)
		{throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
		"ERROR : SetProportions : volcurves vector size <> slices size");}	

		for (int i=0;i<itsSlices.size();i++)
		{itsSlices[i].SetVolCurve((volcurves)[i]);}
	}		

	inline void SetVolCurves(const vector<ARM_VolCurve*>& volcurves)
	{
		if (itsSlices.size() != volcurves.size())
		{throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
		"ERROR : SetProportions : volcurves vector size <> slices size");}	

		for (int i=0;i<itsSlices.size();i++)
		{itsSlices[i].SetVolCurve((volcurves)[i]);}
	}		

	virtual std::vector<ARM_VolCurve*> GetVolCurves(); 

	virtual double GetCorrelation(const std::string&  issuer1,
								  const std::string&  issuer2,
								  double maturity = CREDIT_DEFAULT_VALUE,
								  double strike = CREDIT_DEFAULT_VALUE,
								  double actualYF = CREDIT_DEFAULT_VALUE) ;

	virtual double GetBeta(const std::string& issuer,
						   double maturity = CREDIT_DEFAULT_VALUE,
						   double strike = CREDIT_DEFAULT_VALUE,
						   double actualYF = CREDIT_DEFAULT_VALUE)
	{
		double correlation = GetCorrelation(issuer,issuer,maturity,strike,actualYF);
		correlation = NEG_SQRT(correlation);

		return (correlation);
	}

	virtual double GetCompositeCorrel(const std::string&  issuer1,
								  const std::string&  issuer2,
								  double maturity = CREDIT_DEFAULT_VALUE,
								  double strike = CREDIT_DEFAULT_VALUE,
								  double actualYF = CREDIT_DEFAULT_VALUE) ;


	virtual ARM_Vector* ComputeBetas(int nbissuers,const std::vector<std::string>& labels,const ARM_Vector&nominals, const ARM_Date &Maturity)
	{
		return NULL;
	}

	virtual ICM_QMatrix<double> ComputeCholeskyMatrix()
	{
		ICM_QMatrix<double> matrix(1,1,0.);
		return matrix;
	}

	inline void SetPricerType(ARM_CLASS_NAME pricerName)
	{
		its_PricerType = pricerName;
	}

	inline ARM_CLASS_NAME GetPricerType(void)
	{
		return its_PricerType;
	}

	~ICM_Smile_Correlation() 
	{
	itsSlices.clear();
	};

	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object* src) ;
	/** 
	{
	    ICM_Smile_Correlation* CorrMatrix = (ICM_Smile_Correlation *) src;

		int initial_size = CorrMatrix->GetSize();

		itsSlices = CorrMatrix->itsSlices;
		its_Already_rescal = CorrMatrix->its_Already_rescal;
		its_never_rescal = CorrMatrix->its_never_rescal;
		its_PricerType = CorrMatrix->its_PricerType;
		its_TermStructureRescaling = CorrMatrix->its_TermStructureRescaling;
		its_Rescaling_Full = CorrMatrix->its_Rescaling_Full;
		its_NormalizeByEL = CorrMatrix->its_NormalizeByEL;
		its_IndexJoinOptim = CorrMatrix->its_IndexJoinOptim;
		itsRescalType = CorrMatrix->itsRescalType;
		its_BaseCorrelShift = CorrMatrix->its_BaseCorrelShift;
	}
	**/ 
	inline void SetTermStructureRescaling(const qTERM_STRUCTURE& value) {its_TermStructureRescaling=value;}
	inline qTERM_STRUCTURE GetTermStructureRescaling() {return its_TermStructureRescaling;}

	// -------------
	//	Copy Method 
	// -------------
	void Copy(const ARM_Object* src)
	{
		ICM_Correlation::Copy(src);
 		BitwiseCopy(src);
	}

	// --------------
	//	Clone Method
	// --------------
	virtual ARM_Object* Clone(void)
	{
     ICM_Smile_Correlation* theClone = new ICM_Smile_Correlation();
     theClone->Copy(this);
     return(theClone);
	}

	virtual void ResetBetas() ;
	/*
	{
		int size = itsSlices.size();
		if (its_never_rescal) return;

		for (int i=0; i<size; i++)
		{
			itsSlices[i].GetSmileStrikeLow().clear();
			itsSlices[i].GetSmileStrikeHigh().clear();
			itsSlices[i].GetMaturities().clear();
		}

		its_Already_rescal = false;
	}
	*/ 
	void View(char* id, FILE* ficOut);

	virtual ICM_Correlation* GenerateShiftCorrel(const std::string& label,
												qSENSITIVITY_TYPE typesensi,
												double epsilon = CREDIT_DEFAULT_VALUE) ;


	virtual ICM_Correlation* GenerateShiftBetas(const std::string& label,
												qSENSITIVITY_TYPE typesensi,
												double epsilon /** = CREDIT_DEFAULT_VALUE **/ )
	{
		ICM_Correlation* Correlation = GenerateShiftCorrel(label,
													typesensi,
													epsilon);
		return (Correlation);
	}

	qCorrel_By_Strike GetForcedStrikeType() { return itsForcedStrikeType; }

	virtual void SetForcedStrikeType(qCorrel_By_Strike type)
	{
		itsForcedStrikeType = type;
	}

	virtual void ComputeStrikesEq(ICM_Pricer* pricer,const qRescalType& rescal /*= qRescal_Std_Maturity*/);
	//Rescaling standart
	virtual void ComputeStrikesEq_standart(ICM_Pricer* pricer);

	//Term Structure Rescaling
	void ComputeStrikesEq_Term(ICM_Pricer* pricer);

	double DiffSensis(double Strike_down);
	double DiffSensisEquity(double Strike_up);
	double DiffSensisEquityTerm(double Strike_up);

	//Term Structure DiffSensis
	double DiffSensisTerms(double Strike_down);

	//Rescaling Equity
	virtual void ComputeStrikesEq_Equity(ICM_Pricer* pricer);
	double DiffSensisEqDown(double Strike_down);
	double DiffSensisEqDownTerm(double Strike_down);
	double DiffSensisEqUp(double Strike_up);
	double DiffSensisEqUpTerm(double Strike_up);

	//digital rescaling
	double DiffSensisEqUp_digit(double Strike_up);
	double DiffSensisEqDown_digit(double Strike_down);
	
	void ComputeStrikesEq_Equity_Term(ICM_Pricer* pricer);

	//Rescaling ELoss
	virtual void ComputeStrikesEq_ELoss(ICM_Pricer* pricer);
	virtual void ComputeStrikesEq_ELoss_Term(ICM_Pricer* pricer);

	//Rescaling Equity combined
	void ComputeStrikesEq_Equity_CombinIDX(ICM_Pricer* pricer);

	//Rescaling ING
	virtual void ComputeStrikes_ING(ICM_Pricer* pricer);

	ICM_Mez* GenerateCdoIndex(ARM_Date& Start,
												 ARM_Date& Maturity,
												 const int& NoIndex,
												 const double& PorS);

	//Equivalent Strike Down
	virtual double GetEqStrikeDown(const std::string& indexname);

	//Equivalent Strike Up
	virtual double GetEqStrikeUp(const std::string& indexname);

	//Equivalent Correl Strike Down
	virtual double GetCorrelStrikeDown(double maturity);

	//Equivalent Correl Strike Up
	virtual double GetCorrelStrikeUp(double maturity);

	//Equivalent Smile Strike Down
	virtual void GetSmileStrikeDown(const std::string& indexname, vector<double>& vMatu, vector<double>& vStrikes);

	//Equivalent Smile Strike Up
	virtual void GetSmileStrikeUp(const std::string& indexname, vector<double>& vMatu, vector<double>& vStrikes);

	void Find(const std::string& info,int& nthLine,int& nthCol,int& index); 

	virtual void SetProportionsInfos(const std::string& indexname,
							 const double& proportion,
							 const double& forcedstrikelow = CREDIT_DEFAULT_VALUE,
							 const double& forcedstrikehigh = CREDIT_DEFAULT_VALUE) ;
	/*
	{
		int i =0;

		i = GetLabelNo(indexname);
		itsSlices[i].SetProportion(proportion);
		
		if ((itsSlices[i].GetSmileStrikeLow().size()>0) && (forcedstrikelow != CREDIT_DEFAULT_VALUE))
			itsSlices[i].SetSmileStrikeLow(forcedstrikelow);
		else if (itsSlices[i].GetSmileStrikeLow().size()==0)
		{
			vector<double> SmileStrikeLow;SmileStrikeLow.push_back(forcedstrikelow);
			itsSlices[i].SetSmileStrikeLow(SmileStrikeLow);
		}


		if ((itsSlices[i].GetSmileStrikeHigh().size()>0) && (forcedstrikehigh != CREDIT_DEFAULT_VALUE))
			itsSlices[i].SetSmileStrikeHigh(forcedstrikehigh);
		else if (itsSlices[i].GetSmileStrikeHigh().size()==0)
		{
			vector<double> SmileStrikeHigh;SmileStrikeHigh.push_back(forcedstrikehigh);
			itsSlices[i].SetSmileStrikeHigh(SmileStrikeHigh);
		}

	}
	*/ 

	inline void SetProportions(ARM_Vector* vector) 
	{
		if (vector->GetSize() != itsSlices.size())
		{throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
		"ERROR : SetProportions : proportion vector size <> slices size");}	

		for (int il=0;il<itsSlices.size();il++)
		{itsSlices[il].SetProportion(vector->Elt(il));}
	}

	inline void GetProportions(ARM_Vector*& vector)
	{
		if (vector) delete vector;
		vector = new ARM_Vector(itsSlices.size(),0.);
		for (int il=0;il<itsSlices.size();il++)
		{vector->Elt(il) = itsSlices[il].GetProportion();}
	}

	virtual int GetIndexSize() { return itsSlices.size();}

	inline vector<_Slice_Correl>& GetSlices() {return itsSlices;}
	inline void SetSlices(vector<_Slice_Correl>& vector) {itsSlices=vector;}

	void SetCorrelations(ICM_Smile_Correlation& correl) ;

	virtual void GetCorrelationTerms(ARM_Vector& vector);

	virtual void SetInterpType(const int& interpType) ;

	virtual int GetInterpType() ; 
	

	inline void KeepOneSlice(int i)
	{
		_Slice_Correl Slice = itsSlices[i];
		itsSlices.clear();
		itsSlices.push_back(Slice);
	}
};


ICM_Smile_Correlation* FixedBaseCorrelation(const ARM_Date& AsOf,
											ARM_Currency* Ccy,
											const double& strikeDown,
											const double& strikeUp,
											const double& CorrelDown,
											const double& CorrelUp,
											const string& IndexName,
											double yfmaturity = 5.);

ICM_Smile_Correlation* FixedBaseCorrelationWithExistingBC(ARM_VolCurve* volBC,	
														  const ARM_Date&	   AsOf,
														  ARM_Currency*		    Ccy,
														  const double&	strikeDown,
														  const double& strikeUp,
														  const double& CorrelDown,
														  const double& CorrelUp,
														  const string& IndexName,
														  const double& yfmaturity);

ICM_Smile_Correlation* FixedBaseCorrelationMult(const ARM_Date& AsOf,
											vector<ARM_Currency*> Ccy,
											const vector<double>& strikeDown,
											const vector<double>& strikeUp,
											const vector<double>& CorrelDown,
											const vector<double>& CorrelUp,
											const vector<string>& IndexName,
											const vector<double>& Proportion);

void FixedBaseCorrelationMult_Vector(const ARM_Date& AsOf,
											const vector<ARM_Currency*>& Ccy,
											const vector<double>& strikeDown,
											const vector<double>& strikeUp,
											const vector<double>& CorrelDown,
											const vector<double>& CorrelUp,
											const vector<string>& IndexName,
											vector<ICM_Smile_Correlation*>& Correlation);

ICM_Smile_Correlation* QuickBaseCorrelation(const ARM_Date& AsOf,
										    ARM_Currency* Ccy,
										    const double&	strikeDown,
										    const double& strikeUp,
										    const double& CorrelDown,
										    const double& CorrelUp,
										    const string& IndexName,
										    const double& yfmaturity,
											ARM_Vector* YearTerms,
											ARM_Vector* Strikes,
											ARM_Matrix* BC);

ICM_Smile_Correlation* FixedBaseCorrelation2M(const ARM_Date& AsOf,
											ARM_Currency* Ccy,
											const double& strikeDown,
											const double& strikeUp,
											const double& CorrelDown,
											const double& CorrelUp,
											const double& Maturity,
											const double& PrevCorrelDown,
											const double& PrevCorrelUp,
											const double& PrevMaturity,
											const string& IndexName);

void UpdateSingleCorrelation(ICM_Smile_Correlation* correl,
							 const double& strike_down,
						     const double& strike_up,
							 const double& yt,
							 const double& value,
							 bool& succeed);

ICM_Smile_Correlation* FixedBaseCorrelationEmpty(const ARM_Date& AsOf,
												ARM_Currency* Ccy,
												ARM_Vector& strikes,
												ARM_Vector& yearterms,
												ARM_Matrix& matrix);



// ----------------------------------------------------
// These classes are used for join optimization 
// ----------------------------------------------------
class _context_index
{
public:
	vector<ICM_Mez*>		m_Callib_Index;
	vector<double>			m_Callib_ELIndex;
	vector<double>			m_Callib_ActiveMaturity;

	_context_index() {reset();}

	void reset()
	{
	m_Callib_Index.clear();
	m_Callib_ELIndex.clear();
	m_Callib_ActiveMaturity.clear();
	}

	_context_index(const _context_index& in)
	{
	m_Callib_Index = in.m_Callib_Index;
	m_Callib_ELIndex = in.m_Callib_ELIndex;
	m_Callib_ActiveMaturity = in.m_Callib_ActiveMaturity;
	}

	~_context_index()
	{

	}
};


class _context
{
public:
	ICM_Mez*				m_Callib_PtfBespoke;
	ARM_Model*				m_Callib_model;
	ICM_Smile_Correlation*  m_Active_correl;
	vector<_Slice_Correl>*  m_Slices;

	vector<_context_index>	m_Callib_Index_Vector;

	double					m_Callib_ELBespoke;
	double					m_Callib_K1_PtfBespoke;
	double					m_Callib_K2_PtfBespoke;
	
	int						m_Callib_MaturityNo;
	ARM_CLASS_NAME			m_PricerType;

	_context(const _context& in)
	{
		m_Callib_PtfBespoke=in.m_Callib_PtfBespoke;
		m_Callib_model=in.m_Callib_model;
		m_Active_correl=in.m_Active_correl;
		m_Slices=in.m_Slices;
		m_Callib_Index_Vector=in.m_Callib_Index_Vector;
		m_Callib_ELBespoke=in.m_Callib_ELBespoke;
		m_Callib_K1_PtfBespoke=in.m_Callib_K1_PtfBespoke;
		m_Callib_K2_PtfBespoke=in.m_Callib_K2_PtfBespoke;
		m_PricerType=in.m_PricerType;
		m_Callib_MaturityNo = in.m_Callib_MaturityNo;
	}

	void reset()
	{
	m_Callib_PtfBespoke=NULL;
	m_Callib_model=NULL;
	m_Active_correl=NULL;
	m_Slices=NULL;
	m_Callib_Index_Vector.clear();

	m_Callib_ELBespoke=1.;
	m_Callib_K1_PtfBespoke=-1.;
	m_Callib_K2_PtfBespoke=-1.;
	m_PricerType=ARM_OBJECT;
	m_Callib_MaturityNo = 0;
	}

	_context()
	{reset();}
};

#endif
