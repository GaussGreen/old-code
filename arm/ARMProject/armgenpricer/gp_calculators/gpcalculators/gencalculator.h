
#ifndef _INGPCALCULATORS_GENCALCULATOR_H
#define _INGPCALCULATORS_GENCALCULATOR_H

///gpbase
#include "gpbase/port.h"
#include "gpbase/rootobject.h"

/// gpinfra
#include "gpinfra/gramfunctorargdict.h"
#include "gpinfra/mktdatamanagerrep.h"

/// gpcalib
#include "gpcalib/typedef.h"

///kernel
#include <ccy/currency.h>
#include <util/refvalue.h>
//#include <inst/forex.h>
////#include <crv/fixingsched.h>

class ARM_SwapLeg;
CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_PricingModel;
class ARM_GenSecurity;
class ARM_CalibMethod;
class ARM_DateStripCombiner;
struct ARM_MarketData_ManagerRep;
class ARM_Warning;

/// this is the base class for all the calculators based on the generic pricer
/// Each calculator shares:
///		a) a deal description encapsulated in the generic security
///		b) calibration methods with its associated portfolios.. (autocalibration)
///		c) a propser model with its numerical method
///	Each calculator knows
///		a) where to find the market data
///		b) where to find the pricing model for the calibration model to work

/// This class shares all the common method to all calculators!

class ARM_GenCalculator : public ARM_RootObject
{
public:
	enum DegeneratedCalculator
	{
		No = 0,
		CRAFromCCSO,
		SwaptionFromCCSO,
	};

	/// constructor, copy constructor, assignment operator and destructor
	ARM_GenCalculator(const ARM_Date& asOfDate=ARM_Date());
	ARM_GenCalculator( const ARM_MarketData_ManagerRep& mktDataManager );
	ARM_GenCalculator( const ARM_GenCalculator& rhs );
	ARM_GenCalculator& operator=( const ARM_GenCalculator& rhs );
	virtual ~ARM_GenCalculator();

	/// accessors
    inline const ARM_GenSecurityPtr GetGenSecurity() const              {   return itsGenSecurity;      }
    inline ARM_GenSecurityPtr GetGenSecurity()                          {   return itsGenSecurity;      }   // for non const functions
	inline const ARM_CalibMethodPtr GetCalibMethod() const              {   return itsCalibMethod;      }
	inline ARM_CalibMethodPtr GetCalibMethod()                          {   return itsCalibMethod;      }  // for non const functions
    inline const ARM_PricingModelPtr GetPricingModel() const            {   return itsPricingModel;     }
	inline ARM_PricingModelPtr GetPricingModel()                        {   return itsPricingModel;     } // for non const functions
	inline const ARM_MarketData_ManagerRepPtr GetMktDataManager() const {   return itsMktDataManager;   }
    inline ARM_MarketData_ManagerRepPtr GetMktDataManager()             {   return itsMktDataManager;   }  // for non const functions
    inline const ARM_StringVector& GetKeys() const                      {   return itsMDMKeys;          }
    inline void SetKeys(const ARM_StringVector& keys)                   {   itsMDMKeys = keys;          }
	inline const ARM_Currency* const GetCurrencyUnit() const            {   return itsCurrency;         }
	inline ARM_Currency* GetCurrencyUnit()                              {   return itsCurrency;         }
    void SetCurrencyUnit( ARM_Currency* ccy );
    inline ARM_Currency GetDomesticCcy() const                          {	return itsDomesticCcy;	    }
    inline void SetDomesticCcy(ARM_Currency& DomesticCcy)               {   itsDomesticCcy = DomesticCcy;	    }
    inline ARM_Currency GetForeignCcy() const                           {	return itsForeignCcy;	    }
    inline void SetForeignCcy(ARM_Currency& ForeignCcy)                 {	itsForeignCcy = ForeignCcy;	        }
    inline ARM_Currency GetFundingCcy() const                           {	return(itsFundingCcy);	    }
    inline void SetFundingCcy(ARM_Currency& FundingCcy)                 {	itsFundingCcy = FundingCcy;	        }
    inline ARM_Currency GetBasisCcy() const                             {   return itsBasisCcy;         }
    inline void SetBasisCcy(ARM_Currency& BasisCcy)                     {   itsBasisCcy = BasisCcy;             }
    inline int GetDegeneratedCalculator() const							{   return itsDegeneratedCalculator;         }
    inline void SetDegeneratedCalculator(int degeneratedCalc)           {   itsDegeneratedCalculator = degeneratedCalc; }

    /// to test if basis
    inline bool IsBasis() const {	return false/*( string(itsFundingCcy.GetCcyName()) != string(itsCurrency->GetCcyName()) )*/;}

	ARM_Date GetDateFromTimeAndAsof( double d, double asOfDate ) const { return ARM_Date( d+asOfDate); }
	double GetTimeFromDateAndAsof( const ARM_Date& d, double asOfDate) const { return d.GetJulian()-asOfDate; }

	inline bool IsFullTerminated() const {return itsFullTerminated;} 

	/// common functions
    void ValidateGenSecurity( const ARM_GenSecurityPtr& genSecurity );
    void SetGenSecurity( const ARM_GenSecurityPtr& genSecurity );
    void SetGenSecurityWithoutValidation( const ARM_GenSecurityPtr& genSecurity ) { itsGenSecurity = genSecurity; }
	void SetCalibMethod( const ARM_CalibMethodPtr& CalibMethod );
	void SetPricingModel( const ARM_PricingModelPtr& pricingModel );
    void SetMktDataManager(const ARM_MarketData_ManagerRepPtr& mktDataManager);

    inline void SetFullterminated(bool fullTerminated){itsFullTerminated = fullTerminated;} 

   	inline ARM_MultiTypeDict& GetPricingData() const { return itsPricingData; }

	/// General method to feed a line of a deal description
	void FillRowInfo( const ARM_RowInfo& rowInfo,
		CC_NS(std,vector)<string>::iterator& textIter, 
		CC_NS(std,vector)<ARM_GP_VALUE_TYPE>::iterator& formatIter );

	/// To intialize GenCalculator  and update GenSec 
	/// with a given Mkt Data Manager 
	 void Initialize(ARM_MarketData_ManagerRep* mktDataManager);
	 void UpDateWithSubMktDataManager(const ARM_MarketData_ManagerRep& mktDataManager);
     void UpDateMktDataManagerOnly(const ARM_MarketData_ManagerRep& mktDataManager);
	 void InitializeMktDataManagerOnly(const ARM_MarketData_ManagerRep& mktDataManager);
    
	/// Fct to build the generic security of the calculator
	void CreateAndSetDealDescription(const string& payModelName="", const ARM_StringVector& PricedColumns = ARM_StringVector(0), ARM_CstManagerPtr cstManager = ARM_CstManagerPtr(NULL), bool fixBoundary=false, bool otherPayoffs=false);
	void CreateAndSetCustomDealDescription(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0),const string& payModelName="", const ARM_StringVector& PricedColumns = ARM_StringVector(0), ARM_CstManagerPtr cstManager = ARM_CstManagerPtr(NULL), bool fixBoundary=false, bool otherPayoffs=false);
	
    /// Fct to update internal datsa depending on the market data manager
    void Update();
	void Update(const ARM_MarketData_ManagerRep *mktDataMgr);
	
	/// function with timer
	void CreateAndSetCalibrationAndTimeIt();
	void CreateAndSetModelAndTimeIt();
	double PriceAndTimeIt();
	void CalibrateAndTimeIt();
    void UpdateModelAndTimeIt();
	void UpdateCalibrationAndTimeIt(bool isUpdateStrike=true);
	void CreateAndSetDealDescriptionAndTimeIt(const string& payModelName="", const ARM_StringVector& PricedColumns = ARM_StringVector(0), ARM_CstManagerPtr cstManager = ARM_CstManagerPtr(NULL), bool fixBoundary = false, bool otherPayoffs=false);
	void CreateAndSetCustomDealDescriptionAndTimeIt(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0),const string& payModelName="", const ARM_StringVector& PricedColumns = ARM_StringVector(0), ARM_CstManagerPtr cstManager = ARM_CstManagerPtr(NULL), bool fixBoundary = false, bool otherPayoffs=false);
	void CheckDataAndTimeIt();
	void CheckMktDataAndTimeIt();

    virtual ARM_GenSecurity* GetSubGenSecurity(int indexType) const;

	/// ----------------------------------------------------------
	/// ------------- pure virtual functions
	/// pure virtual methods to force implementation
	virtual ARM_RowInfo ColumnNames() const = 0;
	virtual ARM_RowInfo MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const = 0;
	virtual ARM_DateStripCombiner DatesStructure() const = 0;
	virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const = 0;

	/// pricing function
	virtual void ComputePricingData() const= 0;
	virtual void CreateAndSetModel() = 0;
	virtual void CreateAndSetCalibration() = 0;
    virtual void UpdateModel() = 0;
    virtual void UpdateCalibration(bool isUpdateStrike=true) = 0;
	virtual void Calibrate() = 0;
    virtual double Price() = 0;
    virtual void CheckData() = 0; /// Check internal data consistency
	virtual void CheckMktData() = 0; // Check mkt data consistency
	virtual ARM_Vector* ComputeAll() { return NULL;}; /// Underlying computing
	

    /// Initialisation & updating dedicated to Summit interface
    void Init(const vector< ARM_Object* >& marketDatas, bool toCalibrate = true);
    void Update(const vector< ARM_Object* >& marketDatas);
	
	
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual ARM_CLASS_NAME GetRootName() { return ARM_GENCALCULATOR; }

	// ONLY FOR SUMMIT PARSING
	inline void SetPorS(int PorS) { itsPorS = (PorS == -1) ? -1 : 1; }
	inline int  GetPorS()         { return itsPorS; }

	/// WARNING
	inline ARM_Warning* GetWarning() const { return itsWarning; }

	/// Dates Strip
	inline virtual ARM_DateStripPtr GetOriginalFundDateStrip() const  = 0;
	inline virtual ARM_DateStripPtr GetRefFundDateStrip() const = 0;
	inline virtual ARM_DateStripPtr GetRefDateStrip() const = 0;

	/// Vector
	inline virtual std::vector<double> GetvCpnNominal() const  = 0;
	inline virtual std::vector<double> GetvFundNominal() const = 0;
	inline virtual std::vector<double> GetvFundSpread() const = 0;

	/// the discount curve
	inline virtual ARM_ZeroCurve* GetDomesticZeroCurve() const  = 0;
	inline virtual ARM_ZeroCurve* GetForeignZeroCurve() const  = 0;
	inline virtual ARM_ZeroCurve* GetDomesticDiscountZeroCurve() const  = 0;
	inline virtual ARM_ZeroCurve* GetForeignDiscountZeroCurve() const  = 0;
	inline virtual ARM_Forex* GetForex() const  = 0;
	inline virtual ARM_SwapLeg* GetFundingLeg() const  { return NULL;};

	std::vector<double> ComputeDomesticBasis() const ;

    // As Of payment flag management

    inline void SetDiscPricingMode(int discPricingMode)  { itsAccountingPricingFlag = discPricingMode;  }
    inline int GetDiscPricingMode() {  return(itsAccountingPricingFlag);  }
   inline  ARM_FixingSched* GetFixings() { return(itsLiborFxFixings); }

    inline void SetFixings(ARM_FixingSched* theFixings)
    {
		delete itsLiborFxFixings;
		itsLiborFxFixings = NULL;
		//itsLiborFxFixings = theFixings ? (ARM_FixingSched *) theFixings->Clone() : NULL;
    }

private:

    // How to price an As of cash flow flag

    int itsAccountingPricingFlag; // Default: ARM_DISC_PRICING_METH else
                                  //          ARM_DISC_ACCOUNTING_METH

    // Fixings management
    ARM_FixingSched*                itsLiborFxFixings;
    
    ARM_GenSecurityPtr				itsGenSecurity;
	ARM_CalibMethodPtr              itsCalibMethod;
	ARM_PricingModelPtr             itsPricingModel;
    ARM_MarketData_ManagerRepPtr    itsMktDataManager;
	mutable ARM_MultiTypeDict		itsPricingData;

	bool itsFullTerminated;
	bool itsUpdate;

	int itsDegeneratedCalculator;

    /// Market data selector keys 
    ARM_StringVector				itsMDMKeys;

    /// coupon, funding, basis, domestic & foreign currencies
    /// (not known from the MdM because may not be already built)
    ARM_Currency*  itsCurrency;
    ARM_Currency   itsFundingCcy;
    ARM_Currency   itsBasisCcy;
    ARM_Currency   itsDomesticCcy;
    ARM_Currency   itsForeignCcy;

    void CopyNoCleanUp(const ARM_GenCalculator& rhs);

	typedef void (ARM_GenCalculator::*GenCalcVoidToVoidFunc)();
	typedef double (ARM_GenCalculator::*GenCalcDbleToVoidFunc)();
	void CallFuncAndTimeIt( const GenCalcVoidToVoidFunc& func, const string& funcName );
	double CallFuncAndTimeIt( const GenCalcDbleToVoidFunc& func, const string& funcName );

	// ONLY FOR SUMMIT PARSING : Purchase or Sale flag
	int itsPorS;

	ARM_Warning* itsWarning;

    /// Initialisation to 0 of all columns of the deal description that could be priced
    virtual void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const = 0;

};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

