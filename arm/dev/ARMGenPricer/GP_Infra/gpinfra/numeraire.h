/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file numeraire.h
 *
 *  \brief 
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_NUMERAIRE_H
#define _INGPINFRA_NUMERAIRE_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"

#include "gpbase/rootobject.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h" /// for accessor to the ARM_GP_Vector
#include "gpbase/assignop.h"

#include "timeinfo.h"

#include <map>
using std::map;


CC_BEGIN_NAMESPACE( ARM ) /// macro for namespace ... define namespace only if supported

struct ARM_ZeroCurveFunctor;

//////////////////////////////////////////////
/// \class ARM_Numeraire
/// \brief
/// ARM_Numeraire class specifies the numeraire
/// of a diffusion depending of its type :
///
/// Cash : numeraire = 1.0
///
/// For TerminalZc types a list of
/// N maturities is used (Ti) i=0 to N-1
/// to define a numeraire set Zc(t,Tn(t)) where
/// Tn(t)-1 <= t < Tn(t) :
///
/// TerminalZc : numeraire = Zc(t,TN-1)
//////////////////////////////////////////////
class ARM_Numeraire : public ARM_RootObject
{
public:
    enum NumeraireType
    {
        Cash=0,
        TerminalZc,			/// standard terminal zero coupon with respect to the underlying!
		TerminalEventZc,	/// terminal zero coupon measure with respect to the last event date!
		RollingPayment,		/// rolling on the payment date
		RollingEvent,		/// rolling on the event date
		PowerTerminalZc,
        RollingCash,        /// rolling w.r.t. the numerical method schedule (<=> discrete cash)
        Unknown
    };

	/// text corresponding to the type enum
	static const string NumeraireTypeTable[];

	typedef map<double,ARM_VectorPtr> NumDiscountMap;

/// for easy access
protected:
    /// Maturities (Ti) of the Zc numeraires
    ARM_GP_Vector*     itsTimes;

    /// Maturity index w.r.t the above schedule
    int itsMatIdx;  // n(t)

	/// Spot values of the numeraires : Numeraire(0,Ti)
    double itsValues0;

	mutable NumDiscountMap itsNumDiscountMap;

    void CopyNoCleanUp(const ARM_Numeraire& rhs);
    void CleanUp();

public:
	/// constructor, copy constructor, operator=, destructor
	ARM_Numeraire();
	ARM_Numeraire( const ARM_Numeraire& rhs );
	virtual ~ARM_Numeraire();
    ARM_Numeraire& operator=(const ARM_Numeraire& rhs);

    /// Accessors
    double GetMaturity() const {return (*itsTimes)[itsMatIdx];}
    double GetValue0() const {return itsValues0;}

    /// ------------ pricing part
	/// init part
    virtual void Init(const ARM_ZeroCurveFunctor& discFunctor,const string& modelName,const ARM_TimeInfoPtrVector& timeInfos, const ARM_VectorPtr& numeraireTimes = ARM_VectorPtr(NULL)) = 0;
    virtual void Reset( int pricingDir);
	virtual void ResetLoop();
	void ResetNumDiscountMap();
	virtual void FinalizeInduction(	double evalTime, const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) {};

	// Move the numeraire time
	virtual void MoveNumeraireFwd();
	virtual void MoveNumeraireBckwd();
	
	// Update the discount information (really use for Cash & Rolling numeraire)
	virtual void Update(const ARM_PricingModel& model, const ARM_PricingStatesPtr& states, size_t timeIdx ) = 0;
	virtual void DefaultUpdate(const ARM_PricingModel& model, const ARM_PricingStatesPtr& states) {};

	/// processing for payoff
	virtual void ProcessPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const= 0;
	virtual void ProcessUnPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const= 0;
	virtual bool NeedLocalDiscount() const = 0;
    virtual NumeraireType GetType() const = 0;

	/// Default processing for payoff
	virtual void DefaultProcessPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime,
		const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const;
	virtual void DefaultProcessUnPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const;

	/// ARM_Object part
    virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual string NumeraireTypeName() const = 0;
};



/// Fwd numeraire
class ARM_NumeraireFwd : public ARM_Numeraire
{
public:
	/// constructor, copy constructor, operator=, destructor
	ARM_NumeraireFwd() : ARM_Numeraire() {}
	ARM_NumeraireFwd( const ARM_NumeraireFwd& rhs ) : ARM_Numeraire( rhs) {}
	virtual ~ARM_NumeraireFwd() {};
    ARM_NumeraireFwd& operator=(const ARM_NumeraireFwd& rhs)
	{
		if( this!= &rhs )
			ARM_Numeraire::operator=(rhs);
		return *this;
	}

	virtual void Update(const ARM_PricingModel& model, const ARM_PricingStatesPtr& states, size_t timeIdx)
	{
		DefaultUpdate(model,states);
	}

	/// process payoff and unprocess payoff
	virtual void ProcessPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const
	{
		DefaultProcessPaidPayoffs(
			payModelName, 
			payoffs,
			evalTime, 
			states,
			model);
	}
	virtual void ProcessUnPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const
	{
		DefaultProcessUnPaidPayoffs(
			payModelName,
			payoffs,
			evalTime, 
			states,
			model );
	}

	virtual bool NeedLocalDiscount() const { return false; }

};


class ARM_NumeraireTerminalZc : public ARM_NumeraireFwd
{
public:
	/// constructor, copy constructor, operator=, destructor
	ARM_NumeraireTerminalZc() : ARM_NumeraireFwd() {}
	ARM_NumeraireTerminalZc( const ARM_NumeraireFwd& rhs ) : ARM_NumeraireFwd( rhs) {}
	virtual ~ARM_NumeraireTerminalZc() {};
    ARM_NumeraireTerminalZc& operator=(const ARM_NumeraireTerminalZc& rhs)
	{
		if( this!= &rhs )
			ARM_NumeraireFwd::operator=(rhs);
		return *this;
	}

	/// initialisation
    virtual void Init(const ARM_ZeroCurveFunctor& discFunctor,const string& modelName,const ARM_TimeInfoPtrVector& timeInfos, const ARM_VectorPtr& numeraireTimes = ARM_VectorPtr(NULL));
    virtual NumeraireType GetType() const { return ARM_Numeraire::TerminalZc; }

	/// Standard ARM Support
	virtual ARM_Object* Clone() const { return new ARM_NumeraireTerminalZc(*this); }
	virtual string NumeraireTypeName() const { return "TerminalZc"; }
};

class ARM_NumerairePowerTerminalZc : public ARM_NumeraireFwd
{
public:
	/// constructor, copy constructor, operator=, destructor
	ARM_NumerairePowerTerminalZc( double Alpha = 1 ) : itsAlpha(Alpha), ARM_NumeraireFwd() {}
	ARM_NumerairePowerTerminalZc( const ARM_NumeraireFwd& rhs ) : ARM_NumeraireFwd( rhs) {}
	virtual ~ARM_NumerairePowerTerminalZc() {};
    ARM_NumerairePowerTerminalZc& operator=(const ARM_NumerairePowerTerminalZc& rhs)
	{
		if( this!= &rhs )
			ARM_NumeraireFwd::operator=(rhs);
		return *this;
	}

	virtual void ProcessPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const;
	virtual void ProcessUnPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const;

	/// initialisation
    virtual void Init(const ARM_ZeroCurveFunctor& discFunctor,const string& modelName,const ARM_TimeInfoPtrVector& timeInfos, const ARM_VectorPtr& numeraireTimes = ARM_VectorPtr(NULL));
    virtual NumeraireType GetType() const { return ARM_Numeraire::PowerTerminalZc; }

	/// Standard ARM Support
	virtual ARM_Object* Clone() const { return new ARM_NumeraireTerminalZc(*this); }
	virtual string NumeraireTypeName() const { return "PowerTerminalZc"; }

	/// Accessors

	inline double GetAlpha() { return itsAlpha; }
	inline void SetAlpha( double Alpha ) { itsAlpha = Alpha; }

private: 
	double itsAlpha;

};

class ARM_NumeraireTerminalEventZc : public ARM_NumeraireFwd
{
public:
	/// constructor, copy constructor, operator=, destructor
	ARM_NumeraireTerminalEventZc() : ARM_NumeraireFwd() {}
	ARM_NumeraireTerminalEventZc( const ARM_NumeraireFwd& rhs ) : ARM_NumeraireFwd( rhs) {}
	virtual ~ARM_NumeraireTerminalEventZc() {};
    ARM_NumeraireTerminalEventZc& operator=(const ARM_NumeraireTerminalEventZc& rhs)
	{
		if( this!= &rhs )
			ARM_NumeraireFwd::operator=(rhs);
		return *this;
	}

	/// initialisation
    virtual void Init(const ARM_ZeroCurveFunctor& discFunctor,const string& modelName,const ARM_TimeInfoPtrVector& timeInfos, const ARM_VectorPtr& numeraireTimes = ARM_VectorPtr(NULL));
    virtual NumeraireType GetType() const { return ARM_Numeraire::TerminalEventZc; }

	/// Standard ARM Support
	virtual ARM_Object* Clone() const { return new ARM_NumeraireTerminalEventZc(*this); }
	virtual string NumeraireTypeName() const { return "TerminalEventZc"; }
};


class ARM_NumeraireCash : public ARM_Numeraire
{
public:
	/// constructor, copy constructor, operator=, destructor
	ARM_NumeraireCash( bool renormalisationIsOn = false ) : ARM_Numeraire(), itsDiscountMap(), itsRenormalisationIsOn( renormalisationIsOn ) {}
	ARM_NumeraireCash( const ARM_NumeraireCash& rhs ) : ARM_Numeraire( rhs), itsDiscountMap( rhs.itsDiscountMap), itsRenormalisationIsOn( rhs.itsRenormalisationIsOn) {}
	virtual ~ARM_NumeraireCash() {};
    ASSIGN_OPERATOR(ARM_NumeraireCash)

	bool GetRenormalisation() const { return itsRenormalisationIsOn;}
	void SetRenormalisation(bool renormalisation) { itsRenormalisationIsOn=renormalisation;}

	/// initialisation
    virtual void Init(const ARM_ZeroCurveFunctor& discFunctor,const string& modelName,const ARM_TimeInfoPtrVector& timeInfos, const ARM_VectorPtr& numeraireTimes = ARM_VectorPtr(NULL));
    virtual void Reset( int pricingDir );
	virtual void FinalizeInduction(	double evalTime, const ARM_PricingStatesPtr& states, const ARM_PricingModel& model );

    virtual NumeraireType GetType() const { return ARM_Numeraire::Cash; }
	virtual void Update(const ARM_PricingModel& model, const ARM_PricingStatesPtr& states, size_t timeIdx);

	/// process payoff and unprocess payoff
	virtual void ProcessPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const;
	virtual void ProcessUnPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const;
	virtual bool NeedLocalDiscount() const { return true; }

	/// accessors
	ARM_VectorPtr GetDiscount(double time) const;
	bool IsDiscountEmpty() const { return itsDiscountMap.size() == 0; }

	/// Standard ARM Support
	virtual ARM_Object* Clone() const { return new ARM_NumeraireCash(*this); }
	virtual string NumeraireTypeName() const { return "Cash"; }

protected:
	typedef map<double,ARM_GP_VectorPtr> DiscountMap;
	DiscountMap itsDiscountMap;
private:
	bool itsRenormalisationIsOn;
};


// Rolling cash numeraire
class ARM_NumeraireRollingCash : public ARM_NumeraireCash
{
public:
	/// constructor, copy constructor, operator=, destructor
	ARM_NumeraireRollingCash(bool renormalisationIsOn = false) : ARM_NumeraireCash(renormalisationIsOn) {}
	ARM_NumeraireRollingCash( const ARM_NumeraireRollingCash& rhs ) : ARM_NumeraireCash( rhs) {}
	virtual ~ARM_NumeraireRollingCash() {};
    ARM_NumeraireRollingCash& operator=(const ARM_NumeraireRollingCash& rhs)
	{
		if( this!= &rhs )
			ARM_NumeraireCash::operator=(rhs);
		return *this;
	}

	void Update(const ARM_PricingModel& model,const ARM_PricingStatesPtr& states,size_t timeIdx );

	/// initialisation
    virtual NumeraireType GetType() const { return ARM_Numeraire::RollingCash; }

	virtual bool NeedLocalDiscount() const { return false; }

	/// Standard ARM Support
	virtual ARM_Object* Clone() const { return new ARM_NumeraireRollingCash(*this); }
	virtual string NumeraireTypeName() const { return "RollingCash"; }
};

/// Rolling numeraire
class ARM_NumeraireRolling : public ARM_Numeraire
{
public:
	/// constructor, copy constructor, operator=, destructor
	ARM_NumeraireRolling() : ARM_Numeraire() {}
	ARM_NumeraireRolling( const ARM_NumeraireFwd& rhs ) : ARM_Numeraire( rhs) {}
	virtual ~ARM_NumeraireRolling() {};
    ARM_NumeraireRolling& operator=(const ARM_NumeraireRolling& rhs)
	{
		if( this!= &rhs )
			ARM_Numeraire::operator=(rhs);
		return *this;
	}

	virtual void MoveNumeraireFwd();
	virtual void MoveNumeraireBckwd();

	virtual void Update(const ARM_PricingModel& model, const ARM_PricingStatesPtr& states, size_t timeIdx);
    virtual void Reset( int pricingDir );
	virtual void ResetLoop();

	/// process payoff and unprocess payoff
	virtual void ProcessPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const;

	virtual void ProcessUnPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const;

	virtual bool NeedLocalDiscount() const { return true; }

protected:
	// Rolling discount factor
	ARM_VectorPtrVector itsRollingDiscount;
};

// Rolling payment numeraire
class ARM_NumeraireRollingPayment : public ARM_NumeraireRolling
{
public:
	/// constructor, copy constructor, operator=, destructor
	ARM_NumeraireRollingPayment() : ARM_NumeraireRolling() {}
	ARM_NumeraireRollingPayment( const ARM_NumeraireRollingPayment& rhs ) : ARM_NumeraireRolling( rhs) {}
	virtual ~ARM_NumeraireRollingPayment() {};
    ARM_NumeraireRollingPayment& operator=(const ARM_NumeraireRollingPayment& rhs)
	{
		if( this!= &rhs )
			ARM_NumeraireRolling::operator=(rhs);
		return *this;
	}

	/// initialisation
    virtual void Init(const ARM_ZeroCurveFunctor& discFunctor,const string& modelName,const ARM_TimeInfoPtrVector& timeInfos, const ARM_VectorPtr& numeraireTimes = ARM_VectorPtr(NULL));
    virtual NumeraireType GetType() const { return ARM_Numeraire::RollingPayment; }

	/// Standard ARM Support
	virtual ARM_Object* Clone() const { return new ARM_NumeraireRollingPayment(*this); }
	virtual string NumeraireTypeName() const { return "RollingPayment"; }
};

// Rolling event numeraire
class ARM_NumeraireRollingEvent : public ARM_NumeraireRolling
{
public:
	/// constructor, copy constructor, operator=, destructor
	ARM_NumeraireRollingEvent() : ARM_NumeraireRolling() {}
	ARM_NumeraireRollingEvent( const ARM_NumeraireRollingEvent& rhs ) : ARM_NumeraireRolling( rhs) {}
	virtual ~ARM_NumeraireRollingEvent() {};
    ARM_NumeraireRollingEvent& operator=(const ARM_NumeraireRollingEvent& rhs)
	{
		if( this!= &rhs )
			ARM_NumeraireRolling::operator=(rhs);
		return *this;
	}

	/// initialisation
    virtual void Init(const ARM_ZeroCurveFunctor& discFunctor,const string& modelName,const ARM_TimeInfoPtrVector& timeInfos, const ARM_VectorPtr& numeraireTimes = ARM_VectorPtr(NULL));
    virtual NumeraireType GetType() const { return ARM_Numeraire::RollingEvent; }

	/// Standard ARM Support
	virtual ARM_Object* Clone() const { return new ARM_NumeraireRollingEvent(*this); }
	virtual string NumeraireTypeName() const { return "RollingEvent"; }
};


CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
