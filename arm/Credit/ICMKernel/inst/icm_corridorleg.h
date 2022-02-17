
//////////////////////////////////////////////////////////////////////
// ICM_CorridorLeg.h: interface for the ICM_CorridorLeg class.
//
//////////////////////////////////////////////////////////////////////
/*!
	\file ICM_CorridorLeg.h
	\brief Definition of an ICM_CorridorLeg
*/

#if !defined(AFX_ICM_CorridorLeg_H_)
#define AFX_ICM_CorridorLeg_H_

#include "ICMKernel\inst\icm_leg.h"
//#include "ICMKernel\glob\icm_index_fixing_date.h"


class ICM_CorridorLeg : public ICM_Leg  
{
private:
	// INPUT
	ARM_ReferenceValue* itsSpreads; //cloned
	ARM_ReferenceValue* itsRefValue_kup; // cloned
	ARM_ReferenceValue* itsRefValue_kdw; // cloned
	double				itsLeverageFloatIdx;
	int					itsResetfreq; 

	// OUTPUT
	// nbF : Nb Flows in the corridor leg, beginning at the start date of the corridor
	// nbFA : nb Flows in the corridor leg starting at the asOf date of pricing
	ARM_Vector		itsSwapStart;		// nbF size
	ARM_Vector		itsSwapEnd;			// nbF size
	vector<ARM_Vector>		itsSubStart; // nbF  X nbObservPeriod size
	vector<ARM_Vector>		itsSubEnd;	// nbF  X nbObservPeriod size
	vector<ARM_Vector>		itsSubPayDate; // nbF  X nbObservPeriod size
	vector<ARM_Vector>		itsSubReset;	// nbF  X nbObservPeriod size
	vector<ARM_Vector>		itsSubFwdStart; // nbF  X nbObservPeriod size
	vector<ARM_Vector>		itsSubFwdEnd; // nbF  X nbObservPeriod size
	// some vector can be empty.
	vector<ARM_Vector>		itsResetProbas; // nbF  X nbObservPeriod size
	vector<ARM_Vector>		itsVolK1;		// nbF  X nbObservPeriod size
	vector<ARM_Vector>		itsVolK2;		// nbF  X nbObservPeriod size
	vector<ARM_Vector>		itsFwdSpreadAdj;// nbF  X nbObservPeriod size
	vector<ARM_Vector>		itsResetProbas_Float;	// nbF  X nbObservPeriod size
	vector<ARM_Vector>		itsFwdSpreadAdj_Float;	// nbF  X nbObservPeriod size

public:

	ICM_CorridorLeg():ICM_Leg()
	{Init();}		

	ICM_CorridorLeg(const ARM_Date& startDate, 
            const ARM_Date& endDate, 
			const ARM_Date* refDate,
			const ARM_Date* FstCpnEffDate,
			ARM_ReferenceValue* Spreads,
			ARM_IRIndex* floatingIdx,
			const double LeverageFloatIdx,
			ICM_Credit_Index* creditIdx,
			ARM_ReferenceValue* RefValue_kup,
			ARM_ReferenceValue* RefValue_kdw,
			const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /*= qACCRUED_SETTLED */,
			const int& AccruedDayCount /*= KACTUAL_365*/,
			const int& freq /*= K_ANNUAL*/, 
			const int& Resetfreq /*= K_ANNUAL*/, 
            const int& dayCount/*= K30_360*/, 
            const int& payTiming /*= K_ARREARS*/,
            const int& intRule /*= K_ADJUSTED*/,
            const int& stubRule /*= K_SHORTSTART*/,
            const std::string& discountCcy /* ARM_Currency* discountCcy = ARM_DEFAULT_CURRENCY*/ ,
            const std::string& payCalName /*= NULL*/);

	ICM_CorridorLeg(const ARM_Date& startDate, 
            const ARM_Date& endDate, 
			const ARM_Date* refDate,
			const ARM_Date* FstCpnEffDate,
			ARM_ReferenceValue* Spreads,
			const double Notional,
			const long recieveOrPay,
			ARM_IRIndex* floatingIdx,
			const double LeverageFloatIdx,
			ICM_Credit_Index* creditIdx,
			ARM_ReferenceValue* RefValue_kup,
			ARM_ReferenceValue* RefValue_kdw,
			const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /*= qACCRUED_SETTLED */,
			const int& AccruedDayCount /*= KACTUAL_365*/,
			const int& freq /*= K_ANNUAL*/, 
			const int& Resetfreq /*= K_ANNUAL*/, 
            const int& dayCount/*= K30_360*/, 
            const int& payTiming /*= K_ARREARS*/,
            const int& intRule /*= K_ADJUSTED*/,
            const int& stubRule /*= K_SHORTSTART*/,
            const std::string& discountCcy /* ARM_Currency* discountCcy = ARM_DEFAULT_CURRENCY*/ ,
            const std::string& payCalName /*= NULL*/);

	~ICM_CorridorLeg() 
	{
		if (itsSpreads) delete itsSpreads;
		if (itsRefValue_kup) delete itsRefValue_kup;
		if (itsRefValue_kdw) delete itsRefValue_kdw;

		itsSpreads=NULL;
		itsRefValue_kup=NULL;
		itsRefValue_kdw=NULL;
	}

	void Init(void)
	{
		SetName(ICM_LEG);
		itsSpreads=NULL;
		itsRefValue_kup=NULL;
		itsRefValue_kdw=NULL;
		itsLeverageFloatIdx=1.;
		itsResetfreq=K_QUARTERLY; 
	}

	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg);

	ARM_Object* Clone(void);

	void View(char* id = NULL, FILE* fOut = NULL);

	ARM_ReferenceValue* GetSpreads() {return itsSpreads;}
	ARM_ReferenceValue* GetRefValue_kup() {return itsRefValue_kup;}
	ARM_ReferenceValue* GetRefValue_kdw() {return itsRefValue_kdw;}
	double				GetLeverageFloatIdx() {return itsLeverageFloatIdx;}
	int					GetResetfreq() {return itsResetfreq; }
	
	void SetRefValue_kup(ARM_ReferenceValue* value) 
	{
		if (itsRefValue_kup) delete itsRefValue_kup;
		itsRefValue_kup=(ARM_ReferenceValue*)value->Clone();
	}

	void SetRefValue_kdw(ARM_ReferenceValue* value) 
	{
		if (itsRefValue_kdw) delete itsRefValue_kdw;
		itsRefValue_kdw=(ARM_ReferenceValue*)value->Clone();
	}

	void SetLeverageFloatIdx(const double value) {itsLeverageFloatIdx=value;}
	void SetResetfreq(const int value) {itsResetfreq=value; }
	
	//value must be un %
	void SetSpreads(ARM_ReferenceValue* value) 
	{
		if (itsSpreads) delete itsSpreads;
		itsSpreads=(ARM_ReferenceValue*)value->Clone();
	}

	ARM_Vector* GenResetSched() 
	{
		ARM_Vector* sched_ = NULL;
		ARM_Vector* sched = new ARM_Vector();
		ICM_Security* security = GetCreditInfos();
		int	NbFlows = security->GetAccStartDates().GetSize();

		for (int i=0;i<NbFlows;i++)
		{
		ARM_Date startdate = security->GetAccStartDates().Elt(i);
		ARM_Date enddate = security->GetAccEndDates().Elt(i);

		sched_= CptStartDates(startdate,enddate,itsResetfreq,K_MOD_FOLLOWING,
                            K_LONGSTART,K_ADJUSTED,GetCurrencyUnit()->GetCcyName(),/*adjFirstdate*/1);

		for (int j=0;j<sched_->GetSize();j++)
			{sched->push_back(sched_->Elt(j));}

		if (sched_) delete sched_;
		}

		sched->push_back(GetEndDateNA().GetJulian());

		return sched;
	}

	void GenCorridorSchedule();
	void ExtractBounds(double startdate,double enddate,int& begin,int& end);

	vector<ARM_Vector>& GetResetProbas() {return itsResetProbas;}
	void SetResetProbas(vector<ARM_Vector>& value) {itsResetProbas=value;}

	ARM_Vector&					GetSwapStart()  {return	itsSwapStart;}
	ARM_Vector&					GetSwapEnd() {return itsSwapEnd;}
	const vector<ARM_Vector>&	GetSubStart()  {return	itsSubStart;}
	const vector<ARM_Vector>&	GetSubEnd()  {return itsSubEnd;}
	const vector<ARM_Vector>&	GetSubResetSchedule() {return itsSubReset;}
	const vector<ARM_Vector>&	GetSubFwdStart()  {return	itsSubFwdStart;}
	const vector<ARM_Vector>&	GetSubFwdEnd() {return itsSubFwdEnd;}
	const vector<ARM_Vector>&	GetSubPayDate()  {return itsSubPayDate;}

	vector<ARM_Vector>&	GetVolK1() {return	itsVolK1;}
	vector<ARM_Vector>&	GetVolK2() {return	itsVolK2;}
	vector<ARM_Vector>&	GetFwdSpreadAdj() {return itsFwdSpreadAdj;}
	vector<ARM_Vector>&	GetFwdSpreadAdj_Float() {return itsFwdSpreadAdj_Float;}
	vector<ARM_Vector>&	GetResetProbas_Float() {return itsResetProbas_Float;}

	void SetSwapStart(const ARM_Vector& value) {itsSwapStart=value;}
	void SetSwapEnd(const ARM_Vector& value) {itsSwapEnd=value;}
	void SetSubStart(const vector<ARM_Vector>& value) {itsSubFwdStart=value;}
	void SetSubEnd(const vector<ARM_Vector>& value) {itsSubFwdEnd=value;}
	void SetSubResetSchedule(const vector<ARM_Vector>& value) {itsSubReset=value;}
	void SetSubFwdStart(const vector<ARM_Vector>& value) {itsSubFwdStart=value;}
	void SetSubFwdEnd(const vector<ARM_Vector>& value) {itsSubFwdEnd=value;}
	void SetSubPayDate(const vector<ARM_Vector>& value) {itsSubPayDate=value;}

	void SetVolK1(const vector<ARM_Vector>& value) {itsVolK1=value;}
	void SetVolK2(const vector<ARM_Vector>& value) {itsVolK2=value;}
	void SetFwdSpreadAdj(const vector<ARM_Vector>& value) {itsFwdSpreadAdj=value;}
	void SetFwdSpreadAdj_Float(const vector<ARM_Vector>& value) {itsFwdSpreadAdj_Float=value;}
	void SetResetProbas_Float(const vector<ARM_Vector>& value) {itsResetProbas_Float=value;}

	virtual void CptCoupons(ARM_Object* object,const ARM_Date& asof);

}; 

#endif // !defined(AFX_ICM_CorridorLeg_H_)
