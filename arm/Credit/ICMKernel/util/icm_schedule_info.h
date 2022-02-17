#ifndef __ICM_SCHEDULE_INFO_H__
#define __ICM_SCHEDULE_INFO_H__


#pragma warning (disable : 4786 )

#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\glob\icm_enums.h"
#include "ARMKernel\glob\armglob.h"
#include "ARMKernel\glob\dates.h"

class ICM_Schedule_Info : public ARM_Object 
{
public :
	ARM_Date	itsEffectiveDate;
	ARM_Date	itsMaturityDate;			
	int	itsPayFrequency;
	int itsResetFreq ;
	int	itsDayCount;
	int	itsStubrule;
	int	itsIntRule;
	std::string itsPayCalName;
	int itsPayTiming;
	int itsResetTiming;
	int itsFwdRule;
	bool	itsIncludeMaturity;
	int itsAdj;
	int	itsStartAdj;
	int itsAccDayCount;
	qCDS_ADJ itsCDSAdj;
private :
	ARM_Date*	itspReferenceDate;
	ARM_Date*	itspFirstCpnEffDate;

public :


	ICM_Schedule_Info(const ARM_Date&	EffectiveDate,
							const ARM_Date&	MaturityDate,			
							int	PayFrequency,
							int ResetFreq ,
							int	DayCount,
							int	Stubrule,
							int	IntRule,
							const std::string& PayCalName,
							int PayTiming,
							int ResetTiming,
							int FwdRule,
							bool	IncludeMaturity,
							int Adj,
							int	StartAdj,
							int AccDayCount,
							const ARM_Date*	pReferenceDate, // always copied
							const ARM_Date*	pFirstCpnEffDate, // always copied
							qCDS_ADJ CDSAdj);

// FIXMEFRED: mig.vc8 (28/05/2007 13:23:45): missing return type
	void Set(const ARM_Date&	EffectiveDate,
		const ARM_Date&	MaturityDate,			
		int	PayFrequency,
		int ResetFreq ,
		 int	DayCount,
		 int	Stubrule,
		 int	IntRule,
		const std::string& PayCalName,
		 int PayTiming,
		 int ResetTiming,
		 int FwdRule,
		 bool	IncludeMaturity,
		 int Adj,
		 int	StartAdj,
		 int AccDayCount,
		const ARM_Date*	pReferenceDate,
		const ARM_Date*	pFirstCpnEffDate,
		 qCDS_ADJ CDSAdj);

	void Init(void);
	void Copy(const ARM_Object* srcInfo);
	void BitwiseCopy(const ARM_Object* srcInfo);
	virtual ARM_Object* Clone(void);
	void View(char* id = NULL, FILE* ficOut = NULL);
	~ICM_Schedule_Info();

	void SetpReferenceDate(const ARM_Date* pReferenceDate) {
		if (itspReferenceDate) delete itspReferenceDate;
		itspReferenceDate = (ARM_Date*)(pReferenceDate)->Clone();
	}
	void SetpFirstCpnEffDate(const ARM_Date* pFirstCpnEffDate){
		if (itspFirstCpnEffDate) delete itspReferenceDate;
		itspFirstCpnEffDate = (ARM_Date*)(pFirstCpnEffDate)->Clone();
	}

	const ARM_Date* GetpReferenceDate() const {
		return itspReferenceDate;
	}

	const ARM_Date* GetpFirstCpnEffDate() const {
		return itspFirstCpnEffDate;
	}


private :
	ICM_Schedule_Info() { Init() ;};
	ICM_Schedule_Info(const ICM_Schedule_Info& aICM_Schedule_Info)  { Copy((const ARM_Object*)&aICM_Schedule_Info) ;};
};



#endif
