
//////////////////////////////////////////////////////////////////////
// ICM_GEN.h: interface for the ICM_GenCF class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_ICM_GEN_H_)
#define _ICM_GEN_H_


#include "ARMKernel\inst\security.h"
#include "ICMKernel\util\icm_matrix.h"

/*********************************************************************************/
/*! \class  ICM_GenCF icm_gen.h "icm_gen.h"
 *  \author 
 *	\version 1.0
 *	\date   June 2003
 *	\brief  Pricing <B> Generic Cash Flow </B> */
/***********************************************************************************/


class ICM_GenCF : public ARM_Security  
{
private: 

	ICM_Matrix<ARM_Vector>* itsMatrix;

public:

	void Init();


	ICM_GenCF()
	{
		Init();
	}	

/* ***************************************************************************************************************** */
/*!	\fn ICM_Cds(ARM_Date EffectiveDate,
				ARM_Date ScheduleTerminateDate,
				double FixedRate,
				ARM_Date FirstPeriodReferenceDate ,
				int Frequency = K_QUARTERLY,
				int DayCountBasis = KACTUAL_360,
				double FixedPayerAmount = 100., 
				qPAYMENT_PREMIUM_LEG AccruedOnDefault = qACCRUED_SETTLED,
				qPAYS_ON_DEFAULT AmortizationOnDefault = qPays_Recovery,
				qPAYS_ON_DEFAULT InterestOnDefault = qPays_Recovery,
				ARM_Currency *Ccy = ARM_DEFAULT_CURRENCY, 
				double FloatingPayerAmount = 0.0,
				int stubrule = 3)

	\brief Constructor of an <B> ISDA CDS </B>
	\param EffectiveDate Effective Date as in ISDA definition
	\param ScheduleTerminationDate Maturity Date Unadjusted
	\param FixedRate Premiun in unit of the leg
	\param FirstPeriodReferenceDate of the fixed leg if <B> 1 </B> computed from TerminationDate. 
	\param Ccy Currency of the leg: The default currency is a function of the site 
	\param Frequency Fixed Leg Frequency default value is K_QUATERLY
	\param DayCountBasis Fixed Leg Day Count Basis default value is KACTUAL_360
	\param FixedPayerAmount Fixed payer amount default value is 100
	\param AccruedOnDefault Fixed Payer accrued in case of default
	\param Ccy Currency of the Fixed and Floating leg
	\param FloatingPayerAmount: Amount of the floating leg: If <B> 0 </B> it should be equal to the fixed leg amount
	\note This constructor will make a FixedLeg and Floating Leg 
	<ul>

	<B>itsFixedLeg uses:</B> 
			<li> EffectiveDate </li>
			<li> ScheduleTerminationDate </li>
			<li> FirstPeriodReferenceDate </li>
			<li> Fixed Rate </li>
			<li> Frequency  </li> 
			<li> DayCountBasis </li>
			<li> FixedPayer Amount </li>
			<li> Accrued on Default </li>
			<li> Ccy </li>
	</ul>
	<ul>
	<B>itsFloatingLeg uses:</B> 
			<li> EffectiveDate </li>
 			<li> ScheduleTerminationDate </li>
			<li> Ccy </li>
			<li> FloatingPayerAmount</li>

	</ul>
*/	
/* ***************************************************************************************************************** */

	ICM_GenCF (ICM_Matrix<ARM_Vector>* matrice);

	void Set(ICM_Matrix<ARM_Vector>* matrice);

	virtual ~ICM_GenCF()
	{
	if (itsMatrix)
		delete itsMatrix;
	itsMatrix = NULL;
	}	

	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg);

	ARM_Object* Clone(void);

	void View(char* id, FILE* ficOut);

	// Methodes Set and Get ----------------------------

	void SetMatrix(ICM_Matrix<ARM_Vector>* matrix)
	{
		if (itsMatrix)
			delete itsMatrix;
		itsMatrix = matrix;
	}

	ICM_Matrix<ARM_Vector>* GetMatrix(void)	{ return itsMatrix;	}

};

#endif // !defined(AFX_ICM_CDS_H__0FAFBE8F_B970_476B_AE0F_A498B958D585__INCLUDED_)
