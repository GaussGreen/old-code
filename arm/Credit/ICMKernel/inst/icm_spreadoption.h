#ifndef _ICM_SPREADOPTION_H
#define _ICM_SPREADOPTION_H

#include "ARMKernel\inst\security.h"
#include "ICMKernel/util/icm_ptr.h"
#include "icm_cds.h"
#include <memory>

//	
//		Inheritance from ARM_Security is not really powerful
//		but required .. 
//
class ICM_SpreadOption : public ARM_Security
{
public:
	static const std::string TYPE()  ; 
	static const std::string ACCELERATED(); 
	static const std::string EXSTYLE(); 
	static const std::string KO(); 
	static const std::string UNDMATUSTYLE(); 
private:
	typedef aggreg_ptr<ICM_Cds>	ICM_SecurityPtr ;    
private:
	double					itsStrike;			// option strike 
	bool					itsIsCall;			// K_CALL, K_PUT
	bool					itsIsKO;			// otherwise not KO
	bool					itsIsAccelerated;	// otherwise not accelerated
	ARM_Vector				itsExerciseDates;	
	int						itsExerciseStyle;	// K_EUROPEAN|K_AMERICAN|K_BERMUDAN
	ICM_SecurityPtr			itsSecurity;		// 0.1, Aggregation
	qUnderlying_Maturity_Style	itsUnderlyingMaturityStyle; 
	int						itsExerciseFrequency;	
public:
	//	--	Standard C++ 
	ICM_SpreadOption() ;
	ICM_SpreadOption& operator=(const ICM_SpreadOption&); 
	virtual ~ICM_SpreadOption() ; 
	virtual ARM_Object* Clone()const ; 
	void View(std::ostream& o) const ; 
	//	ARM_Object
	virtual ARM_Object* Clone(); 
	virtual void View(char* id, FILE* ficOut) ;
	virtual void Copy(const ARM_Object* src) { ICMTHROW(ERR_INVALID_OPERATION,"Not Implemented"); }
    virtual ARM_CLASS_NAME GetRootName() ; 

	//	--	Product Description 
	//
	std::string			getProductTypeKey() const; 
	// 
	double				getStrike() const						{	return itsStrike; }
	void				setStrike(double d)						{	itsStrike=d; }
	const ARM_Date 		getExpiryDate() const					{	return itsExerciseDates.Elt(itsExerciseDates.GetNumLines()-1) ; }
	const ARM_Vector&	getExerciseDates() const				{ return itsExerciseDates; }
	void				setExerciseDates(const ARM_Vector&ref)	{ itsExerciseDates=ref ; }
	int					getExerciseStyle() const				{ return itsExerciseStyle; }
	void				setExerciseStyle(int i) 				{ itsExerciseStyle=i; }
	int					getExerciseFrequency() 					{ return itsExerciseFrequency; }
	void				setExerciseFrequency(int i) 			{ itsExerciseFrequency=i; }
	const ICM_Cds&		getSecurity() const						{ return *itsSecurity; }
	void				setSecurity(const ICM_Cds& ref)			{ itsSecurity.adopt(dynamic_cast<ICM_Cds*>(ref.Clone())); }
	// 
	bool				IsCall() const							{	return itsIsCall; }
	void				setIsCall(bool b)						{	itsIsCall=b; } 
	bool				IsKO() const							{	return itsIsKO; }
	void				setIsKO(bool b)							{	itsIsKO=b; } 
	bool				IsAccelerated() const					{	return itsIsAccelerated; }
	void				setIsAccelerated(bool b)				{	itsIsAccelerated=b; } 
	qUnderlying_Maturity_Style	underlyingMaturityStyle()	const	{	return itsUnderlyingMaturityStyle ;} 
	void				setUnderlyingMaturityStyle(qUnderlying_Maturity_Style e)	 {itsUnderlyingMaturityStyle=e ;} 
	//	-- Services
	// 
	static ICM_SpreadOption* build(
		double itsStrike,
		bool isCall, 
		qKoStyle , 
		qAccelerationStyle,
		const ARM_Vector& exerciseDates,
		int exerciseStyle,
		int exerciseFrequency,
		qUnderlying_Maturity_Style	matuStyle,
		const ICM_Cds& ref
		); 
	//	--	Not Implemented
protected:
	void invariant() const ; 
private: 
	virtual ARM_Object& operator = (const ARM_Object& obj) { ICMTHROW(ERR_INVALID_OPERATION,"Not Implemented"); }
};

//	-------------------------------------------------------------------------------------
inline 
ICM_SpreadOption::~ICM_SpreadOption() 
{}
//	-------------------------------------------------------------------------------------

inline const std::string ICM_SpreadOption::TYPE()								{ return "TYPE" ; } 
inline const std::string ICM_SpreadOption::ACCELERATED()						{ return "ACC" ; } 
inline const std::string ICM_SpreadOption::EXSTYLE()							{ return "EXERC" ; } 
inline const std::string ICM_SpreadOption::KO()									{ return "KO" ; } 
inline const std::string ICM_SpreadOption::UNDMATUSTYLE()						{ return "UNDMATU" ; } 

#endif // _ICM_SPREADOPTION_H
