#include "ARMKernel/glob/firsttoinc.h" 
#include "ICMKernel/util/icm_utils.h"
#include "icm_spreadoption.h"


//	-------------------------------------------------------------------------------------------
ICM_SpreadOption::ICM_SpreadOption()
{
	SetName(ICM_SPREADOPTION);
}
//	-------------------------------------------------------------------------------------------
ICM_SpreadOption& 
ICM_SpreadOption::operator=(const ICM_SpreadOption&ref)
{
	if (this!=&ref)
	{
		this->~ICM_SpreadOption(); 
		new(this)ICM_SpreadOption(ref); 
	}
	return *this; 
}
//	-------------------------------------------------------------------------------------------
//	virtual 
ARM_Object* 
ICM_SpreadOption::Clone()const
{
	return new ICM_SpreadOption(*this); 
}
//	-------------------------------------------------------------------------------------------
//	virtual 
ARM_Object* 
ICM_SpreadOption::Clone()
{
	return unconst(*this).Clone(); 
}
//	-------------------------------------------------------------------------------------------
//	virtual 
ARM_CLASS_NAME 
ICM_SpreadOption::GetRootName()
{
	return ARM_SECURITY ;
}
//	-------------------------------------------------------------------------------------------
//	virtual
void 
ICM_SpreadOption::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char  fOutName[200];

    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

       (void) unlink(fOutName);

       fOut = fopen(fOutName, "w"); 
    }
    else
    {
       fOut = ficOut;
    } 

	std::stringstream sstr; 
	View(sstr); 
	// fprintf(fOut, "\t\t\t ----------------- ICM SpreadOption Informations ----------------- \n");
	fprintf(fOut,"%s",sstr.str().c_str()); 
	
    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}
//	-------------------------------------------------------------------------------------------
//	
void 
ICM_SpreadOption::View(std::ostream& o) const 
{
	o<<"----------------- BEGIN ICM SpreadOption Informations -----------------"<<std::endl; 
	o<<"Strike: "<< itsStrike <<std::endl; 
	o<<"IsCall: "<< itsIsCall<<std::endl; 
	o<<"IsKO: "<< itsIsKO<<std::endl; 
	o<<"IsAccelerated: "<< itsIsAccelerated<<std::endl; 
	o<<"ExerciseDates: "<<itsExerciseDates<<std::endl; 
	o<<"ExerciseStyle: "<<itsExerciseStyle<<std::endl; 
	o<<"ExerciseFrequency: "<<itsExerciseFrequency<<std::endl; 
	o<<"UnderlyingMatuStyle: "<<itsUnderlyingMaturityStyle<<std::endl; 
	o<<"----------------- END ICM SpreadOption Informations -----------------"<<std::endl; 
}

//	-------------------------------------------------------------------------------------------
// static 
ICM_SpreadOption* 
ICM_SpreadOption::build(
		double itsStrike,
		bool isCall, 
		qKoStyle koStyle, 
		qAccelerationStyle accStyle,
		const ARM_Vector& exerciseDates,
		int exerciseStyle,
		int exerciseFrequency,
		qUnderlying_Maturity_Style matuStyle,
		const ICM_Cds& ref
		)
{
	std::auto_ptr<ICM_SpreadOption> item (new ICM_SpreadOption) ; 
	if (!item.get()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Memory Issue allocating ICM_SpreadOption"); 
	item->setStrike(itsStrike); 
	item->setIsCall(isCall); 
	item->setIsKO(koStyle==qKO); 
	item->setIsAccelerated(accStyle==qACC); 
	item->setExerciseDates(exerciseDates); 
	item->setExerciseStyle(exerciseStyle); 
	item->setExerciseFrequency(exerciseFrequency); 
	item->setSecurity(ref); 
	item->setUnderlyingMaturityStyle(matuStyle);
	return item.release(); 
}
//	------------------------------------------------------------------------------------------------
std::string			
ICM_SpreadOption::getProductTypeKey() const
{
	std::stringstream sstr ;
	sstr<<"["<<TYPE()<<":ICM_SpreadOption]"; 
	sstr<<"["<<ACCELERATED()<<":"<<IsAccelerated()<<"]" ;
	sstr<<"["<<EXSTYLE()<<":"<<getExerciseStyle()<<"]" ;
	sstr<<"["<<UNDMATUSTYLE()<<":"<<ICM_EnumsCnv::toString(underlyingMaturityStyle())<<"]" ;
	sstr<<"["<<KO()<<":"<<IsKO()<<"]" ;
	return sstr.str(); 
}