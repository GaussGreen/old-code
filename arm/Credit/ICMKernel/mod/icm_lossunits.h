
#ifndef _ICM_LOSS_UNITS
#define _ICM_LOSS_UNITS

#include <map>
#include "ARMKernel/glob/dates.h"
#include "ICMKernel/util/icm_vector.h"

class ICM_LossUnits : public ARM_Object 
{
public:
	typedef struct lossunit_data_t
	{
		bool   IsNull;
		double LU; 
		ARM_Vector LossRates; 
		ARM_Vector LSLossRates; 
		std::vector<int> CollatRank;  
	} ;
	typedef std::map<ARM_Date,lossunit_data_t> cont_t  ; 
private:
	cont_t itsLossUnits; 
public:
	ICM_LossUnits() ; 
	ICM_LossUnits(const ICM_LossUnits&ref); 
	ICM_LossUnits& operator=(const ICM_LossUnits&ref); 
	virtual ~ICM_LossUnits(); 
public:
	unsigned int size() const { return itsLossUnits.size(); }
	void insert(const ARM_Date&,const lossunit_data_t&item) ;
	double getLossUnit(const ARM_Date&date) const; 
	bool IsLossUnitNull(const ARM_Date&date) const; 
	const ARM_Vector& getLossRates(const ARM_Date&date) const; 
	const ARM_Vector& getLSLossRates(const ARM_Date&date) const; 
	const std::vector<int>& getCollatRank(const ARM_Date&date) const; 
	const ARM_Date& getDate(unsigned int i) const; 
	void getDates(ICM_Vector<ARM_Date>&ret) const; 

public:
	virtual void View(char* id = NULL, FILE* ficOut = NULL) ;
} ;


#endif // _ICM_LOSS_UNITS