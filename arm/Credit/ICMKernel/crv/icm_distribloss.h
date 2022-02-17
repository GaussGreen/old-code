
#ifndef _ICM_DISTRIBLOSS_H_ 
#define _ICM_DISTRIBLOSS_H_ 

#include "ARMKernel\glob\armglob.h"
// #include "ICMKernel\util\icm_macro.h"
#include "ICMKernel/util/icm_utils.h"

#include <map>

class	ICM_DistribLoss	: public ARM_Object 
{
typedef std::map<ORDER_DOUBLE_WITH_ORDER,double> map_t; 
static const double prec; 
bool itsUseStdEL;

public:

	ICM_DistribLoss() ; 
	ICM_DistribLoss(const ICM_DistribLoss&ref); 
	ICM_DistribLoss& operator=(const ICM_DistribLoss&ref); 
	virtual ~ICM_DistribLoss(); 

	void Init();
	//	
	virtual ARM_Object* Clone() ; 
	ARM_Object* Clone() const; 
	virtual void View(char* id = NULL, FILE* ficOut = NULL); 
	// not redefined virtual ARM_CLASS_NAME GetRootName(void) ;
	// not redefined virtual void Print() ; 
	// 
	size_t	size()	 const			{	return itsExpectedLosses.size() ; } 
	bool 	empty()	 const			{	return itsExpectedLosses.empty() ; } 
	void	clear()					{	itsExpectedLosses.clear();itsExpectedLosses_disc.clear(); }
	double& operator[](double yf)	{	return (itsUseStdEL?
											itsExpectedLosses[ORDER_DOUBLE_WITH_ORDER(yf,prec)]:
											itsExpectedLosses_disc[ORDER_DOUBLE_WITH_ORDER(yf,prec)]); }
	// 
	double InterpolEL(double yearterm) const  ;
	double InterpolEL_TS(double yearterm) const  ;
	//
	friend std::ostream& operator<<(std::ostream&,const ICM_DistribLoss&); 

	map_t itsExpectedLosses; 
	map_t itsExpectedLosses_disc; //El computes for jumps times
	//vector<double> itsJumpsPoints;

	void SetUseStdEL(const bool& status) {itsUseStdEL = status;}
	void UnifyELs(map_t& ExpectedLosses_disc_unify) const; //copy itsExpectedLosses to  itsExpectedLosses_disc avoiding jump points

	//void SetJumpsPoints(vector<double>& JumpsPoints) {itsJumpsPoints = JumpsPoints;}
} ;

inline void ICM_DistribLoss::Init()
{itsUseStdEL = true;
 //itsJumpsPoints.clear();	
}

//	-------------------------------------------------------------------------
inline 
ICM_DistribLoss::ICM_DistribLoss()
{Init();}
//	--------------------------------------------------------------------------
inline 
ICM_DistribLoss::ICM_DistribLoss(const ICM_DistribLoss&ref) : ARM_Object(ref)
{itsUseStdEL = ref.itsUseStdEL;
itsExpectedLosses=ref.itsExpectedLosses;
itsExpectedLosses_disc=ref.itsExpectedLosses_disc;	
//itsJumpsPoints=ref.itsJumpsPoints;
}
//	--------------------------------------------------------------------------
inline ICM_DistribLoss& 
ICM_DistribLoss::operator=(const ICM_DistribLoss&ref)
{
	if (this!=&ref)
	{
		this->~ICM_DistribLoss(); 
		new(this)ICM_DistribLoss(ref); 
	}
	return *this; 
}
//	--------------------------------------------------------------------------
inline 
ICM_DistribLoss::~ICM_DistribLoss()
{}
//	--------------------------------------------------------------------------
//	virtual 
inline ARM_Object* 
ICM_DistribLoss::Clone()
{
	return new ICM_DistribLoss(*this); 
}
//	--------------------------------------------------------------------------
inline ARM_Object* 
ICM_DistribLoss::Clone() const
{
	return new ICM_DistribLoss(*this);
	//return unconst(this)->Clone(); 
}
//	--------------------------------------------------------------------------
inline 	double 
ICM_DistribLoss::InterpolEL(double yearterm) const 
{
	if (itsExpectedLosses_disc.empty()==false)
		return InterpolEL_TS(yearterm);

	// copied from ICM_Pricer_Security. 
	map_t::const_iterator it = (itsUseStdEL ? itsExpectedLosses.find(yearterm) : itsExpectedLosses_disc.find(yearterm)); 

	if (it!=((itsUseStdEL)?itsExpectedLosses.end():itsExpectedLosses_disc.end())) return it->second; 
	
	int size = (itsUseStdEL?itsExpectedLosses.size():itsExpectedLosses_disc.size());
	//we do linear interpolation if map is not containing exact value
	std::vector<double> values(size);
	std::vector<double> yearterms(size);
	int i=0;
	for (it=(itsUseStdEL?itsExpectedLosses.begin():itsExpectedLosses_disc.begin());
	it!=(itsUseStdEL?itsExpectedLosses.end():itsExpectedLosses_disc.end());++it)
	{
		yearterms[i]=it->first.m_value;
		values[i]=it->second;i++;
	}

	return LinearVectorInterpol(yearterms,values,yearterm);
} 

inline void ICM_DistribLoss::UnifyELs(map_t& ExpectedLosses_disc_unify) const
{
	ExpectedLosses_disc_unify.clear();
	map_t EL = itsExpectedLosses; 

	if (itsExpectedLosses_disc.empty()) {ExpectedLosses_disc_unify = itsExpectedLosses;return;}

	map_t::const_iterator it; 
	for (it=itsExpectedLosses_disc.begin();it!=itsExpectedLosses_disc.end();++it)
	{EL[it->first.m_value]=it->second;}

	ExpectedLosses_disc_unify = EL;
}

#endif // 