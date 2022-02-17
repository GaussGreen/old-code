/****************************************************************************/
/*      Declaration file for structure of data.                             */
/****************************************************************************/

#ifndef	_STRIKE_H_
#define	_STRIKE_H_


#include "kplatform.h"   
#include "General/General.h"
#include "bastypes.h"        // TCouponDates
#include "cdate.h"           // TDate TDateInterval
#include "kdate.h"   
#include "string"

#define LINEAR_INTERP        (long)'L'
#define STAIRCASE_INTERP     (long)'S'

#define HY_CALL               (long)'C'
#define HY_PUT               (long)'P'




class StrikesClass:public std::map<KDate,double>, public CM::Object
{
	long   m_interpMethod;

public:
	StrikesClass(){}
	StrikesClass(std::string strikeSpec);
	StrikesClass(TDate* dates, double* strikes, int numPts, std::string interpType);

	double get_strike(KDate currDate) const;
	KDate  iDate(int i) const;
	double iRate(int i) const;
	KVector(KDate) get_datelist() const;

};

class OptionContext:public CM::Object
{
	long m_longShort;   //1,0,-1
//	long m_callput;
	bool m_isamerican;

	const StrikesClass*  m_strikes;
	const StrikesClass*  m_soft_strikes;

public:
	OptionContext(){m_strikes=0;m_soft_strikes=0;};
	OptionContext(std::string optSpec, const StrikesClass* strikes,const StrikesClass* soft_strikes = 0);
	OptionContext(long longShort,bool isAmerican, const StrikesClass* strikes,const StrikesClass* soft_strikes = 0);
	~OptionContext(){delete m_strikes;delete m_soft_strikes;};

	KDate  iDate(int i) const;
	double iRate(int i) const;
	KDate  iSDate(int i) const;
	double iSRate(int i) const;
	int    strikesSize(){return m_strikes->size();}
	int    softStrikesSize(){return (m_soft_strikes == NULL)? 0:m_soft_strikes->size();}

	long get_option_direction(){return m_longShort;}
	bool get_option_isamerican(){return m_isamerican;}
	double payoff(double y, KDate currDate, double *veag);
	double payoff(double y, double softy, KDate currDate, double *vega);
	KVector(KDate) get_datelist() const;

};


#endif

