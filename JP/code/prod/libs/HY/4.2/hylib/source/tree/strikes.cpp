/****************************************************************************/
/*                                 */
/****************************************************************************/



#include "strikes.h"  
#include "ldate.h"           /* GtoDayCountFraction */
#include "date_sup.h"         /* GtoDaySubtract */
#include "datelist.h"         /* GtoNewBusDateList */  
#include "macros.h"           /* MAX, MIN */    

extern std::ofstream out;

StrikesClass::StrikesClass(std::string strikeSpec)
{
	static char routine[] = "StrikesClass";
	char      spaceString[80]=" _,";     // stirng seperates tokens
	char      *token;

	if(strikeSpec.empty())
	{
		goto done;
	}

	{
		std::string strikeSpecCopy = strikeSpec;
		char*  strikeSpecPtr = &strikeSpecCopy[0];
		// get strike start date
		token = strtok(strikeSpecPtr, spaceString);
		if(token == NULL)
		{
			throw KException(" Need to have strike start date.", routine);
		}
		char strikeStartDateC[32];
		(void)strcpy(strikeStartDateC,token);
		
		// get strike end date
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException(" Need to have strike end date.", routine);
		}
		char strikeEndDateC[32];
		(void)strcpy(strikeEndDateC,token);
		
   
		// get initial strike type
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have initial strike.", routine);
		}

		double strikeStart = atof(token);

		// get end strike type
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have initial strike.", routine);
		}

		double strikeEnd = atof(token);
		

		// get strike step
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have strike step.", routine);
		}
		double strikeStep;
		strikeStep = atof(token);

		// get step frequency
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have step frequency.", routine);
		}
		double strikeStepFreq;
		strikeStepFreq = atol(token);

		
		// get interp methond
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have step frequency.", routine);
		}
		char interpMethodC[32];
		(void)strcpy(interpMethodC,token);
		m_interpMethod = (long)interpMethodC[0];
	
		//finished reding string

		TDate strikeStartDate,strikeEndDate;

		if(GtoStringToDate(strikeStartDateC,
			   &strikeStartDate) == FAILURE)
		{
			throw KException("Strike start date format is wrong.", routine);
		}

		if(GtoStringToDate(strikeEndDateC,
			   &strikeEndDate) == FAILURE)
		{
			throw KException("Strike end date format is wrong.", routine);
		}

		TDateInterval stepInterval;
		/* determining the step interval */
		if(strikeStepFreq == 0)
		{
			if(GtoDateSubtract(strikeEndDate,
							   strikeStartDate,
							   &stepInterval) == FAILURE)
			{
				throw KException("Failed to creat stepInterval.", routine);
			}
		}
		else if(strikeStepFreq == 365)
		{
			if(GtoMakeDateInterval(1,
							      'D',
								   &stepInterval) == FAILURE)
			{
				throw KException("Failed to creat stepInterval.", routine);
			}
		}
		else
		{
			if(GtoFreqAndTypeToInterval(strikeStepFreq,
										'N',
										&stepInterval) == FAILURE)
			{
				throw KException("Failed to creat stepInterval.", routine);
			}
		}


		/* get strike date list */
		TDateList   *strikeDateList=NULL;
		strikeDateList = GtoNewBusDateList(	strikeStartDate,
											strikeEndDate,
											&stepInterval,
											TRUE,   //end stub
											GTO_BAD_DAY_NONE,
											(char*) "NONE");
		if(strikeDateList == NULL)
		{
			throw KException("Failed to creat strikeDateList.", routine);
		}
 
		long nStrikes = strikeDateList->fNumItems;  
		
		double tempStrike;
		for(int i=0; i<nStrikes; i++)
		{
			tempStrike = strikeStart+i*strikeStep;
			if(strikeStep < 0)
			{
				tempStrike = MAX(tempStrike,strikeEnd)/100.0;
			}
			else
			{
				tempStrike = MIN(tempStrike,strikeEnd)/100.0;
			}

			insert(value_type(KDate(strikeDateList->fArray[i]),tempStrike));
//			out<<GtoFormatDate(strikeDateList->fArray[i])<<", "<< tempStrike<<std::endl;
		}

		GtoFreeDateList(strikeDateList);
	}

done:
	;
}

StrikesClass::StrikesClass(TDate* dates, double* strikes, int numPts, std::string interpType)
{
	static char routine[] = "StrikesClass";
	int i;

	for(i=0; i<numPts;i++)
	{
		insert(value_type(dates[i],strikes[i]));
	}

	m_interpMethod = (long) *(interpType.begin());
}



double StrikesClass::get_strike(KDate currDate) const
{
	double strike;
	const_iterator currStrikeIter,nextStrikeIter;

	nextStrikeIter = upper_bound(currDate);

//	if(nextStrikeIter == begin())
//	{
//		strike = nextStrikeIter->second;
//	}
//	else
//	{

		currStrikeIter = nextStrikeIter;
		currStrikeIter--;

//		if(m_interpMethod == STAIRCASE_INTERP)
		if(m_interpMethod == STAIRCASE_INTERP || nextStrikeIter == end())		//HY3.3.2v
		{
			strike = currStrikeIter->second;
		}
		else
		{
			strike = currStrikeIter->second +(nextStrikeIter->second-currStrikeIter->second)/
				(nextStrikeIter->first - currStrikeIter->first)*(currDate-currStrikeIter->first);
		}
//	}

	return strike;
}

KVector(KDate) StrikesClass::get_datelist() const 
{
	KVector(KDate) dates;
	const_iterator iter;

	for (iter = begin(); iter != end(); iter++)
	{
		dates.push_back( iter->first);
	}
	return dates;
}

KDate  StrikesClass::iDate(int i) const
{
//	const_iterator iter = begin();
	const_iterator iter = begin();
	for(int j=0; j<i; j++)
	{
		iter++;
	}
	return iter->first;
}

double  StrikesClass::iRate(int i) const
{
	const_iterator iter = begin();
	for(int j=0; j<i; j++)
	{
		iter++;
	}
	return iter->second;
}

KVector(KDate) OptionContext::get_datelist() const 
{
	KVector(KDate) dates =m_strikes->get_datelist();
	return dates;
}



OptionContext::OptionContext(std::string optSpec, const StrikesClass* strikes,
							 const StrikesClass* soft_strikes)
{
	static char routine[] = "OptionContext";
	char      spaceString[80]=" _,";     // stirng seperates tokens
	char      *token;

	if(optSpec.empty())
	{
		goto done;
	}

	{
		std::string optSpecCopy = optSpec;
		char*  optSpecPtr = &optSpecCopy[0];
		// get long short flag
		token = strtok(optSpecPtr, spaceString);
		if(token == NULL)
		{
			throw KException(" Need to have long short flag.", routine);
		}
		m_longShort = atol(token);
		
		// get option type
//		token = strtok(NULL, spaceString);
//		if(token == NULL)
//		{
//			throw KException(" Need to have option type.", routine);
//		}
//		char optionTypeC[32];
//		(void)strcpy(optionTypeC,token);
//		m_callput = (long)optionTypeC[0];

		// get american
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException(" Need to have american/Euro.", routine);
		}
		char isAmericanC[32];
		(void)strcpy(isAmericanC,token);

		if( toupper( isAmericanC[0] ) == 'E')
		{
			m_isamerican = false;
		}
		else
		{
			m_isamerican = true;
		}
	

        m_strikes = CM::Raw2ConstSharedPointer(const_cast<StrikesClass*>(strikes));
        m_soft_strikes = CM::Raw2ConstSharedPointer(const_cast<StrikesClass*>(soft_strikes));

	}

done:
		;
		
}


OptionContext::OptionContext(long longShort,bool isAmerican,  const StrikesClass* strikes,
							 const StrikesClass* soft_strikes)
{
	m_longShort = longShort;
	m_isamerican = isAmerican;
	m_strikes = CM::Raw2ConstSharedPointer(const_cast<StrikesClass*>(strikes) );
	m_soft_strikes = CM::Raw2ConstSharedPointer(const_cast<StrikesClass*>(soft_strikes));
}

KDate  OptionContext::iDate(int i) const
{
	return m_strikes->iDate(i);
}


double  OptionContext::iRate(int i) const
{
	return m_strikes->iRate(i);
}

KDate  OptionContext::iSDate(int i) const
{
	return m_soft_strikes->iDate(i);
}


double  OptionContext::iSRate(int i) const
{
	return m_soft_strikes->iRate(i);
}

double OptionContext::payoff(double y, double softy, KDate currDate, double *vega)
{
	double soft_strike;
	double returnValue;
	if(m_soft_strikes == NULL)
	{
		returnValue = payoff(y,currDate, vega);
	}
	else
	{
		soft_strike = 100*m_soft_strikes->get_strike(currDate);
	
		if(m_longShort == -1)
		{
			if(softy > soft_strike)
			{
				returnValue = payoff(y,currDate,vega);
			}
			else
			{
				returnValue = 0.0;
			}
		}
		else
		{
			if(softy < soft_strike)
			{
				returnValue = payoff(y,currDate, vega);
			}
			else
			{
				returnValue = 0.0;
			}
		}
	}


	return returnValue;

}

double OptionContext::payoff(double y, KDate currDate, double *vega)
{
	StrikesClass::const_iterator iter, iter_lower;
	double returnValue;
	double strike;
	if(m_longShort == 0)
	{
		returnValue = 0.0;
	}
	else
	{
		if(m_isamerican == FALSE)
		{
			iter = m_strikes->find(currDate);
			//if it is not a strike date
			if(iter == m_strikes->end())
			{
				returnValue = 0.0;
			}
			else
			{
				strike = m_strikes->get_strike(currDate);
				
				if(m_longShort == -1)
				{					
					returnValue = MAX(y-strike,0.0);
					if( (y-strike) > 0.0)
					{
						*vega = 0.0;
					}
				}
				else
				{
					returnValue = MAX(strike - y,0.0);
						if( (strike - y) > 0.0)
					{
						*vega = 0.0;
					}

				}
	//			out << strike<<std::endl;
			}
		}
		else
		{
			iter_lower = m_strikes->lower_bound(currDate);
			iter = m_strikes->upper_bound(currDate);
	//		out<<GtoFormatDate(currDate)<<std::endl;
			// if there is only one strike, american always - no longer true, HY3.3.2v
			// if there are multiple strikes, the period before first strike date is European -- regardless multiple or not, HY3.3.2v
			
//			if(iter == m_strikes->end())					//strikeEndDate = option maturity 
//			{
//				returnValue = 0.0;
//			}
//			else if(iter == m_strikes->begin() && m_strikes->size() != 1) 
//			{
//				returnValue = 0.0;
//			}

			if(iter_lower == m_strikes->end() || iter == m_strikes->begin())		//HY3.3.2v
			{
				returnValue = 0.0;
			}
			else
			{
				strike = m_strikes->get_strike(currDate);
				
				if(m_longShort == -1)
				{
					returnValue = MAX(y-strike,0.0);
					if( (y-strike) > 0.0)
					{
						*vega = 0.0;
					}
				}
				else
				{
					returnValue = MAX(strike - y,0.0);
					if( (strike - y) > 0.0)
					{
						*vega = 0.0;
					}
				}
			}
		}
	}

	return returnValue * m_longShort;
}
