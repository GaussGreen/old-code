/****************************************************************************/
/*                                 */
/****************************************************************************/



#include "ddmap.h"       
#include "fltrate.h"         /* GtoForwardRateSafe */
#include "ldate.h"           /* GtoDayCountFraction */
#include "date_sup.h"         /* GtoDaySubtract */
#include "datelist.h"         /* GtoNewBusDateList */
#include "kratecurve.h"       

DDMap::DDMap(std::string spec, KDate endDate)
{

	static char routine[] = "DDMap::DDMap";
	char      spaceString[80]=" _,";     /* stirng seperates tokens */
	char      *token;

	if(spec.empty())
	{
		goto done;
	}

	{
		std::string specCopy = spec;
		char*  specPtr = &specCopy[0];
		/* get first date */
		token = strtok(specPtr, spaceString);
		if(token == NULL)
		{
			throw KException(" Need to have first date.", routine);
		}
		char firstDate[32];
		(void)strcpy(firstDate,token);

		/* get frequency */
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have frequency.", routine);
  
		}
		long frequency;
		frequency = atol(token);

		/* get rate */
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have rate.\n", routine);
		}

		double numrate = atof(token);

		/* get day count convention */
		long dayCount;
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			dayCount = GTO_ACT_365;
		}
		else
		{
		
			if(GtoStringToDayCountConv(token,&dayCount) == FAILURE)
			{
				throw KException("Wrong Accrual day count convention.", routine);
			}
		}

		TDate firstTDate;
		if(GtoStringToDate(firstDate,
			   &firstTDate) == FAILURE)
		{
			throw KException(" not a valid issue date.", routine);
		}
		

		TDate endTDate;
		if(endDate == 0)
		{
			endTDate = firstTDate + 10*365;
		}
		else
		{
			endTDate = endDate;
		}
		
		TDateInterval payInterval;
		/* determining the payment interval */
		if(frequency == 0)
		{
			if(GtoDateSubtract(endTDate,
							   firstTDate,
							   &payInterval) == FAILURE)
			{
				throw KException("Failed to creat payInterval.", routine);
			}
		}
		else if(frequency == 365)
		{
			if(GtoMakeDateInterval(1,
							  'D',
								&payInterval) == FAILURE)
			{
				throw KException("Failed to creat payInterval.", routine);
			}
		}
		else
		{
			if(GtoFreqAndTypeToInterval(frequency,
										'N',
										&payInterval) == FAILURE)
			{
				throw KException("Failed to creat payInterval.", routine);
			}
		}
		
		double principal = 0.0;
		TDate idate= firstTDate;
		TDate nextDate;
		double yearFract, amount;
		do
		{
			
			GtoDtFwdAny(idate, &payInterval,&nextDate);
			if(GtoDayCountFraction(idate,nextDate,dayCount, &yearFract) == FAILURE)
			{
				throw KException("day fract failed.\n");	
			}

			amount = yearFract * numrate;
			insert(value_type(KDate(nextDate),amount));

			idate = nextDate;
    
		}while( idate <= endTDate);

	}
done:
	;
}


DDMap::DDMap(int numPoints, KDate *dates, double *rates)
{

	static char routine[] = "DDMap::DDMap";

	for (int i = 0; i<numPoints; i++)
	{
		insert(value_type(dates[i],rates[i]));
	}
}


double DDMap::get_amount(KDate currDate) const
{

	const_iterator iter;
	double amount;

	iter = find(currDate);
	if(iter == end())
	{
		amount = 0.0;
	}
	else
	{
		amount = iter->second;
	}

	return amount;
}

double DDMap::get_amount(KDate date1, KDate date2) const
{
	const_iterator iter1, iter2, iter;
	double amount =0.0;

	iter1 = lower_bound(date1);
	iter2 = lower_bound(date2);


	iter = iter1;
	while(iter != iter2)
	{
		amount += get_amount(iter->first);
		iter++;
	}

	return amount;
}

double DDMap::get_average_amount(KDate date1, KDate date2) const
{
	const_iterator iter1, iter2, iter, iterp;
	double average =0.0, totalYear, fractYear;

	iter1 = lower_bound(date1);
	iter2 = lower_bound(date2);

	if(iter1 != iter2)
	{
		iter = iter1;
		iterp = iter;
		iterp--;
		totalYear = GetYearDiff(iterp->first,iter->first);
		fractYear = GetYearDiff(date1,iter->first);
		average = fractYear/totalYear * get_amount(iter->first);

		iter++;
		while(iter != iter2)
		{
			average += get_amount(iter->first);
			iter++;
		}

		// take case the last odd amount
		iterp = iter;
		iterp--;
		totalYear = GetYearDiff(iterp->first,iter->first);
		fractYear = GetYearDiff(iterp->first,date2);
		average += fractYear/totalYear * get_amount(iter->first);
		average = average/(date2-date1)*365;
	}

		

		return average;
}


DDMap DDMap::get_subMap(KDate date1, KDate date2) const
{
	DDMap subMap;
	const_iterator iter1, iter2, iter;


	iter1 = lower_bound(date1);
	iter2 = lower_bound(date2);
	iter = iter1;
	while(iter != iter2)
	{
		subMap.insert(value_type(iter->first,iter->second));
		iter++;
	}

	return subMap;
}


KVector(KDate) DDMap::get_datelist(KDate date1, KDate date2) const 
{
	KVector(KDate) dates;
	DDMap::const_iterator iter1, iter2, iter;

	if(date1 == 0)
	{
		iter1 = begin();
	}
	else
	{
		iter1 = lower_bound(date1);
	}

	if(date2 == 0)
	{
		iter2 = end();
	}
	else
	{
		iter2 = lower_bound(date2);
	}
	

	for (iter = iter1; iter != iter2; iter++)
	{
		dates.push_back( iter->first);
	}
	return dates;
}
	
