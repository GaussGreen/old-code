#include "icm_distribloss.h" 

#include "ICMKernel\util\icm_macro.h"
#include <deque>
// 
const double ICM_DistribLoss::prec(1e-4); 

//	--------------------------------------------------------------------------
//	friend 
std::ostream& 
operator<<(std::ostream&o,const ICM_DistribLoss&ref) 
{
	o<<"\tDistribLoss ___________________________"<<std::endl ;
	o<<"\tYearTerm \tValue "<<std::endl ;

	ICM_DistribLoss::map_t ExpectedLosses_disc_unify;
	ref.UnifyELs(ExpectedLosses_disc_unify);	

	ICM_DistribLoss::map_t::const_iterator it = ref.itsExpectedLosses.begin(); 


	while(it!=ref.itsExpectedLosses.end()) 
	{
		o<<"\t"<<it->first.m_value<<"\t\t"<<ref.InterpolEL(it->first.m_value)<<std::endl ; 
		++it; 
	}
	o<<std::endl ;
	return o; 
}
//	--------------------------------------------------------------------------
// virtual 
void 
ICM_DistribLoss::View(char*id, FILE* ficOut)
{
	FILE* fOut;
	char  fOutName[200];

	if ( ficOut == NULL )
	{
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w"); 
	}
	else	fOut = ficOut;
	std::stringstream sstr ;
	sstr<<*this; 
	fprintf(fOut,"%s\n",sstr.str().c_str()); 
	if ( ficOut == NULL )fclose(fOut);
}


double 
ICM_DistribLoss::InterpolEL_TS(double yearterm) const 
{
	map_t ExpectedLosses = itsExpectedLosses; 
	map_t ExpectedLosses_disc = itsExpectedLosses_disc;
	map_t ExpectedLosses_disc_unify;
	UnifyELs(ExpectedLosses_disc_unify);

	// copied from ICM_Pricer_Security. 
	map_t::const_iterator it; 
	map_t::const_iterator it_disc; 
	map_t::const_iterator it_unify; 

	std::vector<double> yearterms;
	std::vector<double> values;
	std::vector<double> values_disc;
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ jumps;

	bool jump;

	for (it=ExpectedLosses.begin(),it_unify=ExpectedLosses_disc_unify.begin();
		it!=ExpectedLosses.end(),it_unify!=ExpectedLosses_disc_unify.end();
		++it,++it_unify)
	{
		yearterms.push_back(it->first.m_value);
		values.push_back(it->second);
		values_disc.push_back(it_unify->second);

		jump = false;
		
		if (ExpectedLosses_disc.find(it->first.m_value)!=ExpectedLosses_disc.end())
		{jump=true;}

		jumps.push_back(jump);
	}

	int i=0,j=0;
	if (ExpectedLosses.find(yearterm) == ExpectedLosses.end()) //non trouvé, on insert les points interpolés 
	{
		double valueinter = LinearVectorInterpol(yearterms,values,yearterm);
		ExpectedLosses[yearterm]=valueinter;

		valueinter = LinearVectorInterpol(yearterms,values_disc,yearterm);
		ExpectedLosses_disc_unify[yearterm]=valueinter;

		// on regenere les vecteurs impactés
		yearterms.clear();
		values.clear();
		values_disc.clear();
		jumps.clear();

		for (it=ExpectedLosses.begin(),it_unify=ExpectedLosses_disc_unify.begin();
		it!=ExpectedLosses.end(),it_unify!=ExpectedLosses_disc_unify.end();
		++it,++it_unify)
		{
		yearterms.push_back(it->first.m_value);
		values.push_back(it->second);
		values_disc.push_back(it_unify->second);

		jump = false;
		
		if (ExpectedLosses_disc.find(it->first.m_value)!=ExpectedLosses_disc.end())
		{jump=true;}

		jumps.push_back(jump);
		}
	}

	double CumDiffEl = 0.;
	i=0;
	
	while (i<yearterms.size()-1)
	{
		if (((yearterm>yearterms[i+1])||CHECK_EQUAL(yearterm,yearterms[i+1])) 
			&& (jumps[i]==false))
			CumDiffEl+=values[i+1]-values[i];
		else if (((yearterm>yearterms[i+1])||CHECK_EQUAL(yearterm,yearterms[i+1])) 
			&& (jumps[i]==true))
			CumDiffEl+=values[i+1]-values_disc[i];
		else if ((yearterm<yearterms[i+1])||!CHECK_EQUAL(yearterm,yearterms[i+1])) 
			{break;}
		i++;
	}

	return (CumDiffEl);
} 