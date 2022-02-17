/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file integrals.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_Integrals_H
#define _GP_CF_Integrals_H


#include <glob/firsttoinc.h>
#include "gpbase/port.h"



CC_BEGIN_NAMESPACE(ARM)


//////////////////////////////////////////////////////////////////////////
///
///				Oscilatory Integrals Functor
///
//////////////////////////////////////////////////////////////////////////




class OscillatoryIntegral
{
public: 
	enum {AUTOADAPTATIF,MANUAL,CONTROLED,NEWADAPTATIF, NEWMULTADAPTATIF};
	OscillatoryIntegral(vector<double> blist, int ptnb )	// constructor of the manual mode
		: boundaryList(blist),legendrePtnb(ptnb) ,runningMode(MANUAL){Stage_Nb=boundaryList.size();}
	OscillatoryIntegral(int ptnb1,int ptnbN,int NbStage,int NbOscill,double speed, double prec)	// constructor of the controled mode
		: legendrePtnb(ptnbN),legendrePtnb_FirstStage(ptnb1),
			Stage_Nb(NbStage),OscillationPerStage_Nb(NbOscill),oscillationspeed(speed),
			SpecifiedPrecision(prec),runningMode(CONTROLED){}
	OscillatoryIntegral(int autonb,double startpt, int ptnb )				// constructor of the auto adpatatif mode
		: autoAdaptatifComponentnb(autonb),legendrePtnb(ptnb),startpoint(startpt) ,runningMode(AUTOADAPTATIF){}
	OscillatoryIntegral(double lowerbound, int nbpt, double startsearch = 10.) 
		: legendrePtnb(nbpt), startpoint(lowerbound), startSearch(startsearch), runningMode(NEWADAPTATIF){};
	OscillatoryIntegral(double lowerbound, int nbpt, int nbstage, double searchinc)
		: legendrePtnb(nbpt), startpoint(lowerbound), startSearch(searchinc), Stage_Nb(nbstage), runningMode(NEWMULTADAPTATIF){};

	virtual double mapper(double k)=0;
	virtual double inversemapper(double S)=0;
	virtual double mapperjacobian(double k)=0;
	virtual double integrand(double k)=0;
	virtual double oscillatorycomponant(double k)=0; /// this is w(k) where 
				/// we assume that the integrand is like Re[exp(u(k) + i w(k) )]
				/// it is used only in the autoadaptatif case
	double value();
	vector<double> get_boundaryList()
	{ 
		return boundaryList;
	}
	int get_legendrePtnb()
	{ 
		return legendrePtnb;
	}
	int get_autoAdaptatifComponentnb()
	{ 
		return autoAdaptatifComponentnb;
	}

private:
	vector<double> boundaryList;
	int legendrePtnb;
	int autoAdaptatifComponentnb;
	int runningMode;
	double startpoint; /// only for the autoadaptatif case, starting value of the integral

	/// for the controled case
	int legendrePtnb_FirstStage;
	int Stage_Nb;
	int OscillationPerStage_Nb;
	double oscillationspeed;
	double SpecifiedPrecision;
	double startSearch;

	double getUpperBound();
	ARM_GP_Vector getBounds();
};


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
