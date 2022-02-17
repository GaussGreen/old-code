/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 */
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/ostringstream.h"
#include "gpbase/stringmanip.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/utilityport.h"


#include "gpclosedforms/tarnProxy.h"

#include "gpcalib/densityfunctors.h"
#include "gpnumlib/ran2.h"
#include "gpnumlib/gaussiananalytics.h"



CC_BEGIN_NAMESPACE(ARM)


////////////////////////////////////////////////////
///	Class  : ARM_TarnProxy
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_TarnProxy::ARM_TarnProxy(	const ARM_GP_Vector& resetDates,
					const ARM_GP_Vector& df,
					const ARM_GP_Vector& levprec,
					const ARM_GP_Vector& lev,
					const ARM_GP_Vector& fix,
					const ARM_GP_Vector& cap,
					const ARM_GP_Vector& floor,
					const ARM_GP_Vector& fees,
					const ARM_GP_Vector& dcf,
					double target,
					bool globalcap,
					bool globalfloor)
:	itsResetDates(resetDates),
	itsDiscountFactors(df),
	itsLevPrec(levprec),
	itsLev(lev),
	itsAdd(fix),
	itsCap(cap),
	itsFloor(floor),
	itsFees(fees),
	itsDcf(dcf),
	itsTarget(target),
	itsGlobalCap(globalcap),
	itsGlobalFloor(globalfloor),
	itsGrid(0)
{
	if (itsResetDates.size()!=itsDiscountFactors.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_TarnProxy::constructor: pb with size");
	if (itsLevPrec.size()!=itsDiscountFactors.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_TarnProxy::constructor: pb with size");
	if (itsLevPrec.size()!=itsLev.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_TarnProxy::constructor: pb with size");
	if (itsAdd.size()!=itsLev.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_TarnProxy::constructor: pb with size");
	if (itsAdd.size()!=itsCap.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_TarnProxy::constructor: pb with size");
	if (itsFloor.size()!=itsCap.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_TarnProxy::constructor: pb with size");
	if (itsFloor.size()!=itsFees.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_TarnProxy::constructor: pb with size");
	if (itsDcf.size()!=itsFees.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_TarnProxy::constructor: pb with size");

	int nb			= 400;
	double stdDev	= 5.;
	itsGrid = ARM_GP_VectorPtr(new ARM_GP_Vector(2*nb+1));
	for(int i=0;i<2*nb+1;i++)
		(*itsGrid)[i]=-stdDev+i*stdDev/nb;
	
	
}

void ARM_TarnProxy::Build(const ARM_GP_Vector& fwds,const vector< ARM_DensityFunctor* >& densityFunctors)
{
	if (itsResetDates.size()!=fwds.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_TarnProxy::Build: pb with size");
	if (itsResetDates.size()!=densityFunctors.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_TarnProxy::Build: pb with size");

	int n = densityFunctors.size();
	int i;
	itsFunc.resize(n);
	ARM_GP_VectorPtr proba=ARM_GP_VectorPtr(new ARM_GP_Vector(itsGrid->size()));
	for (i =0;i<itsGrid->size();i++)
		(*proba)[i] = ARM_GaussianAnalytics::cdfNormal((*itsGrid)[i]);
	for (i =0;i<n;i++)
		itsFunc[i] = densityFunctors[i]->Quantile(proba,fwds[i],itsResetDates[i]);
}

double ARM_TarnProxy::Index(int i,double sim) const
{
	int nb			= 400;
	double stdDev	= 5.;
	int j;
	if (sim<=-stdDev)
		return itsFunc[i]->Elt(0);
	else if (sim>=stdDev)
		return itsFunc[i]->Elt(2*nb);
	else
	{
		double storageDx = stdDev/nb;
		//standard interpol
		j = (int)ceil((sim - (-stdDev)) / storageDx) ;
		double nextNumMethState = (*itsGrid)[j];
		double prevNumMethState = (*itsGrid)[j-1];
		
		return (  (sim - prevNumMethState) * itsFunc[i]->Elt(j)
						  + (nextNumMethState - sim) * itsFunc[i]->Elt(j-1) ) 
						  / (nextNumMethState - prevNumMethState);
	}
}

void ARM_TarnProxy::Price(double correl, int simul)
{
	if (correl<0)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_TarnProxy::Price: correl positive only");
	double beta =0;
	int n = itsResetDates.size();
	if (n>0)
		beta=-log(correl)/fabs(itsResetDates[0]-itsResetDates[1]);

	ARM_GP_T_TriangularMatrix<double> cov(n);

	int i,j,k;
	for (i=0;i<n;i++)
		for (j=0;j<=i;j++)
			cov(i,j)=exp(-beta*fabs(itsResetDates[i]-itsResetDates[j]));

	cov.CholeskyDecompose();
	ARM_GP_Vector vecIndep(n);
	ARM_GP_Vector vecReal(n);

	itsStdCpn.resize(n,0.);
	itsStdFund.resize(n,0.);
	itsRbCpn.resize(n,0.);
	itsRbFund.resize(n,0.);
	itsGC.resize(n,0.);
	itsGF.resize(n,0.);
	itsProba.resize(n,0.);

	double index,cpn,sumcpn,target;
	bool wasAlive;
	bool isAlive;

	ARM_RandUniform_NRRan2 gen(-156);

	for (k=0;k<simul;k++)
	{
		sumcpn = 0.;
		cpn = 0.;
		gen.draw(vecIndep);
		target = itsTarget;
		wasAlive=true;
		isAlive=true;

		for (i=0;i<n;i++)
			vecIndep[i]=ARM_GaussianAnalytics::cdfNormal_Inv(vecIndep[i]);
		for (i=0;i<n;i++)
		{
			vecReal[i]=0.;
			for (j=0;j<=i;j++)
				vecReal[i]+=cov(i,j)*vecIndep[j];

			index		= Index(i,vecReal[i]);
			cpn			= CC_Min<double>(CC_Max<double>(itsLevPrec[i]*cpn+itsLev[i]*index+itsAdd[i],itsFloor[i]),itsCap[i]);
			itsStdCpn[i]  += cpn*itsDcf[i]*itsDiscountFactors[i];
			itsStdFund[i] += itsFees[i]*itsDiscountFactors[i];

			if (!isAlive)
			{
				itsRbCpn[i]  += cpn*itsDcf[i]*itsDiscountFactors[i];
				itsRbFund[i] += itsFees[i]*itsDiscountFactors[i];
			}

			sumcpn     += itsDcf[i]*cpn;
			isAlive		= sumcpn<itsTarget;

			if (!isAlive && wasAlive)
			{
				itsGC[i]  += (itsDcf[i]*cpn-target)*itsDiscountFactors[i];
				itsProba[i]  += 1.;
			}

			wasAlive	= isAlive;
			itsGF[i]	+= 0.;
			target		= CC_Max<double>(target-itsDcf[i]*cpn,0.);
		}
		itsGF[n-1] += target*itsDiscountFactors[n-1];
	}

	itsStdCpnLeg = 0.;
	itsStdFundLeg = 0.;
	itsRbCpnLeg = 0.;
	itsRbFundLeg = 0.;
	itsGCLeg = 0.;
	itsGFLeg = 0.;
	for (i=0;i<n;i++)
	{
		itsStdCpn[i]	/= simul;
		itsStdFund[i]	/= simul;
		itsRbCpn[i]		/= simul;
		itsRbFund[i]	/= simul;
		itsGC[i]		/= simul;
		itsGF[i]		/= simul;
		itsProba[i]		/= simul;
		itsStdCpnLeg	+= itsStdCpn[i];
		itsStdFundLeg	+= itsStdFund[i];
		itsRbCpnLeg		+= itsRbCpn[i];
		itsRbFundLeg	+= itsRbFund[i];
		itsGCLeg		+= itsGC[i];
		itsGFLeg		+= itsGF[i];
	}
}

double ARM_TarnProxy::GetPrice(int nbPayoff) const
{
	switch( nbPayoff )
    {
        case 1:
			return itsStdFundLeg;

        case 2:
			return itsStdCpnLeg;

        case 3:
			return itsRbFundLeg;

        case 4:
			return itsRbCpnLeg;

		case 5:
			return itsGCLeg;

		case 6:
			return itsGFLeg;

        default :
			double res = (itsStdFundLeg-itsRbFundLeg)-(itsStdCpnLeg-itsRbCpnLeg);
			if (itsGlobalCap)
				res+=itsGCLeg;
			if (itsGlobalFloor)
				res-=itsGFLeg;
			return res;
    }
}

////////////////////////////////////////////////////
///	Class   : ARM_TarnProxy
///	Routines: toString
///	Returns : void
///	Action  : toString
////////////////////////////////////////////////////
string ARM_TarnProxy::toString(const string& indent,const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "ARM_TarnProxy\n";
    os << indent << "----------------------------\n";
	os << indent << "Std FundLeg : \t" << CC_NS(std,fixed) << itsStdFundLeg << " \n";
	os << indent << "Std CpnLeg : \t" << CC_NS(std,fixed) << itsStdCpnLeg << " \n";
	os << indent << "Rbt FundLeg : \t" << CC_NS(std,fixed) << itsRbFundLeg << " \n";
	os << indent << "Rbt CpnLeg : \t" << CC_NS(std,fixed) << itsRbCpnLeg << " \n";
	os << indent << "GlobalCap : \t" << CC_NS(std,fixed) << itsGCLeg << " \n";
	os << indent << "GlobalFloor : \t" << CC_NS(std,fixed) << itsGFLeg << " \n";
	os << "\n\n";

	os << indent << "Proba   \t";
	os << indent << "FundLeg \t";
	os << indent << "CpnLeg  \t";
	os << indent << "RbtFund \t";
	os << indent << "RbtCpn  \t";
	os << indent << "GlbCap  \t";
	os << indent << "GlbFloor\n";
	
    int n = itsResetDates.size();
	for (int k=0;k<n;k++)
	{
		os << indent << CC_NS(std,fixed) << itsProba[k] << "\t";
		os << indent << CC_NS(std,fixed) << itsStdFund[k] << "\t";
		os << indent << CC_NS(std,fixed) << itsStdCpn[k] << "\t";
		os << indent << CC_NS(std,fixed) << itsRbFund[k] << "\t";
		os << indent << CC_NS(std,fixed) << itsRbCpn[k] << "\t";
		os << indent << CC_NS(std,fixed) << itsGC[k] << "\t";
		os << indent << CC_NS(std,fixed) << itsGF[k] << "\n";
	}

	return os.str();	
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/