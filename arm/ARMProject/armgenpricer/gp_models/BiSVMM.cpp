
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/bisvmm.h"
#include "gpmodels/bgmsv1f.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/comparisonfunctor.h"
#include "gpbase/datestrip.h"
#include "gpbase/vectormanip.h"
#include "gpbase/stringmanip.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/discretisationscheme.h"

#include "gpclosedforms/normal.h"
#include "gpnumlib/gaussiananalytics.h"

CC_BEGIN_NAMESPACE( ARM )

ARM_BiSVMM::ARM_BiSVMM(const ARM_ModelNameMap& modelNameMap, const ARM_CurveMatrix& correlationMatrix) :
ARM_MultiAssetsModel(&modelNameMap, &correlationMatrix), itsNbEffectiveReset(0), itsEigenValues(0)
{
	if(modelNameMap.size() != 2)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : must have two models to build bi-svmm");
	}

	itsModel1 = dynamic_cast<ARM_BGMSV1F*>(&*(modelNameMap)[0]->Model() );
	itsModel2 = dynamic_cast<ARM_BGMSV1F*>(&*(modelNameMap)[1]->Model() );

	if(itsModel1 == NULL || itsModel2 == NULL)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : model should be a BGMSV1F");

	if(!GetCorrelMatrix())
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : please provide a correlation");
	}

	if(	GetCorrelMatrix()->rows() != 4 || GetCorrelMatrix()->cols() != 4)
	{
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : correlation should be 4x4!");
	}

	if(itsModel1->GetProxyStatus() || itsModel2->GetProxyStatus())
	{
		itsModel1->SetProxyStatus(true);
		itsModel2->SetProxyStatus(true);
	}

	if(itsModel1->GetResetTimes() != itsModel2->GetResetTimes())
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : two models should have same reset times");
	}

	itsResetTimes = itsModel1->GetResetTimes();
}

ARM_BiSVMM::ARM_BiSVMM(const ARM_BiSVMM& rhs) : ARM_MultiAssetsModel(rhs), itsNbEffectiveReset(rhs.itsNbEffectiveReset), itsResetTimes(rhs.itsResetTimes)
{
	DuplicateCloneablePointorVectorInPlace<std::vector<double>>( rhs.itsEigenValues, itsEigenValues );

	itsModel1 = dynamic_cast<ARM_BGMSV1F*>(&*(*GetModelMap())[0]->Model() );
	itsModel2 = dynamic_cast<ARM_BGMSV1F*>(&*(*GetModelMap())[1]->Model() );	
}

ARM_BiSVMM::~ARM_BiSVMM()
{
	DeletePointorVector<std::vector<double>>( itsEigenValues );
}

void ARM_BiSVMM::ModelStateLocalVariancesAndStdDev(const std::vector<double>& timeSteps)
{
	ARM_MatrixVector modStateLocalVars;
	ARM_MatrixVector auxLocalCov;

	if(itsModel1->GetProxyStatus())
	{
		double lastTimeStep = timeSteps[timeSteps.size() - 1];
		itsNbEffectiveReset = 0;
		while(!(itsResetTimes[itsNbEffectiveReset] > lastTimeStep))
		{
			itsNbEffectiveReset ++;
			if(itsNbEffectiveReset == itsResetTimes.size()) break;
		}
	}
	else
		itsNbEffectiveReset = itsResetTimes.size();

	itsModel1->SetNbEffectiveReset(itsNbEffectiveReset);
	itsModel2->SetNbEffectiveReset(itsNbEffectiveReset);

	ModelStateLocalVariances( timeSteps, auxLocalCov);
	//DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( auxLocalCov, itsModelLocalRealVar);
	
	int totalFwds		= itsNbEffectiveReset;
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;
	
	double fromTime		= timeSteps[0];
	double toTime;

	int i, j, k, startTimePos = 0;

	modStateLocalVars.resize((timeStepsSize-1)*(modelNb+1));

	DeletePointorVector<std::vector<double>>( itsEigenValues );
	itsEigenValues.resize((timeStepsSize-1)*(modelNb+1));

	int factorsNb = 2*itsNbEffectiveReset + 2;
	int nbFwds, idx = 0, currfactorsNb = factorsNb;
	int maxfactorsNb = 0;
	double minratio = ((ARM_ModelParamsBGMSV1F*) itsModel1->GetModelParams())->GetMinRatio();
	bool resizefactors = minratio > 1. || fabs(minratio - 1.) < K_DOUBLE_TOL ? false : true;

	// on commence par calculer toutes les matrices, et les vecteurs propres
	std::vector<double> * eigenValues = new std::vector<double>[timeStepsSize - 1];
	ARM_GP_Matrix ** ACPMatrix = new ARM_GP_Matrix * [timeStepsSize - 1];
	std::vector<double> factorsize(timeStepsSize-1, factorsNb);

	for(i = 0; i < timeStepsSize - 1; i++)
	{
		toTime  = timeSteps[i+1];

		/// get the first bigger time
		while(idx < itsResetTimes.size()
			&& itsResetTimes[idx] < toTime )
			++idx;

		nbFwds = totalFwds - idx;

		int esize = nbFwds * 2 + 2;

		eigenValues[i].resize(esize,0.);

		ACPMatrix[i] = ACPTransformation(auxLocalCov[i],eigenValues[i]);

		// calcul de la dimension effective
		if(resizefactors)
		{
			// variance total
			double varTotal = 0.;
			for(k = 0; k < esize; k++) 
			{
				if(eigenValues[i][k] > 0.) varTotal += eigenValues[i][k];
			}

			double varExpliq = 0.;
			for(k = 0; k < esize; k++)
			{
				if(varExpliq > minratio * varTotal) break;
				if(eigenValues[i][k] > 0.)
				{
					varExpliq += eigenValues[i][k];
				}
			}

			factorsize[i]= k;
			maxfactorsNb = maxfactorsNb > k ? maxfactorsNb : k;
		}
		else
		{
			maxfactorsNb = factorsNb;
		}
	}

	((ARM_ModelParamsBGMSV1F*) itsModel1->GetModelParams())->SetFactorCount(maxfactorsNb);
	((ARM_ModelParamsBGMSV1F*) itsModel2->GetModelParams())->SetFactorCount(0);
	
	factorsNb = FactorCount();

	idx = 0;
	for(i = 0; i < timeStepsSize - 1 ;++i)
	{
		/// initalize the toTime
		toTime  = timeSteps[i+1];

		while(idx < itsResetTimes.size()
			&& itsResetTimes[idx] < toTime )
			++idx;

		nbFwds = totalFwds - idx;

		int esize = nbFwds * 2 + 2;

		currfactorsNb = factorsize[i] > esize ? esize : factorsize[i];

		/// initialize everything
		modStateLocalVars[i] = new ARM_GP_Matrix(esize, currfactorsNb);

		itsEigenValues[i] = new std::vector<double>(currfactorsNb);
		
		int maxRank = -1;
		for(k=0;k<currfactorsNb;++k)
		{
			if(eigenValues[i][k] < K_NEW_DOUBLE_TOL)
			{
				eigenValues[i][k] = 0.;
				if (maxRank == -1)
					maxRank = k;
			}
		}

		if (maxRank == -1) maxRank = currfactorsNb;
		size_t effectiveRank = ( maxRank < currfactorsNb ? maxRank : currfactorsNb);

		//rescaling
		
		for (j = 0;j < esize; j++)
		{
			double sum=0.;
			for(k=0;k<effectiveRank-1;++k){
				sum+=(*ACPMatrix[i])(j,k)*(*ACPMatrix[i])(j,k)*eigenValues[i][k];
			}
			double sgn = ((*ACPMatrix[i])(j,k)>0.?1.:-1.);
			double res = (*auxLocalCov[i])(j,j)-sum;
			if (res>K_NEW_DOUBLE_TOL)
				(*ACPMatrix[i])(j,k)=sgn*sqrt(res/eigenValues[i][k]); 
		}
		
		for(k = 0; k < currfactorsNb; ++k)
		{
			(*itsEigenValues[i])(k) = eigenValues[i][k];
			for(j = 0; j < esize; ++j)
			{	
				(*modStateLocalVars[i])(j,k) =	(*ACPMatrix[i])(j,k);
			}
		}

		fromTime= toTime;
	}

	delete [] eigenValues;
	for(i = 0; i < timeStepsSize - 1; i++) delete [] ACPMatrix[i];
	delete [] ACPMatrix;

	/// set the result
	SetModelStateLocalVars(modStateLocalVars);
}

void ARM_BiSVMM::ModelStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances) const
{
	int totalFwds		= itsNbEffectiveReset;
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;
	
	double fromTime		= timeSteps[0];
	double toTime;


#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ModelStateLocalVariances: localVariances.size() != offsetIndex" );
#endif

#if defined(__GP_STRICT_VALIDATION)
	if( fromTime != 0.0 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ModelStateLocalVariances: first time step != 0" );
#endif

    localVariances.resize((timeStepsSize-1)*(modelNb+1));
	
	int i,j,k,jj,kk, idx = 0;

	int nbFwds;

	for(i = 0;i < timeStepsSize - 1 ; ++i)
	{
		/// initalize the toTime
		toTime  = timeSteps[i+1];

		/// initialize everything
		
		/// get the first bigger time
		while(idx< itsResetTimes.size()
			&& itsResetTimes[idx] < toTime )
			++idx;

		nbFwds = totalFwds - idx;

		if(nbFwds <= 0)
		{
			localVariances[i] = new ARM_GP_Matrix(1,1);
			(*localVariances[i])(0,0) = 0.;
			continue;
		}

		localVariances[i] = new ARM_GP_Matrix(2 * nbFwds + 2, 2 * nbFwds + 2);	

		double rhoj, rhok, rhojk, fac;
		
		// les correls sur le premier modèle
		for(j = 0; j < nbFwds; j++)
		{
			rhoj	= ((ARM_ModelParamsBGMSV1F*) itsModel1->GetModelParams())->GetRho(j + idx);

			(*localVariances[i])(j,j) = 1.;
			(*localVariances[i])(j,nbFwds) = 0.;

			for(k = j+1; k < nbFwds; k++)
			{
				rhok	= ((ARM_ModelParamsBGMSV1F*) itsModel1->GetModelParams())->GetRho(k + idx);
				rhojk	= ((ARM_ModelParamsBGMSV1F*) itsModel1->GetModelParams())->RateRateCorrel(toTime, itsResetTimes[k + idx], k + idx, itsResetTimes[j + idx], j + idx);

				//fac		= (rhoj * rhok + sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok)));
				fac		= sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok));
				
				if(fabs(fac) < K_DOUBLE_TOL)
				{
					rhok = rhok > 0.9999 - K_DOUBLE_TOL ? rhok - 0.0001 : rhok + 0.0001;
					//fac	= (rhoj * rhok + sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok)));
					fac	= sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok));
				}

				(*localVariances[i])(j,k) = (*localVariances[i])(k,j) = (rhojk - rhoj*rhok) / fac;
			}
		}
		(*localVariances[i])(nbFwds,nbFwds) = 1.;

		// les correls sur le second modèle
		for(j = 0, jj = nbFwds + 1; j < nbFwds; j++, jj++)
		{
			rhoj	= ((ARM_ModelParamsBGMSV1F*) itsModel2->GetModelParams())->GetRho(j + idx);

			(*localVariances[i])(jj,jj) = 1.;
			(*localVariances[i])(jj,2*nbFwds+1) = 0.;

			for(k = j+1, kk = j+1+nbFwds+1; k < nbFwds; k++, kk++)
			{
				rhok	= ((ARM_ModelParamsBGMSV1F*) itsModel2->GetModelParams())->GetRho(k + idx);
				rhojk	= ((ARM_ModelParamsBGMSV1F*) itsModel2->GetModelParams())->RateRateCorrel(toTime, itsResetTimes[k + idx], k + idx, itsResetTimes[j + idx], j + idx);

				//fac		= (rhoj * rhok + sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok)));
				fac		= sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok));
				
				if(fabs(fac) < K_DOUBLE_TOL)
				{
					rhok = rhok > 0.9999 - K_DOUBLE_TOL ? rhok - 0.0001 : rhok + 0.0001;
					//fac	= (rhoj * rhok + sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok)));
					fac	= sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok));
				}

				(*localVariances[i])(jj,kk) = (*localVariances[i])(kk,jj) = (rhojk - rhoj*rhok) / fac;
			}
		}
		(*localVariances[i])(2*nbFwds+1,2*nbFwds+1) = 1.;

		ARM_GP_Matrix correlMatrix = GetCorrelMatrix()->Interpolate(timeSteps[i]);

		double rho12jk, rho1j, rho21k, rho2k, rho12j, rho12;

		rho12	= correlMatrix(2,3); // correl entre les deux process de vols
		rho12j	= correlMatrix(0,3); // correl entre le taux 1 et la vol 2
		rho21k	= correlMatrix(1,2); // correl entre le taux 2 et la vol 1
		rho12jk	= correlMatrix(0,1); // correl entre le taux 1 et le taux 2

		// les correls croisées
		for(j = 0; j < nbFwds; j++)
		{
			rho1j = ((ARM_ModelParamsBGMSV1F*) itsModel1->GetModelParams())->GetRho(j + idx);

			for(k = 0; k < nbFwds; k++)
			{
				rho2k = ((ARM_ModelParamsBGMSV1F*) itsModel2->GetModelParams())->GetRho(k + idx);

				fac	= sqrt((1. - rho1j*rho1j)*(1. - rho2k*rho2k));

				if(fabs(fac) < K_DOUBLE_TOL)
				{
					rho2k = rho2k > 0.9999 - K_DOUBLE_TOL ? rho2k - 0.0001 : rho2k + 0.0001;
					fac = sqrt((1. - rho1j*rho1j)*(1. - rho2k*rho2k));
				}

				// correl entre les des taux
				(*localVariances[i])(j, k + nbFwds + 1) = (*localVariances[i])(k + nbFwds + 1,j) = 
					(rho12jk - rho1j*rho21k - rho2k*rho12j + rho1j*rho2k*rho12) / fac;

				// correl entre taux 2 et vol 1
				fac = fabs(fabs(rho2k) - 1.) < K_DOUBLE_TOL ? 0.0000001 : sqrt(1. - rho2k*rho2k);

				(*localVariances[i])(nbFwds, k + nbFwds + 1) = (*localVariances[i])(k + nbFwds + 1, nbFwds) = 
					(rho21k - rho2k*rho12) / fac;
				
			}

			// correl entre taux 1 et vol 2
			fac = fabs(fabs(rho1j) - 1.) < K_DOUBLE_TOL ? 0.0000001 : sqrt(1. - rho1j*rho1j);

			(*localVariances[i])(j, nbFwds + nbFwds + 1) = (*localVariances[i])(nbFwds + nbFwds + 1,j) = 
				(rho12j - rho1j*rho12) / fac;
		}

		// correl vol 1 vol 2
		(*localVariances[i])(nbFwds, nbFwds + nbFwds + 1) = (*localVariances[i])(nbFwds + nbFwds + 1, nbFwds) = rho12;

		fromTime = toTime;
	}
}

void ARM_BiSVMM::NumMethodStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const
{
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	
    localVariances.resize((timeStepsSize-1)*(modelNb+1));
	
	for(int i = 0; i < timeStepsSize-1 ; ++i)
	{
		int factorsNb = itsEigenValues[i]->size();

		/// initialize everything
		localVariances[i] = new ARM_GP_Matrix(factorsNb,factorsNb,0.);

		for (int k = 0; k < factorsNb; k++)
			(*localVariances[i])(k,k) = (*itsEigenValues[i])[k];
	}
}

void ARM_BiSVMM::NumMethodStateGlobalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& globalVariances) const
{
	int totalFwds		= itsNbEffectiveReset;
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( globalVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BiSVMM::NumMethodStateGlobalVariances: globalVariances.size() != offsetIndex" );
#endif
    globalVariances.resize(timeStepsSize*(modelNb+1));
	for (int i=0;i<timeStepsSize;i++)
	{
		globalVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(2*totalFwds+2,1.0);
	}
}

void ARM_BiSVMM::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states, int timeIndex) const
{
	int fwdsNb		= itsNbEffectiveReset;
	int iLast		= fwdsNb-1;
	double currTime = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	
	int iFirst  = 0;
	while( iFirst < fwdsNb 	&& itsResetTimes[iFirst] < nextTime)
	{
		++iFirst;
	}
	
	int aliveFwds	= iLast - iFirst + 1;

	const ARM_MatrixVector& modelLocalStdev	= GetModelStateLocalVars();

	int eigensNb	= (*modelLocalStdev[timeIndex]).cols();
	int statesNb	= states->size();
	int nb			= aliveFwds + 1;

	itsModel1->SetModelNb(0);
	itsModel2->SetModelNb(nb);

	// corrélation des gaussiennes pour chaque modèle
	ARM_GP_Matrix x1(statesNb, nb), x2(statesNb, nb);

	for(int n = 0; n < statesNb; n++)
	{
		for(int i = 0; i < nb; i++)
		{
			x1(n,i) = x2(n,i) = 0.;
			for(int j = 0; j < eigensNb; j++) 
			{
				x1(n,i) += states->GetNumMethodState(n,j) * (*modelLocalStdev[timeIndex])(i,j);
				x2(n,i) += states->GetNumMethodState(n,j) * (*modelLocalStdev[timeIndex])(i + nb,j);
			}
		}
	}

	// simulation dans chaque modèle
	itsModel1->MCFromToNextTime(states, timeIndex, x1);
	itsModel2->MCFromToNextTime(states, timeIndex, x2);
}

ARM_VectorPtr ARM_BiSVMM::ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const
{
	return itsModel1->ComputeNumeraireTimes(timeInfos);
}

string ARM_BiSVMM::toString(const string& indent, const string& nextIndent) const
{
 	const ARM_ModelNameMap& modelMap = *GetModelMap();

    CC_Ostringstream os;

    os << "\n\n";
    os << indent << "Bi BGM with heston stochastic volatility \n";
    os << indent << "----------------------------------------------\n\n";

	os << "Correlation matrix\n";
	os << GetCorrelMatrix()->toString(indent,nextIndent);

	os << endl;

    os << indent << "\n\n------> First SVBGM <------\n";
    os << itsModel1->toString(indent,nextIndent);

    os << indent << "\n\n------> Second SVBGM <------\n";
    os << itsModel2->toString(indent,nextIndent);

	return os.str();
}

CC_END_NAMESPACE()

