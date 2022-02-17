
#include "stdafx.h"
#include "threemoments.h"
#include "MlEqObjects.h"



/****************************************************************
**	class  : MlEqThreeMoment 
**	Routine: initialize 
**	Returns: void
**	Action : 
****************************************************************/




void MlEqThreeMoment::initialize(const std::vector<MlEqAssetHandle>  assets,
								 CVector& basketWeights,
								 CVector& asianWeights,
								 vector < DATE >& asianDates,
								 MlEqCorrelationMatrix& correlation,
								 MlEqConstDateHandle hDate,
								 vector < MlEqStrikeHandle >& pStrike)
{

	m_numberAssets	= assets.size();
	m_hDate = hDate;
	m_covariance.resize(m_numberAssets,m_numberAssets,asianDates.size() );
	m_fwds.resize(m_numberAssets,asianDates.size());



//	if ( pStrike.size() != m_numberAssets ){
//		throw("array of input strikes must coincide with number of assets");
//	}


	int i,j,n;
	double cov,volSi,volSj;

	for ( n = 0 ; n < asianDates.size(); n++ )
	{
		for ( i = 0 ; i < m_numberAssets; i++ )
		{
			volSi = assets[i]->GetNaturalVolatility(*pStrike[i],m_hDate->GetDate(),asianDates[n]) ;

			m_fwds[i][n] =  assets[i]->GetQuantoForward(m_hDate->GetDate(), asianDates[n], false) ;

			MlEqAssetHandle	FxNatural_i		= assets[i]->GetFxNaturalAsset();
			MlEqAssetHandle	FxPay_i			= assets[i]->GetFxPayAsset();
			MlEqAssetHandle	FxComposite_i	= assets[i]->GetFxCompositeAsset();

			for ( j = 0 ; j < m_numberAssets; j++ )
			{
				volSj = assets[j]->GetNaturalVolatility(*pStrike[j],m_hDate->GetDate(),asianDates[n]);

				cov =  volSi*volSj*correlation.GetCorrelation(assets[i]->GetName(),assets[j]->GetName());
				cov += assets[i]->GetStockCompositeFxVariance(*pStrike[i],m_hDate->GetDate(),asianDates[n]) ;
				cov -= assets[i]->GetStockNaturalFxVariance(*pStrike[i],m_hDate->GetDate(),asianDates[n]) ;
				cov += assets[j]->GetStockCompositeFxVariance(*pStrike[j],m_hDate->GetDate(),asianDates[n]) ;
				cov -= assets[j]->GetStockNaturalFxVariance(*pStrike[j],m_hDate->GetDate(),asianDates[n]) ;


				MlEqAssetHandle	FxNatural_j		= assets[j]->GetFxNaturalAsset();
				MlEqAssetHandle	FxPay_j			= assets[j]->GetFxPayAsset();
				MlEqAssetHandle	FxComposite_j	= assets[j]->GetFxCompositeAsset();

				cov += FxComposite_i->GetNaturalVolatility(*pStrike[i],m_hDate->GetDate(),asianDates[n])*
					   FxComposite_j->GetNaturalVolatility(*pStrike[j],m_hDate->GetDate(),asianDates[n])*
					   correlation.GetCorrelation(FxComposite_i->GetName(),FxComposite_j->GetName());

				cov -= FxComposite_i->GetNaturalVolatility(*pStrike[i],m_hDate->GetDate(),asianDates[n])*
					   FxNatural_j->GetNaturalVolatility(*pStrike[j],m_hDate->GetDate(),asianDates[n])*
					   correlation.GetCorrelation(FxComposite_i->GetName(),FxNatural_j->GetName());

				cov -= FxNatural_i->GetNaturalVolatility(*pStrike[i],m_hDate->GetDate(),asianDates[n])*
					   FxComposite_j->GetNaturalVolatility(*pStrike[j],m_hDate->GetDate(),asianDates[n])*
					   correlation.GetCorrelation(FxNatural_i->GetName(),FxComposite_j->GetName());

				cov += FxNatural_i->GetNaturalVolatility(*pStrike[i],m_hDate->GetDate(),asianDates[n])*
					   FxNatural_j->GetNaturalVolatility(*pStrike[j],m_hDate->GetDate(),asianDates[n])*
					   correlation.GetCorrelation(FxNatural_i->GetName(),FxNatural_j->GetName());

				m_covariance[i][j][n] = cov;

			}
		}
	}




	m_basketWeights = basketWeights;

	initialize(
				m_covariance,
				m_fwds,
				asianWeights,
				m_basketWeights,
				asianDates,
				hDate,
				pStrike);


}

/****************************************************************
**	class  : MlEqThreeMoment 
**	Routine: initialize 
**	Returns: void
**	Action : 
****************************************************************/




void MlEqThreeMoment::initialize(MlEqAsset& asset,
								 CVector& asianWeights,
								 vector < DATE >& asianDates,
								 MlEqCorrelationMatrix& correlation,
								 MlEqConstDateHandle hDate,
								 vector < MlEqStrikeHandle >& pStrike)
{
	const std::vector<MlEqAssetHandle>  assets = asset.GetAssets();

	const CVector bWeights = asset.GetWeights();

	CVector basketWeights = bWeights;

	initialize(assets,
				basketWeights,
				asianWeights,
				asianDates,
				correlation,
				hDate,
				pStrike);

}


void MlEqThreeMoment::initialize(
								 cTensor& covariance,
								 CMatrix& forwards,
								 CVector& asianWeights,
								 CVector& basketWeights,
								 vector < DATE >& asianDates,
								 MlEqConstDateHandle hDate,
								 vector < MlEqStrikeHandle >& pStrike)
{


	m_covariance	= covariance;
	m_fwds			= forwards;
	m_basketWeights = basketWeights;
	m_numberAssets	= basketWeights.getsize();
	m_asianWeights	= asianWeights;
	m_asianDates	= asianDates;
	m_hDate = hDate;


//	if ( pStrike.size() != m_numberAssets+1 )
//	{
//		throw("array of input strikes must coincide with number of assets+1 ");
//	}


	m_firstFutureIndex = 0;
	for (int i = 0 ; i < m_asianDates.size(); i++ )
	{
		if ( m_asianDates[i] <= m_hDate->GetDate() )
		{
			m_firstFutureIndex++;
		}
	}

    m_basketForward = 0;
	for (int ai = m_firstFutureIndex ; ai < m_asianWeights.getsize(); ai++ )
	{
		for ( int bi = 0 ; bi < m_numberAssets; bi++ )
		{
m_basketForward  += m_basketWeights[bi]*m_asianWeights[ai]*m_fwds[bi][ai];

		}
	}


	// start calculating effective strike: past settings get moved into strike

	m_effectiveStrikeAdjust = 0.0;
	for (int bi = 0; bi < m_numberAssets; bi++ )
	{
		for (int ai = 0 ; ai < m_firstFutureIndex; ai++ )
		{
				m_effectiveStrikeAdjust -= m_basketWeights[bi]*asianWeights[ai];
			
		}
	}



}



/****************************************************************
**	Class  : MlEqThreeMoment 
**	Routine: operator= 
**	Returns: MlEqThreeMoment&
**	Action : assignment operator
****************************************************************/

MlEqThreeMoment& MlEqThreeMoment::operator=( const MlEqThreeMoment& rhs)
{
	if(this == &rhs) // self-assignment check
		return *this;


	m_asianWeights	= rhs.m_asianWeights;
	m_fwds			= rhs.m_fwds;
	m_covariance	= rhs.m_covariance;
	m_numberAssets	= rhs.m_numberAssets;
	m_basketWeights	= rhs.m_basketWeights;

	m_hDate = rhs.m_hDate;
	m_basketForward = rhs.m_basketForward;
	return *this;
}



/****************************************************************
**	Class  : MlEqThreeMoment 
**	Routine: ~MlEqThreeMoment
**	Returns: nothing
**	Action : 
****************************************************************/

MlEqThreeMoment::~MlEqThreeMoment()
{
}


/****************************************************************
**	Class  : MlEqThreeMoment 
**	Routine: e
**	Returns: e
**	Action : 
****************************************************************/

double MlEqThreeMoment::e(int iasset, int jasset,
						  int iDate, int jDate) 
{
	int tau = (iDate < jDate) ? iDate : jDate;
	double dt = m_hDate->GetYearFraction(m_asianDates[tau]);
	double res =  exp(m_covariance[iasset][jasset][tau]* dt);
	return res;
}


/****************************************************************
**	Class  : MlEqThreeMoment 
**	Routine: calculateReducedMoments
**	Returns: vector of reduced moments
**	Action : 
****************************************************************/

void MlEqThreeMoment::calculateReducedMoments(CVector& moments) 
{

	double first_moment  = 0.0;
	double second_moment = 0.0;
	double third_moment  = 0.0;

	int asset_i = 0;
	for ( int bi = 0; bi < m_numberAssets; bi++ )
	{
		for ( int ai = m_firstFutureIndex; ai < m_asianWeights.getsize(); ai++ )
		{
			first_moment += m_basketWeights[bi] * m_asianWeights[ai]*m_fwds[bi][ai];
		}

		asset_i++;
	}

	int a1,a2,a3,b1,b2,b3;

	//  pre-calculating the e(b1, b2, a1, a2)

	GVector < cTensor > cache;
	cache.resize(m_numberAssets);

	for ( b1=0; b1 < m_numberAssets; b1++ )
	{

cache[b1].resize(m_numberAssets,m_asianWeights.getsize(),m_asianWeights.getsize());

		for ( b2=0; b2 < m_numberAssets; b2++ )
		{
			for ( a1=m_firstFutureIndex; a1 < m_asianWeights.getsize();  a1++ )
			{
				for ( a2=m_firstFutureIndex; a2 < m_asianWeights.getsize(); a2++ )
				{
					cache[b1][b2][a1][a2] = e(b1, b2, a1, a2);          
				}
			}
		}
	}

	// New method, it is almost twice faster than the original one
	second_moment = 0.0;
	for ( b1=0 ; b1 <  m_numberAssets; b1++ )
	{
		for ( b2=0;  b2 < m_numberAssets;  b2++ )
		{
			for ( a1=m_firstFutureIndex;  a1 < m_asianWeights.getsize(); a1++ )
			{
				for ( a2=m_firstFutureIndex; a2 < a1; a2++ )
				{
					second_moment += m_basketWeights[b1] * m_basketWeights[b2] *
						m_asianWeights[a1] * m_asianWeights[a2]*
						m_fwds[b1][a1]*m_fwds[b2][a2]
						* (cache[b1][b2][a1][a2] - 1.0);
				}
			}
		}
	}
	second_moment *= 2.0;
	for ( b1=0; b1 < m_numberAssets; b1++ )
	{
		for ( b2=0; b2 < m_numberAssets; b2++, b2++ )
		{
			for ( a1=m_firstFutureIndex; a1 < m_asianWeights.getsize(); a1++ )
			{
				a2 = a1;
				second_moment += m_basketWeights[b1]* m_basketWeights[b2] *
					m_asianWeights[a1] * m_asianWeights[a2]*
					m_fwds[b1][a1]*m_fwds[b2][a2]*
					(cache[b1][b2][a1][a1] - 1.0);
			}
		}
	}
	

	third_moment = 0;
	for ( b1=0; b1 < m_numberAssets; b1++ )
	{	
		for ( b2=0; b2 <m_numberAssets; b2++ )
		{	
			for (b3=0; b3 < m_numberAssets; b3++ )
			{	
				
				for (a1=0, a1 = m_firstFutureIndex; a1 < m_asianWeights.getsize(); a1++)
				{	
					for (a2=m_firstFutureIndex; a2 < m_asianWeights.getsize(); a2++ )
					{	
						for (a3=m_firstFutureIndex; a3 < m_asianWeights.getsize(); a3++)
						{
							third_moment += m_basketWeights[b1] * m_basketWeights[b2] * m_basketWeights[b3]
								* m_asianWeights[a1] * m_asianWeights[a2] * m_asianWeights[a3]
								* (cache[b1][b2][a1][a2] * cache[b1][b3][a1][a3] * cache[b2][b3][a2][a3]
								- cache[b1][b2][a1][a2] - cache[b1][b3][a1][a3] - cache[b2][b3][a2][a3] + 2.0)
								* m_fwds[b1][a1]*m_fwds[b2][a2]*m_fwds[b3][a3];
						}
					}
				}
			}
		}
	}

	moments.resize(3);
	moments[0] = first_moment;
	moments[1] = second_moment;
	moments[2] = third_moment;
}

/****************************************************************
**	Class  : MlEqThreeMoment 
**	Routine: computeThirdMoments
**	Returns: double
**	Action : 
****************************************************************/

double MlEqThreeMoment::computeThirdMoments() 
{
	int bi1, bi2, bi3;
	int ai1, ai2, ai3;
	double third_moments, temp1, temp2;
	third_moments = temp1 = temp2 = 0.0;

	// t1
	temp1 = 0.0;
	for ( bi1 = 0; bi1 < m_numberAssets; bi1++ )
	{	
		for ( bi2 = 0; bi2 < m_numberAssets; bi2++ )
		{	
			for ( bi3 = 0; bi3 < m_numberAssets; bi3++ )
			{	
				
				for ( ai1= m_firstFutureIndex; ai1 < m_asianWeights.getsize(); ai1++ )
				{	
					for ( ai2= m_firstFutureIndex; ai2 < ai1; ai2++ )
					{	
						for ( ai3= m_firstFutureIndex; ai3 < ai2; ai3++ )
						{
							temp1 += m_basketWeights[bi1] * m_basketWeights[bi2] * m_basketWeights[bi3]
								* m_asianWeights[ai1] * m_asianWeights[ai2] * m_asianWeights[ai3]
								* (e(bi1,bi2,ai1,ai2) * e(bi1,bi3,ai1,ai3) * e(bi2,bi3,ai2,ai3)
								- e(bi1,bi2,ai1,ai2) - e(bi1,bi3,ai1,ai3) - e(bi2,bi3,ai2,ai3) + 2.0)
								*m_fwds[bi1][ai1]* m_fwds[bi2][ai2] * m_fwds[bi3][ai3];
						}
					}
				}
			}
		}
	}
	third_moments += 6.0 * temp1;

	// t2
	temp1 = 0.0;
	for ( bi1 = 0; bi1 < m_numberAssets; bi1++ )
	{	
		for ( bi2 = 0; bi2 < m_numberAssets; bi2++ )
		{	
			for ( bi3 = 0; bi3 < m_numberAssets; bi3++ )
			{	
				
				for ( ai1= m_firstFutureIndex; ai1 < m_asianWeights.getsize(); ai1++ )
				{	
					for ( ai2= m_firstFutureIndex; ai2 < ai1; ai2++ )
					{	
						ai3 = ai2;
						temp1 += m_basketWeights[bi1] * m_basketWeights[bi2] * m_basketWeights[bi3]
							* m_asianWeights[ai1]* m_asianWeights[ai2]* m_asianWeights[ai3]
							* (e(bi1,bi2,ai1,ai2) * e(bi1,bi3,ai1,ai3) * e(bi2,bi3,ai2,ai3)
							- e(bi1,bi2,ai1,ai2) - e(bi1,bi3,ai1,ai3) - e(bi2,bi3,ai2,ai3) + 2.0)
							* m_fwds[bi1][ai1] *m_fwds[bi2][ai2] * m_fwds[bi3][ai3];
					}
				}
			}
		}
	}
	third_moments += 6.0 * temp1;
	
	// t3
	temp1 = 0.0;
	for ( bi1 = 0; bi1 < m_numberAssets; bi1++ )
	{	
		for ( bi2 = 0; bi2 < m_numberAssets; bi2++ )
		{	
			for ( bi3 = 0; bi3 < bi1; bi3++ )
			{	
				
				for ( ai2= m_firstFutureIndex; ai2 < m_asianWeights.getsize(); ai2++ )
				{	
					for ( ai3= m_firstFutureIndex; ai3 < ai2; ai3++ )
					{
						ai1 = ai2;
						temp1 += m_basketWeights[bi1] * m_basketWeights[bi2] * m_basketWeights[bi3]
							* m_asianWeights[ai1] * m_asianWeights[ai2] * m_asianWeights[ai3]
							* (e(bi1,bi2,ai1,ai2) * e(bi1,bi3,ai1,ai3) * e(bi2,bi3,ai2,ai3)
							- e(bi1,bi2,ai1,ai2) - e(bi1,bi3,ai1,ai3) - e(bi2,bi3,ai2,ai3) + 2.0)
							* m_fwds[bi1][ai1]*m_fwds[bi2][ai2]* m_fwds[bi3][ai3];
					}
				}
			}
		}
	}
	third_moments += 6.0 * temp1;

	// t4
	temp1 = 0.0;
	for ( bi1 = 0; bi1 < m_numberAssets; bi1++ )
	{	
		bi2 = bi1;
		for ( bi3 = 0; bi3 < m_numberAssets; bi3++ )
		{	
				
			for ( ai2= m_firstFutureIndex; ai2 < m_asianWeights.getsize(); ai2++ )
			{	
				ai1 = ai2;
				for ( ai3= m_firstFutureIndex; ai3 < ai2; ai3++ )
				{
					temp1 += m_basketWeights[bi1]* m_basketWeights[bi2] * m_basketWeights[bi3]
						* m_asianWeights[ai1] * m_asianWeights[ai2] * m_asianWeights[ai3]
						* (e(bi1,bi2,ai1,ai2) * e(bi1,bi3,ai1,ai3) * e(bi2,bi3,ai2,ai3)
						- e(bi1,bi2,ai1,ai2) - e(bi1,bi3,ai1,ai3) - e(bi2,bi3,ai2,ai3) + 2.0)
						* m_fwds[bi1][ai1]* m_fwds[bi2][ai2]
						* m_fwds[bi3][ai3];
				}
			}
		}
	}
	third_moments += 3.0 * temp1;

	// t5
	temp1 = 0.0;
	for ( bi1 = 0; bi1 < m_numberAssets; bi1++ )
	{	
		for ( bi2 = 0; bi2 < m_numberAssets; bi2++ )
		{	
			bi3 = bi2;
				
			for ( ai1= m_firstFutureIndex; ai1 < m_asianWeights.getsize(); ai1++ )
			{	
				for ( ai2= m_firstFutureIndex; ai2 < ai1; ai2++ )
				{	
					ai3 = ai2;
					temp1 += m_basketWeights[bi1] * m_basketWeights[bi2] * m_basketWeights[bi3]
						* m_asianWeights[ai1] * m_asianWeights[ai2] * m_asianWeights[ai3]
						* (e(bi1,bi2,ai1,ai2) * e(bi1,bi3,ai1,ai3) * e(bi2,bi3,ai2,ai3)
						- e(bi1,bi2,ai1,ai2) - e(bi1,bi3,ai1,ai3) - e(bi2,bi3,ai2,ai3) + 2.0)
						* m_fwds[bi1][ai1]* m_fwds[bi2][ai2]
						* m_fwds[bi3][ai3];
				}
			}
		}
	}
	third_moments += 3.0 * temp1;

	// t6
	temp1 = 0.0;
	for ( bi1 = 0; bi1 < m_numberAssets; bi1++ )
	{	
		for ( bi2 = 0; bi2 < bi1; bi2++ )
		{	
			for ( bi3 = 0; bi3 < bi2; bi3++ )
			{	
				
				for ( ai1= m_firstFutureIndex; ai1 < m_asianWeights.getsize(); ai1++ )
				{	
					ai2 = ai1;
					ai3 = ai1;
					temp1 += m_basketWeights[bi1]* m_basketWeights[bi2] * m_basketWeights[bi3]
						* m_asianWeights[ai1] * m_asianWeights[ai2] * m_asianWeights[ai3]
						* (e(bi1,bi2,ai1,ai2) * e(bi1,bi3,ai1,ai3) * e(bi2,bi3,ai2,ai3)
						- e(bi1,bi2,ai1,ai2) - e(bi1,bi3,ai1,ai3) - e(bi2,bi3,ai2,ai3) + 2.0)
						* m_fwds[bi1][ai1] * m_fwds[bi2][ai2]* m_fwds[bi3][ai3];
				}
			}
		}
	}
	third_moments += 6.0 * temp1;

	// t7
	temp1 = 0.0;
	for ( bi1 = 0; bi1 < m_numberAssets; bi1++ )
	{	
		for ( bi3 = 0; bi3 < bi1; bi3++ )
		{	
			bi2 = bi1;
				
			for ( ai1= m_firstFutureIndex; ai1 < m_asianWeights.getsize(); ai1++ )
			{	
				ai2 = ai1;
				ai3 = ai1;
				temp1 += m_basketWeights[bi1] * m_basketWeights[bi2] * m_basketWeights[bi3]
					* m_asianWeights[ai1]* m_asianWeights[ai2] * m_asianWeights[ai3]
					* (e(bi1,bi2,ai1,ai2) * e(bi1,bi3,ai1,ai3) * e(bi2,bi3,ai2,ai3)
					- e(bi1,bi2,ai1,ai2) - e(bi1,bi3,ai1,ai3) - e(bi2,bi3,ai2,ai3) + 2.0)
					* m_fwds[bi1][ai1] * m_fwds[bi2][ai2]
					* m_fwds[bi3][ai3];
			}
		}
	}
	third_moments += 3.0 * temp1;

	// t8
	temp1 = 0.0;
	for ( bi1 = 0; bi1 < m_numberAssets; bi1++ )
	{	
		for ( bi3 = m_numberAssets; bi3 < bi1; bi3++ )
		{	
			bi2 = bi3;
				
			for ( ai1= m_firstFutureIndex; ai1 < m_asianWeights.getsize(); ai1++ )
			{	
				ai2 = ai1;
				ai3 = ai1;
				temp1 += m_basketWeights[bi1] * m_basketWeights[bi2] * m_basketWeights[bi3]
					* m_asianWeights[ai1] * m_asianWeights[ai2] * m_asianWeights[ai3]
					* (e(bi1,bi2,ai1,ai2) * e(bi1,bi3,ai1,ai3) * e(bi2,bi3,ai2,ai3)
					- e(bi1,bi2,ai1,ai2) - e(bi1,bi3,ai1,ai3) - e(bi2,bi3,ai2,ai3) + 2.0)
					* m_fwds[bi1][ai1]* m_fwds[bi2][ai2]
					* m_fwds[bi3][ai3];
			}
		}
	}
	third_moments += 3.0 * temp1;

	// t9
	temp1 = 0.0;
	for ( bi1 = 0; bi1 < m_numberAssets; bi1++ )
	{	
        bi2 = bi1;
        bi3 = bi1;
				
		for ( ai1= m_firstFutureIndex; ai1 < m_asianWeights.getsize(); ai1++ )
		{	
			ai2 = ai1;
			ai3 = ai1;
			temp1 += m_basketWeights[bi1]* m_basketWeights[bi2] * m_basketWeights[bi3]
				* m_asianWeights[ai1]* m_asianWeights[ai2]* m_asianWeights[ai3]
				* (e(bi1,bi2,ai1,ai2) * e(bi1,bi3,ai1,ai3) * e(bi2,bi3,ai2,ai3)
				- e(bi1,bi2,ai1,ai2) - e(bi1,bi3,ai1,ai3) - e(bi2,bi3,ai2,ai3) + 2.0)
				* m_fwds[bi1][ai1]* m_fwds[bi2][ai2]
				* m_fwds[bi3][ai3];

		}
	}
	third_moments += temp1;
	return third_moments;
}

/****************************************************************
**	Class  : MlEqThreeMoment 
**	Routine: calculate
**	Returns: vector of three moments
**	Action : 
****************************************************************/

CVector MlEqThreeMoment::calculate() 
{

	CVector reducedMoments;
	CVector moments(5);
	double roots[3];

	int n = m_asianDates.size()-1;
	double effectiveMaturity = m_hDate->GetYearFraction(m_asianDates[n]);


	calculateReducedMoments(reducedMoments);
			
	MlEqMaths::math_cubic_roots( 1.0, 3.0, 0.0,
						 -4.0 - pow( reducedMoments[2], 2.0 ) / pow( reducedMoments[1], 3.0 ), roots );
		
		//	take first root for now
			
	double nu = roots[0];	
	if(nu == 1.0)
		throw(  "nu == 1.0" );
		
	double Z0 = sqrt( reducedMoments[1] / ( nu - 1.0 ) );
	double shift = reducedMoments[0] - Z0;
			

	moments[0] = m_basketForward;
	moments[1] = sqrt( log( nu ) / effectiveMaturity );
	moments[2] = m_effectiveStrikeAdjust+shift;
	moments[3] = Z0;;
	moments[4] = effectiveMaturity;


	return moments;
}	


/****************************************************************
**	Class  : MlEqThreeMoment 
**	Routine: copy
**	Returns: MlEqThreeMoment*
**	Action : 
****************************************************************/

MlEqThreeMoment *MlEqThreeMoment::copy() const
{
	return new MlEqThreeMoment( *this );
}



/****************************************************************
**	Class  : MlEqThreeMoment 
**	Routine: calculatePrice
**	Returns: CVector
**	Action : 
****************************************************************/

CVector MlEqThreeMoment::calculatePrice(vector < MlEqStrikeHandle > pStrikes )
{

	CVector moments = calculate();
	MlEqStrike strike;
	CVector results(pStrikes.size());

	double fwd,vol,strikeAdj,mat;

	int n = m_asianDates.size()-1;
	mat = m_hDate->GetYearFraction(m_asianDates[n]);
	int cp = 1;

	double discount_factor = 1.0;


	for ( int istrike = 0 ; istrike < pStrikes.size(); istrike++ )
	{
		MlEqStrike::convertStrikes(strike,*pStrikes[istrike]);

		fwd			= moments[3];
		vol			= moments[1];
		strikeAdj	= moments[2];
		mat			= moments[4];

		results[istrike] = Bs(fwd,vol,mat,strike.m_strike-strikeAdj,discount_factor,cp);
	}

	return results;
}


/****************************************************************
**	Class  :  
**	Routine: Calibrate
**	Returns: void
**	Action : calibrates asset 
**           
****************************************************************/


CVector CalculateEuropeanThreeMomentBasket(const std::vector<MlEqAssetHandle>  assets,const CVector& assetWeights,DATE maturityDate, std::vector < vector < MlEqStrikeHandle > >& pVolStrikes)
{

	int nassets = assetWeights.getsize();
	if ( nassets != assets.size() ){
		throw( "dimension of assets weights must coincide with number of assets");
	}

	MlEqVolatilityStructureHandle  pVol =  assets[0]->GetVolatilityStructure();
	MlEqConstDateHandle dateToDouble = pVol->getDateToDouble();

	if ( pVolStrikes.size() != nassets+1) {
		throw("dimension of strikevector must be number of assets +1");
	}

	for ( int i = 0 ; i < nassets+1; i++ )
	{
		if (  pVolStrikes[0].size() != pVolStrikes[i].size() ){
			throw("same number of strikes must be entered per asset");
		}
	}

	int nstrikes = pVolStrikes[0].size() ;

	vector < DATE > asianDates;
	asianDates.resize(1);

	CVector asianWeights(1);
	asianWeights[0]		= 1.0;
	asianDates[0]		= maturityDate;

	CVector res;
	CVector impliedVols(nstrikes);

	double mat = dateToDouble->GetYearFraction(maturityDate);


	double fwd = 0.0;
	for ( int i = 0 ; i < nassets; i++ ){
		fwd += assetWeights[i]*assets[i]->GetQuantoForward(dateToDouble->GetDate(),maturityDate, false);
	}

	vector < MlEqStrikeHandle >  pStrikes(nassets);

	vector < MlEqStrikeHandle > basketStrike(1);

	MlEqCorrelationMatrixHandle	pcorrelation = assets[0]->GetCorrelationMatrix();

	for ( int istrike = 0 ; istrike < pVolStrikes[0].size(); istrike++ )
	{

		MlEqThreeMoment tm;

		for ( int i = 0 ; i < nassets; i++ ){
			pStrikes[i] = pVolStrikes[i+1][istrike];
		}


		CVector aWeights;
		aWeights = assetWeights;

		tm.initialize(
					  assets,
					  aWeights,
					  asianWeights,
					  asianDates,
					  *pcorrelation,
					  dateToDouble,
					  pStrikes);



		MlEqStrike xstrike;

		basketStrike[0] = pVolStrikes[0][istrike];
		res =  tm.calculatePrice(basketStrike);
		
		MlEqStrike::convertStrikes(xstrike,*basketStrike[0]);

		impliedVols[istrike] =  MlEqBSImpliedVol(res[0],fwd,mat,xstrike.m_strike,1.0,1);

	}


	return impliedVols;
}






