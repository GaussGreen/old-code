	
#ifndef THREEMOMENTSH
#define THREEMOMENTSH
	

#include "mleqdate.h"
	
	
class  MlEqThreeMoment
{		
			
	
protected:	

//	MlEqAssetHandle	m_pBasket;
	int m_numberAssets;

	CVector m_asianWeights;//[idate]
	std::vector < DATE >  m_asianDates;
	int m_firstFutureIndex;

	double m_effectiveStrikeAdjust;
	double m_basketForward;

	MlEqConstDateHandle m_hDate;
	
	double e(int iasset, int jsset,
			 int iDate, int jDate);
	

	void calculateReducedMoments(CVector& moments);

	double computeThirdMoments() ;


	cTensor m_covariance;//[iasset][iasset][idate]
	CMatrix m_fwds;//[iasset][idate]
	
	CVector m_basketWeights;//[iasset]

public:

	void initialize(
									 MlEqAsset& asset,
									 CVector& asianWeights,
									 std::vector < DATE >& asianDates,
									 MlEqCorrelationMatrix& correlation,
									 MlEqConstDateHandle hDate,
									 std::vector < MlEqStrikeHandle >& pStrike);



	void initialize(const std::vector<MlEqAssetHandle>  assets,
								 CVector& basketWeights,
								 CVector& asianWeights,
								 std::vector < DATE >& asianDates,
								 MlEqCorrelationMatrix& correlation,
								 MlEqConstDateHandle hDate,
								 std::vector < MlEqStrikeHandle >& pStrike);

	void initialize(
								 cTensor& covariance,
								 CMatrix& forwards,
								 CVector& asianWeights,
								 CVector& basketWeights,
								 std::vector < DATE >& asianDates,
								 MlEqConstDateHandle hDate,
								 std::vector < MlEqStrikeHandle >& pStrike);



//	MlEqThreeMoment( const MlEqThreeMoment& rhs);
	
	MlEqThreeMoment& operator=( const MlEqThreeMoment& rhs );
	MlEqThreeMoment *copy() const;	
	MlEqThreeMoment(){};

	virtual ~MlEqThreeMoment();

	CVector calculatePrice(std::vector < MlEqStrikeHandle > pStrikes );
	CVector calculate();

};

CVector CalculateEuropeanThreeMomentBasket(const std::vector<MlEqAssetHandle>  assets,const CVector& assetWeights,DATE maturityDate, std::vector < std::vector < MlEqStrikeHandle > >& pVolStrikes);

#endif 

