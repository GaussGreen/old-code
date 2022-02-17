//TreeEq.cpp
#include "tree_asset.h"
#include <math.h>
#include "RootFindBrent.h"
#include "macros.h"
#include "gtobf.h"

const Y_JUMPCOEFF = 3.0;
const double JUMP_HORIZON = 0.25;		// JUMP_HORIZON>0

extern std::ofstream out;
//extern std::ofstream out2("c:\\yy2.dat");
//extern std::ofstream out3("c:\\yy3.dat");

KValarray<int>	AssetTree::get_curr_limit(void)const
{
	KValarray<int> lim(1);
	KTimePoint	&tp = get_curr_tp();
	lim[0] = tp[__d1TreeLimit];
	return lim;
}
KValarray<int>	AssetTree::get_next_limit(void)const
{
	KValarray<int> lim(1);
	KTimePoint	&tp = get_next_tp();
	lim[0] = tp[__d1TreeLimit];
	return lim;
}

KValarray<int>	AssetTree::get_max_limit(void)const
{
	KValarray<int> lim(1);
	KTimePoint	&tp = get_curr_tp();
	lim[0] = tp[__d1TreeMaxLimit];
	return lim;
}


void	AssetTree::memInit(void)
{
	m_ts_probNoJump.resize(this);
	m_ts_asset.resize(this);
	m_assetParam.resize(this);
	m_ts_stock.resize(this);
	m_ts_stock_tplus.resize(this);
	m_ts_stock_tplus = get_stock_price(); 
}

void	AssetTree::set_spotAsset(void)
{

	double spotAsset = m_assetToStockMap->inverse(m_spotStock);

/*	SoverFunc soverFunction(m_assetToStockMap,m_spotStock);

	double spotAsset = RootFindBrent(soverFunction,
									 m_spotStock,
									 ROOTLO,
									 ROOTHI); */
	m_spotAsset = spotAsset;

//	m_spotAsset = m_spotAsset*exp(0.3*0.3);
}


void	AssetTree::drift(void)
{
	double	variance = 0;
	KRateCurve	&zr = get_index_zero();

	for(KTimeLine::iterator iter = tl_begin(); iter!=  tl_end(); iter++)
	{
		iter->second[__d1TreeMidNode] = m_spotAsset;
//	variance += iter->second[__1dTreeSpotVol] * iter->second[__1dTreeSpotVol] 
//				*iter->second[__timeTillNextTp];
	}  
}


double   AssetTree::get_aver_vol(KRateCurve volCurve, KDate date1, KDate date2)
{
//	double vol1 = m_volCurve.get_rate(GTO_LINEAR_INTERP,date1);
//	double vol2 = m_volCurve.get_rate(GTO_LINEAR_INTERP,date2);

	double vol1 = volCurve.get_rate(GTO_LINEAR_INTERP,date1);		//HY3.4v
	double vol2 = volCurve.get_rate(GTO_LINEAR_INTERP,date2);		//HY3.4v

	KDate valueDate = get_valueDate();
	double year1 = GetYearDiff(valueDate, date1);
	double year2 = GetYearDiff(valueDate, date2);
	double var = (vol2*vol2*year2-vol1*vol1*year1)/(year2-year1);
	if(var <0.0)
	{
		throw KException(" var < 0. at date: GtoFormatDate(date1).\n");
	}
	return sqrt(var);
}


// Spot Vol interpolation
void	AssetTree::treeSpotVol(void) 
{
	KTimeLine::iterator iter, nIter;

	double minJump = sqrt(Y_JUMPCOEFF*100);		//HY3.4v
	double maxJump = 0.0;
	double currJump;
//	int i=1;
//	for(iter = tl_begin(); iter != tl_end(); iter++)					
	for(iter = tl_valueDate(); iter != tl_end(); iter++)		//HY4.1v					
	{
		nIter = iter;
		nIter++;
		if(nIter == tl_end())
			break;

		KValarray<double>	spotVol(1);

//		spotVol = get_aver_vol(nIter->first, iter->first);
		spotVol = get_aver_vol(m_volCurve, nIter->first, iter->first);		//HY3.4v
		iter->second[__d1TreeSpotVol] = spotVol[0];
		currJump = spotVol[0]*sqrt(Y_JUMPCOEFF*iter->second[__timeTillNextTp]);
		maxJump = Max(maxJump,currJump);
		minJump = Min(minJump,currJump);		//HY3.4v


		//HY3.4v
		if(m_isOneAssetVol == false)
		{
			spotVol = get_aver_vol(m_volCurveShift, nIter->first, iter->first);		
			iter->second[__d1TreeSpotVolShift] = spotVol[0];
			currJump = spotVol[0]*sqrt(Y_JUMPCOEFF*iter->second[__timeTillNextTp]);
			maxJump = Max(maxJump,currJump);
			minJump = Min(minJump,currJump);
		}
		// end HY3.4v

//		out<< "{i="<<i<<", Date="<<GtoFormatDate(iter->first)<<", vol="<<spotVol[0]<<", Dt="
//			<<iter->second[__timeTillNextTp]<<", Jump="<<currJump<<", maxJump="<<maxJump<<"}"<<std::endl;
//		i++;
	}
	/*the last spot vol is irrelevant*/
	nIter = tl_end();
	nIter--;
	iter = nIter;
	iter--;
	nIter->second[__d1TreeSpotVol] = iter->second[__d1TreeSpotVol];
	nIter->second[__d1TreeSpotVolShift] = iter->second[__d1TreeSpotVolShift];	//HY3.4v

	if(minJump*sqrt(3.0) < maxJump*sqrt(2.0))					//HY3.4v
	{
		m_maxJump = maxJump * sqrt(2.0/Y_JUMPCOEFF);			//HY3.4v
	}
	else													//HY3.4v
	{
		m_maxJump = maxJump;
	}

	//HY4.1v, adjust jump size to have nodes on the barrier
	double spotAsset = get_spotAsset();						//in unit of debtPerShare
	double level;
	if(AssetToStockMapping1* pt = dynamic_cast<AssetToStockMapping1*>(m_assetToStockMap))
	{
		level = (*pt).get_lim();
	}
	double temp = log(spotAsset/level);
	if(temp > m_maxJump)
	{
		m_maxJump = temp/((double)floor(temp/m_maxJump));
//		m_maxJump = temp/((double)floor(temp/m_maxJump + 0.5));
	}
	//end HY4.1v
}

void	AssetTree::limits (void)
{
	double	vol;

	KTimeLine::iterator	iter, nIter;

	iter = tl_begin();
	/* limits for the first node*/
	iter->second[__d1TreeLimit] = 1;
	iter->second[__d1TreeMaxLimit] = 1;
	vol = 0;

	while(iter != tl_end())
	{   
		nIter = iter;
		nIter++;

		if(nIter == tl_valueDate())		//HY4.1v
		{
			m_startLimit = iter->second[__d1TreeLimit] + 1;		//HY4.1v
		}

		if(nIter == tl_end())
			break;

		nIter->second[__d1TreeLimit] = iter->second[__d1TreeLimit] + 1 ;
		nIter->second[__d1TreeMaxLimit] = iter->second[__d1TreeLimit] + 1;
		iter++;
	}
	m_endLimit = iter->second[__d1TreeMaxLimit];				//HY4.1v
}



double	AssetTree::get_max_interval (void)
{
	KTimeLine::iterator	iter;

	iter = tl_begin();
	double maxInterval=0;
	
	while(iter != tl_end())
	{   
		maxInterval = MAX(iter->second[__timeTillNextTp],maxInterval);
		iter++;
	}
	return maxInterval;
}


DTimeSlice&	AssetTree::get_asset_value(void) 
{
//	DTimeSlice* asset = new DTimeSlice(this);
	int	dim = get_curr_limit()[0];
	int	i;
	KTimePoint	&tp = get_curr_tp();
//	double     timeTillNextTp =  tp[__timeTillNextTp];

//	out<<"1:"<<timeTillNextTp<<std::endl;

//	     timeTillNextTp =  get_max_interval();
//	out<<"2:"<<timeTillNextTp<<std::endl;
//	double vol = tp[__d1TreeSpotVol];
	double     jump = get_maxJump();

	m_assetProcess->set_amplitude(m_spotAsset);

	for(i=-dim; i<= dim; i++)
	{
		 double temp = (*m_assetProcess)(i*jump);
		  
		m_ts_asset[i]=temp;
	}

//	out<<dim<<"asset\n"<<m_ts_asset<<std::endl;
//	std::ofstream out("c:\\yy1.dat");
	return m_ts_asset;
}


DTimeSlice&	AssetTree::get_stock_price(void) 
{
//	DTimeSlice* stock = new DTimeSlice(this);
	DTimeSlice ts_asset = get_asset_value();
	int	dim = get_curr_limit()[0];
	int	i;

	for(i=-dim; i<= dim; i++)
	{
		double temp1 = ts_asset[i];
		double temp =(*m_assetToStockMap)(ts_asset[i]);
		m_ts_stock[i] =temp;
	}

//	out<<dim<<"\nasset\n"<<ts_asset<<"\nstock\n"<<m_ts_stock<<std::endl;

	return m_ts_stock;
}


DTimeSlice&	AssetTree::get_probNoJump(void) 
{
	DTimeSlice	ts_asset = get_asset_value();
	double		vollim = get_vollim();
	double		jumpProb = get_jumpProb();
	double		d2;
	int			dim = get_curr_limit()[0];
	int			i;
	bool		isLognormalJump = true;			//always use LN distibution
	double		dist, temp1, temp2, p_nd;		//E2C jump variables

	KTimePoint	&tp = get_curr_tp();
	double     timeTillNextTp =  tp[__timeTillNextTp];

	AssetToStockMapping1* ats = dynamic_cast<AssetToStockMapping1*>(get_assetToStockMap());
	double level = ats->get_lim();

	for(i=-dim; i<= dim; i++)
	{
		if(ts_asset[i] <= level)
		{
			m_ts_probNoJump[i] = 0.0;		
		}
		else
		{
			if (vollim == 0.0)
			{
				m_ts_probNoJump[i] = exp(-jumpProb * timeTillNextTp);	
			}
			else
			{
				if(isLognormalJump)
				{
					d2 = log(level/ts_asset[i])/vollim - vollim/2.0;
					m_ts_probNoJump[i] = exp((log(1.0 - GtoNormalCum(d2))/JUMP_HORIZON - jumpProb) * timeTillNextTp);
				}
				else
				{
					dist = ts_asset[i]/level * exp(vollim * vollim);
					temp1 = -vollim/2.0 + log(dist)/vollim;
					temp2 = -vollim/2.0 - log(dist)/vollim;
					p_nd = GtoNormalCum(temp1) - dist * GtoNormalCum(temp2);

					m_ts_probNoJump[i] = exp((log(p_nd)/JUMP_HORIZON - jumpProb) * timeTillNextTp);
				}
			}
		}
	}
	
	return m_ts_probNoJump;
}


// not used now due to change in definition of divident yield
double AssetTree::get_fwd(KDate date)
{
	KDate valueDate = get_valueDate();

	KRateCurve irCurve = get_index_zero();
	double irRate = irCurve.get_rate(date);
	double repoRate = m_repoCurve.get_rate(date);
	double averYield = m_dividentYield->get_average_amount(valueDate, date);
	double effRate = (1.0+irRate)/(1.0+averYield)-1.0;
	DDMap divident = m_dividentYield->get_subMap(valueDate, date);
	double years = GetYearDiff(valueDate, date);
	double adjustFactor = exp(-m_beta*years);
	double fwd = 1.0;
	double growthFactor = pow((1.0+effRate*adjustFactor)*(1.0+averYield)/(1.0+repoRate),years);


	for(DDMap::iterator iter = divident.begin(); iter != divident.end(); iter++)
	{
		double divident =iter->second;
		double disc = irCurve.get_pv(iter->first);
		fwd -= iter->second*irCurve.get_pv(iter->first);
	}

	if(fwd <0.0)
	{
		throw KException("fwd cannot be negative.\n");
	}

	fwd *= growthFactor;

	return fwd;
}


void	AssetTree::initDebtDrift (void)				//HY4.2
{
	double	divident, fwdRateFactor, fwdFactor, repoRate, dt;
//	double  meanRevert;								//commented out, jump version
	KDate	currentDate, nextDate;
	KDate	valueDate = get_valueDate();
//	double	jumpProb = get_vollim();

	KTimeLine::iterator	iter, nIter;

	iter = tl_valueDate();
	iter->second[__cumDebtDrift] = 1.0;

	while(iter != tl_end())
	{   
		nIter = iter;
		nIter++;

		if(nIter == tl_end())
			break;

		currentDate = iter->first;
		nextDate = nIter->first;

		repoRate = m_repoCurve.get_rate(currentDate);
		divident = m_dividentYield->get_amount(currentDate, nextDate);
		fwdRateFactor = get_index_zero().get_pv(currentDate)/get_index_zero().get_pv(nextDate);
		dt =(nextDate-currentDate)/360.0;		//Act/360
//		meanRevert = exp(-m_beta*(currentDate-valueDate)/365.0);

//		fwdFactor = 1.0+ (fwdRateFactor - 1.0 - divident) * (1.0 - meanRevert);
		fwdFactor = 1.0 - repoRate*dt + (fwdRateFactor - 1.0 - divident);												//jump version
//		fwdFactor = 1.0+ (fwdRateFactor - 1.0 - divident + jumpProb * (nextDate-currentDate)/365.0);	//jump version
		nIter->second[__cumDebtDrift] = iter->second[__cumDebtDrift] * fwdFactor;

		iter++;
	}
}


//double AssetTree::get_fwd_factor()
double AssetTree::get_fwd_factor(bool isEquityDrift)			//HY3.4v, isEquityDrift = false by default
{
	KDate  valueDate = get_valueDate();
	KDate  currentDate = get_curr_date();
	KDate  nextDate = get_next_date();

	double fwdFactor;

	double repoRate = m_repoCurve.get_rate(currentDate);
	double divident = m_dividentYield->get_amount(currentDate, nextDate);
//	out<<"Dt: "<<GtoFormatDate(currentDate)<<"Divident: "<<divident<<std::endl;
	double dt =(nextDate-currentDate)/360.0;

	double fwdRateFactor = get_index_zero().get_pv(currentDate)/get_index_zero().get_pv(nextDate);
	double meanRevert = exp(-m_beta*(currentDate-valueDate)/365.0);
	double meanRevertShift = exp(-m_betaShift*(currentDate-valueDate)/365.0);

	if(isEquityDrift)		//HY3.4v
	{
		fwdFactor = 1.0-repoRate*dt+ (fwdRateFactor - 1.0 - divident)*meanRevertShift;		//HY3.4v
	}
	else					//HY3.4v
	{
		fwdFactor = 1.0-repoRate*dt+ (fwdRateFactor - 1.0 - divident) * meanRevert;
	}

//	double  temp = pow(1.0+m_dividentYield, (nextDate-currentDate)/365.0);
//	return fwdPvofRepo/pow(fwdPvofRate*pow(1+m_dividentYield, (nextDate-currentDate)/365.0),meanRevert);
	return fwdFactor;
}

/*
void	AssetTree::buildLattice (void)			//as of HY4.2d
{
	int		dim, i;
	KTimePoint	&tp = get_curr_tp();
	double	   vol = tp[__d1TreeSpotVol];
	double	   volShift = tp[__d1TreeSpotVolShift];		//HY3.4v

	double     timeTillNextTp =  tp[__timeTillNextTp];
	double     maxJump =  get_maxJump();
	double     probCenter = 1.0 - vol*vol*timeTillNextTp/maxJump/maxJump;
	double     probCenterShift = 1.0 - volShift*volShift*timeTillNextTp/maxJump/maxJump;		//HY3.4v

//	out<<"Dt: "<<timeTillNextTp<<"  vol: "<<vol<<"  maxJump: "<<maxJump;

	dim = tp[__d1TreeLimit];
	DTimeSlice  stock = get_stock_price();
//	out<<"\nStockB\n"<<stock<<std::endl;
	
//	int dim2 = get_curr_limit()[0];

//	out<<" dimFwd: "<<dim<<" disc: "<<1/fwdFactor<<"\n"<<std::endl;

	double   fwdFactor = get_fwd_factor();
	double   fwdFactorShift = fwdFactor;					//HY4.2, no beta shift
//	double   fwdFactorShift = get_fwd_factor(true);				//HY3.4v

	for (i = -dim; i <= dim; i++)                   	
	{ 	
		AssetParam	&param = m_assetParam[i];
		param.m_idx[0] = i-1;
		param.m_idx[1] = i;
		param.m_idx[2] = i+1;

		double s0p=stock[i];
		double s0=m_ts_stock_tplus[i];
		double su=m_ts_stock_tplus[i+1];
		double sd=m_ts_stock_tplus[i-1];
		if( IS_ALMOST_ZERO(s0p ))
		{
			param.m_prob[1] = 1.0;   //center
			param.m_prob[2] = 0.0;   // up
			param.m_prob[0] = 0.0;   // down

			param.m_probDvol[1] = 0.0;
			param.m_probDvol[2] = 0.0;
			param.m_probDvol[0] = 0.0;

			//HY3.4v
			param.m_probShift[1] = 1.0;   //center
			param.m_probShift[2] = 0.0;   // up
			param.m_probShift[0] = 0.0;   // down

			param.m_probDvolShift[1] = 0.0;
			param.m_probDvolShift[2] = 0.0;
			param.m_probDvolShift[0] = 0.0;

			//end HY3.4v
		}
		else
		{
			param.m_prob[1] = probCenter;
			param.m_prob[2] = (s0p*fwdFactor-sd)/(su-sd)-(s0-sd)/(su-sd)*param.m_prob[1];
			param.m_prob[0] = (su-s0p*fwdFactor)/(su-sd)+(s0-su)/(su-sd)*param.m_prob[1];

			param.m_probDvol[1] = 2.0*(probCenter-1.0)/(vol*100.0);			//per 1% vol
			param.m_probDvol[2] = -(s0-sd)/(su-sd)*param.m_probDvol[1];
			param.m_probDvol[0] =  (s0-su)/(su-sd)*param.m_probDvol[1];

			//HY3.4v
			//prob with vol shift
			param.m_probShift[1] = probCenterShift;
			param.m_probShift[2] = (s0p*fwdFactorShift-sd)/(su-sd)-(s0-sd)/(su-sd)*param.m_probShift[1];
			param.m_probShift[0] = (su-s0p*fwdFactorShift)/(su-sd)+(s0-su)/(su-sd)*param.m_probShift[1];

			param.m_probDvolShift[1] = 2.0*(probCenterShift-1.0)/(vol*100.0);			//per 1% vol
			param.m_probDvolShift[2] = -(s0-sd)/(su-sd)*param.m_probDvolShift[1];
			param.m_probDvolShift[0] =  (s0-su)/(su-sd)*param.m_probDvolShift[1];

			//end HY3.4v
		}
//			out<<"Prob:"<<param<<std::endl;
	}  
	//save for the next run
	m_ts_stock_tplus = stock;
}
*/


void	AssetTree::buildLattice (void)				//with jump
{
	int		dim, i;
	double	s0p, fwdFactor;
	double	volShift,volShiftFactor, dVolShiftFactorDVolShift;

	KTimePoint	&tp = get_curr_tp();
	double	   vol = tp[__d1TreeSpotVol];

	double     timeTillNextTp =  tp[__timeTillNextTp];
	double     maxJump =  get_maxJump();
//	double		vollim = get_vollim();

	double		jumpSizeFactor = exp(maxJump);
	double		volFactor = exp(vol * vol * timeTillNextTp);
	double		dVolFactorDVol = 2.0 * volFactor * vol * timeTillNextTp;

	if(m_isOneAssetVol == false)
	{
		volShift = tp[__d1TreeSpotVolShift];		//HY3.4v
		volShiftFactor = exp(volShift * volShift * timeTillNextTp);
		dVolShiftFactorDVolShift = 2.0 * volShiftFactor * volShift * timeTillNextTp;
	}
	
	double		temp = (exp(maxJump/2.0) - exp(-maxJump/2.0));

	dim = tp[__d1TreeLimit];
	
	DTimeSlice  asset = get_asset_value();
	DTimeSlice  stock = get_stock_price();
	DTimeSlice  probNoJump = get_probNoJump();
	
	AssetToStockMapping1* ats = dynamic_cast<AssetToStockMapping1*>(get_assetToStockMap());
	double level = ats->get_lim();

	for (i = -dim; i <= dim; i++)                   	
	{ 	
		AssetParam	&param = m_assetParam[i];
		param.m_idx[0] = i-1;
		param.m_idx[1] = i;
		param.m_idx[2] = i+1;

		s0p = stock[i];

		if( IS_ALMOST_ZERO(s0p ))
		{
			param.m_prob[1] = 1.0;   //center
			param.m_prob[2] = 0.0;   // up
			param.m_prob[0] = 0.0;   // down

			param.m_probDvol[1] = 0.0;
			param.m_probDvol[2] = 0.0;
			param.m_probDvol[0] = 0.0;

			//HY3.4v
			param.m_probShift[1] = 1.0;   //center
			param.m_probShift[2] = 0.0;   // up
			param.m_probShift[0] = 0.0;   // down

			param.m_probDvolShift[1] = 0.0;
			param.m_probDvolShift[2] = 0.0;
			param.m_probDvolShift[0] = 0.0;

			//end HY3.4v
		}
		else
		{
			fwdFactor = level/asset[i] + (1.0 - level/asset[i])/probNoJump[i];

//			param.m_prob[1] = 1.0 - ((volFactor*fwdFactor*fwdFactor - 1.0)
//									  -(fwdFactor - 1.0)*(jumpSizeFactor + 1.0/jumpSizeFactor))/temp/temp;

			param.m_prob[2] = ((volFactor*fwdFactor*fwdFactor - 1.0)*(1.0 - 1.0/jumpSizeFactor)-(fwdFactor - 1.0)
				         *(1.0 - 1.0/jumpSizeFactor/jumpSizeFactor))/(jumpSizeFactor - 1.0/jumpSizeFactor)/temp/temp;
			param.m_prob[0] = -1.0*((volFactor*fwdFactor*fwdFactor - 1.0)*(1.0 - jumpSizeFactor)-(fwdFactor - 1.0)
				         *(1.0 - jumpSizeFactor*jumpSizeFactor))/(jumpSizeFactor - 1.0/jumpSizeFactor)/temp/temp;
			param.m_prob[1] = 1.0 - param.m_prob[0] - param.m_prob[2];

			param.m_probDvol[1] = -dVolFactorDVol * fwdFactor * fwdFactor/temp/temp/100.0;	//per 1% vol
			param.m_probDvol[2] = -param.m_probDvol[1] * (1.0 - 1.0/jumpSizeFactor)/(jumpSizeFactor - 1.0/jumpSizeFactor);
			param.m_probDvol[0] = param.m_probDvol[1] * (1.0 - jumpSizeFactor)/(jumpSizeFactor - 1.0/jumpSizeFactor);

			//prob with vol shift
			if(m_isOneAssetVol == false)
			{
				param.m_probShift[2] = ((volShiftFactor*fwdFactor*fwdFactor - 1.0)*(1.0 - 1.0/jumpSizeFactor)-(fwdFactor - 1.0)
							 *(1.0 - 1.0/jumpSizeFactor/jumpSizeFactor))/(jumpSizeFactor - 1.0/jumpSizeFactor)/temp/temp;
				param.m_probShift[0] = -1.0*((volShiftFactor*fwdFactor*fwdFactor - 1.0)*(1.0 - jumpSizeFactor)-(fwdFactor - 1.0)
							 *(1.0 - jumpSizeFactor*jumpSizeFactor))/(jumpSizeFactor - 1.0/jumpSizeFactor)/temp/temp;
				param.m_probShift[1] = 1.0 - param.m_probShift[0] - param.m_probShift[2];

				param.m_probDvolShift[1] = -dVolShiftFactorDVolShift * fwdFactor * fwdFactor/temp/temp/100.0;	//per 1% vol
				param.m_probDvolShift[2] = -param.m_probDvolShift[1] * (1.0 - 1.0/jumpSizeFactor)/(jumpSizeFactor - 1.0/jumpSizeFactor);
				param.m_probDvolShift[0] = param.m_probDvolShift[1] * (1.0 - jumpSizeFactor)/(jumpSizeFactor - 1.0/jumpSizeFactor);
			}
			else
			{
				param.m_probShift[2] = param.m_prob[2];
				param.m_probShift[0] = param.m_prob[0];
				param.m_probShift[1] = param.m_prob[1];

				param.m_probDvolShift[1] = param.m_probDvol[1];
				param.m_probDvolShift[2] = param.m_probDvol[2];
				param.m_probDvolShift[0] = param.m_probDvol[0];
			}
		}
//			out<<"Prob:"<<param<<std::endl;
	}  
}


//void	AssetTree::dev(DTimeSlice &ts)	
void	AssetTree::dev(DTimeSlice &ts, bool isVolShift)	//HY3.4v
{
	int		i, iu, i0, id;
	double	disc;
	KTimePoint	&tp = get_curr_tp();

	int	dim  = get_curr_limit()[0];
	disc = 1./(1 + tp[__fwdRateTillNextTp]);

//	out<<" dim: "<<dim<<" disc: "<<disc<<"\n"<<std::endl;
//	out<<" \ndim: "<<dim<<"\n"<<std::endl;
	for (i = -dim; i <= dim; i++)
	{
		AssetParam	&param = m_assetParam[i];
		iu = param.m_idx[2];
		i0 = param.m_idx[1];
		id = param.m_idx[0];

		if(!isVolShift)														//HY3.4v
		{
			double temp =  (param.m_prob[2] * ts[iu]+ 
							param.m_prob[1] * ts[i0]+ 
							param.m_prob[0] * ts[id] ) * disc;
			_temp[i]=temp;
		}
		else														//HY3.4v
		{
			double temp =  (param.m_probShift[2] * ts[iu]+ 
							param.m_probShift[1] * ts[i0]+ 
							param.m_probShift[0] * ts[id] ) * disc;		//HY3.4v
			_temp[i]=temp;
		}

//		  out<<"A["<<i<<", "<<_temp.size()<<"]\t";
//		  out<<param<<std::endl;
//		  out<<temp1<<"  "<<temp2<<"  "<<temp3<<std::endl;
	}

// out<<"tempdev\n"<<_temp<<std::endl;

	for(i=-dim; i<= dim; i++)
		ts[i] = _temp[i];
//	out<<"A["<<i<<", "<<ts.size()<<"]\t";
	
}



//void	AssetTree::ddv(DTimeSlice &ts1,DTimeSlice &ts2)	
void	AssetTree::ddv(DTimeSlice &ts1, DTimeSlice &ts2, bool isVolShift)	//HY3.4v	
{
	int		i, iu, i0, id;
	double	disc;
	KTimePoint	&tp = get_curr_tp();

	int	dim  = get_curr_limit()[0];
	disc = 1./(1 + tp[__fwdRateTillNextTp]);

//	out<<" dim: "<<dim<<" disc: "<<disc<<"\n"<<std::endl;
//	out<<" \ndim: "<<dim<<"\n"<<std::endl;
	for (i = -dim; i <= dim; i++)
	{
		AssetParam	&param = m_assetParam[i];
		iu = param.m_idx[2];
		i0 = param.m_idx[1];
		id = param.m_idx[0];

		if(!isVolShift)												//HY3.4v
		{
			 double temp =  ( param.m_prob[2] * ts1[iu]+ 
							  param.m_prob[1] * ts1[i0]+ 
							  param.m_prob[0] * ts1[id] +
							  param.m_probDvol[2] * ts2[iu]+ 
							  param.m_probDvol[1] * ts2[i0]+ 
							  param.m_probDvol[0] * ts2[id]) * disc;
			_temp[i]=temp;
		}
		else														//HY3.4v
		{
			 double temp =  ( param.m_probShift[2] * ts1[iu]+ 
							  param.m_probShift[1] * ts1[i0]+ 
							  param.m_probShift[0] * ts1[id] +
							  param.m_probDvolShift[2] * ts2[iu]+ 
							  param.m_probDvolShift[1] * ts2[i0]+ 
							  param.m_probDvolShift[0] * ts2[id]) * disc;		//HY3.4v
			_temp[i]=temp;
		}

//		  out<<"A["<<i<<", "<<_temp.size()<<"]\t";
//		  out<<param<<std::endl;
//		  out<<temp1<<"  "<<temp2<<"  "<<temp3<<std::endl;
	}

// out<<"tempdev\n"<<_temp<<std::endl;

	for(i=-dim; i<= dim; i++)
		ts1[i] = _temp[i];
//	out<<"A["<<i<<", "<<ts.size()<<"]\t";
	
}
