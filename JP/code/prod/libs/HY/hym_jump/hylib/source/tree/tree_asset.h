#ifndef _TREE_ASSET_H_
#define _TREE_ASSET_H_
#include "tree.h"
#include <map>
#include "kmatrix.h"
#include "kprocess.h"
#include "kplatdep.h"
#include "genfunc.h"   /* BaseFunction */
#include "ddmap.h"
#include <iostream>



class AssetParam
{
public:
	KVector(int)	m_idx;
	KVector(double)	m_prob;
	KVector(double) m_probDvol;
	KVector(double)	m_probShift;		//HY3.4v
	KVector(double) m_probDvolShift;	//HY3.4v

//	double			m_stock;
//	double          m_asset;

	AssetParam():m_idx(3), m_prob(3), m_probDvol(3), m_probShift(3), m_probDvolShift(3)	
	{

	}
	~AssetParam(){}

	friend	std::ostream &operator<<(std::ostream &out, const AssetParam&param)
	{
		out<<"Coordinates:"<<"\t"<<param.m_idx<<"\t";
		out<<"Probabilities:"<<"\t"<<param.m_prob<<"\t";
		out<<"Deriv Probabilities vol:"<<"\t"<<param.m_probDvol<<"\t";

		out<<"Probabilities with volShift zero beta:"<<"\t"<<param.m_probShift<<"\t";			//HY3.4v
		out<<"Deriv Probabilities volShift:"<<"\t"<<param.m_probDvolShift<<"\t";	//HY3.4v
		
		out<<"Probabilities with volShift and beta:"<<"\t"<<param.m_probShift<<"\t";			//HY3.4v

		return out;
	}
};


class AssetTree: public	BaseTree
{
	TimeSlice<AssetParam>	m_assetParam;
	BaseFunction*	m_assetToStockMap;
	BaseFunction*	m_assetProcess;
	KRateCurve		m_volCurve;
	KRateCurve		m_volCurveShift;		//HY3.4v
	KRateCurve		m_repoCurve;
	DDMap*			m_dividentYield;
	double          m_spotStock;
	double          m_spotAsset;
	double          m_beta;
	double          m_betaShift;			//HY3.4v
	double          m_maxJump;
	double			m_vollim;					//HY4.1v
	double			m_jumpProb;					//jump version
	int				m_startLimit;				//HY4.1v, limit on the valueDate, set in limits().
	int				m_endLimit;					//HY4.1v, limit on the last time point, set in limits().
	bool			m_isOneAssetVol;

protected:
	DTimeSlice	m_ts_probNoJump;			//jump version
	DTimeSlice  m_ts_asset;
	DTimeSlice  m_ts_stock;
	DTimeSlice  m_ts_stock_tplus;
	void	memInit(void);
	void	limits (void);
	void	treeSpotVol(void);
	void	drift(void);
	void	buildLattice();
	void	initDebtDrift(void);			//HY4.2
public:
	AssetTree(void){}
	AssetTree(double spotStock, KRateCurve zc, KRateCurve volCurve, KRateCurve volCurveShift, BaseFunction* ats, 
		BaseFunction* ap,KRateCurve repoCurve, double beta, double vollim, double jumpProb, DDMap* dividentYield, bool isCVOption)
	{
		m_spotStock = spotStock;
		m_beta = beta;
		m_vollim = vollim;							//HY4.1v
		m_jumpProb = jumpProb;						//jump version
		m_betaShift = (isCVOption) ? 0.0:beta;
		m_repoCurve = repoCurve;
		m_dividentYield = dividentYield;
		set_index_zero(zc);
		set_vol_curve(volCurve);
		set_vol_curveShift(volCurveShift);		// HY3.4v
		set_assetToStockMap(ats);
		set_assetProcess(ap);
		set_spotAsset();

		m_isOneAssetVol = true;
	}
	~AssetTree(void){}

	void   set_spotAsset();
	double get_spotAsset(){return m_spotAsset;}
	double get_spotStock(){return m_spotStock;}
	double get_aver_vol(KRateCurve volCurve, KDate date1, KDate date2);		//HY3.4v
//	double get_aver_vol(KDate date1, KDate date2);

	//set stock and vol information
	void	set_assetToStockMap( BaseFunction* ats) {m_assetToStockMap=ats;}
	void	set_assetProcess( BaseFunction* ap) {m_assetProcess=ap;}
	void	set_vol_curve(const KRateCurve &vol){m_volCurve = vol;}
	void	set_vol_curveShift(const KRateCurve &volShift){m_volCurveShift = volShift;}		// HY3.4v
	void	set_isOneAssetVolFalse(){m_isOneAssetVol = false;}
	double  get_max_interval();

	BaseFunction*	get_assetToStockMap() {return m_assetToStockMap;}		//HY4.1v

	double get_fwd(KDate date);
	double get_fwd_factor(bool isEquityDrift = false);		//HY3.4v

	double get_maxJump(){return m_maxJump;}
	double get_vollim(){return m_vollim;}				//HY4.1v
	double get_jumpProb(){return m_jumpProb;}				//HY4.1v
	int get_startLimit(){return m_startLimit;}				//HY4.1v
	int get_endLimit(){return m_endLimit;}				//HY4.1v

	virtual KValarray<int>	get_curr_limit(void) const;
	virtual KValarray<int>	get_next_limit(void) const;
	virtual KValarray<int>	get_max_limit(void) const;
//	virtual void	dev(DTimeSlice &);
//	virtual void	ddv(DTimeSlice &ts1,DTimeSlice &ts2);
	
	void dev(DTimeSlice &ts1, bool isVolShift = false);						//HY3.4v
	void ddv(DTimeSlice &ts1, DTimeSlice &ts2, bool isVolShift = false);	//HY3.4v

	DTimeSlice&	get_stock_price(void);
	DTimeSlice&	get_asset_value(void);
//	DTimeSlice&	get_asset_value2(void);

	DTimeSlice&	get_probNoJump(void);

	void	print(std::ostream &out)const{out<<m_assetParam;}
};




#endif

