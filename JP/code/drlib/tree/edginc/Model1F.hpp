//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : Model1F.hpp
//
//   Description : One factor model parent class.
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : February 8, 2002
//
//----------------------------------------------------------------------------

#ifndef MODEL1F_HPP
#define MODEL1F_HPP
#include "edginc/Model.hpp"
#include "edginc/TimeLine.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolRequest.hpp"


DRLIB_BEGIN_NAMESPACE


class TREE_DLL Model1F{
public:

    virtual ~Model1F() {}
    Model1F() {useTree=false;}
    bool isTree() const {return useTree;}
    void chooseTree(bool use_tree) {useTree = use_tree;}
    
    // time line
    CTimeLine   TimePts;    
    
    // penultimate smoothing also uses this
    double	TruncationStd;

 	vector<double> PriceEnd;
    vector<double> DeltaEnd;
	vector<double> GammaEnd;

    virtual int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end) = 0;

    virtual TimeMetricConstSP GetTimeMetric() const = 0;

    class TREE_DLL Product1F{
    public:

        Model1F* model1F;

        Product1F() : model1F(0){}
        virtual ~Product1F() {};

        virtual CAssetConstSP     GetAssetRef() = 0;
        virtual YieldCurveConstSP GetDiscCurveRef() = 0;
		virtual bool GetFwdStartLV() = 0;
        virtual DateTime GetFwdStartDateLV() = 0;
    
        virtual void Init(Control*    control) { 
            // should throw an error if it's called. Don't want to make it pure virtual as it need not be implemented
        };
    
        virtual void InitProd() {};
    
        /** make price refinement - control variate */
        virtual double RefinePrice(double basePrice, double discFactor, bool useCtrlVar) = 0;

        virtual string getCcyTreatment() = 0;

        /** extra output requests */
        virtual void recordOutputRequests(Control* control, Results* results, double fairValue) = 0;

        /** returns a vol request for log-normal vol */
        virtual CVolRequestConstSP GetLNRequest() {
            throw ModelException("Model1F::Product1F::GetLNRequest", "no implementation");
            return CVolRequestConstSP();
        }; // perhaps only log-normal needs this
    };


private:
    bool useTree;

};

DRLIB_END_NAMESPACE
#endif
