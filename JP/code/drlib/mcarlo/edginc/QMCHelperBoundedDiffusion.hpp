//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : QMCHelperBoundedDiffusion
//
//   Description : concrete per-asset implementations of IQMCHelperBoundedDiffusion
//
//   Date        : Thu Apr 20 2006
//
//   Author      : Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>
//
//
//----------------------------------------------------------------------------

#ifndef QMCHelperBoundedDiffusion_HPP
#define QMCHelperBoundedDiffusion_HPP

#include "edginc/DateTime.hpp"

#include "edginc/QMCHelperVectorDiffusionDate.hpp"
#include "edginc/IQMCHelperBoundedDiffusion.hpp"

DRLIB_BEGIN_NAMESPACE

/** Generic implementation of IQMCHelperBoundedDiffusion with two dates to track (maxDiffusiondate and maxCurveMaturity) */
class QMCHelperBoundedDiffusion : public virtual IQMCHelperBoundedDiffusion {
    QMCHelperVectorDiffusionDate   maxDiffDate;
    QMCHelperVectorDiffusionDate   maxCurveMat;

public:
    /** Assets and their dependencies:
        IR -> None or FX
        FX -> domIR and forIR
        CR -> IR
        EQ -> IR
        FIXME: Energy, basis, ...
    */
    QMCHelperBoundedDiffusion(QMCHelperBoundedDiffusion * d1 = NULL, QMCHelperBoundedDiffusion * d2 = NULL) {
        add(d1);
        add(d2);
    }

    // Add dependent asset (should be public?)
    void add( QMCHelperBoundedDiffusion * d) {
        if (d != NULL) {
            maxDiffDate.add(& d->maxDiffDate);
            maxCurveMat.add(& d->maxCurveMat);
        }
    }

    /** Returns pointer to internal storage for MaxDiff(usion)Date or NULL */
    virtual DateTime getMaxDiffDate() const { return maxDiffDate.currentMaturity();}

    /** Returns pointer to internal storage for MaxCurveMaturity calculation or NULL */
    virtual DateTime getMaxCurveMat() const { return maxCurveMat.currentMaturity();}

    virtual void updateMaxDiffDate(const DateTime& date) { maxDiffDate.updateMaturity(date);}
    virtual void updateCurveMat   (const DateTime& date) { maxCurveMat.updateMaturity(date);}
    /// true iff updateMaxDiffDate was never called, i.e. diffusion of this asset is not needed
    virtual bool    empty() const {return maxDiffDate.empty();}
};

// SRMRates sets FX outside of constructor
class QMCIRBoundedDiffusion : public QMCHelperBoundedDiffusion
{
    public:

        QMCIRBoundedDiffusion(QMCHelperBoundedDiffusion *fx) :
            QMCHelperBoundedDiffusion(fx)
        {}

        void setFX(QMCHelperBoundedDiffusion *fx)
        {
            if (fx)
                add(fx);
        }
};

typedef QMCHelperBoundedDiffusion QMCCRBoundedDiffusion;
typedef QMCHelperBoundedDiffusion QMCEQBoundedDiffusion;
typedef QMCHelperBoundedDiffusion QMCFXBoundedDiffusion;
typedef QMCHelperBoundedDiffusion QMCEnergyBoundedDiffusion;

DRLIB_END_NAMESPACE

#endif
