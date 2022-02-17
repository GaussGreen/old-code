/**
 * @file CCMDBetaDSpread.cpp
 *
 *
 */

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/ParSpreadPropShift.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/CCMBetaSens.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Sensitivity of CCM beta-sensitivity to spreads
 *
 * These numbers (previously obtained from Kapital) are used by Credit Hybrids
 * in connection with CCM models parameterised by some numbers �_i for each
 * name i in the basket of names.  The aim is to work out the uncertainty in
 * the modelled price arising from
 *
 *    -  uncertainty about the �_i
 *    -  uncertainty about how sensitive the price is to the �_i
 *
 * with a view to putting away reserves proportional to stddev(resulting price
 * uncertainty).  
 *
 * This class computes CCM_D_BETA_D_PAR_SPREAD_RHO_PROP (for each name in the
 * basket).  CCM_D_BETA_BASKET_D_PAR_SPREAD_RHO_PROP is computed by
 * CCMDBasketBetaDSpread.cpp.
 *
 *
 * <H3>My understanding of the model</H3>
 *
 * <BLOCKQUOTE>
 * Notation: capital letters mean denote random variables, and lower case
 * letters denote constants (except � is also a random variable).  x* denotes
 * our best estimates for the RV X.  Var[...] means variance.  X ~ 0 � v means
 * "X has mean zero and variance v".
 * </BLOCKQUOTE>
 *
 * Specifically, we model the error between our best-guess price v* and the
 * true price V as linear in the error in the �_i ...
 * 
 * -     v* = V + A_i (�_i - �*_i)
 * 
 * ... where the error in the �_i is composed of a component M common to all
 * names, plus an independent component for each name E_i:
 * 
 * -     �*_i = �_i + M + E_i   where M ~ 0 � �, and each E_i ~ 0 � 1-�
 * 
 * (It is easy to show that the correlations between �_i and �_j is then � for
 * all i != j.  Note that the overall variance of the � uncertainty is assumed
 * to be 1 (!), but scaling this variance simply scales the final answers
 * proportionately.)
 * 
 * Our A_i's will be estimated by tweaking the �_i's and seeing how V changes.
 * But, they are themselves sensitive to uncertainty in the data, particularly
 * the spreads.  We assume 50% stddev (parallel) uncertainty in the par curve P
 * ...
 * 
 * -     p*_i,t = P_i,t (1 + S)  for all dates t  where S ~ 0 � 50%�
 * 
 * ... and we model the error in the a*_i as linear in S:
 * 
 * -     a*_i = A_i + k_i S
 * 
 * Then we can expand
 * 
 * -     Var[v* - V] = Var[ A_i (�*_i - �_i) ]
 * -                 = Var[ A_i (M + E_i) ]
 * -                 = Var[ \sum_i A_i M ] + Var[ \sum_i A_i E_i ] +
 *                     2 Cov[ \sum_i A_i M, \sum_i A_i E_i ]
 * -                 = Var[ (\sum_i A_i) M ] + \sum_i Var[ A_i E_i ] +
 *                     0  (see ***)
 * -                 = (Mean[ \sum_i A_i ]� + Var[ \sum_i A_i ]) Var[M] +
 *                     Var[ \sum_i (a*_i + k_i S) E_i ]          (1)
 * 
 * Now,
 * 
 * -     (1) = Var[ \sum_i a*_i E_i + \sum_i k_i S E_i ]
 * -         = Var[ \sum_i a*_i E_i ] + Var[ S \sum_i k_i E_i ] +
 *             Cov[ \sum_i a*_i E_i, S \sum_i k_i E_i ]
 * -         = (1 - �) \sum_i a*_i� + 50%� (1 - �) \sum k_i�
 * 
 * (the Cov[...] being zero by the same argument as ***).  If we split k_i into
 * a general level across all names, plus a name-specific deviation ...
 * 
 * -     k_i = b + c_i  where  \sum_i c_i = 0
 * 
 * ... then overall we get:
 * 
 * -     Var[v* - V] = � (\sum_i a*_i)�
 *                   + (1 - �) \sum_i a*_i�
 *                   + 50%� � (nb)�
 *                   + 50%� (1 - �) \sum_i k_i�
 *
 * where n = num. names = \sum_i 1
 * 
 * These correspond respectively to the "ParallelTerm", "IdiosyncraticTerm",
 * "dB/ds Parallel Term" and "dB/ds Idiosyncratic Term" in Gus's explanation.
 * The 50% s.d. assumed for the spreads uncertainty S is taken from the phrase
 * "broadly comparable to a 50% spread shift".  "Standard_dB/ds" is 50% nb.
 * 
 * (***) To see Cov[ M \sum_i A_i, \sum_i A_i E_i ] = 0: write Q = \sum_i A_i,
 * R = \sum_i A_i E_i and note that M is independent of Q and R, and E[R] is
 * zero. Then Cov[ MQ, R ] = E[M QR] - [MQ] E[R] = E[M] E[QR] - E[MQ] 0 = 0
 * E[QR] = 0
 * 
 *
 * <H3>The greeks we propose to output</H3>
 *
 * In terms of partial derivatives, the formula derived above for the variance
 * in the pricing error is a function of
 * 
 * -     a_i = �V / ��_i for all names i
 * -     k_i = �V / ��_i �S for all names i
 * 
 * (recalling that S is parallel shift in spreads).  nb could be obtained as
 * \sum_i k_i, since the c_i cancel out by definition.  However, in Gus's
 * email, b is in fact estimated directly via an additional tweak.  ["dB/ds
 * Parallel Term".  Assuming that "betaShifted1%" means the same 0.01 absolute
 * tweak to all �_i, the c_i drop out.]  So we must add
 * 
 * -     nb = �V / �� �S
 * 
 * , �� being understood to mean "�_i += 0.01 for all names i".
 * 
 * Sensitivity to spreads is clearly intended to be evaluated with respect to
 * parallel shifts (i.e. all dates at once) only: �S in our notation.  The only
 * question is whether �V / ��_i �S is with respect to a shift in _all_ spread
 * curves, or just in that for name i.  Given that Gus's explanation refers
 * simply to "spreads levels" without qualification, it seems likely that the
 * former is intended.  In any case it's the only possibility that makes sense
 * for �V / �� �S (=~ "Standard_dB/ds").
 * 
 * The format in which we propose to report these greeks is
 * <TABLE>
 * <TR><TH>Symbol<TH>Name in Gus's email<TH>Packet name<TH>Output names<TH>Units
 * <TR>
 *   <TD>�V / ��_i = a_i
 *   <TD>ai
 *   <TD>CCM_BETA_SENS
 *   <TD>all names in basket
 *   <TD>natural �
 * <TR>
 *   <TD>�V / ��_i �S = k_i
 *   <TD>dB/ds for name i
 *   <TD>CCM_D_BETA_D_PAR_SPREAD_RHO_PROP
 *   <TD>all names in basket
 *   <TD>natural � � logarithmic spread
 * <TR>
 *   <TD>�V / �� �S = nb
 *   <TD>dB/ds
 *   <TD>Instrument
 *   <TD>CCM_D_BETA_BASKET_D_PAR_SPREAD_RHO_PROP
 *   <TD>natural � � logarithmic spread
 * </TABLE>
 * 
 * All the outputs are (scalar) doubles, not arrays (but the packages generated
 * by the first two greeks contain an output for each name).  They're all
 * expressed in units of absolute changes to �_i, i.e. multiplying by 0.01
 * gives sensitivity to �_i += 0.01, and a relative shift in spreads,
 * e.g. multiplying by 5% gives sensitivity to spreads changing from 20bp to
 * 21bp.
 * 
 *
 * <H3>Details of the tweaking</H3>
 *
 * (Incidentally, I assume that all the tweaks are intended to be "one-sided"?)
 * 
 * In estimating �V / ��_i (i.e. A_i), we are asked to follow a "proportional"
 * tweaking scheme ...
 * 
 * -     �~_j = 1 - (1 � 0.2) � (1 - �*_j) if j = i ("Proportional beta tweak"),
 *       �*_j otherwise
 * 
 * ... which moves �_i relatively closer to or further from unity, and set
 * 
 * -     �V / ��_i "=" (v(�~_i) - v(�*_i)) / (�~_i - �*_i)
 * 
 * To estimate �V / ��_i �S (i.e. k_i), we are asked to set ...
 * 
 * -     p~_j,t = p*_j,t � 0.99  for all names j, dates t         ("spread down 1%")
 * 
 * -     �~_j = 1 - 0.8 � (1 - �*_j) if j = i ("up beta tweak"),
 *              �*_j otherwise
 * 
 * -     �V / ��_i �S "=" ((v(�~, p~) - v(�*, p~)) / (�~_i - �*_i) -
 *                          (v(�~, p*) - v(�*, p*)) / (�~_i - �*_i)
 *                                          / 0.01)
 * 
 * ... for each i.  The numerator is "ai_betaup(s_d) - ai_betaup(s)".
 * 
 * To estimate �V / �� �S (i.e. b), we are asked to set this time
 * 
 * -     p~_j,t = p_j,t � 0.99 for all names j, dates t
 *       (magnitude actually not specified, "s_d" seems to imply down though)
 *    
 * -     �~_j = �*_j + 0.01  for all j ("betaShifted1%")
 * 
 * -     �V / �� �S "=" ((v(�~, p~) - v(�*, p~)) / 0.01 -
 *                        (v(�~, p*) - v(�*, p*)) / 0.01
 *                                    / 0.01)
 * 
 * The numerator is "B(s) - B(s_d)", and v(�~, p) - v(�*, p) is
 * "V(betaShifted1%) - V(betaCurrent)".  Note that the � tweak does not in this
 * case follow the "proportional" scheme: it can't, because of course not all
 * the �_i would be tweaked by the same absolute amount, and that's a
 * requirement if the idea is to estimate nb.  It would, of course, be possible
 * to use e.g. the smallest tweak suggested by the proportional rule for any of
 * the names.  If any of the �_i are >= 0.99, so that the tweak would
 * invalidate the model, the output is set to "untweakable".
 * 
 *
 * <H3>How the greeks can be used to calculate correlation reserves</H3>
 *
 * In the notation for reserves calculation in Gus's email:
 * 
 * ParallelTerm = � � (sum [ CCM_BETA_SENS.i for i in basket ])�
 * IdiosyncraticTerm = (1 - �) � sum [ CCM_BETA_SENS.i� for i in basket ]
 * dB/ds ParallelTerm = � � (CCM_D_BETA_BASKET_D_PAR_SPREAD_RHO_PROPORTIONAL � 50%)�
 * dB/ds IdiosyncraticTerm = (1 - �) � sum [ CCM_D_BETA_D_PAR_SPREAD_RHO_PROPORTIONAL.i� for i in basket ]
 * 
 * The sum of these is the pricing error variance assuming unit overall
 * variance on � (to plug in a different � variance, just scale it).  The
 * square root of that (plus an "SS Floor Bid/Offer" term) is the reserve.
 */

class CCMDBetaDSpread: public Sensitivity {

public:

    static CClassConstSP const TYPE;
    static const string NAME;

    // Just for GenericSensitivityFactory:
    static const double DEFAULT_SHIFT;

    /**
     * @param spreadScaleShift  The relative amount to tweak the spread curves
     *                          (via ParSpreadPropShift)
     * @param betaScaleShift    The relative amount to tweak the betas
     *                          (via CCMBetaSens)
     */

    CCMDBetaDSpread(double spreadScaleShift = DEFAULT_SHIFT,
                    double betaScaleShift = 0.2):
        Sensitivity(TYPE),
        spreadScaleShift(spreadScaleShift),
        betaScaleShift(betaScaleShift)
    {}

    /** identifies the name used for storing associated results in the output*/

    const string& getSensOutputName() const {
        return NAME;
    }

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one (return value: false) */

    bool discreteShift() const {
        return false;
    }

    /** calculates sensitivity - invoked by calculateSens */

    void calculate(TweakGroup*  tweakGroup,
                   CResults*    results) {
        static const string method("CCMDBetaDSpread::calculate");
        try {
            const string& packetName = getPacketName();

            double sprDownShift = getSpreadScaleShift();
            ParSpreadPropShift sprDown(sprDownShift);
            CCMBetaSens vByBeta(getBetaScaleShift());

            SensMgrOpt sensMgr(tweakGroup);

            OutputNameArrayConstSP betaNames(vByBeta.names(tweakGroup));

            if (sensMgr.allNames(&sprDown)->empty() || betaNames->empty()) {
                results->storeNotApplicable(this);
                return;
            }

            if (!getControl()->sensitivityRequested(vByBeta.getClass())) {
                getControl()->removePacketAfterCalc(vByBeta.getPacketName());
            }

            vByBeta.calculateSens(tweakGroup->getModel(),
                                  tweakGroup->getInstrument(),
                                  getControl(), results);

            sensMgr.shift(&sprDown);

            Results downResults;
            try {
                tweakGroup->getModel()->Price(tweakGroup->getInstrument(),
                                              getControl(), &downResults);
                vByBeta.calculateSens(tweakGroup->getModel(),
                                      tweakGroup->getInstrument(),
                                      getControl(), &downResults);

                sensMgr.restore();
            }
            catch (exception& e) {
                sensMgr.restore();
                for (int i = 0; i < betaNames->size(); ++i) {
                    results->storeGreek(IObjectSP(new Untweakable(e)),
                                        packetName, (*betaNames)[i]);
                }
                return;
            }

            for (int i = 0; i < betaNames->size(); ++i) {
                try {
                    double vByBetaBase = results->retrieveScalarGreek(
                        vByBeta.getPacketName(), (*betaNames)[i]);
                    double vByBetaDown = downResults.retrieveScalarGreek(
                        vByBeta.getPacketName(), (*betaNames)[i]);
                    // normally you'd put this test higher up outside of the 
                    // loop but originally the code threw an exception here
                    // so I'm leaving it like that (doesn't really matter)
                    if (Maths::isZero(sprDownShift)){
                        throw ModelException(method,
                                             "spreadScaleShift is zero");
                    }
                    results->storeScalarGreek(
                        (vByBetaDown - vByBetaBase) / sprDownShift,
                        getPacketName(),
                        (*betaNames)[i]);
                }
                catch (exception& e) {
                    results->storeGreek(IObjectSP(new Untweakable(e)),
                                        packetName, (*betaNames)[i]);
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /**
     * The proportion by which the spread curves are tweaked (via ParSpreadPropShift)
     */

    double getSpreadScaleShift() const {
        return spreadScaleShift;
    }

    /**
     * The proportion by which the betas are tweaked (via CCMBetaSens)
     */

    double getBetaScaleShift() const {
        return betaScaleShift;
    }

private:

    static IObject *defaultCCMDBetaDSpread() {
        return new CCMDBetaDSpread();
    }

    /** Invoked when this class is 'loaded' */

    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CCMDBetaDSpread, clazz);
        SUPERCLASS(Sensitivity);
        EMPTY_SHELL_METHOD(defaultCCMDBetaDSpread);
        FIELD(spreadScaleShift, "Relative amount to shift spreads by");
        FIELD(betaScaleShift, "Relative amount to shift beta by");
        FIELD_MAKE_OPTIONAL(betaScaleShift);

        SensitivityFactory::addSens(
            NAME, new GenericSensitivityFactory<CCMDBetaDSpread>(),
            new CCMDBetaDSpread(),
            ParSpreadPropShift::TYPE);
    }

    double spreadScaleShift;
    double betaScaleShift;
};

const string CCMDBetaDSpread::NAME = "CCM_D_BETA_D_PAR_SPREAD_RHO_PROPORTIONAL";

const double CCMDBetaDSpread::DEFAULT_SHIFT = 0.01;

CClassConstSP const CCMDBetaDSpread::TYPE = CClass::registerClassLoadMethod(
    "CCMDBetaDSpread", typeid(CCMDBetaDSpread), CCMDBetaDSpread::load);

/**
 * Included in RiskMgrLib::linkInClasses() to force CCMDBetaDSpread
 * to get linked into the Windows exe.
 */

bool CCMDBetaDSpreadLinkIn() {
    return CCMDBetaDSpread::TYPE != NULL;
}

DRLIB_END_NAMESPACE

