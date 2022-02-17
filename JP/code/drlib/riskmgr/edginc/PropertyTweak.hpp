/**
 * @file PropertyTweak.hpp
 */

#ifndef DRLIB_PropertyTweak_H
#define DRLIB_PropertyTweak_H

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(OutputName)

/**
 * Parameters for a tweak to a RiskProperty of a market object
 *
 * This is the argument to ITweakableWithRespectTo<TAG>::sensShift(),
 * i.e., what market objects get to see when they're implementing a tweak to a
 * property. For instance, Equity::sensShift(const PropertyTweak<Spot> &shift)
 * implements Spot tweaks to Equity.
 *
 * You probably don't want to construct one of these directly: it's really used
 * as a base class for PropertyTweakHypothesis<TAG> so that market objects
 * don't have to known about the latter.
 *
 * See IRiskQuantityFactory for an overview of how these classes fit into
 * the "declarative" sensitivity framework.
 */

template <class TAG>
struct PropertyTweak {

    /**
     * Tag class for the property being tweaked
     *
     * Could be used to specify e.g. proportional vs absolute tweaks, or
     * whatever.  For most properties it's an empty class --- see e.g. Spot.
     */

    typedef TAG Tag;
    typedef smartConstPtr<Tag> TagConstSP;

    /**
     * Info needed in addition to TAG type and coefficient to fully
     * specify the tweak
     *
     * For term-structured properties (like VolPointwise), this is
     * ExpiryWindow.  For scalar properties (like Spot) it's Void.
     */

    typedef typename TAG::Qualifier Qualifier;
    typedef smartConstPtr<Qualifier> QualifierConstSP;

    /**
     * Constructor (but you more often mean to use a PropertyTweakHypothesis)
     */

    PropertyTweak(double coefficient,
                  OutputNameConstSP marketDataName = OutputNameConstSP(),
                  QualifierConstSP qualifier = QualifierConstSP(),
                  TagConstSP tag = TagConstSP()):
        tag(tag),
        marketDataName(marketDataName),
        qualifier(qualifier),
        coefficient(coefficient)
    {}

    /**
     * (Default constructor, for reflection pre-construction.)
     */

    PropertyTweak(): coefficient(0) {}

    /**
     * Tag object for the property being tweaked
     *
     * Could be used to specify e.g. proportional vs absolute tweaks, or
     * whatever.  For most properties the type is an empty class --- see
     * e.g. Spot --- and 'tag' is just a null smartPtr.
     */

    TagConstSP tag;

    /**
     * Name of object being tweaked
     *
     * Ideally we wouldn't need this, since objects being tweaked ought not to
     * care what they're called, but some intra-object delegation mechanisms
     * require it.
     */

    OutputNameConstSP marketDataName;

    /**
     * Info needed in addition to TAG type and coefficient to fully
     * specify the tweak
     *
     * For term-structured properties (like VolPointwise), Qualifier is
     * ExpiryWindow.  For scalar properties (like Spot) it's Void and
     * 'qualifier' is just a null smartPtr.
     */

    QualifierConstSP qualifier;

    /**
     * Amount by which to tweak
     *
     * What this means is up to the tweaked objects' implementations of
     * ITweakableWithRespectTo<TAG>::sensShift().
     *
     * For instance, Spot is tweaked proportionally: see
     * Equity::sensShift(const PropertyTweak<Spot> &), which does
     * <TT>spot *= (1 + coefficient)</TT>.  For other properties,
     * coefficient may be treated as an absolute perturbation, or
     * whatever.
     *
     * (Note that the TweakOutcome::distance() returned from
     * ITweakableWithRespectTo<TAG>::sensShift() is not necessarily
     * equal to coefficient: while Spot tweaks are proportional,
     * they report absolute distance back.)
     */

    double coefficient;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
