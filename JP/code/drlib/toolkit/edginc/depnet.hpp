/**
 * @file depnet.hpp
 */

#ifndef QLIB_depnet_H
#define QLIB_depnet_H

#include <vector>
#include "edginc/refCountPtr.hpp"
#include "edginc/Object.hpp"
#include "edginc/Array.hpp"
#include "edginc/TRACE.hpp"
#include <boost/type_traits.hpp>

DRLIB_BEGIN_NAMESPACE

/**
 * An experimental utility for building Excel-/Athena-like "dependency graphs"
 *
 * Presently it's aimed at dealing systematically with the logic of recomputing
 * and caching the various steps in a pricing computation, when the underlying
 * market data is tweaked.  See MCPathConfigLV for an example.
 */

namespace depnet {

#define P refCountPtr

// 
// =====================
//  IfPOD::initialize()
// =====================
//
// A utility which we use below to initialize a member (Var::_values) explicitly
// if it's a "plain old data" type like int, or leave it to initialize itself
// if it's a class
// 

// Sadly BOOST_HAS_TRIVIAL_CONSTRUCTOR doesn't work properly (nor do
// any of the other Boost type_traits macros, go figure)

template <class T>
struct typeTraits {
  enum { isPOD = 0 };
};

template <> struct typeTraits<bool> { enum { isPOD = 1 }; };
template <> struct typeTraits<char> { enum { isPOD = 1 }; };
template <> struct typeTraits<unsigned char> { enum { isPOD = 1 }; };
template <> struct typeTraits<int> { enum { isPOD = 1 }; };
template <> struct typeTraits<unsigned> { enum { isPOD = 1 }; };
template <> struct typeTraits<long> { enum { isPOD = 1 }; };
template <> struct typeTraits<unsigned long> { enum { isPOD = 1 }; };
template <> struct typeTraits<float> { enum { isPOD = 1 }; };
template <> struct typeTraits<double> { enum { isPOD = 1 }; };
template <> struct typeTraits<long double> { enum { isPOD = 1 }; };
template <class T> struct typeTraits<T *> { enum { isPOD = 1 }; };
template <class T> struct typeTraits<const T *> { enum { isPOD = 1 }; };

template <class T, bool ISPOD = typeTraits<T>::isPOD>
struct ifPOD;

template <class T>
struct ifPOD<T, true> {
    static void initialize(T &x) {
        x = T();
    }
};

template <class T>
struct ifPOD<T, false> {
    static void initialize(T &x) {
    }
};

// 
// =============
//  AbstractVar
// =============
// 

struct AbstractVar {

    const char *name;
    enum State { base = 0, alt = 1 };

private:

    AbstractVar(const AbstractVar &);
    AbstractVar &operator =(const AbstractVar &);

protected:

    bool _locked;
    mutable bool _dirty, _wouldBeDirty;
    mutable State _state;
    mutable bool _valuesDirty[2];
    mutable unsigned _serials[2];
    mutable auto_ptr<ModelException> _exceptions[2];

    struct Ante {
        mutable unsigned serials[2];
        P<AbstractVar> var;

        Ante(const P<AbstractVar> &var = P<AbstractVar>()): var(var) {
            serials[base] = serials[alt] = 0;
        }
    };

    vector<Ante> _antes;
    vector<AbstractVar *> _deps;

    AbstractVar():
        name(""),
        _locked(false),
        _dirty(true),
        _wouldBeDirty(true),
        _state(base)
    {
        _serials[base] = 0;
        _serials[alt] = 2000000000; // TRACEs nicer than 0x80000000
        _valuesDirty[base] = _valuesDirty[alt] = true;
    }

    void watch(const P<AbstractVar> &v) {
        _antes.push_back(Ante(v));
        v->_deps.push_back(this);
    }

    void markDepsDirty() const {
        for (size_t d = 0; d < _deps.size(); ++d)
            _deps[d]->markDirty();
    }

    void markDirty() const {
        if (!_dirty) {
            _wouldBeDirty = true;
            if (!_locked) {
                _dirty = true;
                markDepsDirty();
            }
        }
    }

    void setFlags() const {
        for (size_t a = 0; a < _antes.size(); ++a)
            _antes[a].var->ensureUpToDate();

        _state = base;
        for (size_t a = 0; a < _antes.size(); ++a)
            if (_antes[a].var->_state == alt) _state = alt;

        for (size_t a = 0; a < _antes.size(); ++a) {
            const AbstractVar &ant = *_antes[a].var;
            if (ant._serials[ant._state] != _antes[a].serials[ant._state])
                _valuesDirty[_state] = true;
        }

        if (_state == alt)
            for (size_t a = 0; a < _antes.size(); ++a) {
                const AbstractVar &ant = *_antes[a].var;
                if (ant._valuesDirty[base] ||
                        ant._serials[base] != _antes[a].serials[base])
                    _valuesDirty[base] = true;
            } 
    }

    virtual void propagate() const = 0;

    void ensureUpToDate() const {
        if (_dirty) {
            setFlags();
            if (_valuesDirty[_state]) {
                propagate();
                commitValue();
            }

            _dirty = _wouldBeDirty = false;
        }
    }

    void commitValue() const {
        for (size_t a = 0; a < _antes.size(); ++a) {
            const AbstractVar &ante = *_antes[a].var;
            _antes[a].serials[_state] = ante._serials[ante._state];
        }

        _valuesDirty[_state] = false;
    }

public:

    void lockValue() {
        if (!_locked) {
            ensureUpToDate();
            _locked = true;
        }
    }

    void unlockValue() {
        _locked = false;
        if (_wouldBeDirty) markDirty();
    }

    virtual ~AbstractVar() {
        for (size_t a = 0; a < _antes.size(); ++a) {
            vector<AbstractVar *> &ds = _antes[a].var->_deps;
            for (size_t d = 0; d < ds.size(); ++d)
                if (ds[d] == this) ds.erase(ds.begin() + d);
        }
    }
};

// 
// =====
//  Var
// =====
// 

template <class T>
class Var: public AbstractVar {
    mutable T _values[2];

    bool (*_equals)(const T &, const T &);

    // precondition: _equals != NULL, _state and _valuesDirty up to date

    bool maybeNewValue(const T &v) const {
        if (_state == alt && !_valuesDirty[base] &&
                !_exceptions[base].get() && (*_equals)(v, _values[base])) {
            _state = base;
            return false;
        }
        else if (!_exceptions[_state].get() && (*_equals)(v, _values[_state])) {
            return false;
        }
        else {
            _exceptions[_state].reset();
            _values[_state] = v;
            ++_serials[_state];
            return true;
        }
    }

    virtual void computeValue(T &) const {}

    void computeValueOrException(T &v) const {
        try {
            computeValue(v);
            _exceptions[_state].reset();
        }
        catch (exception &e) {
            _exceptions[_state].reset(new ModelException(e, "Updating cell"));
        }
    }

    // precondition: _state and _valuesDirty up to date

    void propagate() const {
        if (_equals) {
            T v;
            computeValueOrException(v);
            if (!_exceptions[_state].get()) maybeNewValue(v);
        }
        else {
            computeValueOrException(_values[_state]);
            ++_serials[_state];
        }
    }

public:

    typedef T Value;

    Var(const T &value,
        bool (*noChangeTest)(const T &, const T &) = 0):
        _equals(noChangeTest)
    {
        setValue(value, false);
        ifPOD<T>::initialize(_values[alt]);
    }

    Var(): _equals(0) {
        ifPOD<T>::initialize(_values[base]);
        ifPOD<T>::initialize(_values[alt]);
    }

    void setNoChangeTest(bool (*equals)(const T &, const T &)) {
        _equals = equals;
    }

    const T &value() const {
        ensureUpToDate();
        if (_exceptions[_state].get())
            throw *_exceptions[_state];
        return _values[_state];
    }

    void setValue(const T &v, bool altState, bool noteChange = true) {
        State oldState = _state;
        _state = altState ? alt : base;
        bool wasDirty = _dirty;

        if (_dirty) setFlags();

        bool changed;

        if (_equals && noteChange)
            changed = maybeNewValue(v);
        else {
            changed = noteChange;
            _exceptions[_state].reset();
            _values[_state] = v;
            if (noteChange) ++_serials[_state];
        }

        commitValue();
        _dirty = _wouldBeDirty = false;

        if (!wasDirty && (changed || oldState != _state))
            markDepsDirty();
    }

    void pretendChanged() {
        if (_dirty) setFlags();
        ++_serials[_state];
        if (!_dirty) markDepsDirty();
    }
};

// 
// -----------------
//  setNoChangeTest
// -----------------
// 

template <class T>
struct targetEqual {
    static bool equal(const T &a, const T &b) {
        return a == b;
    }
};

template <class T>
struct targetEqual<const T *> {
    static bool equal(const T * const &a, const T * const&b) {
        return a == b || (a && b && targetEqual<T>::equal(*a, *b));
    }
};

template <class T>
struct targetEqual<T *> {
    static bool equal(const T *&a, const T *&b) {
        return a == b || (a && b && targetEqual<T>::equal(*a, *b));
    }
};

template <>
struct targetEqual<IObject> {
    static bool equal(const IObject &a, const IObject &b) {
        return a.equalTo(&b);
    }
};

template <class T>
struct targetEqual<smartPtr<T> > {
    static bool equal(const smartPtr<T> &a, const smartPtr<T> &b) {
        return targetEqual<T *>::equal(a.get(), b.get());
    }
};

template <class T>
struct targetEqual<smartConstPtr<T> > {
    static bool equal(const smartConstPtr<T> &a, const smartConstPtr<T> &b) {
        return targetEqual<T *>::equal(a.get(), b.get());
    }
};

template <class T>
void setNoChangeTest(Var<T> &v) {
    v.setNoChangeTest(&targetEqual<T>::equal);
}

template <class VarT>
const P<VarT> &withNoChangeTest(const P<VarT> &v) {
    setNoChangeTest(*v);
    return v;
}

template <class T>
Var<T> *withNoChangeTest(Var<T> *v) {
    setNoChangeTest(*v);
    return v;
}

// 
// ------
//  root
// ------
// 

template <class T>
P<Var<T> > root(const T &value) {
    return P<Var<T> >(new Var<T>(value));
}

// 
// ===========
//  MappedVar
// ===========
// 

template <class T>
struct derefed {
    typedef T t;
};

template <class T>
struct derefed<const T &> {
    typedef T t;
};

template <class T>
struct ptrOrRef {
    static const T &ref(const T &r) { return r; }
};

template <class T>
struct ptrOrRef<T *> {
    static const T &ref(const T *p) { return *p; }
};

#define R(T) typename derefed<T>::t

// retrieve VALUE from smartPtr<Var<VALUE> >
#define V(SP) typename SP::element_type::Value

template <class T, class F,
          class A1 = void, class A2 = void, class A3 = void, class A4 = void,
          class A5 = void, class A6 = void, class A7 = void, class A8 = void>
struct MappedVar;

// --- Code generator

// import re
// 
// def reps(n, s):
//   m = min(map(int, re.findall('[1-9]', s)))
//   return [ s.replace(`m`, `i`) for i in range(m, n+1) ]
// 
// def sub(n, s):
//   return re.sub('«([^»]*)»',
//                 lambda t: (', ', '')[',' in t.group(1) or ';' in t.group(1)].join(reps(n, t.group(1))),
//                 s).replace('#', `n`)
// 
// for n in range(1, 7):
//   print sub(n, '''
// // 
// // ------------
// //  # argument
// // ------------
// // 
// 
// template <«class A1, class B1, »class T>
// void apply(void (*f)(«B1, »T&)«, const A1 &a1», T &r) {
//     f(«a1, »r);
// }
// 
// template <«class A1, class B1, »class T, class U>
// void apply(U (*f)(«B1»)«, const A1 &a1», T &r) {
//     r = f(«a1»);
// }
// 
// template <«class A1, class B1, »class T, class U>
// void apply(U (B1::*f)(«B2») const, «const A1 &a1, »T &r) {
//     r = (ptrOrRef<A1>::ref(a1).*f)(«a2»);
// }
// 
// template <«class A1, class B1, »class T>
// void apply(void (B1::*f)(«B2, »T &) const, «const A1 &a1, »T &r) {
//     (ptrOrRef<A1>::ref(a1).*f)(«a2, »r);
// }
// 
// template <«class A1, class B1, »class T>
// void apply(void (T::*setter)(«B1»), «const A1 &a1, »T &r) {
//     (r.*setter)(«a1»);
// }
// 
// template <class T, class F«, class A1»>
// struct MappedVar<T, F«, A1»>: Var<T> {
// 
// private:
// 
//     F _f;
//    « P<Var<A1> > _a1;»
// 
//     void computeValue(T &value) const {
//         apply(_f, «_a1->value(), »value);
//     }
// 
// public:
// 
//     MappedVar(F f«, P<Var<A1> > a1»):
//         _f(f)«, _a1(a1)»
//     {
//        « watch(a1);»
//     }
// };
// 
// template <class T«, class A1, class B1»>
// P<MappedVar<T, void (*)(«B1, »T &)«, V(A1)»> > mapped(
//         void (*f)(«B1, »T &)«, const A1 &a1») {
//     return P<MappedVar<T, void (*)(«B1, »T &)«, V(A1)»> >(new MappedVar<T, void (*)(«B1, »T &)«, V(A1)»>(f«, a1»));
// }
// 
// template <class T«, class A1, class B1»>
// P<MappedVar<R(T), T (*)(«B1»)«, V(A1)»> > mapped(
//         T (*f)(«B1»)«, const A1 &a1») {
//     return P<MappedVar<R(T), T (*)(«B1»)«, V(A1)»> >(new MappedVar<R(T), T (*)(«B1»)«, V(A1)»>(f«, a1»));
// }
// 
// template <class T«, class A1, class B1»>
// P<MappedVar<R(T), T (B1::*)(«B2») const«, V(A1)»> > mapped(
//         T (B1::*f)(«B2») const«, const A1 &a1») {
//     return P<MappedVar<R(T), T (B1::*)(«B2») const«, V(A1)»> >(
//         new MappedVar<R(T), T (B1::*)(«B2») const«, V(A1)»>(f«, a1»));
// }
// 
// template <class T«, class A1, class B1»>
// P<MappedVar<R(T), void (B1::*)(«B2, »T &) const«, V(A1)»> > mapped(
//         void (B1::*f)(«B2, »T &) const«, const A1 &a1») {
//     return P<MappedVar<R(T), void (B1::*)(«B2, »T &) const«, V(A1)»> >(
//         new MappedVar<R(T), void (B1::*)(«B2, »T &) const«, V(A1)»>(f«, a1»));
// }
// 
// template <class T«, class A1, class B1»>
// P<MappedVar<T, void (T::*)(«B1»)«, V(A1)»> > setWith(
//         void (T::*setter)(«B1»)«, const A1 &a1») {
//     return P<MappedVar<T, void (T::*)(«B1»)«, V(A1)»> >(
//         new MappedVar<T, void (T::*)(«B1»)«, V(A1)»>(setter«, a1»));
// }''')

// --- start of generated code

// 
// ------------
//  1 argument
// ------------
// 

template <class A1, class B1, class T>
void apply(T B1::* const m, const A1 &a1, T &r) {
    r = ptrOrRef<A1>::ref(a1).*m;
}

template <class A1, class B1, class T>
void apply(void (*f)(B1, T&), const A1 &a1, T &r) {
    f(a1, r);
}

template <class A1, class B1, class T, class U>
void apply(U (*f)(B1), const A1 &a1, T &r) {
    r = f(a1);
}

template <class A1, class B1, class T, class U>
void apply(U (B1::*f)() const, const A1 &a1, T &r) {
    r = (ptrOrRef<A1>::ref(a1).*f)();
}

template <class A1, class B1, class T>
void apply(void (B1::*f)(T &) const, const A1 &a1, T &r) {
    (ptrOrRef<A1>::ref(a1).*f)(r);
}

template <class A1, class B1, class T>
void apply(void (T::*setter)(B1), const A1 &a1, T &r) {
    (r.*setter)(a1);
}

template <class T, class F, class A1>
struct MappedVar<T, F, A1>: Var<T> {

private:

    F _f;
    P<Var<A1> > _a1;

    void computeValue(T &value) const {
        apply(_f, _a1->value(), value);
    }

public:

    MappedVar(F f, P<Var<A1> > a1):
        _f(f), _a1(a1)
    {
        watch(a1);
    }
};

template <class T, class A1, class B1>
P<MappedVar<T, void (*)(B1, T &), V(A1)> > mapped(
        void (*f)(B1, T &), const A1 &a1) {
    return P<MappedVar<T, void (*)(B1, T &), V(A1)> >(new MappedVar<T, void (*)(B1, T &), V(A1)>(f, a1));
}

template <class T, class A1, class B1>
P<MappedVar<R(T), T (*)(B1), V(A1)> > mapped(
        T (*f)(B1), const A1 &a1) {
    return P<MappedVar<R(T), T (*)(B1), V(A1)> >(new MappedVar<R(T), T (*)(B1), V(A1)>(f, a1));
}

template <class T, class A1, class B1>
P<MappedVar<R(T), T (B1::*)() const, V(A1)> > mapped(
        T (B1::*f)() const, const A1 &a1) {
    return P<MappedVar<R(T), T (B1::*)() const, V(A1)> >(
        new MappedVar<R(T), T (B1::*)() const, V(A1)>(f, a1));
}

template <class T, class A1, class B1>
P<MappedVar<R(T), void (B1::*)(T &) const, V(A1)> > mapped(
        void (B1::*f)(T &) const, const A1 &a1) {
    return P<MappedVar<R(T), void (B1::*)(T &) const, V(A1)> >(
        new MappedVar<R(T), void (B1::*)(T &) const, V(A1)>(f, a1));
}

template <class T, class A1, class B1>
P<MappedVar<T, void (T::*)(B1), V(A1)> > setWith(
        void (T::*setter)(B1), const A1 &a1) {
    return P<MappedVar<T, void (T::*)(B1), V(A1)> >(
        new MappedVar<T, void (T::*)(B1), V(A1)>(setter, a1));
}

// 
// ------------
//  2 argument
// ------------
// 

template <class A1, class B1, class A2, class B2, class T>
void apply(void (*f)(B1, B2, T&), const A1 &a1, const A2 &a2, T &r) {
    f(a1, a2, r);
}

template <class A1, class B1, class A2, class B2, class T, class U>
void apply(U (*f)(B1, B2), const A1 &a1, const A2 &a2, T &r) {
    r = f(a1, a2);
}

template <class A1, class B1, class A2, class B2, class T, class U>
void apply(U (B1::*f)(B2) const, const A1 &a1, const A2 &a2, T &r) {
    r = (ptrOrRef<A1>::ref(a1).*f)(a2);
}

template <class A1, class B1, class A2, class B2, class T>
void apply(void (B1::*f)(B2, T &) const, const A1 &a1, const A2 &a2, T &r) {
    (ptrOrRef<A1>::ref(a1).*f)(a2, r);
}

template <class A1, class B1, class A2, class B2, class T>
void apply(void (T::*setter)(B1, B2), const A1 &a1, const A2 &a2, T &r) {
    (r.*setter)(a1, a2);
}

template <class T, class F, class A1, class A2>
struct MappedVar<T, F, A1, A2>: Var<T> {

private:

    F _f;
    P<Var<A1> > _a1; P<Var<A2> > _a2;

    void computeValue(T &value) const {
        apply(_f, _a1->value(), _a2->value(), value);
    }

public:

    MappedVar(F f, P<Var<A1> > a1, P<Var<A2> > a2):
        _f(f), _a1(a1), _a2(a2)
    {
        watch(a1); watch(a2);
    }
};

template <class T, class A1, class B1, class A2, class B2>
P<MappedVar<T, void (*)(B1, B2, T &), V(A1), V(A2)> > mapped(
        void (*f)(B1, B2, T &), const A1 &a1, const A2 &a2) {
    return P<MappedVar<T, void (*)(B1, B2, T &), V(A1), V(A2)> >(new MappedVar<T, void (*)(B1, B2, T &), V(A1), V(A2)>(f, a1, a2));
}

template <class T, class A1, class B1, class A2, class B2>
P<MappedVar<R(T), T (*)(B1, B2), V(A1), V(A2)> > mapped(
        T (*f)(B1, B2), const A1 &a1, const A2 &a2) {
    return P<MappedVar<R(T), T (*)(B1, B2), V(A1), V(A2)> >(new MappedVar<R(T), T (*)(B1, B2), V(A1), V(A2)>(f, a1, a2));
}

template <class T, class A1, class B1, class A2, class B2>
P<MappedVar<R(T), T (B1::*)(B2) const, V(A1), V(A2)> > mapped(
        T (B1::*f)(B2) const, const A1 &a1, const A2 &a2) {
    return P<MappedVar<R(T), T (B1::*)(B2) const, V(A1), V(A2)> >(
        new MappedVar<R(T), T (B1::*)(B2) const, V(A1), V(A2)>(f, a1, a2));
}

template <class T, class A1, class B1, class A2, class B2>
P<MappedVar<R(T), void (B1::*)(B2, T &) const, V(A1), V(A2)> > mapped(
        void (B1::*f)(B2, T &) const, const A1 &a1, const A2 &a2) {
    return P<MappedVar<R(T), void (B1::*)(B2, T &) const, V(A1), V(A2)> >(
        new MappedVar<R(T), void (B1::*)(B2, T &) const, V(A1), V(A2)>(f, a1, a2));
}

template <class T, class A1, class B1, class A2, class B2>
P<MappedVar<T, void (T::*)(B1, B2), V(A1), V(A2)> > setWith(
        void (T::*setter)(B1, B2), const A1 &a1, const A2 &a2) {
    return P<MappedVar<T, void (T::*)(B1, B2), V(A1), V(A2)> >(
        new MappedVar<T, void (T::*)(B1, B2), V(A1), V(A2)>(setter, a1, a2));
}

// 
// ------------
//  3 argument
// ------------
// 

template <class A1, class B1, class A2, class B2, class A3, class B3, class T>
void apply(void (*f)(B1, B2, B3, T&), const A1 &a1, const A2 &a2, const A3 &a3, T &r) {
    f(a1, a2, a3, r);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class T, class U>
void apply(U (*f)(B1, B2, B3), const A1 &a1, const A2 &a2, const A3 &a3, T &r) {
    r = f(a1, a2, a3);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class T, class U>
void apply(U (B1::*f)(B2, B3) const, const A1 &a1, const A2 &a2, const A3 &a3, T &r) {
    r = (ptrOrRef<A1>::ref(a1).*f)(a2, a3);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class T>
void apply(void (B1::*f)(B2, B3, T &) const, const A1 &a1, const A2 &a2, const A3 &a3, T &r) {
    (ptrOrRef<A1>::ref(a1).*f)(a2, a3, r);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class T>
void apply(void (T::*setter)(B1, B2, B3), const A1 &a1, const A2 &a2, const A3 &a3, T &r) {
    (r.*setter)(a1, a2, a3);
}

template <class T, class F, class A1, class A2, class A3>
struct MappedVar<T, F, A1, A2, A3>: Var<T> {

private:

    F _f;
    P<Var<A1> > _a1; P<Var<A2> > _a2; P<Var<A3> > _a3;

    void computeValue(T &value) const {
        apply(_f, _a1->value(), _a2->value(), _a3->value(), value);
    }

public:

    MappedVar(F f, P<Var<A1> > a1, P<Var<A2> > a2, P<Var<A3> > a3):
        _f(f), _a1(a1), _a2(a2), _a3(a3)
    {
        watch(a1); watch(a2); watch(a3);
    }
};

template <class T, class A1, class B1, class A2, class B2, class A3, class B3>
P<MappedVar<T, void (*)(B1, B2, B3, T &), V(A1), V(A2), V(A3)> > mapped(
        void (*f)(B1, B2, B3, T &), const A1 &a1, const A2 &a2, const A3 &a3) {
    return P<MappedVar<T, void (*)(B1, B2, B3, T &), V(A1), V(A2), V(A3)> >(new MappedVar<T, void (*)(B1, B2, B3, T &), V(A1), V(A2), V(A3)>(f, a1, a2, a3));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3>
P<MappedVar<R(T), T (*)(B1, B2, B3), V(A1), V(A2), V(A3)> > mapped(
        T (*f)(B1, B2, B3), const A1 &a1, const A2 &a2, const A3 &a3) {
    return P<MappedVar<R(T), T (*)(B1, B2, B3), V(A1), V(A2), V(A3)> >(new MappedVar<R(T), T (*)(B1, B2, B3), V(A1), V(A2), V(A3)>(f, a1, a2, a3));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3>
P<MappedVar<R(T), T (B1::*)(B2, B3) const, V(A1), V(A2), V(A3)> > mapped(
        T (B1::*f)(B2, B3) const, const A1 &a1, const A2 &a2, const A3 &a3) {
    return P<MappedVar<R(T), T (B1::*)(B2, B3) const, V(A1), V(A2), V(A3)> >(
        new MappedVar<R(T), T (B1::*)(B2, B3) const, V(A1), V(A2), V(A3)>(f, a1, a2, a3));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3>
P<MappedVar<R(T), void (B1::*)(B2, B3, T &) const, V(A1), V(A2), V(A3)> > mapped(
        void (B1::*f)(B2, B3, T &) const, const A1 &a1, const A2 &a2, const A3 &a3) {
    return P<MappedVar<R(T), void (B1::*)(B2, B3, T &) const, V(A1), V(A2), V(A3)> >(
        new MappedVar<R(T), void (B1::*)(B2, B3, T &) const, V(A1), V(A2), V(A3)>(f, a1, a2, a3));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3>
P<MappedVar<T, void (T::*)(B1, B2, B3), V(A1), V(A2), V(A3)> > setWith(
        void (T::*setter)(B1, B2, B3), const A1 &a1, const A2 &a2, const A3 &a3) {
    return P<MappedVar<T, void (T::*)(B1, B2, B3), V(A1), V(A2), V(A3)> >(
        new MappedVar<T, void (T::*)(B1, B2, B3), V(A1), V(A2), V(A3)>(setter, a1, a2, a3));
}

// 
// ------------
//  4 argument
// ------------
// 

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class T>
void apply(void (*f)(B1, B2, B3, B4, T&), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, T &r) {
    f(a1, a2, a3, a4, r);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class T, class U>
void apply(U (*f)(B1, B2, B3, B4), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, T &r) {
    r = f(a1, a2, a3, a4);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class T, class U>
void apply(U (B1::*f)(B2, B3, B4) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, T &r) {
    r = (ptrOrRef<A1>::ref(a1).*f)(a2, a3, a4);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class T>
void apply(void (B1::*f)(B2, B3, B4, T &) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, T &r) {
    (ptrOrRef<A1>::ref(a1).*f)(a2, a3, a4, r);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class T>
void apply(void (T::*setter)(B1, B2, B3, B4), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, T &r) {
    (r.*setter)(a1, a2, a3, a4);
}

template <class T, class F, class A1, class A2, class A3, class A4>
struct MappedVar<T, F, A1, A2, A3, A4>: Var<T> {

private:

    F _f;
    P<Var<A1> > _a1; P<Var<A2> > _a2; P<Var<A3> > _a3; P<Var<A4> > _a4;

    void computeValue(T &value) const {
        apply(_f, _a1->value(), _a2->value(), _a3->value(), _a4->value(), value);
    }

public:

    MappedVar(F f, P<Var<A1> > a1, P<Var<A2> > a2, P<Var<A3> > a3, P<Var<A4> > a4):
        _f(f), _a1(a1), _a2(a2), _a3(a3), _a4(a4)
    {
        watch(a1); watch(a2); watch(a3); watch(a4);
    }
};

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4>
P<MappedVar<T, void (*)(B1, B2, B3, B4, T &), V(A1), V(A2), V(A3), V(A4)> > mapped(
        void (*f)(B1, B2, B3, B4, T &), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) {
    return P<MappedVar<T, void (*)(B1, B2, B3, B4, T &), V(A1), V(A2), V(A3), V(A4)> >(new MappedVar<T, void (*)(B1, B2, B3, B4, T &), V(A1), V(A2), V(A3), V(A4)>(f, a1, a2, a3, a4));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4>
P<MappedVar<R(T), T (*)(B1, B2, B3, B4), V(A1), V(A2), V(A3), V(A4)> > mapped(
        T (*f)(B1, B2, B3, B4), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) {
    return P<MappedVar<R(T), T (*)(B1, B2, B3, B4), V(A1), V(A2), V(A3), V(A4)> >(new MappedVar<R(T), T (*)(B1, B2, B3, B4), V(A1), V(A2), V(A3), V(A4)>(f, a1, a2, a3, a4));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4>
P<MappedVar<R(T), T (B1::*)(B2, B3, B4) const, V(A1), V(A2), V(A3), V(A4)> > mapped(
        T (B1::*f)(B2, B3, B4) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) {
    return P<MappedVar<R(T), T (B1::*)(B2, B3, B4) const, V(A1), V(A2), V(A3), V(A4)> >(
        new MappedVar<R(T), T (B1::*)(B2, B3, B4) const, V(A1), V(A2), V(A3), V(A4)>(f, a1, a2, a3, a4));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4>
P<MappedVar<R(T), void (B1::*)(B2, B3, B4, T &) const, V(A1), V(A2), V(A3), V(A4)> > mapped(
        void (B1::*f)(B2, B3, B4, T &) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) {
    return P<MappedVar<R(T), void (B1::*)(B2, B3, B4, T &) const, V(A1), V(A2), V(A3), V(A4)> >(
        new MappedVar<R(T), void (B1::*)(B2, B3, B4, T &) const, V(A1), V(A2), V(A3), V(A4)>(f, a1, a2, a3, a4));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4>
P<MappedVar<T, void (T::*)(B1, B2, B3, B4), V(A1), V(A2), V(A3), V(A4)> > setWith(
        void (T::*setter)(B1, B2, B3, B4), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) {
    return P<MappedVar<T, void (T::*)(B1, B2, B3, B4), V(A1), V(A2), V(A3), V(A4)> >(
        new MappedVar<T, void (T::*)(B1, B2, B3, B4), V(A1), V(A2), V(A3), V(A4)>(setter, a1, a2, a3, a4));
}

// 
// ------------
//  5 argument
// ------------
// 

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class T>
void apply(void (*f)(B1, B2, B3, B4, B5, T&), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, T &r) {
    f(a1, a2, a3, a4, a5, r);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class T, class U>
void apply(U (*f)(B1, B2, B3, B4, B5), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, T &r) {
    r = f(a1, a2, a3, a4, a5);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class T, class U>
void apply(U (B1::*f)(B2, B3, B4, B5) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, T &r) {
    r = (ptrOrRef<A1>::ref(a1).*f)(a2, a3, a4, a5);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class T>
void apply(void (B1::*f)(B2, B3, B4, B5, T &) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, T &r) {
    (ptrOrRef<A1>::ref(a1).*f)(a2, a3, a4, a5, r);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class T>
void apply(void (T::*setter)(B1, B2, B3, B4, B5), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, T &r) {
    (r.*setter)(a1, a2, a3, a4, a5);
}

template <class T, class F, class A1, class A2, class A3, class A4, class A5>
struct MappedVar<T, F, A1, A2, A3, A4, A5>: Var<T> {

private:

    F _f;
    P<Var<A1> > _a1; P<Var<A2> > _a2; P<Var<A3> > _a3; P<Var<A4> > _a4; P<Var<A5> > _a5;

    void computeValue(T &value) const {
        apply(_f, _a1->value(), _a2->value(), _a3->value(), _a4->value(), _a5->value(), value);
    }

public:

    MappedVar(F f, P<Var<A1> > a1, P<Var<A2> > a2, P<Var<A3> > a3, P<Var<A4> > a4, P<Var<A5> > a5):
        _f(f), _a1(a1), _a2(a2), _a3(a3), _a4(a4), _a5(a5)
    {
        watch(a1); watch(a2); watch(a3); watch(a4); watch(a5);
    }
};

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5>
P<MappedVar<T, void (*)(B1, B2, B3, B4, B5, T &), V(A1), V(A2), V(A3), V(A4), V(A5)> > mapped(
        void (*f)(B1, B2, B3, B4, B5, T &), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5) {
    return P<MappedVar<T, void (*)(B1, B2, B3, B4, B5, T &), V(A1), V(A2), V(A3), V(A4), V(A5)> >(new MappedVar<T, void (*)(B1, B2, B3, B4, B5, T &), V(A1), V(A2), V(A3), V(A4), V(A5)>(f, a1, a2, a3, a4, a5));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5>
P<MappedVar<R(T), T (*)(B1, B2, B3, B4, B5), V(A1), V(A2), V(A3), V(A4), V(A5)> > mapped(
        T (*f)(B1, B2, B3, B4, B5), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5) {
    return P<MappedVar<R(T), T (*)(B1, B2, B3, B4, B5), V(A1), V(A2), V(A3), V(A4), V(A5)> >(new MappedVar<R(T), T (*)(B1, B2, B3, B4, B5), V(A1), V(A2), V(A3), V(A4), V(A5)>(f, a1, a2, a3, a4, a5));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5>
P<MappedVar<R(T), T (B1::*)(B2, B3, B4, B5) const, V(A1), V(A2), V(A3), V(A4), V(A5)> > mapped(
        T (B1::*f)(B2, B3, B4, B5) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5) {
    return P<MappedVar<R(T), T (B1::*)(B2, B3, B4, B5) const, V(A1), V(A2), V(A3), V(A4), V(A5)> >(
        new MappedVar<R(T), T (B1::*)(B2, B3, B4, B5) const, V(A1), V(A2), V(A3), V(A4), V(A5)>(f, a1, a2, a3, a4, a5));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5>
P<MappedVar<R(T), void (B1::*)(B2, B3, B4, B5, T &) const, V(A1), V(A2), V(A3), V(A4), V(A5)> > mapped(
        void (B1::*f)(B2, B3, B4, B5, T &) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5) {
    return P<MappedVar<R(T), void (B1::*)(B2, B3, B4, B5, T &) const, V(A1), V(A2), V(A3), V(A4), V(A5)> >(
        new MappedVar<R(T), void (B1::*)(B2, B3, B4, B5, T &) const, V(A1), V(A2), V(A3), V(A4), V(A5)>(f, a1, a2, a3, a4, a5));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5>
P<MappedVar<T, void (T::*)(B1, B2, B3, B4, B5), V(A1), V(A2), V(A3), V(A4), V(A5)> > setWith(
        void (T::*setter)(B1, B2, B3, B4, B5), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5) {
    return P<MappedVar<T, void (T::*)(B1, B2, B3, B4, B5), V(A1), V(A2), V(A3), V(A4), V(A5)> >(
        new MappedVar<T, void (T::*)(B1, B2, B3, B4, B5), V(A1), V(A2), V(A3), V(A4), V(A5)>(setter, a1, a2, a3, a4, a5));
}

// 
// ------------
//  6 argument
// ------------
// 

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class T>
void apply(void (*f)(B1, B2, B3, B4, B5, B6, T&), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, T &r) {
    f(a1, a2, a3, a4, a5, a6, r);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class T, class U>
void apply(U (*f)(B1, B2, B3, B4, B5, B6), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, T &r) {
    r = f(a1, a2, a3, a4, a5, a6);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class T, class U>
void apply(U (B1::*f)(B2, B3, B4, B5, B6) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, T &r) {
    r = (ptrOrRef<A1>::ref(a1).*f)(a2, a3, a4, a5, a6);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class T>
void apply(void (B1::*f)(B2, B3, B4, B5, B6, T &) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, T &r) {
    (ptrOrRef<A1>::ref(a1).*f)(a2, a3, a4, a5, a6, r);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class T>
void apply(void (T::*setter)(B1, B2, B3, B4, B5, B6), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, T &r) {
    (r.*setter)(a1, a2, a3, a4, a5, a6);
}

template <class T, class F, class A1, class A2, class A3, class A4, class A5, class A6>
struct MappedVar<T, F, A1, A2, A3, A4, A5, A6>: Var<T> {

private:

    F _f;
    P<Var<A1> > _a1; P<Var<A2> > _a2; P<Var<A3> > _a3; P<Var<A4> > _a4; P<Var<A5> > _a5; P<Var<A6> > _a6;

    void computeValue(T &value) const {
        apply(_f, _a1->value(), _a2->value(), _a3->value(), _a4->value(), _a5->value(), _a6->value(), value);
    }

public:

    MappedVar(F f, P<Var<A1> > a1, P<Var<A2> > a2, P<Var<A3> > a3, P<Var<A4> > a4, P<Var<A5> > a5, P<Var<A6> > a6):
        _f(f), _a1(a1), _a2(a2), _a3(a3), _a4(a4), _a5(a5), _a6(a6)
    {
        watch(a1); watch(a2); watch(a3); watch(a4); watch(a5); watch(a6);
    }
};

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6>
P<MappedVar<T, void (*)(B1, B2, B3, B4, B5, B6, T &), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)> > mapped(
        void (*f)(B1, B2, B3, B4, B5, B6, T &), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6) {
    return P<MappedVar<T, void (*)(B1, B2, B3, B4, B5, B6, T &), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)> >(new MappedVar<T, void (*)(B1, B2, B3, B4, B5, B6, T &), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)>(f, a1, a2, a3, a4, a5, a6));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6>
P<MappedVar<R(T), T (*)(B1, B2, B3, B4, B5, B6), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)> > mapped(
        T (*f)(B1, B2, B3, B4, B5, B6), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6) {
    return P<MappedVar<R(T), T (*)(B1, B2, B3, B4, B5, B6), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)> >(new MappedVar<R(T), T (*)(B1, B2, B3, B4, B5, B6), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)>(f, a1, a2, a3, a4, a5, a6));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6>
P<MappedVar<R(T), T (B1::*)(B2, B3, B4, B5, B6) const, V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)> > mapped(
        T (B1::*f)(B2, B3, B4, B5, B6) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6) {
    return P<MappedVar<R(T), T (B1::*)(B2, B3, B4, B5, B6) const, V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)> >(
        new MappedVar<R(T), T (B1::*)(B2, B3, B4, B5, B6) const, V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)>(f, a1, a2, a3, a4, a5, a6));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6>
P<MappedVar<R(T), void (B1::*)(B2, B3, B4, B5, B6, T &) const, V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)> > mapped(
        void (B1::*f)(B2, B3, B4, B5, B6, T &) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6) {
    return P<MappedVar<R(T), void (B1::*)(B2, B3, B4, B5, B6, T &) const, V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)> >(
        new MappedVar<R(T), void (B1::*)(B2, B3, B4, B5, B6, T &) const, V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)>(f, a1, a2, a3, a4, a5, a6));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6>
P<MappedVar<T, void (T::*)(B1, B2, B3, B4, B5, B6), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)> > setWith(
        void (T::*setter)(B1, B2, B3, B4, B5, B6), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6) {
    return P<MappedVar<T, void (T::*)(B1, B2, B3, B4, B5, B6), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)> >(
        new MappedVar<T, void (T::*)(B1, B2, B3, B4, B5, B6), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6)>(setter, a1, a2, a3, a4, a5, a6));
}

// 
// ------------
//  7 argument
// ------------
// 

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class A7, class B7, class T>
void apply(void (*f)(B1, B2, B3, B4, B5, B6, B7, T&), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, const A7 &a7, T &r) {
    f(a1, a2, a3, a4, a5, a6, a7, r);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class A7, class B7, class T, class U>
void apply(U (*f)(B1, B2, B3, B4, B5, B6, B7), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, const A7 &a7, T &r) {
    r = f(a1, a2, a3, a4, a5, a6, a7);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class A7, class B7, class T, class U>
void apply(U (B1::*f)(B2, B3, B4, B5, B6, B7) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, const A7 &a7, T &r) {
    r = (ptrOrRef<A1>::ref(a1).*f)(a2, a3, a4, a5, a6, a7);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class A7, class B7, class T>
void apply(void (B1::*f)(B2, B3, B4, B5, B6, B7, T &) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, const A7 &a7, T &r) {
    (ptrOrRef<A1>::ref(a1).*f)(a2, a3, a4, a5, a6, a7, r);
}

template <class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class A7, class B7, class T>
void apply(void (T::*setter)(B1, B2, B3, B4, B5, B6, B7), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, const A7 &a7, T &r) {
    (r.*setter)(a1, a2, a3, a4, a5, a6, a7);
}

template <class T, class F, class A1, class A2, class A3, class A4, class A5, class A6, class A7>
struct MappedVar<T, F, A1, A2, A3, A4, A5, A6, A7>: Var<T> {

private:

    F _f;
    P<Var<A1> > _a1; P<Var<A2> > _a2; P<Var<A3> > _a3; P<Var<A4> > _a4; P<Var<A5> > _a5; P<Var<A6> > _a6; P<Var<A7> > _a7;

    void computeValue(T &value) const {
        apply(_f, _a1->value(), _a2->value(), _a3->value(), _a4->value(), _a5->value(), _a6->value(), _a7->value(), value);
    }

public:

    MappedVar(F f, P<Var<A1> > a1, P<Var<A2> > a2, P<Var<A3> > a3, P<Var<A4> > a4, P<Var<A5> > a5, P<Var<A6> > a6, P<Var<A7> > a7):
        _f(f), _a1(a1), _a2(a2), _a3(a3), _a4(a4), _a5(a5), _a6(a6), _a7(a7)
    {
        watch(a1); watch(a2); watch(a3); watch(a4); watch(a5); watch(a6); watch(a7);
    }
};

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class A7, class B7>
P<MappedVar<T, void (*)(B1, B2, B3, B4, B5, B6, B7, T &), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)> > mapped(
        void (*f)(B1, B2, B3, B4, B5, B6, B7, T &), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, const A7 &a7) {
    return P<MappedVar<T, void (*)(B1, B2, B3, B4, B5, B6, B7, T &), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)> >(new MappedVar<T, void (*)(B1, B2, B3, B4, B5, B6, B7, T &), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)>(f, a1, a2, a3, a4, a5, a6, a7));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class A7, class B7>
P<MappedVar<R(T), T (*)(B1, B2, B3, B4, B5, B6, B7), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)> > mapped(
        T (*f)(B1, B2, B3, B4, B5, B6, B7), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, const A7 &a7) {
    return P<MappedVar<R(T), T (*)(B1, B2, B3, B4, B5, B6, B7), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)> >(new MappedVar<R(T), T (*)(B1, B2, B3, B4, B5, B6, B7), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)>(f, a1, a2, a3, a4, a5, a6, a7));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class A7, class B7>
P<MappedVar<R(T), T (B1::*)(B2, B3, B4, B5, B6, B7) const, V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)> > mapped(
        T (B1::*f)(B2, B3, B4, B5, B6, B7) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, const A7 &a7) {
    return P<MappedVar<R(T), T (B1::*)(B2, B3, B4, B5, B6, B7) const, V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)> >(
        new MappedVar<R(T), T (B1::*)(B2, B3, B4, B5, B6, B7) const, V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)>(f, a1, a2, a3, a4, a5, a6, a7));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class A7, class B7>
P<MappedVar<R(T), void (B1::*)(B2, B3, B4, B5, B6, B7, T &) const, V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)> > mapped(
        void (B1::*f)(B2, B3, B4, B5, B6, B7, T &) const, const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, const A7 &a7) {
    return P<MappedVar<R(T), void (B1::*)(B2, B3, B4, B5, B6, B7, T &) const, V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)> >(
        new MappedVar<R(T), void (B1::*)(B2, B3, B4, B5, B6, B7, T &) const, V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)>(f, a1, a2, a3, a4, a5, a6, a7));
}

template <class T, class A1, class B1, class A2, class B2, class A3, class B3, class A4, class B4, class A5, class B5, class A6, class B6, class A7, class B7>
P<MappedVar<T, void (T::*)(B1, B2, B3, B4, B5, B6, B7), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)> > setWith(
        void (T::*setter)(B1, B2, B3, B4, B5, B6, B7), const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4, const A5 &a5, const A6 &a6, const A7 &a7) {
    return P<MappedVar<T, void (T::*)(B1, B2, B3, B4, B5, B6, B7), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)> >(
        new MappedVar<T, void (T::*)(B1, B2, B3, B4, B5, B6, B7), V(A1), V(A2), V(A3), V(A4), V(A5), V(A6), V(A7)>(setter, a1, a2, a3, a4, a5, a6, a7));
}

// --- end of generated code

#undef R

template <class T, class A, class B>
P<MappedVar<T, T B::* const, V(A)> > mapped_member(T B::* const m, const A &a) {
    return P<MappedVar<T, T B::* const, V(A)> > (
        new MappedVar<T, T B::* const, V(A)>(m, a));
}

#undef V

// 
// ===========
//  mappedGet
// ===========
// 

template <class X, class C>
const X &constGet(const array<X, C> &xs, int i) {
    if (i < 0) i = xs.size() + i;
    ASSERT(0 <= i && i < xs.size());
    return xs[i];
}

template <class X, class C>
P<Var<X> > mappedGet(const P<Var<array<X, C> > > &xs,
                     const P<Var<int> > &i) {
    return mapped(&constGet<X, C>, xs, i);
}

template <class X, class C>
P<Var<X> > mappedGet(const P<Var<array<X, C> > > &xs,
                     int i) {
    return mapped(&constGet<X, C>, xs, root(i));
}

// 
// ============
//  GroupedVar
// ============
// 

template <class T>
struct GroupedVar: Var<vector<T> > {

private:

    vector<P<Var<T> > > _xs;

public:

    GroupedVar(const vector<P<Var<T> > > &xs):
        _xs(xs)
    {
        for (size_t a = 0; a < xs.size(); ++a) watch(xs[a]);
    }

    void computeValue(vector<T> &value) const {
        value.resize(_xs.size());
        for (size_t a = 0; a < _xs.size(); ++a)
            value[a] = _xs[a]->value();
    }
};

template <class T>
P<Var<vector<T> > > grouped(const vector<P<Var<T> > > &xs) {
    return P<Var<vector<T> > >(new GroupedVar<T>(xs));
}

#undef P

}

DRLIB_END_NAMESPACE

#endif
