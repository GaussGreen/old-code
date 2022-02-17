/**
 * @file FieldPath.cpp
 */

#include "edginc/config.hpp"
#include "edginc/FieldPath.hpp"
#include "edginc/MarketObject.hpp"

DRLIB_BEGIN_NAMESPACE

static CFieldConstSP fieldOn(IObjectConstSP object, const string& name) {
    if (!object) {
        throw ModelException("Tried to find field \"" + name + "\" on "
                             "a null object");
    }

    for (CClassConstSP c = object->getClass(); !!c;
             c = c->getSuperClass()) {
        CFieldConstSP f = c->hasDeclaredField(name);
        if (f) return f;
    }

    throw ModelException("There is no field \"" + name + "\" on the "
                         "class \"" + object->getClass()->getName() + "\"");
}

FieldPath::FieldPath():
    CObject(TYPE)
{}

FieldPath::FieldPath(CStringArrayConstSP path):
    CObject(TYPE),
    _path(path)
{
    if (!path) {
        throw ModelException("FieldPath::FieldPath()", "'path' is NULL");
    }
}

FieldPathSP FieldPath::SP(CStringArrayConstSP path) {
    return FieldPathSP(new FieldPath(path));
}

FieldPathSP FieldPath::SP(const string& name) {
    return FieldPathSP(new FieldPath(CStringArray::SP(1, name)));
}

FieldPathSP FieldPath::SP() {
    return FieldPathSP(new FieldPath(CStringArray::SP()));
}

FieldPath::~FieldPath() {}

CStringArrayConstSP FieldPath::path() const {
    return _path;
}

string FieldPath::toString(int highlight) const {
    string it;
    for (int g = 0; g < _path->size(); ++g) {
        it.append(".");
        if (g == highlight) it.append(">>");
        it.append((*_path)[g]);
        if (g == highlight) it.append("<<");
    }
    return it;
}

string FieldPath::contextMessage(IObjectConstSP root, int highlight) const {
    const MarketObject* m = dynamic_cast<const MarketObject*>(root.get());
    return "field " + toString(highlight) + " on " +
            (!root ? "NULL object" :
                     "object " +
                         (m ? "\"" + m->getName() + "\" " : string()) +
                         "of type " + root->getClass()->getName());
}

IObjectSP FieldPath::get(IObjectSP root, CClassConstSP type) const {
    try {
        IObjectSP obj = _path->empty() ?
                            root : objectChain(root).objects->back();

        if (type) {
            if (!obj) throw ModelException("Field's value is null");
            if (!type->isInstance(obj)) {
                throw ModelException(
                    "Field's value is of type " + obj->getClass()->getName() +
                    " but was expecting " + type->getName());
            }
        }

        return obj;
    }
    catch (exception& e) {
        throw ModelException(e, "FieldPath::get()", "Looking up field value");
    }
}

IObjectConstSP FieldPath::get(IObjectConstSP root, CClassConstSP type) const {
    return get(IObjectSP::constCast(root), type);
}

void FieldPath::ObjectChain::fieldsUpdated() const {
    for (int o = fields.size() - 1; o >= 0; --o) {
        (*objects)[o]->fieldsUpdated(CFieldArray(1, fields[o]));
    }
}

FieldPath::ObjectChain::ObjectChain() {}
FieldPath::ObjectChain::~ObjectChain() {}

FieldPath::ObjectChain FieldPath::objectChain(IObjectSP root) const {
    int f = 0;
    try {
        if (_path->empty()) {
            throw ModelException("Can't obtain meaningful object chain since "
                                 "field path is empty");
        }

        ObjectChain c;

        c.fields.resize(_path->size());
        c.objects.reset(new ObjectArray(_path->size() + 1));

        (*c.objects)[0] = root;

        for (; f < _path->size(); ++f) {
            if (!(*c.objects)[f]) {
                --f;
                throw ModelException("Null pointer in field chain");
            }

            c.fields[f] = fieldOn((*c.objects)[f], (*_path)[f]);
            (*c.objects)[f + 1] = c.fields[f]->get((*c.objects)[f]);

            MarketObjectWrapper* w =
                dynamic_cast<MarketObjectWrapper*>((*c.objects)[f + 1].get());
            if (w) {
                (*c.objects)[f + 1] = w->getMO();
                if (!(*c.objects)[f + 1]) {
                    throw ModelException(
                        "Wrapper value for " + w->getMOType()->getName() +
                        " value is null");
                }
            }
        }
        
        return c;
    }
    catch (exception& e) {
        throw ModelException(e, "FieldPath::objectChain()",
                             "Locating " + contextMessage(root, f));
    }
}

void FieldPath::set(IObjectSP root, IObjectSP value) const {
    MarketObjectWrapper* w = 0;

    try {
        ObjectChain c = objectChain(root);

        w = dynamic_cast<MarketObjectWrapper*>(c.objects->back().get());

        if (w) {
            if (!value) {
                throw ModelException("Attempt to assign NULL value to wrapper");
            }

            if (!MarketObject::TYPE->isInstance(value)) {
                throw ModelException(
                    "Attempt to assign a " + value->getClass()->getName() +
                    ", which is not a MarketObject, to a wrapper");
            }

            w->setObject(MarketObjectSP::dynamicCast(value));
        }
        else {
            c.fields.back()->set((*c.objects)[c.fields.size() - 1], value);
        }
    }
    catch (exception& e) {
        throw ModelException(e, "FieldPath::set",
            string("Assigning to") + (w ? " wrapper" : "") + " field " +
            contextMessage(root, _path->size() - 1));
    }
}

IObject* FieldPath::defaultOne() {
    return new FieldPath();
}

void FieldPath::load(CClassSP& clazz) {
    REGISTER(FieldPath, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultOne);
    FIELD(_path, "path");
}

CClassConstSP const FieldPath::TYPE = CClass::registerClassLoadMethod(
    "FieldPath", typeid(FieldPath), load);

DRLIB_END_NAMESPACE
