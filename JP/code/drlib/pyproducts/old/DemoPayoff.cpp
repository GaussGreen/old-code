#include <boost/python.hpp>
#include <iostream>

using namespace boost::python;

class QVal {
public:
    QVal(int v=0) {
        value = v;
    }
  int getValue() {
    return value;
  }
    static QVal max(const QVal& v1,const QVal& v2) {
        QVal r;
        if (v1.value > v2.value) r.value = v1.value;
        else r.value = v2.value;
        return r;
    }
    friend QVal operator-(const QVal& v1,const QVal& v2);
    friend QVal operator-(const QVal& v1,int v2);
    friend QVal operator-(int v1,const QVal& v2);
    friend bool operator>(const QVal& v1,int v2);
    friend bool operator>(int v1,const QVal& v2);
private:
    int value;
};

QVal operator-(const QVal& v1,const QVal& v2) {
    QVal r;
    r.value = v1.value - v2.value;
    return r;
}

QVal operator-(const QVal& v1,int v2) {
    QVal r;
    r.value = v1.value - v2;
    return r;
}

QVal operator-(int v1,const QVal& v2) {
    QVal r;
    r.value = v1 - v2.value;
    return r;
}

bool operator>(const QVal& v1,int v2) {
    return (v1.value > v2) ? 1 : 0;
}

bool operator>(int v1,const QVal& v2) {
    return (v1 > v2.value) ? 1 : 0;
}

// Export qlib module with QVal class

BOOST_PYTHON_MODULE(qlib) {
    class_<QVal>("QVal")
        .def(init<int>())
        .def("getValue",&QVal::getValue)
        .def("max",&QVal::max)
        .def(self - self)
        .def(self - int())
        .def(int() - self)
        .def(self > int())
        .def(int() > self);
}

int main(int argc,char** argv) {
  try {
        // Initialize python interpreter
    Py_Initialize();
        // Initialize qlib module defined in BOOST_PYTHON_MODULE(qlib)
        initqlib();
        // Read demo_trade.py
    object module(handle<>(PyImport_ImportModule("demo_trade")));
        // Get payoff function handle
        object payoff = module.attr("payoff");
        // Create payoff function argument instances
        QVal qv1(30);
        QVal qv2(20);
        // Call payoff function with object instances
        // Return another QVal object instance
    QVal qv3 = extract<QVal>(payoff(qv1,qv2));    
        // Print value of return object instance
    std::cout << qv3.getValue() << std::endl;
  } catch (error_already_set) {
    PyErr_Print();
  }
}
