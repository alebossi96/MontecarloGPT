#include "python3.8/Python.h"
#include <numpy/arrayobject.h>
#include "montecarlo.hpp"
#include <vector>
#include<cstdlib>
#include <string>
#include<iostream>

template<typename M>
M fill_el(PyObject * pyObj, std::string name, bool *is_set)
    {
    PyObject *py_mu;
    if(PyObject_HasAttrString(pyObj, name.c_str()))
        {
        py_mu = PyObject_GetAttrString(pyObj, name.c_str());
        M to_fill;
        if (std::is_same<M, double>::value)
            PyArg_Parse(py_mu, "d", &to_fill);
        if (std::is_same<M, bool>::value)
            PyArg_Parse(py_mu, "i", &to_fill);
        *is_set = true;
        return to_fill;
        }
    *is_set = false;
    return -1;
    }
Detector set_detector(PyObject * pyObj)
    {
    bool is_set;
    double x,y,z, r;
    x = fill_el<double>(pyObj, "x",&is_set);
    y = fill_el<double>(pyObj, "y",&is_set);
    z = fill_el<double>(pyObj, "z",&is_set);
    r = fill_el<double>(pyObj, "r",&is_set);
    std::cout <<x<<" "<<y<<" "<<z<<std::endl;
    Detector detector(Vector(x,y,z),r);
    return detector;
    }
extern "C"
    {
    PyObject *mc(PyObject *self, PyObject *args, PyObject *keywds)
        {
        Py_Initialize();//cosa fanno?
        _import_array();
        double g,mu_s;
        PyObject *obj_det;
        static char *kwlist[] = { "g", "mu_s", "detector", NULL};
        PyArg_ParseTupleAndKeywords(args, keywds, "ddO", kwlist, &g, &mu_s, &obj_det);
        std::cout<<"g = "<<g<<" mu_s = "<< mu_s<<std::endl;
        Detector detector =  set_detector(obj_det);
        std::vector<int> dtof = simulate(g,mu_s, detector);
        int * results;
        std::size_t length = dtof.size();
        results =(int *)malloc(length * sizeof(int));
        npy_intp dims[1];
        dims[0] = length;
        for(std::size_t i = 0; i< length; ++i)
            {
            results[i] = dtof[i];
            }
        PyObject* obj = PyArray_SimpleNewFromData(1, dims, NPY_INT, (void*)results);
        std::cout<<g <<"  "<<mu_s<<std::endl;
        return obj;
        }
    PyObject *angle_distribution(PyObject *self, PyObject *args, PyObject *keywds)
        {
        Py_Initialize();//cosa fanno?
        _import_array();
        double g;
        static char *kwlist[] = { "g", NULL};
        PyArg_ParseTupleAndKeywords(args, keywds, "d", kwlist, &g);
        std::cout<<"g = "<<g<<std::endl;
        
        std::array<double, SIZE_LIST_ANGLE> angle_distribution = inverse_transform_sampling(henyey_greenstein_F, g);
        double * results;
        std::size_t length = angle_distribution.size();
        results =(double *)malloc(length * sizeof(double));
        npy_intp dims[1];
        dims[0] = length;
        for(std::size_t i = 0; i< length; ++i)
            {
            results[i] = angle_distribution[i];
            }
        PyObject* obj = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (void*)results);
        return obj;
        }
    static PyMethodDef PyFunctionMethods[] = 
        {
        {"mc", (PyCFunction)(void(*)(void)) mc, METH_VARARGS| METH_KEYWORDS, "docstring"},
        {"angle_distribution", (PyCFunction)(void(*)(void)) angle_distribution, METH_VARARGS| METH_KEYWORDS, "docstring"},
        {NULL, NULL}
        };
    static struct PyModuleDef functionmodule =
        {
        PyModuleDef_HEAD_INIT, "mc", "montecarlo", -1, PyFunctionMethods
        };

    PyMODINIT_FUNC PyInit_montecarlomodule(void)
        {
        return PyModule_Create(&functionmodule);
        }

    }
