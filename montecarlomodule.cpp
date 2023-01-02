#include "python3.8/Python.h"
#include <numpy/arrayobject.h>
#include "montecarlo.hpp"
#include <vector>
#include<cstdlib>
#include <string>
#include<iostream>
extern "C"
    {
    PyObject *mc(PyObject *self, PyObject *args, PyObject *keywds)
        {
        Py_Initialize();//cosa fanno?
        _import_array();
        double g,mu_s;
        static char *kwlist[] = { "g", "mu_s", NULL};
        PyArg_ParseTupleAndKeywords(args, keywds, "dd", kwlist, &g, &mu_s);
        std::vector<int> dtof = simulate(g,mu_s);
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
    
    static PyMethodDef PyFunctionMethods[] = 
        {
        {"mc", (PyCFunction)(void(*)(void)) mc, METH_VARARGS| METH_KEYWORDS, "docstring"},
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
