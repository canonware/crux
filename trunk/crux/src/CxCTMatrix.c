//==============================================================================
//
// <Copyright = jasone>
// <License>
//
//==============================================================================
//
// Version: Crux <Version = Crux>
//
//==============================================================================

#include "../include/_cruxmodule.h"

PyObject *
CxCTMatrixXXXFoo(CxtFastaParserObject *self, PyObject *args)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");

    return rVal;
}

PyObject *
CxCTMatrixXXX(PyObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");

    return rVal;
}

static PyMethodDef CxpCTMatrixMethods[] =
{
    {
	"XXXA",
	(PyCFunction) CxCTMatrixXXXFoo,
	METH_VARARGS,
	"XXXA"
    },
    {
	"XXXB",
	(PyCFunction) CxCTMatrixXXX,
	METH_NOARGS,
	"XXXB"
    },
    {NULL, NULL}
};

void
CxCTMatrixInit(void)
{
    PyMethodDef *def;
    PyObject *module, *moduleDict, *classDict, *className, *cTMatrixClass;
    PyObject *func, *method;

    module = Py_InitModule("CTMatrix", CxpCTMatrixMethods);
    moduleDict = PyModule_GetDict(module);
    classDict = PyDict_New();
    className = PyString_FromString("CTMatrix");
    cTMatrixClass = PyClass_New(NULL, classDict, className);
    PyDict_SetItemString(moduleDict, "CTMatrix", cTMatrixClass);
    Py_DECREF(classDict);
    Py_DECREF(className);
    Py_DECREF(cTMatrixClass);

    for (def = CxpCTMatrixMethods; def->ml_name != NULL; def++)
    {
	func = PyCFunction_New(def, NULL);
	method = PyMethod_New(func, NULL, cTMatrixClass);
	PyDict_SetItemString(classDict, def->ml_name, method);
	Py_DECREF(func);
	Py_DECREF(method);
    }
}
