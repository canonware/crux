/******************************************************************************
 *
 * <Copyright = jasone>
 * <License>
 *
 ******************************************************************************
 *
 * Version: Crux <Version = crux>
 *
 ******************************************************************************/

#include "../include/_cruxmodule.h"

static PyMethodDef cruxFuncs[] =
{
    {NULL}
};

PyObject *CxgException;

void
init_crux(void)
{
    PyObject *m;

    /* Initialize the exception handling machinery. */
    CxXepInit();

    /* Do processor-specific initialization, if necessary. */
#ifdef CxmCpuInit
    CxmCpuInit();
#endif

    m = Py_InitModule3("_crux", cruxFuncs, "crux extensions");

    CxgException = PyErr_NewException("_crux.Exception", PyExc_Exception, NULL);
    Py_INCREF(CxgException);
    PyModule_AddObject(m, "Exception", CxgException);

    CxTreeInit();
    CxNodeInit();
    CxEdgeInit();
    CxRingInit();
}
