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
#ifndef HAVE_ASPRINTF
#include "asprintf_c"
#endif

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
    CxFastaParserInit();
    CxDistMatrixInit();
}

void
CxError(PyObject *exception, const char *format, ...)
{
    va_list ap;
    char *str;

    va_start(ap, format);
    vasprintf(&str, format, ap);
    va_end(ap);

    PyErr_SetString(exception, str);
    free(str);
}
