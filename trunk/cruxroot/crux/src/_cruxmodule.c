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

static PyObject *
hello(PyObject *self)
{
    PyObject *retval;

    xep_begin();
    xep_try
    {
	retval = Py_BuildValue("s", "hello");
    }
    xep_catch(CW_CRUXX_OOM)
    {
	xep_handled();
	retval = PyErr_NoMemory();
    }
    xep_end();

    return retval;
}

static PyMethodDef crux_funcs[] =
{
    {NULL}
};

void
init_crux(void)
{
    /* Initialize the exception handling machinery. */
    xep_init();

    /* Do processor-specific initialization, if necessary. */
#ifdef CRUX_CPU_INIT
    CRUX_CPU_INIT();
#endif

    Py_InitModule3("_crux", crux_funcs, "crux extensions");
    crux_tree_init();
    crux_node_init();
    crux_edge_init();
}
