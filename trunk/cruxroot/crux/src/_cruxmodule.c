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

void
init_crux(void)
{
    /* Initialize the exception handling machinery. */
    CxXepInit();

    /* Do processor-specific initialization, if necessary. */
#ifdef CxmCpuInit
    CxmCpuInit();
#endif

    Py_InitModule3("_crux", cruxFuncs, "crux extensions");
    CxTreeInit();
    CxNodeInit();
    CxEdgeInit();
    CxRingInit();
}
