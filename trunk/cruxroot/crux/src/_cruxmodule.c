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
