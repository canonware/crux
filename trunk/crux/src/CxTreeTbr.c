//==============================================================================
//
// <Copyright = jasone>
// <License>
//
//==============================================================================
//
// Version: Crux <Version = crux>
//
//==============================================================================

#include "../include/_cruxmodule.h"

bool
CxTreeTbr(CxtTreeObject *self, CxtEdgeObject *aBisect,
	  CxtEdgeObject *aReconnectA, CxtEdgeObject *aReconnectB)
{
    bool rVal;

//     CxTrTbr(self->tr, aBisect, aReconnectA, aReconnectB);

    rVal = false;
    //RETURN:
    return rVal;
}

PyObject *
CxTreeTbrPargs(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtEdgeObject *bisect, *reconnectA, *reconnectB;

    if (PyArg_ParseTuple(args, "(O!O!O!)",
			 &CxtEdge, &bisect,
			 &CxtEdge, &reconnectA,
			 &CxtEdge, &reconnectB)
	== 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    if (CxTreeTbr(self, bisect, reconnectA, reconnectB))
    {
	rVal = PyErr_NoMemory();
	goto RETURN;
    }

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
}

bool
CxTreeTbrNneighborsGet(CxtTreeObject *self, unsigned *rNneighbors)
{
    bool rVal;

    *rNneighbors = CxTrTbrNneighborsGet(self->tr);

    rVal = false;
    //RETURN:
    return rVal;
}

PyObject *
CxTreeTbrNneighborsGetPargs(CxtTreeObject *self)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    unsigned nneighbors;

    if (CxTreeTbrNneighborsGet(self, &nneighbors))
    {
	rVal = PyErr_NoMemory();
	goto RETURN;
    }

    rVal = Py_BuildValue("I", nneighbors);
    RETURN:
    return rVal;
}

bool
CxTreeTbrNeighborGet(CxtTreeObject *self, unsigned aNeighbor,
		     CxtEdgeObject **rBisect, CxtEdgeObject **rReconnectA,
		     CxtEdgeObject **rReconnectB)
{
    bool rVal;

//     CxTrTbrNeighborGet(self->tr, neighbor,
// 		       &bisect, &reconnectA, &reconnectB);

    rVal = false;
    //RETURN:
    return rVal;
}

PyObject *
CxTreeTbrNeighborGetPargs(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    unsigned nneighbors, neighbor;
    CxtEdgeObject *bisect, *reconnectA, *reconnectB;

    if (PyArg_ParseTuple(args, "I", &neighbor) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    if (CxTreeTbrNneighborsGet(self, &nneighbors))
    {
	rVal = PyErr_NoMemory();
	goto RETURN;
    }

    if (neighbor >= nneighbors)
    {
	CxError(CxgTreeValueError,
		"neighbor: %u is out of range [0..%u]",
		neighbor, nneighbors);
	rVal = NULL;
	goto RETURN;
    }

    if (CxTreeTbrNeighborGet(self, neighbor, &bisect, &reconnectA, &reconnectB))
    {
	rVal = PyErr_NoMemory();
	goto RETURN;
    }

    rVal = Py_BuildValue("(O!O!O!)",
			 &CxtEdge, bisect,
			 &CxtEdge, reconnectA,
			 &CxtEdge, reconnectB);
    RETURN:
    return rVal;
}
