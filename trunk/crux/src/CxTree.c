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

#define CxmTree_c
#include "../include/_cruxmodule.h"

#include <math.h>

//#define CxmTreeGCVerbose
#ifdef CxmTreeGCVerbose
#undef Py_INCREF
#define Py_INCREF(op)							\
	fprintf(stderr, "%s:%d:%s(): INCREF(%p) --> %d\n",		\
		__FILE__, __LINE__, __func__, op,			\
		((op)->ob_refcnt) + 1);					\
	(_Py_INC_REFTOTAL  _Py_REF_DEBUG_COMMA				\
	(op)->ob_refcnt++)

#undef Py_DECREF
#define Py_DECREF(op)							\
	fprintf(stderr, "%s:%d:%s(): DECREF(%p) --> %d\n",		\
		__FILE__, __LINE__, __func__, op,			\
		((op)->ob_refcnt) - 1);					\
	if (_Py_DEC_REFTOTAL  _Py_REF_DEBUG_COMMA			\
	   --(op)->ob_refcnt != 0)					\
		_Py_CHECK_REFCNT(op)					\
	else								\
		_Py_Dealloc((PyObject *)(op))
#endif

//==============================================================================
// Begin Tree.

static PyObject *
CxpTreeNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtTreeObject *self;

    self = (CxtTreeObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	self->tr = CxTrNew();
	CxTrAuxSet(self->tr, self);
	self->seq = 1; // 0 is skipped so that it can mean "uninitialized".

	self->treeAux = NULL;
	self->nTreeAux = 0;
	self->nodeAux = NULL;
	self->nNodeAux = 0;
	self->edgeAux = NULL;
	self->nEdgeAux = 0;
	self->ringAux = NULL;
	self->nRingAux = 0;

	self->GcCleared = false;
	rVal = (PyObject *) self;
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	rVal = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, rVal, rVal->ob_refcnt);
#endif
    return rVal;
}

static int
CxpTreeTraverse(CxtTreeObject *self, visitproc visit, void *arg)
{
    int rVal;
    CxtTrNode base;

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    if (self->GcCleared == false)
    {
	base = CxTrBaseGet(self->tr);
	if (base != CxmTrNodeNone)
	{
	    CxtNodeObject *node;

	    node = (CxtNodeObject *) CxTrNodeAuxGet(self->tr, base);
	    if (visit((PyObject *) node, arg) < 0)
	    {
		rVal = -1;
		goto RETURN;
	    }
	}
    }

    rVal = 0;
    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif
    return rVal;
}

static int
CxpTreeClear(CxtTreeObject *self)
{
    CxtTrNode base;

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    if (self->GcCleared == false)
    {
	base = CxTrBaseGet(self->tr);
	if (base != CxmTrNodeNone)
	{
	    CxtNodeObject *node;

	    node = (CxtNodeObject *) CxTrNodeAuxGet(self->tr, base);
	    CxTrBaseSet(self->tr, CxmTrNodeNone);
	    Py_DECREF(node);
	}

	self->GcCleared = true;
    }

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif
    return 0;
}

static void
CxpTreeDelete(CxtTreeObject *self)
{
    unsigned i;

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    CxpTreeClear(self);

    // Clean up aux.
    if (self->aux != NULL)
    {
	for (i = 0; i < self->nAux; i++)
	{
	    if (self->treeAux[i].cleanupTree != NULL)
	    {
		self->treeAux[i].cleanupTree(self, self->aux[i]);
	    }
	}

	free(self->aux);
    }

    // Clean up aux control data structures.
    for (i = 0; i < self->nRingAux; i++)
    {
	if (self->ringAux[i].cleanupFinal != NULL)
	{
	    self->ringAux[i].cleanupFinal(self, self->ringAux[i].data);
	}
    }
    if (self->ringAux != NULL)
    {
	free(self->ringAux);
    }

    for (i = 0; i < self->nEdgeAux; i++)
    {
	if (self->edgeAux[i].cleanupFinal != NULL)
	{
	    self->edgeAux[i].cleanupFinal(self, self->edgeAux[i].data);
	}
    }
    if (self->edgeAux != NULL)
    {
	free(self->edgeAux);
    }

    for (i = 0; i < self->nNodeAux; i++)
    {
	if (self->nodeAux[i].cleanupFinal != NULL)
	{
	    self->nodeAux[i].cleanupFinal(self, self->nodeAux[i].data);
	}
    }
    if (self->nodeAux != NULL)
    {
	free(self->nodeAux);
    }

    for (i = 0; i < self->nTreeAux; i++)
    {
	if (self->treeAux[i].cleanupFinal != NULL)
	{
	    self->treeAux[i].cleanupFinal(self, self->treeAux[i].data);
	}
    }
    if (self->treeAux != NULL)
    {
	free(self->treeAux);
    }

    CxTrDelete(self->tr);
    self->ob_type->tp_free((PyObject*) self);

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p\n",
 	    __FILE__, __LINE__, __func__, self);
#endif
}

static PyObject *CxpTreeNewCode;

CxtTreeObject *
CxTreeNew(void)
{
    CxtTreeObject *rVal;
    PyObject *globals, *locals, *obj;

    globals = PyEval_GetGlobals();
    if (globals == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }
    locals = Py_BuildValue("{}");
    if (locals == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }

    obj = PyEval_EvalCode((PyCodeObject *) CxpTreeNewCode,
			  globals,
			  locals);
    if (obj == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }
    Py_DECREF(obj);

    rVal = (CxtTreeObject *) PyDict_GetItemString(locals, "tree");
    if (rVal == NULL)
    {
	goto RETURN;
    }
    Py_INCREF(rVal);

    RETURN:
    Py_DECREF(locals);
    return rVal;
}

unsigned
CxTreeNtaxaGet(CxtTreeObject *self)
{
    return CxTrNtaxaGet(self->tr);
}

PyObject *
CxTreeNtaxaGetPargs(CxtTreeObject *self)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	rVal = Py_BuildValue("i", CxTrNtaxaGet(self->tr));
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	rVal = PyErr_NoMemory();
    }
    CxmXepEnd();

    return rVal;
}

PyObject *
CxTreeNedgesCget(CxtTreeObject *self)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	rVal = Py_BuildValue("i", CxTrNedgesGet(self->tr));
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	rVal = PyErr_NoMemory();
    }
    CxmXepEnd();

    return rVal;
}

CxtNodeObject *
CxTreeBaseGet(CxtTreeObject *self)
{
    CxtNodeObject *rVal;
    CxtTrNode base;

    base = CxTrBaseGet(self->tr);
    if (base == CxmTrNodeNone)
    {
	rVal = NULL;
    }
    else
    {
	rVal = (CxtNodeObject *) CxTrNodeAuxGet(self->tr, base);
    }

    return rVal;
}

PyObject *
CxTreeBaseGetPargs(CxtTreeObject *self)
{
    PyObject *rVal;
    CxtTrNode base;

    base = CxTrBaseGet(self->tr);
    if (base == CxmTrNodeNone)
    {
	Py_INCREF(Py_None);
	rVal = Py_None;
    }
    else
    {
	rVal = (PyObject *) CxTrNodeAuxGet(self->tr, base);
	Py_INCREF(rVal);
    }

    return rVal;
}

void
CxTreeBaseSet(CxtTreeObject *self, CxtNodeObject *aNode)
{
    CxtNodeObject *oldNode;
    CxtTrNode oldTrNode;

    // Decref if clobbering an already-set base (but wait until after the
    // new base is set so that the nodes/edges/rings are always reachable.
    oldTrNode = CxTrBaseGet(self->tr);
    if (oldTrNode != CxmTrNodeNone)
    {
	oldNode = (CxtNodeObject *) CxTrNodeAuxGet(self->tr, oldTrNode);
	CxTrBaseSet(self->tr, CxmTrNodeNone);
    }
    else
    {
	oldNode = NULL;
    }

    if (aNode != NULL)
    {
	Py_INCREF(aNode);
	CxTrBaseSet(self->tr, aNode->node);
	self->seq++;
    }

    if (oldNode != NULL)
    {
	Py_DECREF(oldNode);
    }
}

PyObject *
CxTreeBaseSetPargs(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal;
    CxtNodeObject *node;

    node = NULL;
    if (PyArg_ParseTuple(args, "|O!", &CxtNode, &node) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }
    if (node != NULL && node->tree != self)
    {
	CxError(CxgTreeValueError, "Node does not belong to this tree");
	rVal = NULL;
	goto RETURN;
    }

    CxTreeBaseSet(self, node);

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
}

static bool
CxpTreeIterateNodeCallbackDefault(CxtNodeObject *aNode,
				  CxtTreeIteratorStage aStage, void *aContext)
{
    return false;
}

static bool
CxpTreeIterateEdgeCallbackDefault(CxtEdgeObject *aEdge,
				  CxtTreeIteratorStage aStage, void *aContext)
{
    return false;
}

static bool
CxpTreeIterateRingCallbackDefault(CxtRingObject *aRing,
				  CxtTreeIteratorStage aStage, void *aContext)
{
    return false;
}

static bool
CxpTreeIterateRecurse(CxtRingObject *aRing,
		      CxtTreeIterateNodeCallback *aNodeCallback,
		      CxtTreeIterateEdgeCallback *aEdgeCallback,
		      CxtTreeIterateRingCallback *aRingCallback,
		      void *aContext)
{
    bool rVal;
    CxtEdgeObject *edge;
    CxtRingObject *ringOther, *ringCur;
    CxtNodeObject *node;

    edge = CxRingEdge(aRing);
    ringOther = CxRingOther(aRing);
    node = CxRingNode(ringOther);

    if (aRingCallback(aRing, CxTreeIteratorStagePre, aContext)
	|| aEdgeCallback(edge, CxTreeIteratorStagePre, aContext)
	|| aRingCallback(ringOther, CxTreeIteratorStagePre, aContext)
	|| aNodeCallback(node, CxTreeIteratorStagePre, aContext))
    {
	rVal = true;
	goto RETURN;
    }

    ringCur = CxRingNext(ringOther);
    if (ringCur != ringOther)
    {
	while (true)
	{
	    if (CxpTreeIterateRecurse(ringCur, aNodeCallback, aEdgeCallback,
				      aRingCallback, aContext))
	    {
		rVal = true;
		goto RETURN;
	    }

	    ringCur = CxRingNext(ringCur);
	    if (ringCur == ringOther)
	    {
		break;
	    }

	    if (aNodeCallback(node, CxTreeIteratorStageIn, aContext))
	    {
		rVal = true;
		goto RETURN;
	    }
	}
    }

    if (aNodeCallback(node, CxTreeIteratorStagePost, aContext)
	|| aRingCallback(ringOther, CxTreeIteratorStagePost, aContext)
	|| aEdgeCallback(edge, CxTreeIteratorStagePost, aContext)
	|| aRingCallback(aRing, CxTreeIteratorStagePost, aContext))
    {
	rVal = true;
	goto RETURN;
    }

    rVal = false;
    RETURN:
    return rVal;
}

bool
CxTreeIterate(CxtTreeObject *aTree,
	      CxtTreeIterateNodeCallback *aNodeCallback,
	      CxtTreeIterateEdgeCallback *aEdgeCallback,
	      CxtTreeIterateRingCallback *aRingCallback,
	      void *aContext)
{
    bool rVal;
    CxtNodeObject *base;
    CxtTreeIterateNodeCallback *nodeCallback;
    CxtTreeIterateEdgeCallback *edgeCallback;
    CxtTreeIterateRingCallback *ringCallback;

    // Set callback functions pointers.
    nodeCallback = (aNodeCallback != NULL)
	? aNodeCallback
	: CxpTreeIterateNodeCallbackDefault;
    edgeCallback = (aEdgeCallback != NULL)
	? aEdgeCallback
	: CxpTreeIterateEdgeCallbackDefault;
    ringCallback = (aRingCallback != NULL)
	? aRingCallback
	: CxpTreeIterateRingCallbackDefault;

    if ((base = CxTreeBaseGet(aTree)) != NULL)
    {
	CxtRingObject *ringStart;

	if (nodeCallback(base, CxTreeIteratorStagePre, aContext))
	{
	    rVal = true;
	    goto RETURN;
	}

	if ((ringStart = CxNodeRing(base)) != NULL)
	{
	    CxtRingObject *ringCur;

	    ringCur = ringStart;
	    do
	    {
		if (CxpTreeIterateRecurse(ringCur, nodeCallback, edgeCallback,
					  ringCallback, aContext))
		{
		    rVal = true;
		    goto RETURN;
		}

		ringCur = CxRingNext(ringCur);
	    } while (ringCur != ringStart);
	}

	if (nodeCallback(base, CxTreeIteratorStagePost, aContext))
	{
	    rVal = true;
	    goto RETURN;
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

bool
CxTreeAuxRegister(CxtTreeObject *self, const char *aKey, void *aData,
		  CxtTreeAuxCleanup *aCleanupFinal,
		  CxtTreeAuxCleanup *aCleanupTree, unsigned *rInd)
{
    bool rVal;

    CxmCheckPtr(aKey);
#ifdef CxmDebug
    {
	unsigned i;

	for (i = 0; i < self->nTreeAux; i++)
	{
	    CxmAssert(self->treeAux[i].key != aKey);
	}
    }
#endif

    // Allocate space for aux registration.
    if (self->treeAux == NULL)
    {
	self->treeAux = (CxtTreeAux *) malloc(sizeof(CxtTreeAux));
	if (self->treeAux == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}
    }
    else
    {
	CxtTreeAux *tTreeAux;

	tTreeAux = (CxtTreeAux *) realloc(self->treeAux,
					  sizeof(CxtTreeAux)
					  * (self->nTreeAux + 1));
	if (self->treeAux == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}

	self->treeAux = tTreeAux;
    }

    self->treeAux[self->nTreeAux].key = aKey;
    self->treeAux[self->nTreeAux].data = aData;
    self->treeAux[self->nTreeAux].cleanupFinal = aCleanupFinal;

    *rInd = self->nTreeAux;
    self->nTreeAux++;

    rVal = false;
    RETURN:
    return rVal;
}

bool
CxTreeAuxSearch(CxtTreeObject *self, const char *aKey, unsigned *rInd)
{
    bool rVal;
    unsigned i;

    for (i = 0; i < self->nTreeAux; i++)
    {
	if (strcmp(self->treeAux[i].key, aKey) == 0)
	{
	    *rInd = i;
	    rVal = false;
	    goto RETURN;
	}
    }

    rVal = true;
    RETURN:
    return rVal;
}

bool
CxTreeNodeAuxRegister(CxtTreeObject *self, const char *aKey, void *aData,
		      CxtTreeAuxCleanup *aCleanupFinal,
		      CxtNodeAuxCleanup *aCleanupNode, unsigned *rInd)
{
    bool rVal;

    CxmCheckPtr(aKey);
#ifdef CxmDebug
    {
	unsigned i;

	for (i = 0; i < self->nNodeAux; i++)
	{
	    CxmAssert(self->nodeAux[i].key != aKey);
	}
    }
#endif

    // Allocate space for aux registration.
    if (self->nodeAux == NULL)
    {
	self->nodeAux = (CxtNodeAux *) malloc(sizeof(CxtNodeAux));
	if (self->nodeAux == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}
    }
    else
    {
	CxtNodeAux *tNodeAux;

	tNodeAux = (CxtNodeAux *) realloc(self->nodeAux,
					  sizeof(CxtNodeAux)
					  * (self->nNodeAux + 1));
	if (self->nodeAux == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}

	self->nodeAux = tNodeAux;
    }

    self->nodeAux[self->nNodeAux].key = aKey;
    self->nodeAux[self->nNodeAux].data = aData;
    self->nodeAux[self->nNodeAux].cleanupFinal = aCleanupFinal;
    self->nodeAux[self->nNodeAux].cleanupNode = aCleanupNode;

    *rInd = self->nNodeAux;
    self->nNodeAux++;

    rVal = false;
    RETURN:
    return rVal;
}

bool
CxTreeNodeAuxSearch(CxtTreeObject *self, const char *aKey, unsigned *rInd)
{
    bool rVal;
    unsigned i;

    for (i = 0; i < self->nNodeAux; i++)
    {
	if (strcmp(self->nodeAux[i].key, aKey) == 0)
	{
	    *rInd = i;
	    rVal = false;
	    goto RETURN;
	}
    }

    rVal = true;
    RETURN:
    return rVal;
}

bool
CxTreeEdgeAuxRegister(CxtTreeObject *self, const char *aKey, void *aData,
		      CxtTreeAuxCleanup *aCleanupFinal,
		      CxtEdgeAuxCleanup *aCleanupEdge, unsigned *rInd)
{
    bool rVal;

    CxmCheckPtr(aKey);
#ifdef CxmDebug
    {
	unsigned i;

	for (i = 0; i < self->nEdgeAux; i++)
	{
	    CxmAssert(self->edgeAux[i].key != aKey);
	}
    }
#endif

    // Allocate space for aux registration.
    if (self->edgeAux == NULL)
    {
	self->edgeAux = (CxtEdgeAux *) malloc(sizeof(CxtEdgeAux));
	if (self->edgeAux == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}
    }
    else
    {
	CxtEdgeAux *tEdgeAux;

	tEdgeAux = (CxtEdgeAux *) realloc(self->edgeAux,
					  sizeof(CxtEdgeAux)
					  * (self->nEdgeAux + 1));
	if (self->edgeAux == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}

	self->edgeAux = tEdgeAux;
    }

    self->edgeAux[self->nEdgeAux].key = aKey;
    self->edgeAux[self->nEdgeAux].data = aData;
    self->edgeAux[self->nEdgeAux].cleanupFinal = aCleanupFinal;
    self->edgeAux[self->nEdgeAux].cleanupEdge = aCleanupEdge;

    *rInd = self->nEdgeAux;
    self->nEdgeAux++;

    rVal = false;
    RETURN:
    return rVal;
}

bool
CxTreeEdgeAuxSearch(CxtTreeObject *self, const char *aKey, unsigned *rInd)
{
    bool rVal;
    unsigned i;

    for (i = 0; i < self->nEdgeAux; i++)
    {
	if (strcmp(self->edgeAux[i].key, aKey) == 0)
	{
	    *rInd = i;
	    rVal = false;
	    goto RETURN;
	}
    }

    rVal = true;
    RETURN:
    return rVal;
}

bool
CxTreeRingAuxRegister(CxtTreeObject *self, const char *aKey, void *aData,
		      CxtTreeAuxCleanup *aCleanupFinal,
		      CxtRingAuxCleanup *aCleanupRing, unsigned *rInd)
{
    bool rVal;

    CxmCheckPtr(aKey);
#ifdef CxmDebug
    {
	unsigned i;

	for (i = 0; i < self->nRingAux; i++)
	{
	    CxmAssert(self->ringAux[i].key != aKey);
	}
    }
#endif

    // Allocate space for aux registration.
    if (self->ringAux == NULL)
    {
	self->ringAux = (CxtRingAux *) malloc(sizeof(CxtRingAux));
	if (self->ringAux == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}
    }
    else
    {
	CxtRingAux *tRingAux;

	tRingAux = (CxtRingAux *) realloc(self->ringAux,
					  sizeof(CxtRingAux)
					  * (self->nRingAux + 1));
	if (self->ringAux == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}

	self->ringAux = tRingAux;
    }

    self->ringAux[self->nRingAux].key = aKey;
    self->ringAux[self->nRingAux].data = aData;
    self->ringAux[self->nRingAux].cleanupFinal = aCleanupFinal;
    self->ringAux[self->nRingAux].cleanupRing = aCleanupRing;

    *rInd = self->nRingAux;
    self->nRingAux++;

    rVal = false;
    RETURN:
    return rVal;
}

bool
CxTreeRingAuxSearch(CxtTreeObject *self, const char *aKey, unsigned *rInd)
{
    bool rVal;
    unsigned i;

    for (i = 0; i < self->nRingAux; i++)
    {
	if (strcmp(self->ringAux[i].key, aKey) == 0)
	{
	    *rInd = i;
	    rVal = false;
	    goto RETURN;
	}
    }

    rVal = true;
    RETURN:
    return rVal;
}

bool
CxTreeAuxSet(CxtTreeObject *self, unsigned aInd, void *aAux)
{
    bool rVal;

    CxmAssert(aInd < self->nTreeAux);

    if (self->nAux <= aInd)
    {
	// Allocate space for aux vector.
	if (self->aux == NULL)
	{
	    self->aux = (void **) calloc(aInd + 1, sizeof(void *));
	    if (self->aux == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }
	}
	else
	{
	    void **tAux;
	    unsigned i;

	    tAux = (void **) realloc(self->aux, (aInd + 1) * sizeof(void *));
	    if (tAux == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }

	    self->aux = tAux;
	    for (i = self->nAux; i < aInd + 1; i++)
	    {
		self->aux[i] = NULL;
	    }
	}
	self->nAux = aInd + 1;
    }

    self->aux[aInd] = aAux;

    rVal = false;
    RETURN:
    return rVal;
}

static PyMethodDef CxpTreeMethods[] =
{
    {
	"ntaxaGet",
	(PyCFunction) CxTreeNtaxaGetPargs,
	METH_NOARGS,
	"ntaxaGet"
    },
    {
	"nedgesGet",
	(PyCFunction) CxTreeNedgesCget,
	METH_NOARGS,
	"nedgesGet"
    },
    {
	"baseGet",
	(PyCFunction) CxTreeBaseGetPargs,
	METH_NOARGS,
	"baseGet"
    },
    {
	"baseSet",
	(PyCFunction) CxTreeBaseSetPargs,
	METH_VARARGS,
	"baseSet"
    },
    {
	"canonize",
	(PyCFunction) CxTreeCanonize,
	METH_NOARGS,
	"canonize"
    },
    {
	"collapse",
	(PyCFunction) CxTreeCollapse,
	METH_NOARGS,
	"collapse"
    },
    {
	"_rfSequence",
	(PyCFunction) CxTreeRfSequence,
	METH_VARARGS,
	"_rfSequence"
    },
    {
	"_rfPair",
	(PyCFunction) CxTreeRfPair,
	METH_VARARGS,
	"_rfPair"
    },
    {
	"tbr",
	(PyCFunction) CxTreeTbrPargs,
	METH_VARARGS,
	"tbr"
    },
    {
	"tbrNneighborsGet",
	(PyCFunction) CxTreeTbrNneighborsGetPargs,
	METH_NOARGS,
	"tbrNneighborsGet"
    },
    {
	"tbrNeighborGet",
	(PyCFunction) CxTreeTbrNeighborGetPargs,
	METH_VARARGS,
	"tbrNeighborGet"
    },
    {
	"_mpPrepare",
	(PyCFunction) CxTreeMpPrepare,
	METH_VARARGS,
	"_mpPrepare"
    },
    {
	"mpFinish",
	(PyCFunction) CxTreeMpFinish,
	METH_NOARGS,
	"mpFinish"
    },
    {
	"mp",
	(PyCFunction) CxTreeMp,
	METH_NOARGS,
	"mp"
    },
    {
	"tbrBestNeighborsMp",
	(PyCFunction) CxTreeTbrBestNeighborsMp,
	METH_VARARGS,
	"tbrBestNeighborsMp"
    },
    {
	"tbrBetterNeighborsMp",
	(PyCFunction) CxTreeTbrBetterNeighborsMp,
	METH_VARARGS,
	"tbrBetterNeighborsMp"
    },
    {
	"tbrAllNeighborsMp",
	(PyCFunction) CxTreeTbrAllNeighborsMp,
	METH_NOARGS,
	"tbrAllNeighborsMp"
    },
    {
	"nheldGet",
	(PyCFunction) CxTreeNheldGet,
	METH_NOARGS,
	"nheldGet"
    },
    {
	"heldGet",
	(PyCFunction) CxTreeheldGet,
	METH_VARARGS,
	"heldGet"
    },
    {NULL, NULL}
};

PyTypeObject CxtTree =
{
    PyObject_HEAD_INIT(NULL)
    0,			// int ob_size
    "C_Tree.C_Tree",	// char *tp_name
    sizeof(CxtTreeObject),	// int tp_basicsize
    0,			// int tp_itemsize
    (destructor) CxpTreeDelete,	// destructor tp_dealloc
    0,			// printfunc tp_print
    0,			// getattrfunc tp_getattr
    0,			// setattrfunc tp_setattr
    0,			// cmpfunc tp_compare
    0,			// reprfunc tp_repr
    0,			// PyNumberMethods *tp_as_number
    0,			// PySequenceMethods *tp_as_sequence
    0,			// PyMappingMethods *tp_as_mapping
    0,			// hashfunc tp_hash
    0,			// ternaryfunc tp_call
    0,			// reprfunc tp_str
    PyObject_GenericGetAttr,	// getattrofunc tp_getattro
    0,			// setattrofunc tp_setattro
    0,			// PyBufferProcs *tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, // long tp_flags
    "Tree(): Create the C portion of a tree.",	// char *tp_doc
    (traverseproc) CxpTreeTraverse,	// traverseproc tp_traverse
    (inquiry) CxpTreeClear,	// inquiry tp_clear
    0,			// richcmpfunc tp_richcompare
    0,			// long tp_weaklistoffset
    0,			// getiterfunc tp_iter
    0,			// iternextfunc tp_iternext
    CxpTreeMethods,	// struct PyMethodDef *tp_methods
    0,			// struct PyMemberDef *tp_members
    0,			// struct PyGetSetDef *tp_getset
    0,			// struct _typeobject *tp_base
    0,			// PyObject *tp_dict
    0,			// descrgetfunc tp_descr_get
    0,			// descrsetfunc tp_descr_set
    0,			// long tp_dictoffset
    0,			// initproc tp_init
    0,			// allocfunc tp_alloc
    CxpTreeNew,		// newfunc tp_new
    _PyObject_Del,	// freefunc tp_free
    0			// inquiry tp_is_gc
};

static PyMethodDef CxpTreeFuncs[] =
{
    {NULL}
};

PyObject *CxgTreeException;
PyObject *CxgTreeValueError;
PyObject *CxgTreeTypeError;

void
CxTreeInit(void)
{
    PyObject *m;

    // Create new type.
    if (PyType_Ready(&CxtTree) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_Tree", CxpTreeFuncs, "Tree extensions");
    Py_INCREF(&CxtTree);
    PyModule_AddObject(m, "C_Tree", (PyObject *) &CxtTree);

    // Create exception objects.
    // Exception.
    CxgTreeException = PyErr_NewException("C_Tree.Exception", CxgException,
					  NULL);
    Py_INCREF(CxgTreeException);
    PyModule_AddObject(m, "Exception", CxgTreeException);

    // ValueError.
    CxgTreeValueError = PyErr_NewException("C_Tree.ValueError",
					   CxgTreeException,
					   NULL);
    Py_INCREF(CxgTreeValueError);
    PyModule_AddObject(m, "ValueError", CxgTreeValueError);

    // TypeError.
    CxgTreeTypeError = PyErr_NewException("C_Tree.TypeError",
					  CxgTreeException,
					  NULL);
    Py_INCREF(CxgTreeTypeError);
    PyModule_AddObject(m, "TypeError", CxgTreeTypeError);

    // Pre-compile Python code that is used for creating a tree.
    CxpTreeNewCode = Py_CompileString("\
import crux.Tree\n\
tree = crux.Tree.Tree()\n\
",
				      "<string>",
				      Py_file_input);
    CxmCheckPtr(CxpTreeNewCode);
}

// End Tree.
//==============================================================================
// Begin Node.

static PyObject *
CxpNodeNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtNodeObject *self;
    CxtTreeObject *tree;

    if (PyArg_ParseTuple(args, "O!", &CxtTree, &tree) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    self = (CxtNodeObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	Py_INCREF(tree);
	self->tree = tree;
	self->node = CxTrNodeNew(tree->tr);
	CxTrNodeAuxSet(tree->tr, self->node, self);
	self->aux = NULL;
	self->nAux = 0;

	self->GcCleared = false;
	rVal = (PyObject *) self;
    }
    CxmXepCatch(CxmXepOOM)
    {
	Py_DECREF(tree);
	CxmXepHandled();
	rVal = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, rVal, rVal->ob_refcnt);
#endif
    return rVal;
}

static int
CxpNodeTraverse(CxtNodeObject *self, visitproc visit, void *arg)
{
    int rVal;
    CxtTrRing trRing, trCurRing;
    CxtRingObject *ring;

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    if (self->GcCleared == false)
    {
	if (visit((PyObject *) self->tree, arg) < 0)
	{
	    rVal = -1;
	    goto RETURN;
	}

	// Report all rings.  It is not good enough to simply report one, since
	// Python's mark/sweep GC apparently keeps track of how many times each
	// object is visited.
	trRing = CxTrNodeRingGet(self->tree->tr, self->node);
	if (trRing != CxmTrRingNone)
	{
	    trCurRing = trRing;
	    do
	    {
		ring = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr,
							trCurRing);
		if (visit((PyObject *) ring, arg) < 0)
		{
		    rVal = -1;
		    goto RETURN;
		}

		trCurRing = CxTrRingNextGet(self->tree->tr, trCurRing);
	    } while (trCurRing != trRing);
	}
    }

    rVal = 0;
    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif
    return rVal;
}

static int
CxpNodeClear(CxtNodeObject *self)
{
    CxtTrRing trRing, trCurRing;
    CxtRingObject *ring;
    CxtTrNode base;

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    if (self->GcCleared == false)
    {
	// Detach from rings.
	trRing = CxTrNodeRingGet(self->tree->tr, self->node);
	if (trRing != CxmTrRingNone)
	{
	    trCurRing = trRing;
	    do
	    {
		ring = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr,
							trCurRing);

		// Get next ring before detaching.
		trRing = trCurRing;
		trCurRing = CxTrRingNextGet(self->tree->tr, trCurRing);

		CxEdgeDetach(ring->edge);
	    } while (trCurRing != trRing);
	}

	// Detach from tree if tree base.
	base = CxTrBaseGet(self->tree->tr);
	if (base == self->node)
	{
	    CxtNodeObject *node;

	    node = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, base);
	    CxTrBaseSet(self->tree->tr, CxmTrNodeNone);
	    Py_DECREF(node);
	}

	// Delete node.
	CxTrNodeDelete(self->tree->tr, self->node);

	// Drop reference to tree.
	Py_DECREF(self->tree);

	self->GcCleared = true;
    }

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif
    return 0;
}

static void
CxpNodeDelete(CxtNodeObject *self)
{
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    CxpNodeClear(self);

    if (self->aux != NULL)
    {
	unsigned i;

	for (i = 0; i < self->nAux; i++)
	{
	    if (self->tree->nodeAux[i].cleanupNode != NULL)
	    {
		self->tree->nodeAux[i].cleanupNode(self, self->aux[i]);
	    }
	}

	free(self->aux);
    }

    self->ob_type->tp_free((PyObject*) self);

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p\n",
 	    __FILE__, __LINE__, __func__, self);
#endif
}

static PyObject *CxpNodeNewCode;

CxtNodeObject *
CxNodeNew(CxtTreeObject *aTree)
{
    CxtNodeObject *rVal;
    PyObject *globals, *locals, *obj;

    globals = PyEval_GetGlobals();
    if (globals == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }
    locals = Py_BuildValue("{sO}", "tree", (PyObject *) aTree);
    if (locals == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }

    obj = PyEval_EvalCode((PyCodeObject *) CxpNodeNewCode,
			  globals,
			  locals);
    if (obj == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }
    Py_DECREF(obj);

    rVal = (CxtNodeObject *) PyDict_GetItemString(locals, "node");
    if (rVal == NULL)
    {
	goto RETURN;
    }
    Py_INCREF(rVal);

    RETURN:
    Py_DECREF(locals);
    return rVal;
}

PyObject *
CxNodeTree(CxtNodeObject *self)
{
    return Py_BuildValue("O", self->tree);
}

uint32_t
CxNodeTaxonNumGet(CxtNodeObject *self)
{
    return CxTrNodeTaxonNumGet(self->tree->tr, self->node);
}

PyObject *
CxNodeTaxonNumGetPargs(CxtNodeObject *self)
{
    PyObject *rVal;
    uint32_t taxonNum;

    taxonNum = CxTrNodeTaxonNumGet(self->tree->tr, self->node);

    if (taxonNum == CxmTrNodeTaxonNone)
    {
	Py_INCREF(Py_None);
	rVal = Py_None;
    }
    else
    {
	rVal = Py_BuildValue("i", taxonNum);
    }

    return rVal;
}

void
CxNodeTaxonNumSet(CxtNodeObject *self, uint32_t aTaxonNum)
{
    CxTrNodeTaxonNumSet(self->tree->tr, self->node, aTaxonNum);
    self->tree->seq++;
}

PyObject *
CxNodeTaxonNumSetPargs(CxtNodeObject *self, PyObject *args)
{
    PyObject *rVal;
    uint32_t taxonNum;

    taxonNum = CxmTrNodeTaxonNone;
    if (PyArg_ParseTuple(args, "|i", &taxonNum) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    CxNodeTaxonNumSet(self, taxonNum);

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
}

CxtRingObject *
CxNodeRing(CxtNodeObject *self)
{
    CxtRingObject *rVal;
    CxtTrRing trRing;

    trRing = CxTrNodeRingGet(self->tree->tr, self->node);
    if (trRing != CxmTrRingNone)
    {
	rVal = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRing);
    }
    else
    {
	rVal = NULL;
    }

    return rVal;
}

PyObject *
CxNodeRingPargs(CxtNodeObject *self)
{
    PyObject *rVal;
    CxtRingObject *ring;
    CxtTrRing trRing;

    trRing = CxTrNodeRingGet(self->tree->tr, self->node);
    if (trRing != CxmTrRingNone)
    {
	ring = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRing);
	Py_INCREF(ring);
	rVal = (PyObject *) ring;
    }
    else
    {
	Py_INCREF(Py_None);
	rVal = Py_None;
    }

    return rVal;
}

void
CxNodeRingSet(CxtNodeObject *self, CxtRingObject *aRing)
{
    CxmAssert(CxRingNode(aRing) == self);

    CxTrNodeRingSet(self->tree->tr, self->node, aRing->ring);
}

unsigned
CxNodeDegree(CxtNodeObject *self)
{
    return CxTrNodeDegree(self->tree->tr, self->node);
}

PyObject *
CxNodeDegreePargs(CxtNodeObject *self)
{
    return Py_BuildValue("i", CxTrNodeDegree(self->tree->tr, self->node));
}

bool
CxNodeAuxSet(CxtNodeObject *self, unsigned aInd, void *aAux)
{
    bool rVal;

    CxmAssert(aInd < self->tree->nNodeAux);

    if (self->nAux <= aInd)
    {
	// Allocate space for aux vector.
	if (self->aux == NULL)
	{
	    self->aux = (void **) calloc(aInd + 1, sizeof(void *));
	    if (self->aux == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }
	}
	else
	{
	    void **tAux;
	    unsigned i;

	    tAux = (void **) realloc(self->aux, (aInd + 1) * sizeof(void *));
	    if (tAux == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }

	    self->aux = tAux;
	    for (i = self->nAux; i < aInd + 1; i++)
	    {
		self->aux[i] = NULL;
	    }
	}
	self->nAux = aInd + 1;
    }

    self->aux[aInd] = aAux;

    rVal = false;
    RETURN:
    return rVal;
}

static PyMethodDef CxpNodeMethods[] =
{
    {
	"tree",
	(PyCFunction) CxNodeTree,
	METH_NOARGS,
	"tree"
    },
    {
	"taxonNumGet",
	(PyCFunction) CxNodeTaxonNumGetPargs,
	METH_NOARGS,
	"taxonNumGet"
    },
    {
	"taxonNumSet",
	(PyCFunction) CxNodeTaxonNumSetPargs,
	METH_VARARGS,
	"taxonNumSet"
    },
    {
	"ring",
	(PyCFunction) CxNodeRingPargs,
	METH_NOARGS,
	"ring"
    },
    {
	"degree",
	(PyCFunction) CxNodeDegreePargs,
	METH_NOARGS,
	"degree"
    },
    {NULL, NULL}
};

PyTypeObject CxtNode =
{
    PyObject_HEAD_INIT(NULL)
    0,			// int ob_size
    "C_Node.C_Node",	// char *tp_name
    sizeof(CxtNodeObject),	// int tp_basicsize
    0,			// int tp_itemsize
    (destructor) CxpNodeDelete,	// destructor tp_dealloc
    0,			// printfunc tp_print
    0,			// getattrfunc tp_getattr
    0,			// setattrfunc tp_setattr
    0,			// cmpfunc tp_compare
    0,			// reprfunc tp_repr
    0,			// PyNumberMethods *tp_as_number
    0,			// PySequenceMethods *tp_as_sequence
    0,			// PyMappingMethods *tp_as_mapping
    0,			// hashfunc tp_hash
    0,			// ternaryfunc tp_call
    0,			// reprfunc tp_str
    PyObject_GenericGetAttr,	// getattrofunc tp_getattro
    0,			// setattrofunc tp_setattro
    0,			// PyBufferProcs *tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, // long tp_flags
    "Node(): Create the C portion of a node.",	// char *tp_doc
    (traverseproc) CxpNodeTraverse,	// traverseproc tp_traverse
    (inquiry) CxpNodeClear,	// inquiry tp_clear
    0,			// richcmpfunc tp_richcompare
    0,			// long tp_weaklistoffset
    0,			// getiterfunc tp_iter
    0,			// iternextfunc tp_iternext
    CxpNodeMethods,	// struct PyMethodDef *tp_methods
    0,			// struct PyMemberDef *tp_members
    0,			// struct PyGetSetDef *tp_getset
    0,			// struct _typeobject *tp_base
    0,			// PyObject *tp_dict
    0,			// descrgetfunc tp_descr_get
    0,			// descrsetfunc tp_descr_set
    0,			// long tp_dictoffset
    0,			// initproc tp_init
    0,			// allocfunc tp_alloc
    CxpNodeNew,		// newfunc tp_new
    _PyObject_Del,	// freefunc tp_free
    0			// inquiry tp_is_gc
};

static PyMethodDef CxpNodeFuncs[] =
{
    {NULL}
};

void
CxNodeInit(void)
{
    PyObject *m;

    // Create new type.
    if (PyType_Ready(&CxtNode) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_Node", CxpNodeFuncs, "Node extensions");
    Py_INCREF(&CxtNode);
    PyModule_AddObject(m, "C_Node", (PyObject *) &CxtNode);

    // Pre-compile Python code that is used for creating a node.
    CxpNodeNewCode = Py_CompileString("\
import crux.Node\n\
node = crux.Node.Node(tree)\n\
",
				      "<string>",
				      Py_file_input);
    CxmCheckPtr(CxpNodeNewCode);
}

// End Node.
//==============================================================================
// Begin Edge.

static PyObject *
CxpEdgeNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtEdgeObject *self;
    CxtTreeObject *tree;
    volatile uint32_t tryStage = 0;

    if (PyArg_ParseTuple(args, "O!", &CxtTree, &tree) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    self = (CxtEdgeObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	// Avoid traversing this object in GC until it is fully constructed, and
	// its rings are fully constructed.
	PyObject_GC_UnTrack(self);

	Py_INCREF(tree);
	self->tree = tree;
	self->edge = CxTrEdgeNew(tree->tr);
	CxTrEdgeAuxSet(tree->tr, self->edge, self);
	self->aux = NULL;
	self->nAux = 0;
	tryStage = 1;

	// Create associated ring objects.
	self->ringA = CxRingNew(self, 0);
	tryStage = 2;

	// Avoid traversing ringA until after ringB has been constructed.
	PyObject_GC_UnTrack(self->ringA);

	self->ringB = CxRingNew(self, 1);

	self->GcDetached = false;
	self->GcCleared = false;
	rVal = (PyObject *) self;

	// It's okay to traverse the edge and rings now.
	PyObject_GC_Track(self->ringA);
	PyObject_GC_Track(self);
    }
    CxmXepCatch(CxmXepOOM)
    {
	switch (tryStage)
	{
	    case 2:
	    {
		Py_DECREF(self->ringA);
	    }
	    case 1:
	    {
		CxTrEdgeDelete(self->tree->tr, self->edge);
		type->tp_free((PyObject *) self);
		Py_DECREF(tree);
	    }
	    case 0:
	    {
		break;
	    }
	    default:
	    {
		CxmNotReached();
	    }
	}
	CxmXepHandled();
	rVal = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, rVal, rVal->ob_refcnt);
#endif
    return rVal;
}

static int
CxpEdgeTraverse(CxtEdgeObject *self, visitproc visit, void *arg)
{
    int rVal;

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    if (self->GcCleared == false)
    {
	if (visit((PyObject *) self->tree, arg) < 0)
	{
	    rVal = -1;
	    goto RETURN;
	}

	if (visit((PyObject *) self->ringA, arg) < 0)
	{
	    rVal = -1;
	    goto RETURN;
	}

	if (visit((PyObject *) self->ringB, arg) < 0)
	{
	    rVal = -1;
	    goto RETURN;
	}
    }

    rVal = 0;
    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif
    return rVal;
}

static int
CxpEdgeClear(CxtEdgeObject *self)
{
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    if (self->GcCleared == false)
    {
	// Detach from nodes, if not already done.
	if (self->GcDetached == false)
	{
	    self->GcDetached = true;
	    self->ringA->GcDetached = true;
	    self->ringB->GcDetached = true;

	    if (CxTrRingNodeGet(self->tree->tr, self->ringA->ring)
		!= CxmTrNodeNone)
	    {
		CxEdgeDetach(self);
	    }
	}

	// Drop references to rings.
	Py_DECREF(self->ringA);
	Py_DECREF(self->ringB);

	self->GcCleared = true;
    }

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif
    return 0;
}

static void
CxpEdgeDelete(CxtEdgeObject *self)
{
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    CxpEdgeClear(self);

    if (self->aux != NULL)
    {
	unsigned i;

	for (i = 0; i < self->nAux; i++)
	{
	    if (self->tree->edgeAux[i].cleanupEdge != NULL)
	    {
		self->tree->edgeAux[i].cleanupEdge(self, self->aux[i]);
	    }
	}

	free(self->aux);
    }

    // Delete the CxTrEdge (and associated CxTrRing objects).  The CxtRingObject
    // clear/delete code held onto a reference to this CxtEdgeOjbect long enough
    // to ensure that both CxtRingObject's associated with this CxtEdgeObject
    // have already been deallocated.
    CxTrEdgeDelete(self->tree->tr, self->edge);

    // Drop reference to tree.
    Py_DECREF(self->tree);

    self->ob_type->tp_free((PyObject*) self);

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p\n",
 	    __FILE__, __LINE__, __func__, self);
#endif
}

static PyObject *CxpEdgeNewCode;

CxtEdgeObject *
CxEdgeNew(CxtTreeObject *aTree)
{
    CxtEdgeObject *rVal;
    PyObject *globals, *locals, *obj;

    globals = PyEval_GetGlobals();
    if (globals == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }
    locals = Py_BuildValue("{sO}", "tree", (PyObject *) aTree);
    if (locals == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }

    obj = PyEval_EvalCode((PyCodeObject *) CxpEdgeNewCode,
			  globals,
			  locals);
    if (obj == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }
    Py_DECREF(obj);

    rVal = (CxtEdgeObject *) PyDict_GetItemString(locals, "edge");
    if (rVal == NULL)
    {
	goto RETURN;
    }
    Py_INCREF(rVal);

    RETURN:
    Py_DECREF(locals);
    return rVal;
}

PyObject *
CxEdgeTree(CxtEdgeObject *self)
{
    return Py_BuildValue("O", self->tree);
}

void
CxEdgeRingsGet(CxtEdgeObject *self, CxtRingObject **rRingA,
	       CxtRingObject **rRingB)
{
    *rRingA = self->ringA;
    *rRingB = self->ringB;
}

PyObject *
CxEdgeRingsGetPargs(CxtEdgeObject *self)
{
    return Py_BuildValue("(OO)", self->ringA, self->ringB);
}

double
CxEdgeLengthGet(CxtEdgeObject *self)
{
    return CxTrEdgeLengthGet(self->tree->tr, self->edge);
}

PyObject *
CxEdgeLengthGetPargs(CxtEdgeObject *self)
{
    return Py_BuildValue("d", CxTrEdgeLengthGet(self->tree->tr, self->edge));
}

void
CxEdgeLengthSet(CxtEdgeObject *self, double aLength)
{
    CxTrEdgeLengthSet(self->tree->tr, self->edge, aLength);
    self->tree->seq++;
}

PyObject *
CxEdgeLengthSetPargs(CxtEdgeObject *self, PyObject *args)
{
    PyObject *rVal;
    double length;

    length = 0.0;
    if (PyArg_ParseTuple(args, "|d", &length) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    CxEdgeLengthSet(self, length);

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
}

void
CxEdgeAttach(CxtEdgeObject *self, CxtNodeObject *aNodeA,
	     CxtNodeObject *aNodeB)
{
    // Rings refer to nodes.
    Py_INCREF(aNodeA);
    Py_INCREF(aNodeB);

    // Nodes refer to rings.  In actuality, a node only refers to one element of
    // a ring, but doing all the "correct" reference management for the ring
    // would be hard (and slow), so instead, pretend that a node refers to each
    // ring element.
    Py_INCREF(self->ringA);
    Py_INCREF(self->ringB);

    // Attach.
    //
    // CxpTreeCanonize() assumes that attaching inserts the edge at the tail
    // of the ring; make sure to keep that code in sync with this function.
    CxTrEdgeAttach(self->tree->tr, self->edge, aNodeA->node, aNodeB->node);

    self->tree->seq++;
}

PyObject *
CxEdgeAttachPargs(CxtEdgeObject *self, PyObject *args)
{
    PyObject *rVal;
    CxtNodeObject *nodeA, *nodeB;

    if (PyArg_ParseTuple(args, "O!O!", &CxtNode, &nodeA,
			 &CxtNode, &nodeB) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }
    // Make sure that the edge is currently detached.
    if (CxTrRingNodeGet(self->tree->tr, self->ringA->ring) != CxmTrNodeNone)
    {
	CxError(CxgEdgeValueError, "Edge is already attached");
	rVal = NULL;
	goto RETURN;
    }
    if (nodeA->tree != self->tree || nodeB->tree != self->tree)
    {
	CxError(CxgEdgeValueError, "Node does not belong to this tree");
	rVal = NULL;
	goto RETURN;
    }

    CxEdgeAttach(self, nodeA, nodeB);

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
}

bool
CxEdgeDetach(CxtEdgeObject *self)
{
    bool rVal;
    CxtTrNode trNodeA, trNodeB;
    CxtNodeObject *nodeA, *nodeB;
    CxtTrRing trRingA, trRingB;

    // Make sure that the edge is currently attached.
    if (CxTrRingNodeGet(self->tree->tr, self->ringA->ring) == CxmTrNodeNone)
    {
	rVal = true;
	goto RETURN;
    }

    trRingA = CxTrEdgeRingGet(self->tree->tr, self->edge, 0);
    CxmAssert(trRingA != CxmTrRingNone);
    trRingB = CxTrEdgeRingGet(self->tree->tr, self->edge, 1);

    trNodeA = CxTrRingNodeGet(self->tree->tr, trRingA);
    nodeA = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, trNodeA);

    trNodeB = CxTrRingNodeGet(self->tree->tr, trRingB);
    nodeB = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, trNodeB);

    // Detach.
    CxTrEdgeDetach(self->tree->tr, self->edge);

    // Rings refer to nodes.
    Py_DECREF(nodeA);
    Py_DECREF(nodeB);

    // Nodes refer to rings.
    Py_DECREF(self->ringA);
    Py_DECREF(self->ringB);

    self->tree->seq++;

    rVal = false;
    RETURN:
    return rVal;
}

PyObject *
CxEdgeDetachPargs(CxtEdgeObject *self)
{
    PyObject *rVal;

    if (CxEdgeDetach(self))
    {
	CxError(CxgEdgeValueError, "Edge is already detached");
	rVal = NULL;
    }
    else
    {
	Py_INCREF(Py_None);
	rVal = Py_None;
    }

    return rVal;
}

bool
CxEdgeAuxSet(CxtEdgeObject *self, unsigned aInd, void *aAux)
{
    bool rVal;

    CxmAssert(aInd < self->tree->nEdgeAux);

    if (self->nAux <= aInd)
    {
	// Allocate space for aux vector.
	if (self->aux == NULL)
	{
	    self->aux = (void **) calloc(aInd + 1, sizeof(void *));
	    if (self->aux == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }
	}
	else
	{
	    void **tAux;
	    unsigned i;

	    tAux = (void **) realloc(self->aux, (aInd + 1) * sizeof(void *));
	    if (tAux == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }

	    self->aux = tAux;
	    for (i = self->nAux; i < aInd + 1; i++)
	    {
		self->aux[i] = NULL;
	    }
	}
	self->nAux = aInd + 1;
    }

    self->aux[aInd] = aAux;

    rVal = false;
    RETURN:
    return rVal;
}

static PyMethodDef CxpEdgeMethods[] =
{
    {
	"tree",
	(PyCFunction) CxEdgeTree,
	METH_NOARGS,
	"tree"
    },
    {
	"rings",
	(PyCFunction) CxEdgeRingsGetPargs,
	METH_NOARGS,
	"rings"
    },
    {
	"lengthGet",
	(PyCFunction) CxEdgeLengthGetPargs,
	METH_NOARGS,
	"lengthGet"
    },
    {
	"lengthSet",
	(PyCFunction) CxEdgeLengthSetPargs,
	METH_VARARGS,
	"lengthSet"
    },
    {
	"attach",
	(PyCFunction) CxEdgeAttachPargs,
	METH_VARARGS,
	"attach"
    },
    {
	"detach",
	(PyCFunction) CxEdgeDetachPargs,
	METH_NOARGS,
	"detach"
    },
    {NULL, NULL}
};

PyTypeObject CxtEdge =
{
    PyObject_HEAD_INIT(NULL)
    0,			// int ob_size
    "C_Edge.C_Edge",	// char *tp_name
    sizeof(CxtEdgeObject),	// int tp_basicsize
    0,			// int tp_itemsize
    (destructor) CxpEdgeDelete,	// destructor tp_dealloc
    0,			// printfunc tp_print
    0,			// getattrfunc tp_getattr
    0,			// setattrfunc tp_setattr
    0,			// cmpfunc tp_compare
    0,			// reprfunc tp_repr
    0,			// PyNumberMethods *tp_as_number
    0,			// PySequenceMethods *tp_as_sequence
    0,			// PyMappingMethods *tp_as_mapping
    0,			// hashfunc tp_hash
    0,			// ternaryfunc tp_call
    0,			// reprfunc tp_str
    PyObject_GenericGetAttr,	// getattrofunc tp_getattro
    0,			// setattrofunc tp_setattro
    0,			// PyBufferProcs *tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, // long tp_flags
    "Edge(): Create the C portion of an edge.",	// char *tp_doc
    (traverseproc) CxpEdgeTraverse,	// traverseproc tp_traverse
    (inquiry) CxpEdgeClear,	// inquiry tp_clear
    0,			// richcmpfunc tp_richcompare
    0,			// long tp_weaklistoffset
    0,			// getiterfunc tp_iter
    0,			// iternextfunc tp_iternext
    CxpEdgeMethods,	// struct PyMethodDef *tp_methods
    0,			// struct PyMemberDef *tp_members
    0,			// struct PyGetSetDef *tp_getset
    0,			// struct _typeobject *tp_base
    0,			// PyObject *tp_dict
    0,			// descrgetfunc tp_descr_get
    0,			// descrsetfunc tp_descr_set
    0,			// long tp_dictoffset
    0,			// initproc tp_init
    0,			// allocfunc tp_alloc
    CxpEdgeNew,		// newfunc tp_new
    _PyObject_Del,	// freefunc tp_free
    0			// inquiry tp_is_gc
};

static PyMethodDef CxpEdgeFuncs[] =
{
    {NULL}
};

PyObject *CxgEdgeException;
PyObject *CxgEdgeValueError;

void
CxEdgeInit(void)
{
    PyObject *m;

    // Create new type.
    if (PyType_Ready(&CxtEdge) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_Edge", CxpEdgeFuncs, "Edge extensions");
    Py_INCREF(&CxtEdge);
    PyModule_AddObject(m, "C_Edge", (PyObject *) &CxtEdge);

    // Create exception objects.
    CxgEdgeException = PyErr_NewException("C_Edge.Exception", CxgException,
					  NULL);
    Py_INCREF(CxgEdgeException);
    PyModule_AddObject(m, "Exception", CxgEdgeException);

    CxgEdgeValueError = PyErr_NewException("C_Edge.ValueError",
					   CxgEdgeException,
					   NULL);
    Py_INCREF(CxgEdgeValueError);
    PyModule_AddObject(m, "ValueError", CxgEdgeValueError);

    // Pre-compile Python code that is used for creating a wrapped edge.
    CxpEdgeNewCode = Py_CompileString("\
import crux.Edge\n\
edge = crux.Edge.Edge(tree)\n\
",
				      "<string>",
				      Py_file_input);
    CxmCheckPtr(CxpEdgeNewCode);
}

// End Edge.
//==============================================================================
// Begin Ring.

static PyObject *
CxpRingNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtRingObject *self;
    CxtEdgeObject *edge;
    uint32_t end;

    if (PyArg_ParseTuple(args, "O!i", &CxtEdge, &edge, &end) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    self = (CxtRingObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }

    if (end != 0 && end != 1)
    {
	type->tp_free((PyObject *) self);

	CxError(CxgRingValueError, "End must be 0 or 1");
	rVal = NULL;
	goto RETURN;
    }

    Py_INCREF(edge->tree);
    Py_INCREF(edge);

    self->tree = edge->tree;
    self->edge = edge;

    Py_INCREF(self);
    Py_INCREF(self);

    self->ring = CxTrEdgeRingGet(edge->tree->tr, edge->edge, end);
    if (CxTrRingAuxGet(self->tree->tr, self->ring) != NULL)
    {
	Py_DECREF(self);
	Py_DECREF(self);
	Py_DECREF(edge);
	Py_DECREF(edge->tree);
	type->tp_free((PyObject *) self);

	CxError(CxgRingValueError, "Internal error (aux should be NULL)");
	rVal = NULL;
	goto RETURN;
    }
    CxTrRingAuxSet(self->tree->tr, self->ring, self);
    self->aux = NULL;
    self->nAux = 0;

    self->GcDetached = false;
    self->GcCleared = false;

    rVal = (PyObject *) self;

    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, rVal, rVal->ob_refcnt);
#endif
    return rVal;
}

static int
CxpRingTraverse(CxtRingObject *self, visitproc visit, void *arg)
{
    int rVal;
    CxtTrNode trNode;
    CxtTrRing trRing;
    CxtRingObject *ring;

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    if (self->GcCleared == false)
    {
	if (visit((PyObject *) self->tree, arg) < 0)
	{
	    rVal = -1;
	    goto RETURN;
	}

	trNode = CxTrRingNodeGet(self->tree->tr, self->ring);
	if (trNode != CxmTrNodeNone)
	{
	    CxtNodeObject *node;

	    node = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, trNode);
	    if (visit((PyObject *) node, arg) < 0)
	    {
		rVal = -1;
		goto RETURN;
	    }
	}

	// Report edge.
	if (visit((PyObject *) self->edge, arg) < 0)
	{
	    rVal = -1;
	    goto RETURN;
	}

	// Report next ring element.
	trRing = CxTrRingNextGet(self->tree->tr, self->ring);
	ring = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRing);
	if (visit((PyObject *) ring, arg) < 0)
	{
	    rVal = -1;
	    goto RETURN;
	}

	// Report previous ring element.
	trRing = CxTrRingPrevGet(self->tree->tr, self->ring);
	ring = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRing);
	if (visit((PyObject *) ring, arg) < 0)
	{
	    rVal = -1;
	    goto RETURN;
	}
    }

    rVal = 0;
    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif
    return rVal;
}

static int
CxpRingClear(CxtRingObject *self)
{
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    if (self->GcCleared == false)
    {
	// Detach from nodes, if not already done.
	if (self->GcDetached == false)
	{
	    self->edge->GcDetached = true;
	    self->edge->ringA->GcDetached = true;
	    self->edge->ringB->GcDetached = true;

	    if (CxTrRingNodeGet(self->tree->tr, self->ring) != CxmTrNodeNone)
	    {
		CxEdgeDetach(self->edge);
	    }
	}

	// Drop reference to tree.
	Py_DECREF(self->tree);

	// Drop references to self.
	Py_DECREF(self);
	Py_DECREF(self);

	self->GcCleared = true;
    }

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif
    return 0;
}

static void
CxpRingDelete(CxtRingObject *self)
{
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    CxpRingClear(self);

    if (self->aux != NULL)
    {
	unsigned i;

	for (i = 0; i < self->nAux; i++)
	{
	    if (self->tree->ringAux[i].cleanupRing != NULL)
	    {
		self->tree->ringAux[i].cleanupRing(self, self->aux[i]);
	    }
	}

	free(self->aux);
    }

    // Drop the reference to the associated edge, now that there is no chance of
    // accessing it again.
    Py_DECREF(self->edge);

    self->ob_type->tp_free((PyObject*) self);

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p\n",
 	    __FILE__, __LINE__, __func__, self);
#endif
}

static PyObject *CxpRingNewCode;

CxtRingObject *
CxRingNew(CxtEdgeObject *aEdge, uint32_t aEnd)
{
    CxtRingObject *rVal;
    PyObject *globals, *locals, *obj;

    globals = PyEval_GetGlobals();
    if (globals == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }
    locals = Py_BuildValue("{sOsi}",
			   "edge", (PyObject *) aEdge,
			   "end", aEnd);
    if (locals == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }

    obj = PyEval_EvalCode((PyCodeObject *) CxpRingNewCode,
			  globals,
			  locals);
    if (obj == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }
    Py_DECREF(obj);

    rVal = (CxtRingObject *) PyDict_GetItemString(locals, "ring");
    if (rVal == NULL)
    {
	goto RETURN;
    }
    Py_INCREF(rVal);

    RETURN:
    Py_DECREF(locals);
    return rVal;
}

CxtTreeObject *
CxRingTree(CxtRingObject *self)
{
    return self->tree;
}

PyObject *
CxRingTreePargs(CxtRingObject *self)
{
    return Py_BuildValue("O", self->tree);
}

CxtNodeObject *
CxRingNode(CxtRingObject *self)
{
    CxtNodeObject *rVal;
    CxtTrNode node;

    node = CxTrRingNodeGet(self->tree->tr, self->ring);
    if (node != CxmTrNodeNone)
    {
	rVal = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, node);
    }
    else
    {
	rVal = NULL;
    }

    return rVal;
}

PyObject *
CxRingNodePargs(CxtRingObject *self)
{
    PyObject *rVal;
    CxtTrNode node;

    node = CxTrRingNodeGet(self->tree->tr, self->ring);
    if (node != CxmTrNodeNone)
    {
	rVal = (PyObject *) CxTrNodeAuxGet(self->tree->tr, node);
	Py_INCREF(rVal);
    }
    else
    {
	Py_INCREF(Py_None);
	rVal = Py_None;
    }

    return rVal;
}

CxtEdgeObject *
CxRingEdge(CxtRingObject *self)
{
    CxtEdgeObject *rVal;
    CxtTrEdge edge;

    edge = CxTrRingEdgeGet(self->tree->tr, self->ring);
    rVal = (CxtEdgeObject *) CxTrEdgeAuxGet(self->tree->tr, edge);

    return rVal;
}

PyObject *
CxRingEdgePargs(CxtRingObject *self)
{
    PyObject *rVal;
    CxtTrEdge edge;

    edge = CxTrRingEdgeGet(self->tree->tr, self->ring);
    rVal = (PyObject *) CxTrEdgeAuxGet(self->tree->tr, edge);
    Py_INCREF(rVal);

    return rVal;
}

CxtRingObject *
CxRingOther(CxtRingObject *self)
{
    CxtRingObject *rVal;
    CxtTrRing other;

    other = CxTrRingOtherGet(self->tree->tr, self->ring);
    rVal = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, other);

    return rVal;
}

PyObject *
CxRingOtherPargs(CxtRingObject *self)
{
    PyObject *rVal;
    CxtTrRing other;

    other = CxTrRingOtherGet(self->tree->tr, self->ring);
    rVal = (PyObject *) CxTrRingAuxGet(self->tree->tr, other);
    Py_INCREF(rVal);

    return rVal;
}

CxtRingObject *
CxRingNext(CxtRingObject *self)
{
    CxtRingObject *rVal;
    CxtTrRing nextRing;

    nextRing = CxTrRingNextGet(self->tree->tr, self->ring);
    rVal = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, nextRing);

    return rVal;
}

PyObject *
CxRingNextPargs(CxtRingObject *self)
{
    PyObject *rVal;
    CxtTrRing nextRing;

    nextRing = CxTrRingNextGet(self->tree->tr, self->ring);
    rVal = (PyObject *) CxTrRingAuxGet(self->tree->tr, nextRing);
    Py_INCREF(rVal);

    return rVal;
}

CxtRingObject *
CxRingPrev(CxtRingObject *self)
{
    CxtRingObject *rVal;
    CxtTrRing prevRing;

    prevRing = CxTrRingPrevGet(self->tree->tr, self->ring);
    rVal = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, prevRing);

    return rVal;
}

PyObject *
CxRingPrevPargs(CxtRingObject *self)
{
    PyObject *rVal;
    CxtTrRing prevRing;

    prevRing = CxTrRingPrevGet(self->tree->tr, self->ring);
    rVal = (PyObject *) CxTrRingAuxGet(self->tree->tr, prevRing);
    Py_INCREF(rVal);

    return rVal;
}

bool
CxRingAuxSet(CxtRingObject *self, unsigned aInd, void *aAux)
{
    bool rVal;

    CxmAssert(aInd < self->tree->nRingAux);

    if (self->nAux <= aInd)
    {
	// Allocate space for aux vector.
	if (self->aux == NULL)
	{
	    self->aux = (void **) calloc(aInd + 1, sizeof(void *));
	    if (self->aux == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }
	}
	else
	{
	    void **tAux;
	    unsigned i;

	    tAux = (void **) realloc(self->aux, (aInd + 1) * sizeof(void *));
	    if (tAux == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }

	    self->aux = tAux;
	    for (i = self->nAux; i < aInd + 1; i++)
	    {
		self->aux[i] = NULL;
	    }
	}
	self->nAux = aInd + 1;
    }

    self->aux[aInd] = aAux;

    rVal = false;
    RETURN:
    return rVal;
}

static PyMethodDef CxpRingMethods[] =
{
    {
	"tree",
	(PyCFunction) CxRingTreePargs,
	METH_NOARGS,
	"tree"
    },
    {
	"node",
	(PyCFunction) CxRingNodePargs,
	METH_NOARGS,
	"node"
    },
    {
	"edge",
	(PyCFunction) CxRingEdgePargs,
	METH_NOARGS,
	"edge"
    },
    {
	"other",
	(PyCFunction) CxRingOtherPargs,
	METH_NOARGS,
	"other"
    },
    {
	"next",
	(PyCFunction) CxRingNextPargs,
	METH_NOARGS,
	"next"
    },
    {
	"prev",
	(PyCFunction) CxRingPrevPargs,
	METH_NOARGS,
	"prev"
    },
    {NULL, NULL}
};

PyTypeObject CxtRing =
{
    PyObject_HEAD_INIT(NULL)
    0,			// int ob_size
    "C_Ring.C_Ring",	// char *tp_name
    sizeof(CxtRingObject),	// int tp_basicsize
    0,			// int tp_itemsize
    (destructor) CxpRingDelete,	// destructor tp_dealloc
    0,			// printfunc tp_print
    0,			// getattrfunc tp_getattr
    0,			// setattrfunc tp_setattr
    0,			// cmpfunc tp_compare
    0,			// reprfunc tp_repr
    0,			// PyNumberMethods *tp_as_number
    0,			// PySequenceMethods *tp_as_sequence
    0,			// PyMappingMethods *tp_as_mapping
    0,			// hashfunc tp_hash
    0,			// ternaryfunc tp_call
    0,			// reprfunc tp_str
    PyObject_GenericGetAttr,	// getattrofunc tp_getattro
    0,			// setattrofunc tp_setattro
    0,			// PyBufferProcs *tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, // long tp_flags
    "Ring(): Create the C portion of an ring.",	// char *tp_doc
    (traverseproc) CxpRingTraverse,	// traverseproc tp_traverse
    (inquiry) CxpRingClear,	// inquiry tp_clear
    0,			// richcmpfunc tp_richcompare
    0,			// long tp_weaklistoffset
    0,			// getiterfunc tp_iter
    0,			// iternextfunc tp_iternext
    CxpRingMethods,	// struct PyMethodDef *tp_methods
    0,			// struct PyMemberDef *tp_members
    0,			// struct PyGetSetDef *tp_getset
    0,			// struct _typeobject *tp_base
    0,			// PyObject *tp_dict
    0,			// descrgetfunc tp_descr_get
    0,			// descrsetfunc tp_descr_set
    0,			// long tp_dictoffset
    0,			// initproc tp_init
    0,			// allocfunc tp_alloc
    CxpRingNew,		// newfunc tp_new
    _PyObject_Del,	// freefunc tp_free
    0			// inquiry tp_is_gc
};

static PyMethodDef CxpRingFuncs[] =
{
    {NULL}
};

PyObject *CxgRingException;
PyObject *CxgRingValueError;

void
CxRingInit(void)
{
    PyObject *m;

    // Create new type.
    if (PyType_Ready(&CxtRing) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_Ring", CxpRingFuncs, "Ring extensions");
    Py_INCREF(&CxtRing);
    PyModule_AddObject(m, "C_Ring", (PyObject *) &CxtRing);

    // Create exception objects.
    CxgRingException = PyErr_NewException("C_Ring.Exception", CxgException,
					  NULL);
    Py_INCREF(CxgRingException);
    PyModule_AddObject(m, "Exception", CxgRingException);

    CxgRingValueError = PyErr_NewException("C_Ring.ValueError",
					   CxgRingException,
					   NULL);
    Py_INCREF(CxgRingValueError);
    PyModule_AddObject(m, "ValueError", CxgRingValueError);

    // Pre-compile Python code that is used for creating a wrapped ring.
    CxpRingNewCode = Py_CompileString("\
import crux.Ring\n\
ring = crux.Ring.Ring(edge, end)\n\
",
				      "<string>",
				      Py_file_input);
    CxmCheckPtr(CxpRingNewCode);
}

// End Ring.
//==============================================================================
