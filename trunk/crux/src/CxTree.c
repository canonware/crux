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

/******************************************************************************/
/* Begin Tree. */

static PyObject *
CxpTreeNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtTreeObject *self;

    self = (CxtTreeObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	self->tr = CxTrNew();
	CxTrAuxSet(self->tr, self);
	memset(self->aux, 0x0, sizeof(void *) * CxmTreeObjectAuxCount);
	self->GcCleared = false;
	retval = (PyObject *) self;
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	retval = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, retval, retval->ob_refcnt);
#endif
    return retval;
}

static int
CxpTreeTraverse(CxtTreeObject *self, visitproc visit, void *arg)
{
    int retval;
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
		retval = -1;
		goto RETURN;
	    }
	}
    }

    retval = 0;
    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif
    return retval;
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
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    CxpTreeClear(self);
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
    CxtTreeObject *retval;
    PyObject *globals, *locals, *obj;

    globals = PyEval_GetGlobals();
    if (globals == NULL)
    {
	retval = NULL;
	goto RETURN;
    }
    locals = Py_BuildValue("{}");
    if (locals == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    obj = PyEval_EvalCode((PyCodeObject *) CxpTreeNewCode,
			  globals,
			  locals);
    if (obj == NULL)
    {
	retval = NULL;
	goto RETURN;
    }
    Py_DECREF(obj);

    retval = (CxtTreeObject *) PyDict_GetItemString(locals, "tree");
    if (retval == NULL)
    {
	goto RETURN;
    }
    Py_INCREF(retval);

    RETURN:
    Py_DECREF(locals);
    return retval;
}

unsigned
CxTreeNtaxaGet(CxtTreeObject *self)
{
    return CxTrNtaxaGet(self->tr);
}

PyObject *
CxTreeNtaxaGetPargs(CxtTreeObject *self)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	retval = Py_BuildValue("i", CxTrNtaxaGet(self->tr));
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	retval = PyErr_NoMemory();
    }
    CxmXepEnd();

    return retval;
}

PyObject *
CxTreeNedgesCget(CxtTreeObject *self)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;

    CxmXepBegin();
    CxmXepTry
    {
	retval = Py_BuildValue("i", CxTrNedgesGet(self->tr));
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	retval = PyErr_NoMemory();
    }
    CxmXepEnd();

    return retval;
}

CxtNodeObject *
CxTreeBaseGet(CxtTreeObject *self)
{
    CxtNodeObject *retval;
    CxtTrNode base;

    base = CxTrBaseGet(self->tr);
    if (base == CxmTrNodeNone)
    {
	retval = NULL;
    }
    else
    {
	retval = (CxtNodeObject *) CxTrNodeAuxGet(self->tr, base);
    }

    return retval;
}

PyObject *
CxTreeBaseGetPargs(CxtTreeObject *self)
{
    PyObject *retval;
    CxtTrNode base;

    base = CxTrBaseGet(self->tr);
    if (base == CxmTrNodeNone)
    {
	Py_INCREF(Py_None);
	retval = Py_None;
    }
    else
    {
	retval = (PyObject *) CxTrNodeAuxGet(self->tr, base);
	Py_INCREF(retval);
    }

    return retval;
}

void
CxTreeBaseSet(CxtTreeObject *self, CxtNodeObject *aNode)
{
    CxtNodeObject *oldNode;
    CxtTrNode oldTrNode;

    /* Decref if clobbering an already-set base (but wait until after the
     * new base is set so that the nodes/edges/rings are always reachable. */
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
    }

    if (oldNode != NULL)
    {
	Py_DECREF(oldNode);
    }
}

PyObject *
CxTreeBaseSetPargs(CxtTreeObject *self, PyObject *args)
{
    PyObject *retval;
    CxtNodeObject *node;

    node = NULL;
    if (PyArg_ParseTuple(args, "|O!", &CxtNode, &node) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (node != NULL && node->tree != self)
    {
	CxError(CxgTreeValueError, "Node does not belong to this tree");
	retval = NULL;
	goto RETURN;
    }

    CxTreeBaseSet(self, node);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
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
	"_rnj",
	(PyCFunction) CxTreeRnj,
	METH_VARARGS,
	"_rnj"
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
	(PyCFunction) CxTreeTbr,
	METH_VARARGS,
	"tbr"
    },
    {
	"tbrNneighborsGet",
	(PyCFunction) CxTreeTbrNneighborsGet,
	METH_NOARGS,
	"tbrNneighborsGet"
    },
    {
	"tbrNeighborGet",
	(PyCFunction) CxTreeTbrNeighborGet,
	METH_VARARGS,
	"tbrNeighborGet"
    },
    {
	"mpPrepare",
	(PyCFunction) CxTreeMpPrepare,
	METH_VARARGS,
	"mpPrepare"
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
    0,			/* int ob_size */
    "C_Tree.C_Tree",	/* char *tp_name */
    sizeof(CxtTreeObject),	/* int tp_basicsize */
    0,			/* int tp_itemsize */
    (destructor) CxpTreeDelete,	/* destructor tp_dealloc */
    0,			/* printfunc tp_print */
    0,			/* getattrfunc tp_getattr */
    0,			/* setattrfunc tp_setattr */
    0,			/* cmpfunc tp_compare */
    0,			/* reprfunc tp_repr */
    0,			/* PyNumberMethods *tp_as_number */
    0,			/* PySequenceMethods *tp_as_sequence */
    0,			/* PyMappingMethods *tp_as_mapping */
    0,			/* hashfunc tp_hash */
    0,			/* ternaryfunc tp_call */
    0,			/* reprfunc tp_str */
    PyObject_GenericGetAttr,	/* getattrofunc tp_getattro */
    0,			/* setattrofunc tp_setattro */
    0,			/* PyBufferProcs *tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* long tp_flags */
    "Tree(): Create the C portion of a tree.",	/* char *tp_doc */
    (traverseproc) CxpTreeTraverse,	/* traverseproc tp_traverse */
    (inquiry) CxpTreeClear,	/* inquiry tp_clear */
    0,			/* richcmpfunc tp_richcompare */
    0,			/* long tp_weaklistoffset */
    0,			/* getiterfunc tp_iter */
    0,			/* iternextfunc tp_iternext */
    CxpTreeMethods,	/* struct PyMethodDef *tp_methods */
    0,			/* struct PyMemberDef *tp_members */
    0,			/* struct PyGetSetDef *tp_getset */
    0,			/* struct _typeobject *tp_base */
    0,			/* PyObject *tp_dict */
    0,			/* descrgetfunc tp_descr_get */
    0,			/* descrsetfunc tp_descr_set */
    0,			/* long tp_dictoffset */
    0,			/* initproc tp_init */
    0,			/* allocfunc tp_alloc */
    CxpTreeNew,		/* newfunc tp_new */
    _PyObject_Del,	/* freefunc tp_free */
    0			/* inquiry tp_is_gc */
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

    /* Create new type. */
    if (PyType_Ready(&CxtTree) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_Tree", CxpTreeFuncs, "Tree extensions");
    Py_INCREF(&CxtTree);
    PyModule_AddObject(m, "C_Tree", (PyObject *) &CxtTree);

    /* Create exception objects. */
    /* Exception. */
    CxgTreeException = PyErr_NewException("C_Tree.Exception", CxgException,
					  NULL);
    Py_INCREF(CxgTreeException);
    PyModule_AddObject(m, "Exception", CxgTreeException);

    /* ValueError. */
    CxgTreeValueError = PyErr_NewException("C_Tree.ValueError",
					   CxgTreeException,
					   NULL);
    Py_INCREF(CxgTreeValueError);
    PyModule_AddObject(m, "ValueError", CxgTreeValueError);

    /* TypeError. */
    CxgTreeTypeError = PyErr_NewException("C_Tree.TypeError",
					  CxgTreeException,
					  NULL);
    Py_INCREF(CxgTreeTypeError);
    PyModule_AddObject(m, "TypeError", CxgTreeTypeError);

    /* Pre-compile Python code that is used for creating a tree. */
    CxpTreeNewCode = Py_CompileString("\
import crux.Tree\n\
tree = crux.Tree.Tree()\n\
",
				      "<string>",
				      Py_file_input);
    CxmCheckPtr(CxpTreeNewCode);
}

/* End Tree. */
/******************************************************************************/
/* Begin Node. */

static PyObject *
CxpNodeNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtNodeObject *self;
    CxtTreeObject *tree;

    if (PyArg_ParseTuple(args, "O!", &CxtTree, &tree) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    self = (CxtNodeObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	Py_INCREF(tree);
	self->tree = tree;
	self->node = CxTrNodeNew(tree->tr);
	CxTrNodeAuxSet(tree->tr, self->node, self);
	memset(self->aux, 0x0, sizeof(void *) * CxmNodeObjectAuxCount);
	self->GcCleared = false;

	retval = (PyObject *) self;
    }
    CxmXepCatch(CxmXepOOM)
    {
	Py_DECREF(tree);
	CxmXepHandled();
	retval = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, retval, retval->ob_refcnt);
#endif
    return retval;
}

static int
CxpNodeTraverse(CxtNodeObject *self, visitproc visit, void *arg)
{
    int retval;
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
	    retval = -1;
	    goto RETURN;
	}

	/* Report all rings.  It is not good enough to simply report one, since
	 * Python's mark/sweep GC apparently keeps track of how many times each
	 * object is visited. */
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
		    retval = -1;
		    goto RETURN;
		}

		trCurRing = CxTrRingNextGet(self->tree->tr, trCurRing);
	    } while (trCurRing != trRing);
	}
    }

    retval = 0;
    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif
    return retval;
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
	/* Detach from rings. */
	trRing = CxTrNodeRingGet(self->tree->tr, self->node);
	if (trRing != CxmTrRingNone)
	{
	    PyObject *obj;

	    trCurRing = trRing;
	    do
	    {
		ring = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr,
							trCurRing);

		/* Get next ring before detaching. */
		trRing = trCurRing;
		trCurRing = CxTrRingNextGet(self->tree->tr, trCurRing);

		obj = CxEdgeDetach(ring->edge);
		CxmCheckPtr(obj);
		Py_DECREF(obj);
	    } while (trCurRing != trRing);
	}

	/* Detach from tree if tree base. */
	base = CxTrBaseGet(self->tree->tr);
	if (base == self->node)
	{
	    CxtNodeObject *node;

	    node = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, base);
	    CxTrBaseSet(self->tree->tr, CxmTrNodeNone);
	    Py_DECREF(node);
	}

	/* Delete node. */
	CxTrNodeDelete(self->tree->tr, self->node);

	/* Drop reference to tree. */
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
    CxtNodeObject *retval;
    PyObject *globals, *locals, *obj;

    globals = PyEval_GetGlobals();
    if (globals == NULL)
    {
	retval = NULL;
	goto RETURN;
    }
    locals = Py_BuildValue("{sO}", "tree", (PyObject *) aTree);
    if (locals == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    obj = PyEval_EvalCode((PyCodeObject *) CxpNodeNewCode,
			  globals,
			  locals);
    if (obj == NULL)
    {
	retval = NULL;
	goto RETURN;
    }
    Py_DECREF(obj);

    retval = (CxtNodeObject *) PyDict_GetItemString(locals, "node");
    if (retval == NULL)
    {
	goto RETURN;
    }
    Py_INCREF(retval);

    RETURN:
    Py_DECREF(locals);
    return retval;
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
    PyObject *retval;
    uint32_t taxonNum;

    taxonNum = CxTrNodeTaxonNumGet(self->tree->tr, self->node);

    if (taxonNum == CxmTrNodeTaxonNone)
    {
	Py_INCREF(Py_None);
	retval = Py_None;
    }
    else
    {
	retval = Py_BuildValue("i", taxonNum);
    }

    return retval;
}

void
CxNodeTaxonNumSet(CxtNodeObject *self, uint32_t aTaxonNum)
{
    CxTrNodeTaxonNumSet(self->tree->tr, self->node, aTaxonNum);
}

PyObject *
CxNodeTaxonNumSetPargs(CxtNodeObject *self, PyObject *args)
{
    PyObject *retval;
    uint32_t taxonNum;

    taxonNum = CxmTrNodeTaxonNone;
    if (PyArg_ParseTuple(args, "|i", &taxonNum) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    CxNodeTaxonNumSet(self, taxonNum);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

CxtRingObject *
CxNodeRing(CxtNodeObject *self)
{
    CxtRingObject *retval;
    CxtTrRing trRing;

    trRing = CxTrNodeRingGet(self->tree->tr, self->node);
    if (trRing != CxmTrRingNone)
    {
	retval = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRing);
    }
    else
    {
	retval = NULL;
    }

    return retval;
}

PyObject *
CxNodeRingPargs(CxtNodeObject *self)
{
    PyObject *retval;
    CxtRingObject *ring;
    CxtTrRing trRing;

    trRing = CxTrNodeRingGet(self->tree->tr, self->node);
    if (trRing != CxmTrRingNone)
    {
	ring = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRing);
	Py_INCREF(ring);
	retval = (PyObject *) ring;
    }
    else
    {
	Py_INCREF(Py_None);
	retval = Py_None;
    }

    return retval;
}

PyObject *
CxNodeDegree(CxtNodeObject *self)
{
    return Py_BuildValue("i", CxTrNodeDegree(self->tree->tr, self->node));
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
	(PyCFunction) CxNodeDegree,
	METH_NOARGS,
	"degree"
    },
    {NULL, NULL}
};

PyTypeObject CxtNode =
{
    PyObject_HEAD_INIT(NULL)
    0,			/* int ob_size */
    "C_Node.C_Node",	/* char *tp_name */
    sizeof(CxtNodeObject),	/* int tp_basicsize */
    0,			/* int tp_itemsize */
    (destructor) CxpNodeDelete,	/* destructor tp_dealloc */
    0,			/* printfunc tp_print */
    0,			/* getattrfunc tp_getattr */
    0,			/* setattrfunc tp_setattr */
    0,			/* cmpfunc tp_compare */
    0,			/* reprfunc tp_repr */
    0,			/* PyNumberMethods *tp_as_number */
    0,			/* PySequenceMethods *tp_as_sequence */
    0,			/* PyMappingMethods *tp_as_mapping */
    0,			/* hashfunc tp_hash */
    0,			/* ternaryfunc tp_call */
    0,			/* reprfunc tp_str */
    PyObject_GenericGetAttr,	/* getattrofunc tp_getattro */
    0,			/* setattrofunc tp_setattro */
    0,			/* PyBufferProcs *tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* long tp_flags */
    "Node(): Create the C portion of a node.",	/* char *tp_doc */
    (traverseproc) CxpNodeTraverse,	/* traverseproc tp_traverse */
    (inquiry) CxpNodeClear,	/* inquiry tp_clear */
    0,			/* richcmpfunc tp_richcompare */
    0,			/* long tp_weaklistoffset */
    0,			/* getiterfunc tp_iter */
    0,			/* iternextfunc tp_iternext */
    CxpNodeMethods,	/* struct PyMethodDef *tp_methods */
    0,			/* struct PyMemberDef *tp_members */
    0,			/* struct PyGetSetDef *tp_getset */
    0,			/* struct _typeobject *tp_base */
    0,			/* PyObject *tp_dict */
    0,			/* descrgetfunc tp_descr_get */
    0,			/* descrsetfunc tp_descr_set */
    0,			/* long tp_dictoffset */
    0,			/* initproc tp_init */
    0,			/* allocfunc tp_alloc */
    CxpNodeNew,		/* newfunc tp_new */
    _PyObject_Del,	/* freefunc tp_free */
    0			/* inquiry tp_is_gc */
};

static PyMethodDef CxpNodeFuncs[] =
{
    {NULL}
};

void
CxNodeInit(void)
{
    PyObject *m;

    /* Create new type. */
    if (PyType_Ready(&CxtNode) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_Node", CxpNodeFuncs, "Node extensions");
    Py_INCREF(&CxtNode);
    PyModule_AddObject(m, "C_Node", (PyObject *) &CxtNode);

    /* Pre-compile Python code that is used for creating a node. */
    CxpNodeNewCode = Py_CompileString("\
import crux.Node\n\
node = crux.Node.Node(tree)\n\
",
				      "<string>",
				      Py_file_input);
    CxmCheckPtr(CxpNodeNewCode);
}

/* End Node. */
/******************************************************************************/
/* Begin Edge. */

static PyObject *
CxpEdgeNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtEdgeObject *self;
    CxtTreeObject *tree;
    volatile uint32_t tryStage = 0;

    if (PyArg_ParseTuple(args, "O!", &CxtTree, &tree) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    self = (CxtEdgeObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    CxmXepBegin();
    CxmXepTry
    {
	/* Avoid traversing this object in GC until it is fully constructed, and
	 * its rings are fully constructed. */
	PyObject_GC_UnTrack(self);

	Py_INCREF(tree);
	self->tree = tree;
	self->edge = CxTrEdgeNew(tree->tr);
	CxTrEdgeAuxSet(tree->tr, self->edge, self);
	memset(self->aux, 0x0, sizeof(void *) * CxmEdgeObjectAuxCount);
	tryStage = 1;

	/* Create associated ring objects. */
	self->ringA = CxRingNew(self, 0);
	tryStage = 2;

	/* Avoid traversing ringA until after ringB has been constructed. */
	PyObject_GC_UnTrack(self->ringA);

	self->ringB = CxRingNew(self, 1);

	self->GcDetached = false;
	self->GcCleared = false;
	retval = (PyObject *) self;

	/* It's okay to traverse the edge and rings now. */
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
	retval = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, retval, retval->ob_refcnt);
#endif
    return retval;
}

static int
CxpEdgeTraverse(CxtEdgeObject *self, visitproc visit, void *arg)
{
    int retval;

#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Enter: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif

    if (self->GcCleared == false)
    {
	if (visit((PyObject *) self->tree, arg) < 0)
	{
	    retval = -1;
	    goto RETURN;
	}

	if (visit((PyObject *) self->ringA, arg) < 0)
	{
	    retval = -1;
	    goto RETURN;
	}

	if (visit((PyObject *) self->ringB, arg) < 0)
	{
	    retval = -1;
	    goto RETURN;
	}
    }

    retval = 0;
    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif
    return retval;
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
	/* Detach from nodes, if not already done. */
	if (self->GcDetached == false)
	{
	    self->GcDetached = true;
	    self->ringA->GcDetached = true;
	    self->ringB->GcDetached = true;

	    if (CxTrRingNodeGet(self->tree->tr, self->ringA->ring)
		!= CxmTrNodeNone)
	    {
		PyObject *obj = CxEdgeDetach(self);
		CxmCheckPtr(obj);
		Py_DECREF(obj);
	    }
	}

	/* Drop references to rings. */
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

    /* Delete the CxTrEdge (and associated CxTrRing objects).  The
     * CxtRingObject clear/delete code held onto a reference to this
     * CxtEdgeOjbect long enough to ensure that both CxtRingObject's
     * associated with this CxtEdgeObject have already been deallocated. */
    CxTrEdgeDelete(self->tree->tr, self->edge);

    /* Drop reference to tree. */
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
    CxtEdgeObject *retval;
    PyObject *globals, *locals, *obj;

    globals = PyEval_GetGlobals();
    if (globals == NULL)
    {
	retval = NULL;
	goto RETURN;
    }
    locals = Py_BuildValue("{sO}", "tree", (PyObject *) aTree);
    if (locals == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    obj = PyEval_EvalCode((PyCodeObject *) CxpEdgeNewCode,
			  globals,
			  locals);
    if (obj == NULL)
    {
	retval = NULL;
	goto RETURN;
    }
    Py_DECREF(obj);

    retval = (CxtEdgeObject *) PyDict_GetItemString(locals, "edge");
    if (retval == NULL)
    {
	goto RETURN;
    }
    Py_INCREF(retval);

    RETURN:
    Py_DECREF(locals);
    return retval;
}

PyObject *
CxEdgeTree(CxtEdgeObject *self)
{
    return Py_BuildValue("O", self->tree);
}

PyObject *
CxEdgeRingsGet(CxtEdgeObject *self)
{
    return Py_BuildValue("(OO)", self->ringA, self->ringB);
}

PyObject *
CxEdgeLengthGet(CxtEdgeObject *self)
{
    return Py_BuildValue("d", CxTrEdgeLengthGet(self->tree->tr, self->edge));
}

void
CxEdgeLengthSet(CxtEdgeObject *self, double aLength)
{
    CxTrEdgeLengthSet(self->tree->tr, self->edge, aLength);
}

PyObject *
CxEdgeLengthSetPargs(CxtEdgeObject *self, PyObject *args)
{
    PyObject *retval;
    double length;

    length = 0.0;
    if (PyArg_ParseTuple(args, "|d", &length) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    CxEdgeLengthSet(self, length);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

void
CxEdgeAttach(CxtEdgeObject *self, CxtNodeObject *aNodeA,
	     CxtNodeObject *aNodeB)
{
    /* Rings refer to nodes. */
    Py_INCREF(aNodeA);
    Py_INCREF(aNodeB);

    /* Nodes refer to rings.  In actuality, a node only refers to one element of
     * a ring, but doing all the "correct" reference management for the ring
     * would be hard (and slow), so instead, pretend that a node refers to each
     * ring element. */
    Py_INCREF(self->ringA);
    Py_INCREF(self->ringB);

    /* Attach. */
    CxTrEdgeAttach(self->tree->tr, self->edge, aNodeA->node, aNodeB->node);
}

PyObject *
CxEdgeAttachPargs(CxtEdgeObject *self, PyObject *args)
{
    PyObject *retval;
    CxtNodeObject *nodeA, *nodeB;

    if (PyArg_ParseTuple(args, "O!O!", &CxtNode, &nodeA,
			 &CxtNode, &nodeB) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    /* Make sure that the edge is currently detached. */
    if (CxTrRingNodeGet(self->tree->tr, self->ringA->ring) != CxmTrNodeNone)
    {
	CxError(CxgEdgeValueError, "Edge is already attached");
	retval = NULL;
	goto RETURN;
    }
    if (nodeA->tree != self->tree || nodeB->tree != self->tree)
    {
	CxError(CxgEdgeValueError, "Node does not belong to this tree");
	retval = NULL;
	goto RETURN;
    }

    CxEdgeAttach(self, nodeA, nodeB);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

PyObject *
CxEdgeDetach(CxtEdgeObject *self)
{
    PyObject *retval;
    CxtTrNode trNodeA, trNodeB;
    CxtNodeObject *nodeA, *nodeB;
    CxtTrRing trRingA, trRingB;

    /* Make sure that the edge is currently attached. */
    if (CxTrRingNodeGet(self->tree->tr, self->ringA->ring) == CxmTrNodeNone)
    {
	CxError(CxgEdgeValueError, "Edge is already detached");
	retval = NULL;
	goto RETURN;
    }

    trRingA = CxTrEdgeRingGet(self->tree->tr, self->edge, 0);
    CxmAssert(trRingA != CxmTrRingNone);
    trRingB = CxTrEdgeRingGet(self->tree->tr, self->edge, 1);

    trNodeA = CxTrRingNodeGet(self->tree->tr, trRingA);
    nodeA = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, trNodeA);

    trNodeB = CxTrRingNodeGet(self->tree->tr, trRingB);
    nodeB = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, trNodeB);

    /* Detach. */
    CxTrEdgeDetach(self->tree->tr, self->edge);

    /* Rings refer to nodes. */
    Py_DECREF(nodeA);
    Py_DECREF(nodeB);

    /* Nodes refer to rings. */
    Py_DECREF(self->ringA);
    Py_DECREF(self->ringB);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
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
	(PyCFunction) CxEdgeRingsGet,
	METH_NOARGS,
	"rings"
    },
    {
	"lengthGet",
	(PyCFunction) CxEdgeLengthGet,
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
	(PyCFunction) CxEdgeDetach,
	METH_NOARGS,
	"detach"
    },
    {NULL, NULL}
};

PyTypeObject CxtEdge =
{
    PyObject_HEAD_INIT(NULL)
    0,			/* int ob_size */
    "C_Edge.C_Edge",	/* char *tp_name */
    sizeof(CxtEdgeObject),	/* int tp_basicsize */
    0,			/* int tp_itemsize */
    (destructor) CxpEdgeDelete,	/* destructor tp_dealloc */
    0,			/* printfunc tp_print */
    0,			/* getattrfunc tp_getattr */
    0,			/* setattrfunc tp_setattr */
    0,			/* cmpfunc tp_compare */
    0,			/* reprfunc tp_repr */
    0,			/* PyNumberMethods *tp_as_number */
    0,			/* PySequenceMethods *tp_as_sequence */
    0,			/* PyMappingMethods *tp_as_mapping */
    0,			/* hashfunc tp_hash */
    0,			/* ternaryfunc tp_call */
    0,			/* reprfunc tp_str */
    PyObject_GenericGetAttr,	/* getattrofunc tp_getattro */
    0,			/* setattrofunc tp_setattro */
    0,			/* PyBufferProcs *tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* long tp_flags */
    "Edge(): Create the C portion of an edge.",	/* char *tp_doc */
    (traverseproc) CxpEdgeTraverse,	/* traverseproc tp_traverse */
    (inquiry) CxpEdgeClear,	/* inquiry tp_clear */
    0,			/* richcmpfunc tp_richcompare */
    0,			/* long tp_weaklistoffset */
    0,			/* getiterfunc tp_iter */
    0,			/* iternextfunc tp_iternext */
    CxpEdgeMethods,	/* struct PyMethodDef *tp_methods */
    0,			/* struct PyMemberDef *tp_members */
    0,			/* struct PyGetSetDef *tp_getset */
    0,			/* struct _typeobject *tp_base */
    0,			/* PyObject *tp_dict */
    0,			/* descrgetfunc tp_descr_get */
    0,			/* descrsetfunc tp_descr_set */
    0,			/* long tp_dictoffset */
    0,			/* initproc tp_init */
    0,			/* allocfunc tp_alloc */
    CxpEdgeNew,		/* newfunc tp_new */
    _PyObject_Del,	/* freefunc tp_free */
    0			/* inquiry tp_is_gc */
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

    /* Create new type. */
    if (PyType_Ready(&CxtEdge) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_Edge", CxpEdgeFuncs, "Edge extensions");
    Py_INCREF(&CxtEdge);
    PyModule_AddObject(m, "C_Edge", (PyObject *) &CxtEdge);

    /* Create exception objects. */
    CxgEdgeException = PyErr_NewException("C_Edge.Exception", CxgException,
					  NULL);
    Py_INCREF(CxgEdgeException);
    PyModule_AddObject(m, "Exception", CxgEdgeException);

    CxgEdgeValueError = PyErr_NewException("C_Edge.ValueError",
					   CxgEdgeException,
					   NULL);
    Py_INCREF(CxgEdgeValueError);
    PyModule_AddObject(m, "ValueError", CxgEdgeValueError);

    /* Pre-compile Python code that is used for creating a wrapped edge. */
    CxpEdgeNewCode = Py_CompileString("\
import crux.Edge\n\
edge = crux.Edge.Edge(tree)\n\
",
				      "<string>",
				      Py_file_input);
    CxmCheckPtr(CxpEdgeNewCode);
}

/* End Edge. */
/******************************************************************************/
/* Begin Ring. */

static PyObject *
CxpRingNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtRingObject *self;
    CxtEdgeObject *edge;
    uint32_t end;

    if (PyArg_ParseTuple(args, "O!i", &CxtEdge, &edge, &end) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    self = (CxtRingObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    if (end != 0 && end != 1)
    {
	type->tp_free((PyObject *) self);

	CxError(CxgRingValueError, "End must be 0 or 1");
	retval = NULL;
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
	retval = NULL;
	goto RETURN;
    }
    CxTrRingAuxSet(self->tree->tr, self->ring, self);
    memset(self->aux, 0x0, sizeof(void *) * CxmRingObjectAuxCount);

    self->GcDetached = false;
    self->GcCleared = false;

    retval = (PyObject *) self;

    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, retval, retval->ob_refcnt);
#endif
    return retval;
}

static int
CxpRingTraverse(CxtRingObject *self, visitproc visit, void *arg)
{
    int retval;
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
	    retval = -1;
	    goto RETURN;
	}

	trNode = CxTrRingNodeGet(self->tree->tr, self->ring);
	if (trNode != CxmTrNodeNone)
	{
	    CxtNodeObject *node;

	    node = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, trNode);
	    if (visit((PyObject *) node, arg) < 0)
	    {
		retval = -1;
		goto RETURN;
	    }
	}

	/* Report edge. */
	if (visit((PyObject *) self->edge, arg) < 0)
	{
	    retval = -1;
	    goto RETURN;
	}

	/* Report next ring element. */
	trRing = CxTrRingNextGet(self->tree->tr, self->ring);
	ring = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRing);
	if (visit((PyObject *) ring, arg) < 0)
	{
	    retval = -1;
	    goto RETURN;
	}

	/* Report previous ring element. */
	trRing = CxTrRingPrevGet(self->tree->tr, self->ring);
	ring = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRing);
	if (visit((PyObject *) ring, arg) < 0)
	{
	    retval = -1;
	    goto RETURN;
	}
    }

    retval = 0;
    RETURN:
#ifdef CxmTreeGCVerbose
    fprintf(stderr, "%s:%d:%s() Leave: %p (%d)\n",
 	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
#endif
    return retval;
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
	/* Detach from nodes, if not already done. */
	if (self->GcDetached == false)
	{
	    self->edge->GcDetached = true;
	    self->edge->ringA->GcDetached = true;
	    self->edge->ringB->GcDetached = true;

	    if (CxTrRingNodeGet(self->tree->tr, self->ring) != CxmTrNodeNone)
	    {
		PyObject *obj = CxEdgeDetach(self->edge);
		CxmCheckPtr(obj);
		Py_DECREF(obj);
	    }
	}

	/* Drop reference to tree. */
	Py_DECREF(self->tree);

	/* Drop references to self. */
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

    /* Drop the reference to the associated edge, now that there is no chance of
     * accessing it again. */
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
    CxtRingObject *retval;
    PyObject *globals, *locals, *obj;

    globals = PyEval_GetGlobals();
    if (globals == NULL)
    {
	retval = NULL;
	goto RETURN;
    }
    locals = Py_BuildValue("{sOsi}",
			   "edge", (PyObject *) aEdge,
			   "end", aEnd);
    if (locals == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    obj = PyEval_EvalCode((PyCodeObject *) CxpRingNewCode,
			  globals,
			  locals);
    if (obj == NULL)
    {
	retval = NULL;
	goto RETURN;
    }
    Py_DECREF(obj);

    retval = (CxtRingObject *) PyDict_GetItemString(locals, "ring");
    if (retval == NULL)
    {
	goto RETURN;
    }
    Py_INCREF(retval);

    RETURN:
    Py_DECREF(locals);
    return retval;
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
    CxtNodeObject *retval;
    CxtTrNode node;

    node = CxTrRingNodeGet(self->tree->tr, self->ring);
    if (node != CxmTrNodeNone)
    {
	retval = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, node);
    }
    else
    {
	retval = NULL;
    }

    return retval;
}

PyObject *
CxRingNodePargs(CxtRingObject *self)
{
    PyObject *retval;
    CxtTrNode node;

    node = CxTrRingNodeGet(self->tree->tr, self->ring);
    if (node != CxmTrNodeNone)
    {
	retval = (PyObject *) CxTrNodeAuxGet(self->tree->tr, node);
	Py_INCREF(retval);
    }
    else
    {
	Py_INCREF(Py_None);
	retval = Py_None;
    }

    return retval;
}

CxtEdgeObject *
CxRingEdge(CxtRingObject *self)
{
    CxtEdgeObject *retval;
    CxtTrEdge edge;

    edge = CxTrRingEdgeGet(self->tree->tr, self->ring);
    retval = (CxtEdgeObject *) CxTrEdgeAuxGet(self->tree->tr, edge);

    return retval;
}

PyObject *
CxRingEdgePargs(CxtRingObject *self)
{
    PyObject *retval;
    CxtTrEdge edge;

    edge = CxTrRingEdgeGet(self->tree->tr, self->ring);
    retval = (PyObject *) CxTrEdgeAuxGet(self->tree->tr, edge);
    Py_INCREF(retval);

    return retval;
}

CxtRingObject *
CxRingOther(CxtRingObject *self)
{
    CxtRingObject *retval;
    CxtTrRing other;

    other = CxTrRingOtherGet(self->tree->tr, self->ring);
    retval = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, other);

    return retval;
}

PyObject *
CxRingOtherPargs(CxtRingObject *self)
{
    PyObject *retval;
    CxtTrRing other;

    other = CxTrRingOtherGet(self->tree->tr, self->ring);
    retval = (PyObject *) CxTrRingAuxGet(self->tree->tr, other);
    Py_INCREF(retval);

    return retval;
}

CxtRingObject *
CxRingNext(CxtRingObject *self)
{
    CxtRingObject *retval;
    CxtTrRing nextRing;

    nextRing = CxTrRingNextGet(self->tree->tr, self->ring);
    retval = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, nextRing);

    return retval;
}

PyObject *
CxRingNextPargs(CxtRingObject *self)
{
    PyObject *retval;
    CxtTrRing nextRing;

    nextRing = CxTrRingNextGet(self->tree->tr, self->ring);
    retval = (PyObject *) CxTrRingAuxGet(self->tree->tr, nextRing);
    Py_INCREF(retval);

    return retval;
}

CxtRingObject *
CxRingPrev(CxtRingObject *self)
{
    CxtRingObject *retval;
    CxtTrRing prevRing;

    prevRing = CxTrRingPrevGet(self->tree->tr, self->ring);
    retval = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, prevRing);

    return retval;
}

PyObject *
CxRingPrevPargs(CxtRingObject *self)
{
    PyObject *retval;
    CxtTrRing prevRing;

    prevRing = CxTrRingPrevGet(self->tree->tr, self->ring);
    retval = (PyObject *) CxTrRingAuxGet(self->tree->tr, prevRing);
    Py_INCREF(retval);

    return retval;
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
    0,			/* int ob_size */
    "C_Ring.C_Ring",	/* char *tp_name */
    sizeof(CxtRingObject),	/* int tp_basicsize */
    0,			/* int tp_itemsize */
    (destructor) CxpRingDelete,	/* destructor tp_dealloc */
    0,			/* printfunc tp_print */
    0,			/* getattrfunc tp_getattr */
    0,			/* setattrfunc tp_setattr */
    0,			/* cmpfunc tp_compare */
    0,			/* reprfunc tp_repr */
    0,			/* PyNumberMethods *tp_as_number */
    0,			/* PySequenceMethods *tp_as_sequence */
    0,			/* PyMappingMethods *tp_as_mapping */
    0,			/* hashfunc tp_hash */
    0,			/* ternaryfunc tp_call */
    0,			/* reprfunc tp_str */
    PyObject_GenericGetAttr,	/* getattrofunc tp_getattro */
    0,			/* setattrofunc tp_setattro */
    0,			/* PyBufferProcs *tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* long tp_flags */
    "Ring(): Create the C portion of an ring.",	/* char *tp_doc */
    (traverseproc) CxpRingTraverse,	/* traverseproc tp_traverse */
    (inquiry) CxpRingClear,	/* inquiry tp_clear */
    0,			/* richcmpfunc tp_richcompare */
    0,			/* long tp_weaklistoffset */
    0,			/* getiterfunc tp_iter */
    0,			/* iternextfunc tp_iternext */
    CxpRingMethods,	/* struct PyMethodDef *tp_methods */
    0,			/* struct PyMemberDef *tp_members */
    0,			/* struct PyGetSetDef *tp_getset */
    0,			/* struct _typeobject *tp_base */
    0,			/* PyObject *tp_dict */
    0,			/* descrgetfunc tp_descr_get */
    0,			/* descrsetfunc tp_descr_set */
    0,			/* long tp_dictoffset */
    0,			/* initproc tp_init */
    0,			/* allocfunc tp_alloc */
    CxpRingNew,		/* newfunc tp_new */
    _PyObject_Del,	/* freefunc tp_free */
    0			/* inquiry tp_is_gc */
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

    /* Create new type. */
    if (PyType_Ready(&CxtRing) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_Ring", CxpRingFuncs, "Ring extensions");
    Py_INCREF(&CxtRing);
    PyModule_AddObject(m, "C_Ring", (PyObject *) &CxtRing);

    /* Create exception objects. */
    CxgRingException = PyErr_NewException("C_Ring.Exception", CxgException,
					  NULL);
    Py_INCREF(CxgRingException);
    PyModule_AddObject(m, "Exception", CxgRingException);

    CxgRingValueError = PyErr_NewException("C_Ring.ValueError",
					   CxgRingException,
					   NULL);
    Py_INCREF(CxgRingValueError);
    PyModule_AddObject(m, "ValueError", CxgRingValueError);

    /* Pre-compile Python code that is used for creating a wrapped ring. */
    CxpRingNewCode = Py_CompileString("\
import crux.Ring\n\
ring = crux.Ring.Ring(edge, end)\n\
",
				      "<string>",
				      Py_file_input);
    CxmCheckPtr(CxpRingNewCode);
}

/* End Ring. */
/******************************************************************************/
