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

#if (0) // XXX Remove this code.
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

#include <math.h>

static PyTypeObject CxtTree;
static PyTypeObject CxtNode;
static PyTypeObject CxtEdge;
static PyTypeObject CxtRing;

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
	retval = (PyObject *) self;
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	retval = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
    return retval;
}

static int
CxpTreeTraverse(CxtTreeObject *self, visitproc visit, void *arg)
{
    int retval;
    CxtTrNode base;

//    fprintf(stderr, "%s:%d:%s(%p) Enter (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
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

    retval = 0;
    RETURN:
//    fprintf(stderr, "%s:%d:%s(%p) Leave (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    return retval;
}

static int
CxpTreeClear(CxtTreeObject *self)
{
    CxtTrNode base;

//    fprintf(stderr, "%s:%d:%s(%p) Enter (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    base = CxTrBaseGet(self->tr);
    if (base != CxmTrNodeNone)
    {
	CxtNodeObject *node;

	node = (CxtNodeObject *) CxTrNodeAuxGet(self->tr, base);
	CxTrBaseSet(self->tr, CxmTrNodeNone);
	Py_DECREF(node);
    }

//    fprintf(stderr, "%s:%d:%s(%p) Leave (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    return 0;
}

static void
CxpTreeDelete(CxtTreeObject *self)
{
//    fprintf(stderr, "%s:%d:%s(%p) Enter (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    CxpTreeClear(self);
    CxTrDelete(self->tr);
    self->ob_type->tp_free((PyObject*) self);
//    fprintf(stderr, "%s:%d:%s(%p) Leave (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
}

static PyObject *CxpTreeNewCode;

CxtTreeObject *
CxTreeNew(void)
{
    CxtTreeObject *retval;
    PyObject *globals, *locals, *obj;

    globals = PyEval_GetGlobals();
    locals = Py_BuildValue("{}");


    obj = PyEval_EvalCode((PyCodeObject *) CxpTreeNewCode,
			  globals,
			  locals);
    CxmCheckPtr(obj);
    Py_DECREF(obj);

    retval = (CxtTreeObject *) PyDict_GetItemString(locals, "tree");
    CxmCheckPtr(retval);
    Py_INCREF(retval);

    Py_DECREF(locals);

    return retval;
}

PyObject *
CxTreeNtaxaGet(CxtTreeObject *self)
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

PyObject *
CxTreeBaseGet(CxtTreeObject *self)
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
	PyErr_SetString(CxgTreeValueError, "Node does not belong to this tree");
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
	(PyCFunction) CxTreeNtaxaGet,
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
	(PyCFunction) CxTreeBaseGet,
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
	"_nj",
	(PyCFunction) CxTreeNj,
	METH_VARARGS,
	"_nj"
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

static PyTypeObject CxtTree =
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
    CxgTreeException = PyErr_NewException("C_Tree.Exception", CxgException,
					  NULL);
    Py_INCREF(CxgTreeException);
    PyModule_AddObject(m, "Exception", CxgTreeException);

    CxgTreeValueError = PyErr_NewException("C_Tree.ValueError",
					   CxgTreeException,
					   NULL);
    Py_INCREF(CxgTreeValueError);
    PyModule_AddObject(m, "ValueError", CxgTreeValueError);

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

	self->valid = true;

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
    return retval;
}

static int
CxpNodeTraverse(CxtNodeObject *self, visitproc visit, void *arg)
{
    int retval;
    CxtTrRing trRing;
    CxtRingObject *ring;

//    fprintf(stderr, "%s:%d:%s(%p) Enter (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    if (self->valid)
    {
	if (visit((PyObject *) self->tree, arg) < 0)
	{
	    retval = -1;
	    goto RETURN;
	}

	trRing = CxTrNodeRingGet(self->tree->tr, self->node);
	if (trRing != CxmTrRingNone)
	{
	    ring = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRing);
	    if (visit((PyObject *) ring, arg) < 0)
	    {
		retval = -1;
		goto RETURN;
	    }
	}
    }

    retval = 0;
    RETURN:
//    fprintf(stderr, "%s:%d:%s(%p) Leave (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    return retval;
}

static int
CxpNodeClear(CxtNodeObject *self)
{
    CxtTrRing trRing;
    CxtRingObject *ring;
    CxtEdgeObject *edge;
    PyObject *obj;

//    fprintf(stderr, "%s:%d:%s(%p) Enter (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    if (self->valid)
    {
	while (true)
	{
	    trRing = CxTrNodeRingGet(self->tree->tr, self->node);
	    if (trRing == CxmTrRingNone)
	    {
		break;
	    }
	    ring = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRing);
	    edge = (CxtEdgeObject *) CxRingEdge(ring);
	    obj = CxEdgeDetach(edge);
	    Py_DECREF(obj);
	}

	if (CxTrBaseGet(self->tree->tr) == self->node)
	{
	    CxTrBaseSet(self->tree->tr, CxmTrNodeNone);
	    Py_DECREF(self);
	}
	CxTrNodeDelete(self->tree->tr, self->node);
	Py_DECREF(self->tree);

	self->valid = false;
    }

//    fprintf(stderr, "%s:%d:%s(%p) Leave (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    return 0;
}

static void
CxpNodeDelete(CxtNodeObject *self)
{
//    fprintf(stderr, "%s:%d:%s(%p) Enter (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    CxpNodeClear(self);
    self->ob_type->tp_free((PyObject*) self);
//    fprintf(stderr, "%s:%d:%s(%p) Leave (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
}

static PyObject *CxpNodeNewCode;

CxtNodeObject *
CxNodeNew(CxtTreeObject *aTree)
{
    CxtNodeObject *retval;
    PyObject *globals, *locals, *obj;

    globals = PyEval_GetGlobals();
    locals = Py_BuildValue("{sO}", "tree", (PyObject *) aTree);

    obj = PyEval_EvalCode((PyCodeObject *) CxpNodeNewCode,
			  globals,
			  locals);
    CxmCheckPtr(obj);
    Py_DECREF(obj);

    retval = (CxtNodeObject *) PyDict_GetItemString(locals, "node");
    CxmCheckPtr(retval);
    Py_INCREF(retval);

    Py_DECREF(locals);

    return retval;
}

PyObject *
CxNodeTree(CxtNodeObject *self)
{
    return Py_BuildValue("O", self->tree);
}

PyObject *
CxNodeTaxonNumGet(CxtNodeObject *self)
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

PyObject *
CxNodeRing(CxtNodeObject *self)
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
	(PyCFunction) CxNodeTaxonNumGet,
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
	(PyCFunction) CxNodeRing,
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

static PyTypeObject CxtNode =
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
    CxtRingObject *ringA, *ringB;
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
	Py_INCREF(tree);
	self->tree = tree;
	self->edge = CxTrEdgeNew(tree->tr);
	CxTrEdgeAuxSet(tree->tr, self->edge, self);
	tryStage = 1;

	/* Create associated ring objects. */
	ringA = CxRingNew(self, 0);
	//Py_INCREF(ringA);
	tryStage = 2;

	ringB = CxRingNew(self, 1);
	//Py_INCREF(ringB);

	self->valid = true;

	retval = (PyObject *) self;
    }
    CxmXepCatch(CxmXepOOM)
    {
	switch (tryStage)
	{
	    case 2:
	    {
		Py_DECREF(ringA);
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
    return retval;
}

static int
CxpEdgeTraverse(CxtEdgeObject *self, visitproc visit, void *arg)
{
    int retval;
    uint32_t i;
    CxtTrRing trRing;
    CxtRingObject *ring;

//    fprintf(stderr, "%s:%d:%s(%p) Enter (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    if (self->valid)
    {
	if (visit((PyObject *) self->tree, arg) < 0)
	{
	    retval = -1;
	    goto RETURN;
	}

	for (i = 0; i < 2; i++)
	{
	    trRing = CxTrEdgeRingGet(self->tree->tr, self->edge, i);
	    ring = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRing);
	    if (visit((PyObject *) ring, arg) < 0)
	    {
		retval = -1;
		goto RETURN;
	    }
	}
    }

    retval = 0;
    RETURN:
//    fprintf(stderr, "%s:%d:%s(%p) Leave (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    return retval;
}

/* Clear edge and associated rings all at once.  This can be called either when
 * an edge or a ring is cleared. */
static int
CxpEdgeRingClear(CxtEdgeObject *aEdge, CxtRingObject *aRingA,
		 CxtRingObject *aRingB)
{
    CxtTrRing trRing;
    PyObject *obj;

//    fprintf(stderr, "%s:%d:%s(%p, %p, %p) Enter (%d, %d, %d)\n",
//	    __FILE__, __LINE__, __func__, aEdge, aRingA, aRingB,
//	    aEdge->ob_refcnt, aRingA->ob_refcnt, aRingB->ob_refcnt);

    /* Detach from nodes, if necessary. */
    trRing = CxTrEdgeRingGet(aEdge->tree->tr, aEdge->edge, 0);
    if (CxTrRingNodeGet(aEdge->tree->tr, trRing) != CxmTrNodeNone)
    {
	obj = CxEdgeDetach(aEdge);
	Py_DECREF(obj);
    }

    /* Avoid running this block of code twice.  Do this before decref's, in
     * order to avoid recursion. */
    aEdge->valid = false;

    /* Get rid of circular references between edge and rings. */
    Py_DECREF(aEdge);
    Py_DECREF(aEdge);
    Py_DECREF(aRingA);
    Py_DECREF(aRingB);

    /* Get rid of self references of rings. */
    Py_DECREF(aRingA);
    Py_DECREF(aRingA);
    Py_DECREF(aRingB);
    Py_DECREF(aRingB);

    /* Deallocate the CxtTrEdge. */
    CxTrEdgeDelete(aEdge->tree->tr, aEdge->edge);

    /* Get rid of references from edge and rings to tree. */
    Py_DECREF(aEdge->tree);
    Py_DECREF(aEdge->tree);
    Py_DECREF(aEdge->tree);

//    fprintf(stderr, "%s:%d:%s(%p, %p, %p) Leave (%d, %d, %d)\n",
//	    __FILE__, __LINE__, __func__, aEdge, aRingA, aRingB,
//	    aEdge->ob_refcnt, aRingA->ob_refcnt, aRingB->ob_refcnt);
    return 0;
}

static int
CxpEdgeClear(CxtEdgeObject *self)
{
    int retval;

//    fprintf(stderr, "%s:%d:%s(%p) Enter (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    if (self->valid)
    {
	CxtTrRing trRingA, trRingB;
	CxtRingObject *ringA, *ringB;

	trRingA = CxTrEdgeRingGet(self->tree->tr, self->edge, 0);
	ringA = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRingA);

	trRingB = CxTrEdgeRingGet(self->tree->tr, self->edge, 1);
	ringB = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRingB);

	retval = CxpEdgeRingClear(self, ringA, ringB);
    }
    else
    {
	retval = 0;
    }

//    fprintf(stderr, "%s:%d:%s(%p) Leave (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    return retval;
}

static void
CxpEdgeDelete(CxtEdgeObject *self)
{
//    fprintf(stderr, "%s:%d:%s(%p) Enter (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    CxpEdgeClear(self);
    self->ob_type->tp_free((PyObject*) self);
//    fprintf(stderr, "%s:%d:%s(%p) Leave (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
}

static PyObject *CxpEdgeNewCode;

CxtEdgeObject *
CxEdgeNew(CxtTreeObject *aTree)
{
    CxtEdgeObject *retval;
    PyObject *globals, *locals, *obj;

    globals = PyEval_GetGlobals();
    locals = Py_BuildValue("{sO}", "tree", (PyObject *) aTree);

    obj = PyEval_EvalCode((PyCodeObject *) CxpEdgeNewCode,
			  globals,
			  locals);
    CxmCheckPtr(obj);
    Py_DECREF(obj);

    retval = (CxtEdgeObject *) PyDict_GetItemString(locals, "edge");
    CxmCheckPtr(retval);
    Py_INCREF(retval);

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
    CxtTrRing trRingA, trRingB;
    CxtRingObject *ringA, *ringB;

    trRingA = CxTrEdgeRingGet(self->tree->tr, self->edge, 0);
    ringA = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRingA);

    trRingB = CxTrEdgeRingGet(self->tree->tr, self->edge, 1);
    ringB = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRingB);

    return Py_BuildValue("(OO)", ringA, ringB);
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
    if (length < 0.0)
    {
	PyErr_SetString(CxgEdgeValueError, "Length must be non-negative");
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
    CxtTrRing trRingA, trRingB;
    CxtRingObject *ringA, *ringB;

    trRingA = CxTrEdgeRingGet(self->tree->tr, self->edge, 0);
    ringA = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRingA);

    trRingB = CxTrEdgeRingGet(self->tree->tr, self->edge, 1);
    ringB = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRingB);

    /* Rings refer to nodes. */
    Py_INCREF(aNodeA);
    Py_INCREF(aNodeB);
    /* Nodes refer to rings.  In actuality, a node only refers to one element of
     * a ring, but doing all the "correct" reference management for the ring
     * would be hard (and slow), so instead, pretend that a node refers to each
     * ring element.  However, the node only needs to report one element of the
     * ring during object traversal, and the ring elements take care of the rest
     * by reporting the next element in the ring. */
    Py_INCREF(ringA);
    Py_INCREF(ringB);

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
    if (nodeA->tree != self->tree || nodeB->tree != self->tree)
    {
	PyErr_SetString(CxgEdgeValueError, "Node does not belong to this tree");
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
    CxtTrRing trRingA;

    trRingA = CxTrEdgeRingGet(self->tree->tr, self->edge, 0);
    if (trRingA != CxmTrRingNone)
    {
	CxtTrNode trNodeA, trNodeB;
	CxtNodeObject *nodeA, *nodeB;
	CxtTrRing trRingB;
	CxtRingObject *ringA, *ringB;

	ringA = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRingA);

	trRingB = CxTrEdgeRingGet(self->tree->tr, self->edge, 1);
	ringB = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRingB);

	trNodeA = CxTrRingNodeGet(self->tree->tr, trRingA);
	nodeA = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, trNodeA);

	trNodeB = CxTrRingNodeGet(self->tree->tr, trRingB);
	nodeB = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, trNodeB);	

	CxTrEdgeDetach(self->tree->tr, self->edge);
	/* Rings refer to nodes. */
	Py_DECREF(nodeA);
	Py_DECREF(nodeB);
	/* Nodes refer to rings. */
	Py_DECREF(ringA);
	Py_DECREF(ringB);
    }

    Py_INCREF(Py_None);
    return Py_None;
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

static PyTypeObject CxtEdge =
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

	PyErr_SetString(CxgRingValueError, "End must be 0 or 1");
	retval = NULL;
	goto RETURN;
    }
    
    Py_INCREF(edge->tree);
    Py_INCREF(edge);
    /* The ring linkage causes a new ring to ref itself (twice). */
    Py_INCREF(self);
    Py_INCREF(self);
    self->tree = edge->tree;
    self->edge = edge;
    self->ring = CxTrEdgeRingGet(edge->tree->tr, edge->edge, end);
    if (CxTrRingAuxGet(self->tree->tr, self->ring) != NULL)
    {
	Py_DECREF(self);
	Py_DECREF(self);
	Py_DECREF(edge);
	Py_DECREF(edge->tree);
	type->tp_free((PyObject *) self);

	PyErr_SetString(CxgRingValueError,
			"Internal error (aux should be NULL)");
	retval = NULL;
	goto RETURN;
    }
    CxTrRingAuxSet(self->tree->tr, self->ring, self);

    retval = (PyObject *) self;

    RETURN:
    return retval;
}

static int
CxpRingTraverse(CxtRingObject *self, visitproc visit, void *arg)
{
    int retval;
    CxtTrNode trNode;
    CxtTrEdge trEdge;
    CxtTrRing trRing;
    CxtEdgeObject *edge;
    CxtRingObject *ring;

//    fprintf(stderr, "%s:%d:%s(%p) Enter (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    if (self->edge->valid)
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

	trEdge = CxTrRingEdgeGet(self->tree->tr, self->ring);
	edge = (CxtEdgeObject *) CxTrEdgeAuxGet(self->tree->tr, trEdge);
	if (visit((PyObject *) edge, arg) < 0)
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
//    fprintf(stderr, "%s:%d:%s(%p) Leave (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    return retval;
}

static int
CxpRingClear(CxtRingObject *self)
{
    int retval;

//    fprintf(stderr, "%s:%d:%s(%p) Enter (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    if (self->edge->valid)
    {
	CxtTrRing trRingA, trRingB;
	CxtRingObject *ringA, *ringB;

	trRingA = CxTrEdgeRingGet(self->tree->tr, self->edge->edge, 0);
	ringA = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRingA);

	trRingB = CxTrEdgeRingGet(self->tree->tr, self->edge->edge, 1);
	ringB = (CxtRingObject *) CxTrRingAuxGet(self->tree->tr, trRingB);

	retval = CxpEdgeRingClear(self->edge, ringA, ringB);
    }
    else
    {
	retval = 0;
    }

//    fprintf(stderr, "%s:%d:%s(%p) Leave (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    return retval;
}

static void
CxpRingDelete(CxtRingObject *self)
{
//    fprintf(stderr, "%s:%d:%s(%p) Enter (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
    CxpRingClear(self);
    self->ob_type->tp_free((PyObject*) self);
//    fprintf(stderr, "%s:%d:%s(%p) Leave (%d)\n",
//	    __FILE__, __LINE__, __func__, self, self->ob_refcnt);
}

static PyObject *CxpRingNewCode;

CxtRingObject *
CxRingNew(CxtEdgeObject *aEdge, uint32_t aEnd)
{
    CxtRingObject *retval;
    PyObject *globals, *locals, *obj;

    globals = PyEval_GetGlobals();
    locals = Py_BuildValue("{sOsi}",
			   "edge", (PyObject *) aEdge,
			   "end", aEnd);

    obj = PyEval_EvalCode((PyCodeObject *) CxpRingNewCode,
			  globals,
			  locals);
    CxmCheckPtr(obj);
    Py_DECREF(obj);

    retval = (CxtRingObject *) PyDict_GetItemString(locals, "ring");
    CxmCheckPtr(retval);
    Py_INCREF(retval);

    Py_DECREF(locals);

    return retval;
}

PyObject *
CxRingTree(CxtRingObject *self)
{
    return Py_BuildValue("O", self->tree);
}

PyObject *
CxRingNode(CxtRingObject *self)
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

PyObject *
CxRingEdge(CxtRingObject *self)
{
    PyObject *retval;
    CxtTrEdge edge;

    edge = CxTrRingEdgeGet(self->tree->tr, self->ring);
    retval = (PyObject *) CxTrEdgeAuxGet(self->tree->tr, edge);
    Py_INCREF(retval);

    return retval;
}

PyObject *
CxRingOther(CxtRingObject *self)
{
    PyObject *retval;
    CxtTrRing other;

    other = CxTrRingOtherGet(self->tree->tr, self->ring);
    retval = (PyObject *) CxTrRingAuxGet(self->tree->tr, other);
    Py_INCREF(retval);

    return retval;
}

PyObject *
CxRingNext(CxtRingObject *self)
{
    PyObject *retval;
    CxtTrRing nextRing;

    nextRing = CxTrRingNextGet(self->tree->tr, self->ring);
    retval = (PyObject *) CxTrRingAuxGet(self->tree->tr, nextRing);
    Py_INCREF(retval);

    return retval;
}

PyObject *
CxRingPrev(CxtRingObject *self)
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
	(PyCFunction) CxRingTree,
	METH_NOARGS,
	"tree"
    },
    {
	"node",
	(PyCFunction) CxRingNode,
	METH_NOARGS,
	"node"
    },
    {
	"edge",
	(PyCFunction) CxRingEdge,
	METH_NOARGS,
	"edge"
    },
    {
	"other",
	(PyCFunction) CxRingOther,
	METH_NOARGS,
	"other"
    },
    {
	"next",
	(PyCFunction) CxRingNext,
	METH_NOARGS,
	"next"
    },
    {
	"prev",
	(PyCFunction) CxRingPrev,
	METH_NOARGS,
	"prev"
    },
    {NULL, NULL}
};

static PyTypeObject CxtRing =
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
