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

#include <math.h>

static PyTypeObject CxtTree;
static PyTypeObject CxtNode;
static PyTypeObject CxtEdge;

/******************************************************************************/
/* Begin tree. */

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
    return retval;
}

static int
CxpTreeClear(CxtTreeObject *self)
{
    CxtTrNode base;

    base = CxTrBaseGet(self->tr);
    if (base != CxmTrNodeNone)
    {
	CxtNodeObject *node;

	node = (CxtNodeObject *) CxTrNodeAuxGet(self->tr, base);
	CxTrBaseSet(self->tr, CxmTrNodeNone);
	Py_DECREF(node);
    }

    return 0;
}

static void
CxpTreeDelete(CxtTreeObject *self)
{
    CxpTreeClear(self);
    CxTrDelete(self->tr);
    self->ob_type->tp_free((PyObject*) self);
}

static PyObject *CxpTreeNewCode;

CxtTreeObject *
CxTreeNew(void)
{
    CxtTreeObject *retval;
    PyObject *globals, *locals;

    globals = PyEval_GetGlobals();
    locals = Py_BuildValue("{}");

    retval = (CxtTreeObject *)
	PyEval_EvalCode((PyCodeObject *) CxpTreeNewCode, globals, locals);

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
CxTreeBaseSetCargs(CxtTreeObject *self, CxtNodeObject *aNode)
{
    CxtNodeObject *oldNode;
    CxtTrNode oldTrNode;

    /* Decref if clobbering an already-set base. */
    oldTrNode = CxTrBaseGet(self->tr);
    if (oldTrNode != CxmTrNodeNone)
    {
	oldNode = (CxtNodeObject *) CxTrNodeAuxGet(self->tr, oldTrNode);
	CxTrBaseSet(self->tr, CxmTrNodeNone);
	Py_DECREF(oldNode);
    }

    if (aNode != NULL)
    {
	/* Circular reference. */
	Py_INCREF(aNode);
	CxTrBaseSet(self->tr, aNode->node);
    }
}

PyObject *
CxTreeBaseSet(CxtTreeObject *self, PyObject *args)
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
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    CxTreeBaseSetCargs(self, node);

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
	(PyCFunction) CxTreeBaseSet,
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
    "_Tree.Tree",	/* char *tp_name */
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

void
CxTreeInit(void)
{
    PyObject *m;

    if (PyType_Ready(&CxtTree) < 0)
    {
	return;
    }
    m = Py_InitModule3("_Tree", CxpTreeFuncs, "Tree extensions");
    Py_INCREF(&CxtTree);
    PyModule_AddObject(m, "Tree", (PyObject *) &CxtTree);

    /* Pre-compile Python code that is used for creating a tree. */
    CxpTreeNewCode = Py_CompileString("Tree.Tree()",
				      "<string>",
				      Py_eval_input);
}

/* End tree. */
/******************************************************************************/
/* Begin node. */

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
CxpNodeTraverse(CxtNodeObject *self, visitproc visit, void *arg)
{
    int retval;
    CxtTrEdge firstEdge, curEdge;
    uint32_t end;
    CxtEdgeObject *edgeObj;

    if (visit((PyObject *) self->tree, arg) < 0)
    {
	retval = -1;
	goto RETURN;
    }

    CxTrNodeEdgeGet(self->tree->tr, self->node, &firstEdge, &end);
    if (firstEdge != CxmTrEdgeNone)
    {
	curEdge = firstEdge;
	do
	{
	    edgeObj = (CxtEdgeObject *) CxTrEdgeAuxGet(self->tree->tr, curEdge);
	    if (visit((PyObject *) edgeObj, arg) < 0)
	    {
		retval = -1;
		goto RETURN;
	    }

	    CxTrEdgeNextGet(self->tree->tr, curEdge, end, &curEdge, &end);
	} while (curEdge != firstEdge);
    }

    retval = 0;
    RETURN:
    return retval;
}

static int
CxpNodeClear(CxtNodeObject *self)
{
    CxtTrEdge edge;
    CxtEdgeObject *edgeObj;

    while (true)
    {
	CxTrNodeEdgeGet(self->tree->tr, self->node, &edge, NULL);
	if (edge == CxmTrEdgeNone)
	{
	    break;
	}
	edgeObj = CxTrEdgeAuxGet(self->tree->tr, edge);
	CxEdgeDetach(edgeObj);
	Py_DECREF(edgeObj);
    }

    if (CxTrBaseGet(self->tree->tr) == self->node)
    {
	CxTrBaseSet(self->tree->tr, CxmTrNodeNone);
	Py_DECREF(self);
    }
    CxTrNodeDelete(self->tree->tr, self->node);
    Py_DECREF(self->tree);

    return 0;
}

static void
CxpNodeDelete(CxtNodeObject *self)
{
    CxpNodeClear(self);
    self->ob_type->tp_free((PyObject*) self);
}

static PyObject *CxpNodeNewCode;

CxtNodeObject *
CxNodeNew(CxtTreeObject *aTree)
{
    CxtNodeObject *retval;
    PyObject *globals, *locals;

    globals = PyEval_GetGlobals();
    locals = Py_BuildValue("{sO}", "tree", (PyObject *) aTree);

    retval = (CxtNodeObject *)
	PyEval_EvalCode((PyCodeObject *) CxpNodeNewCode,
			globals,
			locals);

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
CxNodeTaxonNumSetCargs(CxtNodeObject *self, uint32_t aTaxonNum)
{
    CxTrNodeTaxonNumSet(self->tree->tr, self->node, aTaxonNum);
}

PyObject *
CxNodeTaxonNumSet(CxtNodeObject *self, PyObject *args)
{
    PyObject *retval;
    uint32_t taxonNum;

    taxonNum = CxmTrNodeTaxonNone;
    if (PyArg_ParseTuple(args, "|i", &taxonNum) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    CxNodeTaxonNumSetCargs(self, taxonNum);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

PyObject *
CxNodeEdge(CxtNodeObject *self)
{
    PyObject *retval;
    CxtEdgeObject *edgeObj;
    CxtTrEdge edge;
    uint32_t end;

    CxTrNodeEdgeGet(self->tree->tr, self->node, &edge, &end);
    if (edge != CxmTrEdgeNone)
    {
	edgeObj = (CxtEdgeObject *) CxTrEdgeAuxGet(self->tree->tr, edge);
	retval = Py_BuildValue("(Oi)", edgeObj, end);
    }
    else
    {
	retval = Py_BuildValue("(ss)", NULL, NULL);
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
	(PyCFunction) CxNodeTaxonNumSet,
	METH_VARARGS,
	"taxonNumSet"
    },
    {
	"edge",
	(PyCFunction) CxNodeEdge,
	METH_NOARGS,
	"edge"
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
    "_Node.Node",	/* char *tp_name */
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

    if (PyType_Ready(&CxtNode) < 0)
    {
	return;
    }
    m = Py_InitModule3("_Node", CxpNodeFuncs, "Node extensions");
    Py_INCREF(&CxtNode);
    PyModule_AddObject(m, "Node", (PyObject *) &CxtNode);

    /* Pre-compile Python code that is used for creating a node. */
    CxpNodeNewCode = Py_CompileString("Node.Node(tree)",
				      "<string>",
				      Py_eval_input);
}

/* End node. */
/******************************************************************************/
/* Begin edge. */

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
CxpEdgeTraverse(CxtEdgeObject *self, visitproc visit, void *arg)
{
    int retval;
    CxtTrNode trNode;

    if (visit((PyObject *) self->tree, arg) < 0)
    {
	retval = -1;
	goto RETURN;
    }

    trNode = CxTrEdgeNodeGet(self->tree->tr, self->edge, 0);
    if (trNode != CxmTrNodeNone)
    {
	CxtNodeObject *node;

	node = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, trNode);
	if (visit((PyObject *) node, arg) < 0)
	{
	    retval = -1;
	    goto RETURN;
	}

	trNode = CxTrEdgeNodeGet(self->tree->tr, self->edge, 1);
	node = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, trNode);
	if (visit((PyObject *) node, arg) < 0)
	{
	    retval = -1;
	    goto RETURN;
	}
    }

    retval = 0;
    RETURN:
    return retval;
}

static int
CxpEgeClear(CxtEdgeObject *self)
{
    CxtTrNode nodeA;

    nodeA = CxTrEdgeNodeGet(self->tree->tr, self->edge, 0);
    if (nodeA != CxmTrNodeNone)
    {
	CxtTrNode nodeB;
	CxtNodeObject *nodeAObj, *nodeBObj;

	nodeAObj = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, nodeA);

	nodeB = CxTrEdgeNodeGet(self->tree->tr, self->edge, 1);
	nodeBObj = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, nodeB);

	CxTrEdgeDetach(self->tree->tr, self->edge);
	Py_DECREF(nodeAObj);
	Py_DECREF(nodeBObj);
    }

    CxTrEdgeDelete(self->tree->tr, self->edge);
    Py_DECREF(self->tree);

    return 0;
}

static void
CxpEdgeDelete(CxtEdgeObject *self)
{
    CxpEgeClear(self);
    self->ob_type->tp_free((PyObject*) self);
}

static PyObject *CxpEdgeNewCode;

CxtEdgeObject *
CxEdgeNew(CxtTreeObject *aTree)
{
    CxtEdgeObject *retval;
    PyObject *globals, *locals;

    globals = PyEval_GetGlobals();
    locals = Py_BuildValue("{sO}", "tree", (PyObject *) aTree);

    retval = (CxtEdgeObject *)
	PyEval_EvalCode((PyCodeObject *) CxpEdgeNewCode,
			globals,
			locals);

    Py_DECREF(locals);

    return retval;
}

PyObject *
CxEdgeTree(CxtEdgeObject *self)
{
    return Py_BuildValue("O", self->tree);
}

PyObject *
CxEdgeNodeCargs(CxtEdgeObject *self, uint32_t aInd)
{
    PyObject *retval;
    CxtTrNode node;

    node = CxTrEdgeNodeGet(self->tree->tr, self->edge, aInd);
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
CxEdgeNode(CxtEdgeObject *self, PyObject *args)
{
    PyObject *retval;
    uint32_t ind;

    if (PyArg_ParseTuple(args, "i", &ind) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (ind != 0 && ind != 1)
    {
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    retval = CxEdgeNodeCargs(self, ind);

    RETURN:
    return retval;
}

void
CxEdgeNextCargs(CxtEdgeObject *self, uint32_t aInd, CxtEdgeObject **rEdge,
		uint32_t *rNextEnd)
{
    CxtTrEdge nextEdge;

    CxTrEdgeNextGet(self->tree->tr, self->edge, aInd, &nextEdge, rNextEnd);
    *rEdge = (CxtEdgeObject *) CxTrEdgeAuxGet(self->tree->tr, nextEdge);
}

PyObject *
CxEdgeNext(CxtEdgeObject *self, PyObject *args)
{
    PyObject *retval;
    CxtEdgeObject *edgeObj;
    uint32_t ind, nextEnd;

    if (PyArg_ParseTuple(args, "i", &ind) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (ind != 0 && ind != 1)
    {
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    CxEdgeNextCargs(self, ind, &edgeObj, &nextEnd);

    retval = Py_BuildValue("(Oi)", edgeObj, nextEnd);
    RETURN:
    return retval;
}

void
CxEdgePrevCargs(CxtEdgeObject *self, uint32_t aInd, CxtEdgeObject **rEdge,
		uint32_t *rPrevEnd)
{
    CxtTrEdge prevEdge;

    CxTrEdgePrevGet(self->tree->tr, self->edge, aInd, &prevEdge, rPrevEnd);
    *rEdge = (CxtEdgeObject *) CxTrEdgeAuxGet(self->tree->tr, prevEdge);
}

PyObject *
CxEdgePrev(CxtEdgeObject *self, PyObject *args)
{
    PyObject *retval;
    CxtEdgeObject *edgeObj;
    uint32_t ind, prevEnd;

    if (PyArg_ParseTuple(args, "i", &ind) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (ind != 0 && ind != 1)
    {
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    CxEdgePrevCargs(self, ind, &edgeObj, &prevEnd);

    retval = Py_BuildValue("(Oi)", edgeObj, prevEnd);
    RETURN:
    return retval;
}

PyObject *
CxEdgeLengthGet(CxtEdgeObject *self)
{
    return Py_BuildValue("d", CxTrEdgeLengthGet(self->tree->tr, self->edge));
}

void
CxEdgeLengthSetCargs(CxtEdgeObject *self, double aLength)
{
    CxTrEdgeLengthSet(self->tree->tr, self->edge, aLength);
}

PyObject *
CxEdgeLengthSet(CxtEdgeObject *self, PyObject *args)
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
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    CxEdgeLengthSetCargs(self, length);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

void
CxEdgeAttachCargs(CxtEdgeObject *self, CxtNodeObject *aNodeA,
		  CxtNodeObject *aNodeB)
{
    Py_INCREF(aNodeA);
    Py_INCREF(aNodeB);
    /* Cyclic references (nodes refer to edge). */
    Py_INCREF(self);
    Py_INCREF(self);
    CxTrEdgeAttach(self->tree->tr, self->edge, aNodeA->node, aNodeB->node);
}

PyObject *
CxEdgeAttach(CxtEdgeObject *self, PyObject *args)
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
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    CxEdgeAttachCargs(self, nodeA, nodeB);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

PyObject *
CxEdgeDetach(CxtEdgeObject *self)
{
    CxtTrNode nodeA;

    nodeA = CxTrEdgeNodeGet(self->tree->tr, self->edge, 0);
    if (nodeA != CxmTrNodeNone)
    {
	CxtTrNode nodeB;
	CxtNodeObject *nodeAObj, *nodeBObj;

	nodeAObj = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, nodeA);

	nodeB = CxTrEdgeNodeGet(self->tree->tr, self->edge, 1);
	nodeBObj = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, nodeB);

	CxTrEdgeDetach(self->tree->tr, self->edge);
	Py_DECREF(nodeAObj);
	Py_DECREF(nodeBObj);
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
	"node",
	(PyCFunction) CxEdgeNode,
	METH_VARARGS,
	"node"
    },
    {
	"next",
	(PyCFunction) CxEdgeNext,
	METH_VARARGS,
	"next"
    },
    {
	"prev",
	(PyCFunction) CxEdgePrev,
	METH_VARARGS,
	"prev"
    },
    {
	"lengthGet",
	(PyCFunction) CxEdgeLengthGet,
	METH_NOARGS,
	"lengthGet"
    },
    {
	"lengthSet",
	(PyCFunction) CxEdgeLengthSet,
	METH_VARARGS,
	"lengthSet"
    },
    {
	"attach",
	(PyCFunction) CxEdgeAttach,
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
    "_Edge.Edge",	/* char *tp_name */
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
    (inquiry) CxpEgeClear,	/* inquiry tp_clear */
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

void
CxEdgeInit(void)
{
    PyObject *m;

    if (PyType_Ready(&CxtEdge) < 0)
    {
	return;
    }
    m = Py_InitModule3("_Edge", CxpEdgeFuncs, "Edge extensions");
    Py_INCREF(&CxtEdge);
    PyModule_AddObject(m, "Edge", (PyObject *) &CxtEdge);

    /* Pre-compile Python code that is used for creating a wrapped edge. */
    CxpEdgeNewCode = Py_CompileString("Edge.Edge(tree)",
				      "<string>",
				      Py_eval_input);
}

/* End edge. */
/******************************************************************************/
