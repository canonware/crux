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

static PyTypeObject Tree_Type;
static PyTypeObject Node_Type;
static PyTypeObject Edge_Type;

/******************************************************************************/
/* Begin tree. */

static PyObject *
tree_p_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
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
tree_p_traverse(CxtTreeObject *self, visitproc visit, void *arg)
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
tree_p_clear(CxtTreeObject *self)
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
tree_p_delete(CxtTreeObject *self)
{
    tree_p_clear(self);
    CxTrDelete(self->tr);
    self->ob_type->tp_free((PyObject*) self);
}

static PyObject *tree_p_new_code;

CxtTreeObject *
CxTreeNew(void)
{
    CxtTreeObject *retval;
    PyObject *globals, *locals;

    globals = PyEval_GetGlobals();
    locals = Py_BuildValue("{}");

    retval = (CxtTreeObject *)
	PyEval_EvalCode((PyCodeObject *) tree_p_new_code, globals, locals);

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
CxTreeBaseSetCargs(CxtTreeObject *self, CxtNodeObject *a_node)
{
    CxtNodeObject *old_node;
    CxtTrNode old_tr_node;

    /* Decref if clobbering an already-set base. */
    old_tr_node = CxTrBaseGet(self->tr);
    if (old_tr_node != CxmTrNodeNone)
    {
	old_node = (CxtNodeObject *) CxTrNodeAuxGet(self->tr, old_tr_node);
	CxTrBaseSet(self->tr, CxmTrNodeNone);
	Py_DECREF(old_node);
    }

    if (a_node != NULL)
    {
	/* Circular reference. */
	Py_INCREF(a_node);
	CxTrBaseSet(self->tr, a_node->node);
    }
}

PyObject *
CxTreeBaseSet(CxtTreeObject *self, PyObject *args)
{
    PyObject *retval;
    CxtNodeObject *node;

    node = NULL;
    if (PyArg_ParseTuple(args, "|O!", &Node_Type, &node) == 0)
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

static PyMethodDef tree_p_methods[] =
{
    {
	"ntaxa_get",
	(PyCFunction) CxTreeNtaxaGet,
	METH_NOARGS,
	"ntaxa_get"
    },
    {
	"nedges_get",
	(PyCFunction) CxTreeNedgesCget,
	METH_NOARGS,
	"nedges_get"
    },
    {
	"base_get",
	(PyCFunction) CxTreeBaseGet,
	METH_NOARGS,
	"base_get"
    },
    {
	"base_set",
	(PyCFunction) CxTreeBaseSet,
	METH_VARARGS,
	"base_set"
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
	"tbr_nneighbors_get",
	(PyCFunction) CxTreeTbrNneighborsGet,
	METH_NOARGS,
	"tbr_nneighbors_get"
    },
    {
	"tbr_neighbor_get",
	(PyCFunction) CxTreeTbrNeighborGet,
	METH_VARARGS,
	"tbr_neighbor_get"
    },
    {
	"mp_prepare",
	(PyCFunction) CxTreeMpPrepare,
	METH_VARARGS,
	"mp_prepare"
    },
    {
	"mp_finish",
	(PyCFunction) CxTreeMpFinish,
	METH_NOARGS,
	"mp_finish"
    },
    {
	"mp",
	(PyCFunction) CxTreeMp,
	METH_NOARGS,
	"mp"
    },
    {
	"tbr_best_neighbors_mp",
	(PyCFunction) CxTreeTbrBestNeighborsMp,
	METH_VARARGS,
	"tbr_best_neighbors_mp"
    },
    {
	"tbr_better_neighbors_mp",
	(PyCFunction) CxTreeTbrBetterNeighborsMp,
	METH_VARARGS,
	"tbr_better_neighbors_mp"
    },
    {
	"tbr_all_neighbors_mp",
	(PyCFunction) CxTreeTbrAllNeighborsMp,
	METH_NOARGS,
	"tbr_all_neighbors_mp"
    },
    {
	"nheld_get",
	(PyCFunction) CxTreeNheldGet,
	METH_NOARGS,
	"nheld_get"
    },
    {
	"held_get",
	(PyCFunction) CxTreeheldGet,
	METH_VARARGS,
	"held_get"
    },
    {NULL, NULL}
};

static PyTypeObject Tree_Type =
{
    PyObject_HEAD_INIT(NULL)
    0,			/* int ob_size */
    "_tree.Tree",	/* char *tp_name */
    sizeof(CxtTreeObject),	/* int tp_basicsize */
    0,			/* int tp_itemsize */
    (destructor) tree_p_delete,	/* destructor tp_dealloc */
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
    (traverseproc) tree_p_traverse,	/* traverseproc tp_traverse */
    (inquiry) tree_p_clear,	/* inquiry tp_clear */
    0,			/* richcmpfunc tp_richcompare */
    0,			/* long tp_weaklistoffset */
    0,			/* getiterfunc tp_iter */
    0,			/* iternextfunc tp_iternext */
    tree_p_methods,	/* struct PyMethodDef *tp_methods */
    0,			/* struct PyMemberDef *tp_members */
    0,			/* struct PyGetSetDef *tp_getset */
    0,			/* struct _typeobject *tp_base */
    0,			/* PyObject *tp_dict */
    0,			/* descrgetfunc tp_descr_get */
    0,			/* descrsetfunc tp_descr_set */
    0,			/* long tp_dictoffset */
    0,			/* initproc tp_init */
    0,			/* allocfunc tp_alloc */
    tree_p_new,		/* newfunc tp_new */
    _PyObject_Del,	/* freefunc tp_free */
    0			/* inquiry tp_is_gc */
};

static PyMethodDef tree_p_funcs[] =
{
    {NULL}
};

void
CxTreeInit(void)
{
    PyObject *m;

    if (PyType_Ready(&Tree_Type) < 0)
    {
	return;
    }
    m = Py_InitModule3("_tree", tree_p_funcs, "tree extensions");
    Py_INCREF(&Tree_Type);
    PyModule_AddObject(m, "Tree", (PyObject *) &Tree_Type);

    /* Pre-compile Python code that is used for creating a tree. */
    tree_p_new_code = Py_CompileString("tree.tree()",
				       "<string>",
				       Py_eval_input);
}

/* End tree. */
/******************************************************************************/
/* Begin node. */

static PyObject *
node_p_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtNodeObject *self;
    CxtTreeObject *tree;

    if (PyArg_ParseTuple(args, "O!", &Tree_Type, &tree) == 0)
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
node_p_traverse(CxtNodeObject *self, visitproc visit, void *arg)
{
    int retval;
    CxtTrEdge first_edge, cur_edge;
    uint32_t end;
    CxtEdgeObject *edge_obj;

    if (visit((PyObject *) self->tree, arg) < 0)
    {
	retval = -1;
	goto RETURN;
    }

    CxTrNodeEdgeGet(self->tree->tr, self->node, &first_edge, &end);
    if (first_edge != CxmTrEdgeNone)
    {
	cur_edge = first_edge;
	do
	{
	    edge_obj = (CxtEdgeObject *) CxTrEdgeAuxGet(self->tree->tr, cur_edge);
	    if (visit((PyObject *) edge_obj, arg) < 0)
	    {
		retval = -1;
		goto RETURN;
	    }

	    CxTrEdgeNextGet(self->tree->tr, cur_edge, end, &cur_edge, &end);
	} while (cur_edge != first_edge);
    }

    retval = 0;
    RETURN:
    return retval;
}

static int
node_p_clear(CxtNodeObject *self)
{
    CxtTrEdge edge;
    CxtEdgeObject *edge_obj;

    while (true)
    {
	CxTrNodeEdgeGet(self->tree->tr, self->node, &edge, NULL);
	if (edge == CxmTrEdgeNone)
	{
	    break;
	}
	edge_obj = CxTrEdgeAuxGet(self->tree->tr, edge);
	CxEdgeDetach(edge_obj);
	Py_DECREF(edge_obj);
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
node_p_delete(CxtNodeObject *self)
{
    node_p_clear(self);
    self->ob_type->tp_free((PyObject*) self);
}

static PyObject *tree_p_CxNodeNew_code;

CxtNodeObject *
CxNodeNew(CxtTreeObject *a_tree)
{
    CxtNodeObject *retval;
    PyObject *globals, *locals;

    globals = PyEval_GetGlobals();
    locals = Py_BuildValue("{sO}", "tree", (PyObject *) a_tree);

    retval = (CxtNodeObject *)
	PyEval_EvalCode((PyCodeObject *) tree_p_CxNodeNew_code,
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
    uint32_t taxon_num;

    taxon_num = CxTrNodeTaxonNumGet(self->tree->tr, self->node);

    if (taxon_num == CxmTrNodeTaxonNone)
    {
	Py_INCREF(Py_None);
	retval = Py_None;
    }
    else
    {
	retval = Py_BuildValue("i", taxon_num);
    }

    return retval;
}

void
CxNodeTaxonNumSetCargs(CxtNodeObject *self, uint32_t a_taxon_num)
{
    CxTrNodeTaxonNumSet(self->tree->tr, self->node, a_taxon_num);
}

PyObject *
CxNodeTaxonNumSet(CxtNodeObject *self, PyObject *args)
{
    PyObject *retval;
    uint32_t taxon_num;

    taxon_num = CxmTrNodeTaxonNone;
    if (PyArg_ParseTuple(args, "|i", &taxon_num) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    CxNodeTaxonNumSetCargs(self, taxon_num);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

PyObject *
CxNodeEdge(CxtNodeObject *self)
{
    PyObject *retval;
    CxtEdgeObject *edge_obj;
    CxtTrEdge edge;
    uint32_t end;

    CxTrNodeEdgeGet(self->tree->tr, self->node, &edge, &end);
    if (edge != CxmTrEdgeNone)
    {
	edge_obj = (CxtEdgeObject *) CxTrEdgeAuxGet(self->tree->tr, edge);
	retval = Py_BuildValue("(Oi)", edge_obj, end);
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

static PyMethodDef node_p_methods[] =
{
    {
	"tree",
	(PyCFunction) CxNodeTree,
	METH_NOARGS,
	"tree"
    },
    {
	"taxon_num_get",
	(PyCFunction) CxNodeTaxonNumGet,
	METH_NOARGS,
	"taxon_num_get"
    },
    {
	"taxon_num_set",
	(PyCFunction) CxNodeTaxonNumSet,
	METH_VARARGS,
	"taxon_num_set"
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

static PyTypeObject Node_Type =
{
    PyObject_HEAD_INIT(NULL)
    0,			/* int ob_size */
    "_node.Node",	/* char *tp_name */
    sizeof(CxtNodeObject),	/* int tp_basicsize */
    0,			/* int tp_itemsize */
    (destructor) node_p_delete,	/* destructor tp_dealloc */
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
    (traverseproc) node_p_traverse,	/* traverseproc tp_traverse */
    (inquiry) node_p_clear,	/* inquiry tp_clear */
    0,			/* richcmpfunc tp_richcompare */
    0,			/* long tp_weaklistoffset */
    0,			/* getiterfunc tp_iter */
    0,			/* iternextfunc tp_iternext */
    node_p_methods,	/* struct PyMethodDef *tp_methods */
    0,			/* struct PyMemberDef *tp_members */
    0,			/* struct PyGetSetDef *tp_getset */
    0,			/* struct _typeobject *tp_base */
    0,			/* PyObject *tp_dict */
    0,			/* descrgetfunc tp_descr_get */
    0,			/* descrsetfunc tp_descr_set */
    0,			/* long tp_dictoffset */
    0,			/* initproc tp_init */
    0,			/* allocfunc tp_alloc */
    node_p_new,		/* newfunc tp_new */
    _PyObject_Del,	/* freefunc tp_free */
    0			/* inquiry tp_is_gc */
};

static PyMethodDef node_p_funcs[] =
{
    {NULL}
};

void
CxNodeInit(void)
{
    PyObject *m;

    if (PyType_Ready(&Node_Type) < 0)
    {
	return;
    }
    m = Py_InitModule3("_node", node_p_funcs, "node extensions");
    Py_INCREF(&Node_Type);
    PyModule_AddObject(m, "Node", (PyObject *) &Node_Type);

    /* Pre-compile Python code that is used for creating a node. */
    tree_p_CxNodeNew_code = Py_CompileString("node.node(tree)",
					    "<string>",
					    Py_eval_input);
}

/* End node. */
/******************************************************************************/
/* Begin edge. */

static PyObject *
edge_p_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtEdgeObject *self;
    CxtTreeObject *tree;

    if (PyArg_ParseTuple(args, "O!", &Tree_Type, &tree) == 0)
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
edge_p_traverse(CxtEdgeObject *self, visitproc visit, void *arg)
{
    int retval;
    CxtTrNode tr_node;

    if (visit((PyObject *) self->tree, arg) < 0)
    {
	retval = -1;
	goto RETURN;
    }

    tr_node = CxTrEdgeNodeGet(self->tree->tr, self->edge, 0);
    if (tr_node != CxmTrNodeNone)
    {
	CxtNodeObject *node;

	node = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, tr_node);
	if (visit((PyObject *) node, arg) < 0)
	{
	    retval = -1;
	    goto RETURN;
	}

	tr_node = CxTrEdgeNodeGet(self->tree->tr, self->edge, 1);
	node = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, tr_node);
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
edge_p_clear(CxtEdgeObject *self)
{
    CxtTrNode node_a;

    node_a = CxTrEdgeNodeGet(self->tree->tr, self->edge, 0);
    if (node_a != CxmTrNodeNone)
    {
	CxtTrNode node_b;
	CxtNodeObject *node_a_obj, *node_b_obj;

	node_a_obj = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, node_a);

	node_b = CxTrEdgeNodeGet(self->tree->tr, self->edge, 1);
	node_b_obj = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, node_b);

	CxTrEdgeDetach(self->tree->tr, self->edge);
	Py_DECREF(node_a_obj);
	Py_DECREF(node_b_obj);
    }

    CxTrEdgeDelete(self->tree->tr, self->edge);
    Py_DECREF(self->tree);

    return 0;
}

static void
edge_p_delete(CxtEdgeObject *self)
{
    edge_p_clear(self);
    self->ob_type->tp_free((PyObject*) self);
}

static PyObject *tree_p_CxEdgeNew_code;

CxtEdgeObject *
CxEdgeNew(CxtTreeObject *a_tree)
{
    CxtEdgeObject *retval;
    PyObject *globals, *locals;

    globals = PyEval_GetGlobals();
    locals = Py_BuildValue("{sO}", "tree", (PyObject *) a_tree);

    retval = (CxtEdgeObject *)
	PyEval_EvalCode((PyCodeObject *) tree_p_CxEdgeNew_code,
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
CxEdgeNodeCargs(CxtEdgeObject *self, uint32_t a_ind)
{
    PyObject *retval;
    CxtTrNode node;

    node = CxTrEdgeNodeGet(self->tree->tr, self->edge, a_ind);
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
CxEdgeNextCargs(CxtEdgeObject *self, uint32_t a_ind, CxtEdgeObject **r_edge,
		uint32_t *r_next_end)
{
    CxtTrEdge next_edge;

    CxTrEdgeNextGet(self->tree->tr, self->edge, a_ind, &next_edge, r_next_end);
    *r_edge = (CxtEdgeObject *) CxTrEdgeAuxGet(self->tree->tr, next_edge);
}

PyObject *
CxEdgeNext(CxtEdgeObject *self, PyObject *args)
{
    PyObject *retval;
    CxtEdgeObject *edge_obj;
    uint32_t ind, next_end;

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

    CxEdgeNextCargs(self, ind, &edge_obj, &next_end);

    retval = Py_BuildValue("(Oi)", edge_obj, next_end);
    RETURN:
    return retval;
}

void
CxEdgePrevCargs(CxtEdgeObject *self, uint32_t a_ind, CxtEdgeObject **r_edge,
		uint32_t *r_prev_end)
{
    CxtTrEdge prev_edge;

    CxTrEdgePrevGet(self->tree->tr, self->edge, a_ind, &prev_edge, r_prev_end);
    *r_edge = (CxtEdgeObject *) CxTrEdgeAuxGet(self->tree->tr, prev_edge);
}

PyObject *
CxEdgePrev(CxtEdgeObject *self, PyObject *args)
{
    PyObject *retval;
    CxtEdgeObject *edge_obj;
    uint32_t ind, prev_end;

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

    CxEdgePrevCargs(self, ind, &edge_obj, &prev_end);

    retval = Py_BuildValue("(Oi)", edge_obj, prev_end);
    RETURN:
    return retval;
}

PyObject *
CxEdgeLengthGet(CxtEdgeObject *self)
{
    return Py_BuildValue("d", CxTrEdgeLengthGet(self->tree->tr, self->edge));
}

void
CxEdgeLengthSetCargs(CxtEdgeObject *self, double a_length)
{
    CxTrEdgeLengthSet(self->tree->tr, self->edge, a_length);
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
CxEdgeAttachCargs(CxtEdgeObject *self, CxtNodeObject *a_node_a, CxtNodeObject *a_node_b)
{
    Py_INCREF(a_node_a);
    Py_INCREF(a_node_b);
    /* Cyclic references (nodes refer to edge). */
    Py_INCREF(self);
    Py_INCREF(self);
    CxTrEdgeAttach(self->tree->tr, self->edge, a_node_a->node, a_node_b->node);
}

PyObject *
CxEdgeAttach(CxtEdgeObject *self, PyObject *args)
{
    PyObject *retval;
    CxtNodeObject *node_a, *node_b;

    if (PyArg_ParseTuple(args, "O!O!", &Node_Type, &node_a,
			 &Node_Type, &node_b) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (node_a->tree != self->tree || node_b->tree != self->tree)
    {
	Py_INCREF(PyExc_ValueError);
	retval = PyExc_ValueError;
	goto RETURN;
    }

    CxEdgeAttachCargs(self, node_a, node_b);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

PyObject *
CxEdgeDetach(CxtEdgeObject *self)
{
    CxtTrNode node_a;

    node_a = CxTrEdgeNodeGet(self->tree->tr, self->edge, 0);
    if (node_a != CxmTrNodeNone)
    {
	CxtTrNode node_b;
	CxtNodeObject *node_a_obj, *node_b_obj;

	node_a_obj = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, node_a);

	node_b = CxTrEdgeNodeGet(self->tree->tr, self->edge, 1);
	node_b_obj = (CxtNodeObject *) CxTrNodeAuxGet(self->tree->tr, node_b);

	CxTrEdgeDetach(self->tree->tr, self->edge);
	Py_DECREF(node_a_obj);
	Py_DECREF(node_b_obj);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef edge_p_methods[] =
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
	"length_get",
	(PyCFunction) CxEdgeLengthGet,
	METH_NOARGS,
	"length_get"
    },
    {
	"length_set",
	(PyCFunction) CxEdgeLengthSet,
	METH_VARARGS,
	"length_set"
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

static PyTypeObject Edge_Type =
{
    PyObject_HEAD_INIT(NULL)
    0,			/* int ob_size */
    "_edge.Edge",	/* char *tp_name */
    sizeof(CxtEdgeObject),	/* int tp_basicsize */
    0,			/* int tp_itemsize */
    (destructor) edge_p_delete,	/* destructor tp_dealloc */
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
    (traverseproc) edge_p_traverse,	/* traverseproc tp_traverse */
    (inquiry) edge_p_clear,	/* inquiry tp_clear */
    0,			/* richcmpfunc tp_richcompare */
    0,			/* long tp_weaklistoffset */
    0,			/* getiterfunc tp_iter */
    0,			/* iternextfunc tp_iternext */
    edge_p_methods,	/* struct PyMethodDef *tp_methods */
    0,			/* struct PyMemberDef *tp_members */
    0,			/* struct PyGetSetDef *tp_getset */
    0,			/* struct _typeobject *tp_base */
    0,			/* PyObject *tp_dict */
    0,			/* descrgetfunc tp_descr_get */
    0,			/* descrsetfunc tp_descr_set */
    0,			/* long tp_dictoffset */
    0,			/* initproc tp_init */
    0,			/* allocfunc tp_alloc */
    edge_p_new,		/* newfunc tp_new */
    _PyObject_Del,	/* freefunc tp_free */
    0			/* inquiry tp_is_gc */
};

static PyMethodDef edge_p_funcs[] =
{
    {NULL}
};

void
CxEdgeInit(void)
{
    PyObject *m;

    if (PyType_Ready(&Edge_Type) < 0)
    {
	return;
    }
    m = Py_InitModule3("_edge", edge_p_funcs, "edge extensions");
    Py_INCREF(&Edge_Type);
    PyModule_AddObject(m, "Edge", (PyObject *) &Edge_Type);

    /* Pre-compile Python code that is used for creating a wrapped edge. */
    tree_p_CxEdgeNew_code = Py_CompileString("edge.edge(tree)",
					    "<string>",
					    Py_eval_input);
}

/* End edge. */
/******************************************************************************/
