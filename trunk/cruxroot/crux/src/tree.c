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
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;
    TreeObject *self;

    self = (TreeObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    xep_begin();
    xep_try
    {
	self->tr = tr_new(NULL, NULL, NULL, NULL);//XXX Remove args to tr_new().
	tr_aux_set(self->tr, self);
	retval = (PyObject *) self;
    }
    xep_catch(CW_CRUXX_OOM)
    {
	xep_handled();
	retval = PyErr_NoMemory();
    }
    xep_end();

    RETURN:
    return retval;
}

static int
tree_p_traverse(TreeObject *self, visitproc visit, void *arg)
{
    int retval;
    cw_tr_node_t base;

    base = tr_base_get(self->tr);
    if (base != CW_TR_NODE_NONE)
    {
	NodeObject *node;

	node = (NodeObject *) tr_node_aux_get(self->tr, base);
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
tree_p_clear(TreeObject *self)
{
    cw_tr_node_t base;

    base = tr_base_get(self->tr);
    if (base != CW_TR_NODE_NONE)
    {
	NodeObject *node;

	node = (NodeObject *) tr_node_aux_get(self->tr, base);
	tr_base_set(self->tr, CW_TR_NODE_NONE);
	Py_DECREF(node);
    }

    return 0;
}

static void
tree_p_delete(TreeObject *self)
{
    tree_p_clear(self);
    tr_delete(self->tr);
    self->ob_type->tp_free((PyObject*) self);
}

PyObject *
tree_ntaxa_get(TreeObject *self)
{
    PyObject *retval
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;

    xep_begin();
    xep_try
    {
	retval = Py_BuildValue("i", tr_ntaxa_get(self->tr));
    }
    xep_catch(CW_CRUXX_OOM)
    {
	xep_handled();
	retval = PyErr_NoMemory();
    }
    xep_end();

    return retval;
}

PyObject *
tree_nedges_get(TreeObject *self)
{
    PyObject *retval
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;

    xep_begin();
    xep_try
    {
	retval = Py_BuildValue("i", tr_nedges_get(self->tr));
    }
    xep_catch(CW_CRUXX_OOM)
    {
	xep_handled();
	retval = PyErr_NoMemory();
    }
    xep_end();

    return retval;
}

PyObject *
tree_base_get(TreeObject *self)
{
    PyObject *retval;
    cw_tr_node_t base;

    base = tr_base_get(self->tr);
    if (base == CW_TR_NODE_NONE)
    {
	Py_INCREF(Py_None);
	retval = Py_None;
    }
    else
    {
	retval = (PyObject *) tr_node_aux_get(self->tr, base);
	Py_INCREF(retval);
    }

    return retval;
}

void
tree_base_set_cargs(TreeObject *self, NodeObject *a_node)
{
    NodeObject *old_node;
    cw_tr_node_t old_tr_node;

    /* Decref if clobbering an already-set base. */
    old_tr_node = tr_base_get(self->tr);
    if (old_tr_node != CW_TR_NODE_NONE)
    {
	old_node = (NodeObject *) tr_node_aux_get(self->tr, old_tr_node);
	tr_base_set(self->tr, CW_TR_NODE_NONE);
	Py_DECREF(old_node);
    }

    if (a_node != NULL)
    {
	/* Circular reference. */
	Py_INCREF(a_node);
	tr_base_set(self->tr, a_node->node);
    }
}

PyObject *
tree_base_set(TreeObject *self, PyObject *args)
{
    PyObject *retval;
    NodeObject *node;

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

    tree_base_set_cargs(self, node);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

static PyMethodDef tree_p_methods[] =
{
    {
	"ntaxa_get",
	(PyCFunction) tree_ntaxa_get,
	METH_NOARGS,
	"ntaxa_get"
    },
    {
	"nedges_get",
	(PyCFunction) tree_nedges_get,
	METH_NOARGS,
	"nedges_get"
    },
    {
	"base_get",
	(PyCFunction) tree_base_get,
	METH_NOARGS,
	"base_get"
    },
    {
	"base_set",
	(PyCFunction) tree_base_set,
	METH_VARARGS,
	"base_set"
    },
    {
	"canonize",
	(PyCFunction) tree_canonize,
	METH_NOARGS,
	"canonize"
    },
    {
	"_nj",
	(PyCFunction) tree__nj,
	METH_VARARGS,
	"_nj"
    },
    {
	"tbr",
	(PyCFunction) tree_tbr,
	METH_VARARGS,
	"tbr"
    },
    {
	"tbr_nneighbors_get",
	(PyCFunction) tree_tbr_nneighbors_get,
	METH_NOARGS,
	"tbr_nneighbors_get"
    },
    {
	"tbr_neighbor_get",
	(PyCFunction) tree_tbr_neighbor_get,
	METH_VARARGS,
	"tbr_neighbor_get"
    },
    {
	"mp_prepare",
	(PyCFunction) tree_mp_prepare,
	METH_VARARGS,
	"mp_prepare"
    },
    {
	"mp_finish",
	(PyCFunction) tree_mp_finish,
	METH_NOARGS,
	"mp_finish"
    },
    {
	"mp",
	(PyCFunction) tree_mp,
	METH_NOARGS,
	"mp"
    },
    {
	"tbr_best_neighbors_mp",
	(PyCFunction) tree_tbr_best_neighbors_mp,
	METH_VARARGS,
	"tbr_best_neighbors_mp"
    },
    {
	"tbr_better_neighbors_mp",
	(PyCFunction) tree_tbr_better_neighbors_mp,
	METH_VARARGS,
	"tbr_better_neighbors_mp"
    },
    {
	"tbr_all_neighbors_mp",
	(PyCFunction) tree_tbr_all_neighbors_mp,
	METH_NOARGS,
	"tbr_all_neighbors_mp"
    },
    {
	"nheld_get",
	(PyCFunction) tree_nheld_get,
	METH_NOARGS,
	"nheld_get"
    },
    {
	"held_get",
	(PyCFunction) tree_held_get,
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
    sizeof(TreeObject),	/* int tp_basicsize */
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

static PyObject *tree_p_new_code;

TreeObject *
tree_new(void)
{
    TreeObject *retval;
    PyObject *globals, *locals;

    globals = PyEval_GetGlobals();
    locals = Py_BuildValue("{}");

    retval = (TreeObject *)
	PyEval_EvalCode((PyCodeObject *) tree_p_new_code, globals, locals);

    Py_DECREF(locals);

    return retval;
}

void
crux_tree_init(void)
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
    tree_p_new_code = Py_CompileString("crux.tree()",
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
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;
    NodeObject *self;
    TreeObject *tree;

    if (PyArg_ParseTuple(args, "O!", &Tree_Type, &tree) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    self = (NodeObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    xep_begin();
    xep_try
    {
	Py_INCREF(tree);
	self->tree = tree;
	self->node = tr_node_new(tree->tr);
	tr_node_aux_set(tree->tr, self->node, self);

	retval = (PyObject *) self;
    }
    xep_catch(CW_CRUXX_OOM)
    {
	xep_handled();
	retval = PyErr_NoMemory();
    }
    xep_end();

    RETURN:
    return retval;
}

static int
node_p_traverse(NodeObject *self, visitproc visit, void *arg)
{
    int retval;
    cw_tr_edge_t first_edge, cur_edge;
    uint32_t end;
    EdgeObject *edge_obj;

    if (visit((PyObject *) self->tree, arg) < 0)
    {
	retval = -1;
	goto RETURN;
    }

    tr_node_edge_get(self->tree->tr, self->node, &first_edge, &end);
    if (first_edge != CW_TR_EDGE_NONE)
    {
	cur_edge = first_edge;
	do
	{
	    edge_obj = (EdgeObject *) tr_edge_aux_get(self->tree->tr, cur_edge);
	    if (visit((PyObject *) edge_obj, arg) < 0)
	    {
		retval = -1;
		goto RETURN;
	    }

	    tr_edge_next_get(self->tree->tr, cur_edge, end, &cur_edge, &end);
	} while (cur_edge != first_edge);
    }

    retval = 0;
    RETURN:
    return retval;
}

static int
node_p_clear(NodeObject *self)
{
    cw_tr_edge_t edge;
    EdgeObject *edge_obj;

    while (true)
    {
	tr_node_edge_get(self->tree->tr, self->node, &edge, NULL);
	if (edge == CW_TR_EDGE_NONE)
	{
	    break;
	}
	edge_obj = tr_edge_aux_get(self->tree->tr, edge);
	edge_detach(edge_obj);
	Py_DECREF(edge_obj);
    }

    if (tr_base_get(self->tree->tr) == self->node)
    {
	tr_base_set(self->tree->tr, CW_TR_NODE_NONE);
	Py_DECREF(self);
    }
    tr_node_delete(self->tree->tr, self->node);
    Py_DECREF(self->tree);

    return 0;
}

static void
node_p_delete(NodeObject *self)
{
    node_p_clear(self);
    self->ob_type->tp_free((PyObject*) self);
}

PyObject *
node_tree(NodeObject *self)
{
    return Py_BuildValue("O", self->tree);
}

PyObject *
node_taxon_num_get(NodeObject *self)
{
    PyObject *retval;
    uint32_t taxon_num;

    taxon_num = tr_node_taxon_num_get(self->tree->tr, self->node);

    if (taxon_num == CW_TR_NODE_TAXON_NONE)
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
node_taxon_num_set_cargs(NodeObject *self, uint32_t a_taxon_num)
{
    tr_node_taxon_num_set(self->tree->tr, self->node, a_taxon_num);
}

PyObject *
node_taxon_num_set(NodeObject *self, PyObject *args)
{
    PyObject *retval;
    uint32_t taxon_num;

    taxon_num = CW_TR_NODE_TAXON_NONE;
    if (PyArg_ParseTuple(args, "|i", &taxon_num) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    node_taxon_num_set_cargs(self, taxon_num);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

PyObject *
node_edge(NodeObject *self)
{
    PyObject *retval;
    EdgeObject *edge_obj;
    cw_tr_edge_t edge;
    uint32_t end;

    tr_node_edge_get(self->tree->tr, self->node, &edge, &end);
    if (edge != CW_TR_EDGE_NONE)
    {
	edge_obj = (EdgeObject *) tr_edge_aux_get(self->tree->tr, edge);
	retval = Py_BuildValue("(Oi)", edge_obj, end);
    }
    else
    {
	retval = Py_BuildValue("(ss)", NULL, NULL);
    }

    return retval;
}

PyObject *
node_degree(NodeObject *self)
{
    return Py_BuildValue("i", tr_node_degree(self->tree->tr, self->node));
}

static PyMethodDef node_p_methods[] =
{
    {
	"tree",
	(PyCFunction) node_tree,
	METH_NOARGS,
	"tree"
    },
    {
	"taxon_num_get",
	(PyCFunction) node_taxon_num_get,
	METH_NOARGS,
	"taxon_num_get"
    },
    {
	"taxon_num_set",
	(PyCFunction) node_taxon_num_set,
	METH_VARARGS,
	"taxon_num_set"
    },
    {
	"edge",
	(PyCFunction) node_edge,
	METH_NOARGS,
	"edge"
    },
    {
	"degree",
	(PyCFunction) node_degree,
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
    sizeof(NodeObject),	/* int tp_basicsize */
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

static PyObject *tree_p_node_new_code;

// XXX Move in file (do tree and edge as well).
NodeObject *
node_new(TreeObject *a_tree)
{
    NodeObject *retval;
    PyObject *globals, *locals;

    globals = PyEval_GetGlobals();
    locals = Py_BuildValue("{sO}", "tree", (PyObject *) a_tree);

    retval = (NodeObject *)
	PyEval_EvalCode((PyCodeObject *) tree_p_node_new_code,
			globals,
			locals);

    Py_DECREF(locals);

    return retval;
}

void
crux_node_init(void)
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
    tree_p_node_new_code = Py_CompileString("crux.node(tree)",
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
#ifdef CW_CC_SILENCE
	= NULL
#endif
	;
    EdgeObject *self;
    TreeObject *tree;

    if (PyArg_ParseTuple(args, "O!", &Tree_Type, &tree) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    self = (EdgeObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    xep_begin();
    xep_try
    {
	Py_INCREF(tree);
	self->tree = tree;
	self->edge = tr_edge_new(tree->tr);
	tr_edge_aux_set(tree->tr, self->edge, self);

	retval = (PyObject *) self;
    }
    xep_catch(CW_CRUXX_OOM)
    {
	xep_handled();
	retval = PyErr_NoMemory();
    }
    xep_end();

    RETURN:
    return retval;
}

static int
edge_p_traverse(EdgeObject *self, visitproc visit, void *arg)
{
    int retval;
    cw_tr_node_t tr_node;

    if (visit((PyObject *) self->tree, arg) < 0)
    {
	retval = -1;
	goto RETURN;
    }

    tr_node = tr_edge_node_get(self->tree->tr, self->edge, 0);
    if (tr_node != CW_TR_NODE_NONE)
    {
	NodeObject *node;

	node = (NodeObject *) tr_node_aux_get(self->tree->tr, tr_node);
	if (visit((PyObject *) node, arg) < 0)
	{
	    retval = -1;
	    goto RETURN;
	}

	tr_node = tr_edge_node_get(self->tree->tr, self->edge, 1);
	node = (NodeObject *) tr_node_aux_get(self->tree->tr, tr_node);
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
edge_p_clear(EdgeObject *self)
{
    cw_tr_node_t node_a;

    node_a = tr_edge_node_get(self->tree->tr, self->edge, 0);
    if (node_a != CW_TR_NODE_NONE)
    {
	cw_tr_node_t node_b;
	NodeObject *node_a_obj, *node_b_obj;

	node_a_obj = (NodeObject *) tr_node_aux_get(self->tree->tr, node_a);

	node_b = tr_edge_node_get(self->tree->tr, self->edge, 1);
	node_b_obj = (NodeObject *) tr_node_aux_get(self->tree->tr, node_b);

	tr_edge_detach(self->tree->tr, self->edge);
	Py_DECREF(node_a_obj);
	Py_DECREF(node_b_obj);
    }

    tr_edge_delete(self->tree->tr, self->edge);
    Py_DECREF(self->tree);

    return 0;
}

static void
edge_p_delete(EdgeObject *self)
{
    edge_p_clear(self);
    self->ob_type->tp_free((PyObject*) self);
}

PyObject *
edge_tree(EdgeObject *self)
{
    return Py_BuildValue("O", self->tree);
}

PyObject *
edge_node_cargs(EdgeObject *self, uint32_t a_ind)
{
    PyObject *retval;
    cw_tr_node_t node;

    node = tr_edge_node_get(self->tree->tr, self->edge, a_ind);
    if (node != CW_TR_NODE_NONE)
    {
	retval = (PyObject *) tr_node_aux_get(self->tree->tr, node);
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
edge_node(EdgeObject *self, PyObject *args)
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

    retval = edge_node_cargs(self, ind);

    RETURN:
    return retval;
}

void
edge_next_cargs(EdgeObject *self, uint32_t a_ind, EdgeObject **r_edge,
		uint32_t *r_next_end)
{
    cw_tr_edge_t next_edge;

    tr_edge_next_get(self->tree->tr, self->edge, a_ind, &next_edge, r_next_end);
    *r_edge = (EdgeObject *) tr_edge_aux_get(self->tree->tr, next_edge);
}

PyObject *
edge_next(EdgeObject *self, PyObject *args)
{
    PyObject *retval;
    EdgeObject *edge_obj;
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

    edge_next_cargs(self, ind, &edge_obj, &next_end);

    retval = Py_BuildValue("(Oi)", edge_obj, next_end);
    RETURN:
    return retval;
}

void
edge_prev_cargs(EdgeObject *self, uint32_t a_ind, EdgeObject **r_edge,
		uint32_t *r_prev_end)
{
    cw_tr_edge_t prev_edge;

    tr_edge_prev_get(self->tree->tr, self->edge, a_ind, &prev_edge, r_prev_end);
    *r_edge = (EdgeObject *) tr_edge_aux_get(self->tree->tr, prev_edge);
}

PyObject *
edge_prev(EdgeObject *self, PyObject *args)
{
    PyObject *retval;
    EdgeObject *edge_obj;
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

    edge_prev_cargs(self, ind, &edge_obj, &prev_end);

    retval = Py_BuildValue("(Oi)", edge_obj, prev_end);
    RETURN:
    return retval;
}

PyObject *
edge_length_get(EdgeObject *self)
{
    return Py_BuildValue("d", tr_edge_length_get(self->tree->tr, self->edge));
}

void
edge_length_set_cargs(EdgeObject *self, double a_length)
{
    tr_edge_length_set(self->tree->tr, self->edge, a_length);
}

PyObject *
edge_length_set(EdgeObject *self, PyObject *args)
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

    edge_length_set_cargs(self, length);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

void
edge_attach_cargs(EdgeObject *self, NodeObject *a_node_a, NodeObject *a_node_b)
{
    Py_INCREF(a_node_a);
    Py_INCREF(a_node_b);
    /* Cyclic references (nodes refer to edge). */
    Py_INCREF(self);
    Py_INCREF(self);
    tr_edge_attach(self->tree->tr, self->edge, a_node_a->node, a_node_b->node);
}

PyObject *
edge_attach(EdgeObject *self, PyObject *args)
{
    PyObject *retval;
    NodeObject *node_a, *node_b;

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

    edge_attach_cargs(self, node_a, node_b);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

PyObject *
edge_detach(EdgeObject *self)
{
    cw_tr_node_t node_a;

    node_a = tr_edge_node_get(self->tree->tr, self->edge, 0);
    if (node_a != CW_TR_NODE_NONE)
    {
	cw_tr_node_t node_b;
	NodeObject *node_a_obj, *node_b_obj;

	node_a_obj = (NodeObject *) tr_node_aux_get(self->tree->tr, node_a);

	node_b = tr_edge_node_get(self->tree->tr, self->edge, 1);
	node_b_obj = (NodeObject *) tr_node_aux_get(self->tree->tr, node_b);

	tr_edge_detach(self->tree->tr, self->edge);
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
	(PyCFunction) edge_tree,
	METH_NOARGS,
	"tree"
    },
    {
	"node",
	(PyCFunction) edge_node,
	METH_VARARGS,
	"node"
    },
    {
	"next",
	(PyCFunction) edge_next,
	METH_VARARGS,
	"next"
    },
    {
	"prev",
	(PyCFunction) edge_prev,
	METH_VARARGS,
	"prev"
    },
    {
	"length_get",
	(PyCFunction) edge_length_get,
	METH_NOARGS,
	"length_get"
    },
    {
	"length_set",
	(PyCFunction) edge_length_set,
	METH_VARARGS,
	"length_set"
    },
    {
	"attach",
	(PyCFunction) edge_attach,
	METH_VARARGS,
	"attach"
    },
    {
	"detach",
	(PyCFunction) edge_detach,
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
    sizeof(EdgeObject),	/* int tp_basicsize */
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

static PyObject *tree_p_edge_new_code;

EdgeObject *
edge_new(TreeObject *a_tree)
{
    EdgeObject *retval;
    PyObject *globals, *locals;

    globals = PyEval_GetGlobals();
    locals = Py_BuildValue("{sO}", "tree", (PyObject *) a_tree);

    retval = (EdgeObject *)
	PyEval_EvalCode((PyCodeObject *) tree_p_edge_new_code,
			globals,
			locals);

    Py_DECREF(locals);

    return retval;
}

void
crux_edge_init(void)
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
    tree_p_edge_new_code = Py_CompileString("crux.edge(tree)",
					    "<string>",
					    Py_eval_input);
}

/* End edge. */
/******************************************************************************/
