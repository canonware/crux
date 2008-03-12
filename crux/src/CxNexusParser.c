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

static PyTypeObject CxtNexusParser;
static PyTypeObject CxtNexusParserBase;

//==============================================================================
// Begin CxNexusParser.

static PyObject *
CxpNexusParserNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtNexusParserObject *self;

    self = (CxtNexusParserObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }

    // XXX Initialize internal variables.

    rVal = (PyObject *) self;
    RETURN:
    return rVal;
}

static int
CxpNexusParserTraverse(CxtNexusParserObject *self, visitproc visit, void *arg)
{
    return 0;
}

static int
CxpNexusParserClear(CxtNexusParserObject *self)
{
    return 0;
}

static void
CxpNexusParserDelete(CxtNexusParserObject *self)
{
    // XXX Clean up internal variables.

    self->ob_type->tp_free((PyObject*) self);
}

PyObject *
CxNexusParserInputName(CxtNexusParserObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserTokenGet(CxtNexusParserObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserTokenUnget(CxtNexusParserObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserTextRange(CxtNexusParserObject *self, PyObject *args)
{
    PyObject *rVal;

    if (PyArg_ParseTuple(args, "XXX") == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    RETURN:
    return rVal;
}

PyObject *
CxNexusParserHandlerSearch(CxtNexusParserObject *self, PyObject *args)
{
    PyObject *rVal;

    if (PyArg_ParseTuple(args, "XXX") == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    RETURN:
    return rVal;
}

PyObject *
CxNexusParserRoot(CxtNexusParserObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserPrint(CxtNexusParserObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

static PyMethodDef CxpNexusParserMethods[] =
{
    {
	"inputName",
	(PyCFunction) CxNexusParserInputName,
	METH_NOARGS,
	"inputName"
    },
    {
	"tokenGet",
	(PyCFunction) CxNexusParserTokenGet,
	METH_NOARGS,
	"tokenGet"
    },
    {
	"tokenUnget",
	(PyCFunction) CxNexusParserTokenUnget,
	METH_NOARGS,
	"tokenUnget"
    },
    {
	"textRange",
	(PyCFunction) CxNexusParserTextRange,
	METH_VARARGS,
	"textRange"
    },
    {
	"handlerSearch",
	(PyCFunction) CxNexusParserHandlerSearch,
	METH_VARARGS,
	"handlerSearch"
    },
    {
	"root",
	(PyCFunction) CxNexusParserRoot,
	METH_NOARGS,
	"root"
    },
    {
	"print",
	(PyCFunction) CxNexusParserPrint,
	METH_NOARGS,
	"print"
    },
    {NULL, NULL}
};

static PyTypeObject CxtNexusParser =
{
    PyObject_HEAD_INIT(NULL)
    0,			// int ob_size
    "C_NexusParser.C_NexusParser",	// char *tp_name
    sizeof(CxtNexusParserObject),	// int tp_basicsize
    0,			// int tp_itemsize
    (destructor) CxpNexusParserDelete,	// destructor tp_dealloc
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
    "NexusParser(): Create the C portion of a NEXUS parser.",	// char *tp_doc
    (traverseproc) CxpNexusParserTraverse,	// traverseproc tp_traverse
    (inquiry) CxpNexusParserClear,	// inquiry tp_clear
    0,			// richcmpfunc tp_richcompare
    0,			// long tp_weaklistoffset
    0,			// getiterfunc tp_iter
    0,			// iternextfunc tp_iternext
    CxpNexusParserMethods,	// struct PyMethodDef *tp_methods
    0,			// struct PyMemberDef *tp_members
    0,			// struct PyGetSetDef *tp_getset
    0,			// struct _typeobject *tp_base
    0,			// PyObject *tp_dict
    0,			// descrgetfunc tp_descr_get
    0,			// descrsetfunc tp_descr_set
    0,			// long tp_dictoffset
    0,			// initproc tp_init
    0,			// allocfunc tp_alloc
    CxpNexusParserNew,		// newfunc tp_new
    _PyObject_Del,	// freefunc tp_free
    0			// inquiry tp_is_gc
};

static PyMethodDef CxpNexusParserFuncs[] =
{
    {NULL}
};

void
CxNexusParserInit(void)
{
    PyObject *m;

    // Create new type.
    if (PyType_Ready(&CxtNexusParser) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_NexusParser", CxpNexusParserFuncs,
		       "NexusParser extensions");
    Py_INCREF(&CxtNexusParser);
    PyModule_AddObject(m, "C_NexusParser", (PyObject *) &CxtNexusParser);
}

// End NexusParser.
//==============================================================================
// Begin NexusParserBase.

static PyObject *
CxpNexusParserBaseNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtNexusParserBaseObject *self;

    self = (CxtNexusParserBaseObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }

    // XXX Initialize internal variables.

    rVal = (PyObject *) self;
    RETURN:
    return rVal;
}

static int
CxpNexusParserBaseTraverse(CxtNexusParserBaseObject *self, visitproc visit,
			   void *arg)
{
    return 0;
}

static int
CxpNexusParserBaseClear(CxtNexusParserBaseObject *self)
{
    return 0;
}

static void
CxpNexusParserBaseDelete(CxtNexusParserBaseObject *self)
{
    // XXX Clean up internal variables.

    self->ob_type->tp_free((PyObject*) self);
}

PyObject *
CxNexusParserBaseInput(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseText(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseOffset(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseLine(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseColumn(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseLengthGet(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseLengthSet(CxtNexusParserBaseObject *self, PyObject *args)
{
    PyObject *rVal;

    if (PyArg_ParseTuple(args, "XXX") == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    RETURN:
    return rVal;
}

PyObject *
CxNexusParserBaseCText(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseFinish(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseFinished(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBasePrint(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseChildAppend(CxtNexusParserBaseObject *self, PyObject *args)
{
    PyObject *rVal;

    if (PyArg_ParseTuple(args, "XXX") == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    RETURN:
    return rVal;
}

PyObject *
CxNexusParserBaseNext(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBasePrev(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseParent(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseIndex(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseNChildren(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseChild(CxtNexusParserBaseObject *self, PyObject *args)
{
    PyObject *rVal;

    if (PyArg_ParseTuple(args, "XXX") == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    RETURN:
    return rVal;
}

PyObject *
CxNexusParserBaseLeft(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseRight(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

PyObject *
CxNexusParserBaseAccessCallback(CxtNexusParserBaseObject *self)
{
    PyObject *rVal;

    CxmError("XXX Not implemented");
    rVal = NULL; // XXX

    return rVal;
}

static PyMethodDef CxpNexusParserBaseMethods[] =
{
    {
	"input",
	(PyCFunction) CxNexusParserBaseInput,
	METH_NOARGS,
	"input"
    },
    {
	"text",
	(PyCFunction) CxNexusParserBaseText,
	METH_NOARGS,
	"text"
    },
    {
	"offset",
	(PyCFunction) CxNexusParserBaseOffset,
	METH_NOARGS,
	"offset"
    },
    {
	"line",
	(PyCFunction) CxNexusParserBaseLine,
	METH_NOARGS,
	"line"
    },
    {
	"column",
	(PyCFunction) CxNexusParserBaseColumn,
	METH_NOARGS,
	"column"
    },
    {
	"lengthGet",
	(PyCFunction) CxNexusParserBaseLengthGet,
	METH_NOARGS,
	"lengthGet"
    },
    {
	"lengthSet",
	(PyCFunction) CxNexusParserBaseLengthSet,
	METH_VARARGS,
	"lengthSet"
    },
    {
	"cText",
	(PyCFunction) CxNexusParserBaseCText,
	METH_NOARGS,
	"cText"
    },
    {
	"finish",
	(PyCFunction) CxNexusParserBaseFinish,
	METH_NOARGS,
	"finish"
    },
    {
	"finished",
	(PyCFunction) CxNexusParserBaseFinished,
	METH_NOARGS,
	"finished"
    },
    {
	"print",
	(PyCFunction) CxNexusParserBasePrint,
	METH_NOARGS,
	"print"
    },
    {
	"childAppend",
	(PyCFunction) CxNexusParserBaseChildAppend,
	METH_VARARGS,
	"childAppend"
    },
    {
	"next",
	(PyCFunction) CxNexusParserBaseNext,
	METH_NOARGS,
	"next"
    },
    {
	"prev",
	(PyCFunction) CxNexusParserBasePrev,
	METH_NOARGS,
	"prev"
    },
    {
	"parent",
	(PyCFunction) CxNexusParserBaseParent,
	METH_NOARGS,
	"parent"
    },
    {
	"index",
	(PyCFunction) CxNexusParserBaseIndex,
	METH_NOARGS,
	"index"
    },
    {
	"nChildren",
	(PyCFunction) CxNexusParserBaseNChildren,
	METH_NOARGS,
	"nChildren"
    },
    {
	"child",
	(PyCFunction) CxNexusParserBaseChild,
	METH_VARARGS,
	"child"
    },
    {
	"left",
	(PyCFunction) CxNexusParserBaseLeft,
	METH_NOARGS,
	"left"
    },
    {
	"right",
	(PyCFunction) CxNexusParserBaseRight,
	METH_NOARGS,
	"right"
    },
    {
	"accessCallback",
	(PyCFunction) CxNexusParserBaseAccessCallback,
	METH_NOARGS,
	"accessCallback"
    },
    {NULL, NULL}
};

static PyTypeObject CxtNexusParserBase =
{
    PyObject_HEAD_INIT(NULL)
    0,			// int ob_size
    "C_NexusParserBase.C_NexusParserBase",	// char *tp_name
    sizeof(CxtNexusParserBaseObject),	// int tp_basicsize
    0,			// int tp_itemsize
    (destructor) CxpNexusParserBaseDelete,	// destructor tp_dealloc
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
    "NexusParserBase(): Create the C portion of a NEXUS parser.",// char *tp_doc
    (traverseproc) CxpNexusParserBaseTraverse,	// traverseproc tp_traverse
    (inquiry) CxpNexusParserBaseClear,	// inquiry tp_clear
    0,			// richcmpfunc tp_richcompare
    0,			// long tp_weaklistoffset
    0,			// getiterfunc tp_iter
    0,			// iternextfunc tp_iternext
    CxpNexusParserBaseMethods,	// struct PyMethodDef *tp_methods
    0,			// struct PyMemberDef *tp_members
    0,			// struct PyGetSetDef *tp_getset
    0,			// struct _typeobject *tp_base
    0,			// PyObject *tp_dict
    0,			// descrgetfunc tp_descr_get
    0,			// descrsetfunc tp_descr_set
    0,			// long tp_dictoffset
    0,			// initproc tp_init
    0,			// allocfunc tp_alloc
    CxpNexusParserBaseNew,		// newfunc tp_new
    _PyObject_Del,	// freefunc tp_free
    0			// inquiry tp_is_gc
};

static PyMethodDef CxpNexusParserBaseFuncs[] =
{
    {NULL}
};

void
CxNexusParserBaseInit(void)
{
    PyObject *m;

    // Create new type.
    if (PyType_Ready(&CxtNexusParserBase) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_NexusParserBase", CxpNexusParserBaseFuncs,
		       "NexusParserBase extensions");
    Py_INCREF(&CxtNexusParserBase);
    PyModule_AddObject(m, "C_NexusParserBase",
		       (PyObject *) &CxtNexusParserBase);
}

// End NexusParserBase.
//==============================================================================
