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

static PyTypeObject CxtDistMatrix;

static PyObject *
CxpDistMatrixNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtDistMatrixObject *self;
    PyObject *input;
    CxDistMatrixInputType inputType;

    input = NULL;
    if (PyArg_ParseTuple(args, "O", &input) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (PyFile_Check(input))
    {
	inputType = CxDistMatrixInputFile;
    }
    else if (PyString_Check(input))
    {
	inputType = CxDistMatrixInputString;
    }
    else if (input == NULL /* XXX TaxonMap */)
    {
	inputType = CxDistMatrixInputTaxonMap;
    }
    else
    {
	PyErr_SetString(CxgDistMatrixTypeError,
			"input: file, string, or TaxonMap expected");
    }

    self = (CxtDistMatrixObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    switch (inputType)
    {
	case CxDistMatrixInputFile:
	case CxDistMatrixInputString:
	{
	    self->buf = NULL;
	    self->bufLen = 0;
	    self->tokenLen = 0;
	    self->line = 0;
	    // XXX Parse.
	    break;
	}
	case CxDistMatrixInputTaxonMap:
	{
	    self->map = input;
	    Py_INCREF(self->map);
	}
	default:
	{
	    CxmNotReached();
	}
    }

    retval = (PyObject *) self;
    RETURN:
    return retval;
}

static int
CxpDistMatrixTraverse(CxtDistMatrixObject *self, visitproc visit, void *arg)
{
    // XXX Traverse map.
    return 0;
}

static int
CxpDistMatrixClear(CxtDistMatrixObject *self)
{
    // XXX Clear map.
    return 0;
}

static void
CxpDistMatrixDelete(CxtDistMatrixObject *self)
{
    CxpDistMatrixClear(self);
    if (self->buf != NULL)
    {
	free(self->buf);
    }
    self->ob_type->tp_free((PyObject*) self);
}

CxmpInline void
CxpDistMatrixAppendC(CxtDistMatrixObject *self, char c)
{
    /* Expand the token buffer, if necessary. */
    if (self->tokenLen + 1 > self->bufLen)
    {
	if (self->buf == NULL)
	{
#define CxmpDistMatrixBufLenStart 128
	    self->buf = (char *) CxmMalloc(CxmpDistMatrixBufLenStart);
	    self->bufLen = CxmpDistMatrixBufLenStart;
#undef CxmpDistMatrixBufLenStart
	}
	else
	{
	    self->buf = (char *) CxmRealloc(self->buf, self->bufLen * 2);
	    self->bufLen *= 2;
	}
    }

    /* Append the character. */
    self->buf[self->tokenLen] = c;
    self->tokenLen++;
}

CxmpInline bool
CxpDistMatrixGetC(CxtDistMatrixObject *self, char *r_c, int *r_line,
		   int *r_column)
{
    bool retval;
    char c;
    int line = self->line;
    int column = self->column;

    if (self->fileInput)
    {
	if (fread(&c, 1, 1, self->i.f.file) == 0)
	{
	    retval = true;
	    goto RETURN;
	}
    }
    else
    {
	c = self->i.s.string[self->offset];
	self->offset++;
	if (c == '\0')
	{
	    retval = true;
	    goto RETURN;
	}
    }

    /* Update line and column info. */
    if (c == '\n')
    {
	self->line++;
	self->column = 0;
    }
    else
    {
	self->column++;
    }

    /* Set returns. */
    *r_c = c;
    *r_line = line;
    *r_column = column;

    retval = false;
    RETURN:
    return retval;
}

static void
CxpSyntaxError(CxtDistMatrixObject *self, const char *format, ...)
{
    va_list ap;
    char *str;

    va_start(ap, format);
    vasprintf(&str, format, ap);
    va_end(ap);

    PyErr_SetString(CxgDistMatrixSyntaxError, str);
    free(str);
}

PyObject *
CxDistMatrixNtaxaGet(CxtDistMatrixObject *self)
{
    return Py_BuildValue("i", self->ntaxa);
}

PyObject *
CxDistMatrixTaxonMapGet(CxtDistMatrixObject *self)
{
    Py_INCREF(self->map);
    return self->map;
}

PyObject *
CxDistMatrixDistanceGet(CxtDistMatrixObject *self, PyObject *args)
{
    PyObject *retval;

    CxmError("XXX Not implemented");

    return retval;
}

PyObject *
CxDistMatrixDistanceSet(CxtDistMatrixObject *self, PyObject *args)
{
    PyObject *retval;

    CxmError("XXX Not implemented");

    return retval;
}

PyObject *
CxDistMatrixPrint(CxtDistMatrixObject *self)
{
    PyObject *retval;

    CxmError("XXX Not implemented");

    return retval;
}

PyObject *
CxDistMatrixPrints(CxtDistMatrixObject *self)
{
    PyObject *retval;

    CxmError("XXX Not implemented");

    return retval;
}

PyObject *
CxDistMatrixLine(CxtDistMatrixObject *self)
{
    return Py_BuildValue("i", self->line);
}

static PyMethodDef CxpDistMatrixMethods[] =
{
    {
	"ntaxaGet",
	(PyCFunction) CxDistMatrixNtaxaGet,
	METH_NOARGS,
	"ntaxaGet"
    },
    {
	"taxonMapGet",
	(PyCFunction) CxDistMatrixTaxonMapGet,
	METH_NOARGS,
	"taxonMapGet"
    },
    {
	"distanceGet",
	(PyCFunction) CxDistMatrixDistanceGet,
	METH_VARARGS,
	"distanceGet"
    },
    {
	"distanceSet",
	(PyCFunction) CxDistMatrixDistanceSet,
	METH_VARARGS,
	"distanceSet"
    },
    {
	"print",
	(PyCFunction) CxDistMatrixPrint,
	METH_VARARGS,
	"print"
    },
    {
	"prints",
	(PyCFunction) CxDistMatrixPrints,
	METH_VARARGS,
	"prints"
    },
    {NULL, NULL}
};

static PyTypeObject CxtDistMatrix =
{
    PyObject_HEAD_INIT(NULL)
    0,			/* int ob_size */
    "C_DistMatrix.C_DistMatrix",	/* char *tp_name */
    sizeof(CxtDistMatrixObject),	/* int tp_basicsize */
    0,			/* int tp_itemsize */
    (destructor) CxpDistMatrixDelete,	/* destructor tp_dealloc */
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
    "DistMatrix(): Create the C portion of a tree.",	/* char *tp_doc */
    (traverseproc) CxpDistMatrixTraverse,	/* traverseproc tp_traverse */
    (inquiry) CxpDistMatrixClear,	/* inquiry tp_clear */
    0,			/* richcmpfunc tp_richcompare */
    0,			/* long tp_weaklistoffset */
    0,			/* getiterfunc tp_iter */
    0,			/* iternextfunc tp_iternext */
    CxpDistMatrixMethods,	/* struct PyMethodDef *tp_methods */
    0,			/* struct PyMemberDef *tp_members */
    0,			/* struct PyGetSetDef *tp_getset */
    0,			/* struct _typeobject *tp_base */
    0,			/* PyObject *tp_dict */
    0,			/* descrgetfunc tp_descr_get */
    0,			/* descrsetfunc tp_descr_set */
    0,			/* long tp_dictoffset */
    0,			/* initproc tp_init */
    0,			/* allocfunc tp_alloc */
    CxpDistMatrixNew,		/* newfunc tp_new */
    _PyObject_Del,	/* freefunc tp_free */
    0			/* inquiry tp_is_gc */
};

static PyMethodDef CxpDistMatrixFuncs[] =
{
    {NULL}
};

PyObject *CxgDistMatrixException;
PyObject *CxgDistMatrixValueError;
PyObject *CxgDistMatrixTypeError;
PyObject *CxgDistMatrixSyntaxError;

void
CxDistMatrixInit(void)
{
    PyObject *m;

    /* Create new type. */
    if (PyType_Ready(&CxtDistMatrix) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_DistMatrix", CxpDistMatrixFuncs,
		       "DistMatrix extensions");
    Py_INCREF(&CxtDistMatrix);
    PyModule_AddObject(m, "C_DistMatrix", (PyObject *) &CxtDistMatrix);

    /* Create exception objects. */
    /* Exception. */
    CxgDistMatrixException = PyErr_NewException("C_DistMatrix.Exception",
						 CxgException,
						 NULL);
    Py_INCREF(CxgDistMatrixException);
    PyModule_AddObject(m, "Exception", CxgDistMatrixException);

    /* ValueError. */
    CxgDistMatrixValueError = PyErr_NewException("C_DistMatrix.ValueError",
						  CxgDistMatrixException,
						  NULL);
    Py_INCREF(CxgDistMatrixValueError);
    PyModule_AddObject(m, "ValueError", CxgDistMatrixValueError);

    /* TypeError. */
    CxgDistMatrixTypeError = PyErr_NewException("C_DistMatrix.TypeError",
						  CxgDistMatrixException,
						  NULL);
    Py_INCREF(CxgDistMatrixTypeError);
    PyModule_AddObject(m, "TypeError", CxgDistMatrixTypeError);

    /* SyntaxError. */
    CxgDistMatrixSyntaxError = PyErr_NewException("C_DistMatrix.SyntaxError",
						   CxgDistMatrixException,
						   NULL);
    Py_INCREF(CxgDistMatrixSyntaxError);
    PyModule_AddObject(m, "SyntaxError", CxgDistMatrixSyntaxError);
}
