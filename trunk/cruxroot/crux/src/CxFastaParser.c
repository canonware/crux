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

static PyTypeObject CxtFastaParser;

static PyObject *
CxpFastaParserNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtFastaParserObject *self;

    self = (CxtFastaParserObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    self->buf = NULL;
    self->bufLen = 0;
    self->tokenLen = 0;
    self->line = 0;

    retval = (PyObject *) self;
    RETURN:
    return retval;
}

static int
CxpFastaParserTraverse(CxtFastaParserObject *self, visitproc visit, void *arg)
{
    return 0;
}

static int
CxpFastaParserClear(CxtFastaParserObject *self)
{
    return 0;
}

static void
CxpFastaParserDelete(CxtFastaParserObject *self)
{
    CxpFastaParserClear(self);
    if (self->buf != NULL)
    {
	free(self->buf);
    }
    self->ob_type->tp_free((PyObject*) self);
}

CxmpInline void
CxpFastaParserAppendC(CxtFastaParserObject *self, char c)
{
    /* Expand the token buffer, if necessary. */
    if (self->tokenLen + 1 > self->bufLen)
    {
	if (self->buf == NULL)
	{
#define CxmpFastaParserBufLenStart 1024
	    self->buf = (char *) CxmMalloc(CxmpFastaParserBufLenStart);
	    self->bufLen = CxmpFastaParserBufLenStart;
#undef CxmpFastaParserBufLenStart
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
CxpFastaParserGetC(CxtFastaParserObject *self, char *rC, int *rLine,
		   int *rColumn)
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
	c = self->i.s.string[self->i.s.offset];
	self->i.s.offset++;
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
    *rC = c;
    *rLine = line;
    *rColumn = column;

    retval = false;
    RETURN:
    return retval;
}

PyObject *
CxFastaParserParse(CxtFastaParserObject *self, PyObject *args)
{
    PyObject *retval;
    PyObject *input;
    char *charType;
    bool dnaChars;

    if (PyArg_ParseTuple(args, "Os", &input, &charType) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    /* Determine input type. */
    if (PyFile_Check(input))
    {
	self->fileInput = true;
	self->i.f.file = PyFile_AsFile(input);
    }
    else if (PyString_Check(input))
    {
	self->fileInput = false;
	self->i.s.string = PyString_AsString(input);
	self->i.s.offset = 0;
    }
    else
    {
	CxError(CxgFastaParserTypeError, "input: file or string expected");
	retval = NULL;
	goto RETURN;
    }

    /* Determine character data type. */
    if (strcmp(charType, "DNA") == 0)
    {
	dnaChars = true;
    }
    else if (strcmp(charType, "protein") == 0)
    {
	dnaChars = false;
    }
    else
    {
	CxError(CxgFastaParserValueError,
		"charType: 'DNA' or 'protein' expected");
	retval = NULL;
	goto RETURN;
    }

    /* Parse. */
    CxmXepBegin();
    CxmXepTry
    {
	char c;
	int line, column;
	enum
	{
	    CxpStateStart,
	    CxpStateLabel,
	    CxpStateComment,
	    CxpStateChars
	} state;
	PyObject *result;

	state = CxpStateStart;
	self->line = 1;
	self->column = 0;
	self->tokenLen = 0;

	while (CxpFastaParserGetC(self, &c, &line, &column) == false)
	{
	    switch (state)
	    {
		case CxpStateStart:
		{
		    if (c == '>' && column == 0)
		    {
			state = CxpStateLabel;
		    }
		    break;
		}
		case CxpStateLabel:
		{
		    switch (c)
		    {
			case '\n':
			{
			    if (self->tokenLen == 0)
			    {
				CxError(CxgFastaParserSyntaxError,
					"At %d:%d (token '%.*s', char '%c'):"
					" Empty label",
					line, column, self->tokenLen, self->buf,
					c);
				goto ERROR;
			    }

			    result = PyEval_CallMethod((PyObject *) self,
						       "labelAccept", "()");
			    Py_DECREF(result);
			    self->tokenLen = 0;
			    state = CxpStateChars;
			    break;
			}
			case ' ': case '\t':
			{
			    if (self->tokenLen == 0)
			    {
				CxError(CxgFastaParserSyntaxError,
					"At %d:%d (token '%.*s', char '%c'):"
					" Empty label",
					line, column, self->tokenLen, self->buf,
					c);
				goto ERROR;
			    }

			    result = PyEval_CallMethod((PyObject *) self,
						       "labelAccept", "()");
			    Py_DECREF(result);
			    self->tokenLen = 0;
			    state = CxpStateComment;
			    break;
			}
			default:
			{
			    CxpFastaParserAppendC(self, c);
			}
		    }
		    break;
		}
		case CxpStateComment:
		{
		    if (c == '\n')
		    {
			if (self->tokenLen > 0)
			{
			    result = PyEval_CallMethod((PyObject *) self,
						       "commentAccept", "()");
			    Py_DECREF(result);
			    self->tokenLen = 0;
			    state = CxpStateChars;
			}
			else
			{
			    state = CxpStateChars;
			}
		    }
		    else
		    {
			CxpFastaParserAppendC(self, c);
		    }
		    break;
		}
		case CxpStateChars:
		{
		    if (c == '>' && column == 0)
		    {
			if (self->tokenLen == 0)
			{
			    CxError(CxgFastaParserSyntaxError,
				    "At %d:%d (token '%.*s', char '%c'):"
				    " Missing character data",
				    line, column, self->tokenLen, self->buf,
				    c);
			    goto ERROR;
			}

			result = PyEval_CallMethod((PyObject *) self,
						   "charsAccept", "()");
			Py_DECREF(result);
			self->tokenLen = 0;
			state = CxpStateLabel;
		    }
		    else
		    {
			if (dnaChars)
			{
			    switch (c)
			    {
				case 'N': case 'n': case 'X': case 'x':
				case 'V': case 'v': case 'H': case 'h':
				case 'M': case 'm': case 'D': case 'd':
				case 'R': case 'r': case 'W': case 'w':
				case 'A': case 'a': case 'B': case 'b':
				case 'S': case 's': case 'Y': case 'y':
				case 'C': case 'c': case 'K': case 'k':
				case 'G': case 'g': case 'T': case 't':
				case '-':
				{
				    CxpFastaParserAppendC(self, c);
				    break;
				}
				case '\r': case '\n': case '\t': case ' ':
				{
				    break;
				}
				default:
				{
				    CxError(CxgFastaParserSyntaxError,
					    "At %d:%d (token '%.*s',"
					    " char '%c'):"
					    " Invalid DNA character data",
					    line, column,
					    self->tokenLen, self->buf, c);
				    goto ERROR;
				}
			    }
			}
			else
			{
			    switch (c)
			    {
				case 'A': case 'a': case 'B': case 'b':
				case 'C': case 'c': case 'D': case 'd':
				case 'E': case 'e': case 'F': case 'f':
				case 'G': case 'g': case 'H': case 'h':
				case 'I': case 'i': case 'K': case 'k':
				case 'L': case 'l': case 'M': case 'm':
				case 'N': case 'n': case 'P': case 'p':
				case 'Q': case 'q': case 'R': case 'r':
				case 'S': case 's': case 'T': case 't':
				case 'U': case 'u': case 'V': case 'v':
				case 'W': case 'w': case 'X': case 'x':
				case 'Y': case 'y': case 'Z': case 'z':
				case '-':
				{
				    CxpFastaParserAppendC(self, c);
				    break;
				}
				case '\r': case '\n': case '\t': case ' ':
				{
				    break;
				}
				default:
				{
				    CxError(CxgFastaParserSyntaxError,
					    "At %d:%d (token '%.*s',"
					    " char '%c'):"
					    " Invalid protein character data",
					    line, column,
					    self->tokenLen, self->buf, c);
				    goto ERROR;
				}
			    }
			}
		    }
		    break;
		}
		default:
		{
		    CxmNotReached();
		}
	    }
	}

	/* Make sure that the input ends with character data. */
	if (state != CxpStateChars)
	{
	    CxError(CxgFastaParserSyntaxError,
		    "At %d:%d (token '%.*s', char '%c'):"
		    " Input ended while reading label",
		    line, column, self->tokenLen, self->buf, c);
	    goto ERROR;
	}

	/* Accept the last token. */
	if (self->tokenLen == 0)
	{
	    CxError(CxgFastaParserSyntaxError,
		    "At %d:%d (token '%.*s', char '%c'):"
		    " Missing character data",
		    line, column, self->tokenLen, self->buf, c);
	    goto ERROR;
	}
	result = PyEval_CallMethod((PyObject *) self, "charsAccept", "()");
	Py_DECREF(result);

	Py_INCREF(Py_None);
	retval = Py_None;
	break;

	ERROR:
	retval = NULL;
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	PyErr_NoMemory();
	retval = NULL;
    }
    CxmXepEnd();

    RETURN:
    return retval;
}

PyObject *
CxFastaParserToken(CxtFastaParserObject *self)
{
    PyObject *retval;

    if (self->buf != NULL)
    {
	retval = PyString_FromStringAndSize(self->buf, self->tokenLen);
    }
    else
    {
	retval = PyString_FromString("");
    }

    return retval;
}

PyObject *
CxFastaParserLine(CxtFastaParserObject *self)
{
    return Py_BuildValue("i", self->line);
}

static PyMethodDef CxpFastaParserMethods[] =
{
    {
	"_parse",
	(PyCFunction) CxFastaParserParse,
	METH_VARARGS,
	"_parse"
    },
    {
	"token",
	(PyCFunction) CxFastaParserToken,
	METH_NOARGS,
	"token"
    },
    {
	"line",
	(PyCFunction) CxFastaParserLine,
	METH_NOARGS,
	"line"
    },
    {NULL, NULL}
};

static PyTypeObject CxtFastaParser =
{
    PyObject_HEAD_INIT(NULL)
    0,			/* int ob_size */
    "C_FastaParser.C_FastaParser",	/* char *tp_name */
    sizeof(CxtFastaParserObject),	/* int tp_basicsize */
    0,			/* int tp_itemsize */
    (destructor) CxpFastaParserDelete,	/* destructor tp_dealloc */
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
    "FastaParser(): Create the C portion of a tree.",	/* char *tp_doc */
    (traverseproc) CxpFastaParserTraverse,	/* traverseproc tp_traverse */
    (inquiry) CxpFastaParserClear,	/* inquiry tp_clear */
    0,			/* richcmpfunc tp_richcompare */
    0,			/* long tp_weaklistoffset */
    0,			/* getiterfunc tp_iter */
    0,			/* iternextfunc tp_iternext */
    CxpFastaParserMethods,	/* struct PyMethodDef *tp_methods */
    0,			/* struct PyMemberDef *tp_members */
    0,			/* struct PyGetSetDef *tp_getset */
    0,			/* struct _typeobject *tp_base */
    0,			/* PyObject *tp_dict */
    0,			/* descrgetfunc tp_descr_get */
    0,			/* descrsetfunc tp_descr_set */
    0,			/* long tp_dictoffset */
    0,			/* initproc tp_init */
    0,			/* allocfunc tp_alloc */
    CxpFastaParserNew,		/* newfunc tp_new */
    _PyObject_Del,	/* freefunc tp_free */
    0			/* inquiry tp_is_gc */
};

static PyMethodDef CxpFastaParserFuncs[] =
{
    {NULL}
};

PyObject *CxgFastaParserException;
PyObject *CxgFastaParserValueError;
PyObject *CxgFastaParserTypeError;
PyObject *CxgFastaParserSyntaxError;

void
CxFastaParserInit(void)
{
    PyObject *m;

    /* Create new type. */
    if (PyType_Ready(&CxtFastaParser) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_FastaParser", CxpFastaParserFuncs,
		       "FastaParser extensions");
    Py_INCREF(&CxtFastaParser);
    PyModule_AddObject(m, "C_FastaParser", (PyObject *) &CxtFastaParser);

    /* Create exception objects. */
    /* Exception. */
    CxgFastaParserException = PyErr_NewException("C_FastaParser.Exception",
						 CxgException,
						 NULL);
    Py_INCREF(CxgFastaParserException);
    PyModule_AddObject(m, "Exception", CxgFastaParserException);

    /* ValueError. */
    CxgFastaParserValueError = PyErr_NewException("C_FastaParser.ValueError",
						  CxgFastaParserException,
						  NULL);
    Py_INCREF(CxgFastaParserValueError);
    PyModule_AddObject(m, "ValueError", CxgFastaParserValueError);

    /* TypeError. */
    CxgFastaParserTypeError = PyErr_NewException("C_FastaParser.TypeError",
						 CxgFastaParserException,
						 NULL);
    Py_INCREF(CxgFastaParserTypeError);
    PyModule_AddObject(m, "TypeError", CxgFastaParserTypeError);

    /* SyntaxError. */
    CxgFastaParserSyntaxError = PyErr_NewException("C_FastaParser.SyntaxError",
						   CxgFastaParserException,
						   NULL);
    Py_INCREF(CxgFastaParserSyntaxError);
    PyModule_AddObject(m, "SyntaxError", CxgFastaParserSyntaxError);
}
