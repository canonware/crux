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

static void
CxpDistMatrixError(PyObject *exception, const char *format, ...)
{
    va_list ap;
    char *str;

    va_start(ap, format);
    vasprintf(&str, format, ap);
    va_end(ap);

    PyErr_SetString(exception, str);
    free(str);
}

CxmInline void
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
CxpDistMatrixGetC(CxtDistMatrixObject *self, char *rC, long *rLine,
		  long *rColumn)
{
    bool retval;
    char c;
    long line = self->line;
    long column = self->column;

    if (self->inputType == CxDistMatrixInputFile)
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

static bool
CxpDistMatrixGetToken(CxtDistMatrixObject *self,
		      CxtDistMatrixTokenType *rTokenType, long *rLine,
		      long *rColumn)
{
    bool retval = false;
    enum
    {
	CxDistMatrixTokenStateStart,
	CxDistMatrixTokenStateInt,
	CxDistMatrixTokenStateDec,
	CxDistMatrixTokenStateExp,
	CxDistMatrixTokenStateExpCont,
	CxDistMatrixTokenStateLabel
    } CxDistMatrixTokenState;
    char c;
    long line, column;

    self->tokenLen = 0;

    CxDistMatrixTokenState = CxDistMatrixTokenStateStart;
    while (true)
    {
	if (CxpDistMatrixGetC(self, &c, &line, &column))
	{
	    /* End of input.  Accept final token. */
	    retval = true;
	    switch (CxDistMatrixTokenState)
	    {
		case CxDistMatrixTokenStateStart:
		{
		    *rTokenType = CxtDistMatrixTokenNone;
		    goto RETURN;
		}
		case CxDistMatrixTokenStateInt:
		{
		    *rTokenType = CxtDistMatrixTokenInt;
		    goto RETURN;
		}
		case CxDistMatrixTokenStateDec:
		case CxDistMatrixTokenStateExp:
		case CxDistMatrixTokenStateExpCont:
		{
		    *rTokenType = CxtDistMatrixTokenDec;
		    goto RETURN;
		}
		case CxDistMatrixTokenStateLabel:
		{
		    *rTokenType = CxtDistMatrixTokenLabel;
		    goto RETURN;
		}
		default:
		{
		    CxmNotReached();
		}
	    }
	}
	else
	{
	    switch (CxDistMatrixTokenState)
	    {
		case CxDistMatrixTokenStateStart:
		{
		    *rLine = line;
		    *rColumn = column;
		    switch (c)
		    {
			case ' ': case '\t': case '\n': case '\r':
			{
			    /* Stay in the start state. */
			    break;
			}
			case '+': case '-':
			case '0': case '1': case '2': case '3': case '4':
			case '5': case '6': case '7': case '8': case '9':
			{
			    CxpDistMatrixAppendC(self, c);
			    CxDistMatrixTokenState = CxDistMatrixTokenStateInt;
			    break;
			}
			default:
			{
			    CxpDistMatrixAppendC(self, c);
			    CxDistMatrixTokenState
				= CxDistMatrixTokenStateLabel;
			}
		    }
		    break;
		}
		case CxDistMatrixTokenStateInt:
		{
		    switch (c)
		    {
			case ' ': case '\t': case '\n': case '\r':
			{
			    /* Accept. */
			    *rTokenType = CxtDistMatrixTokenInt;
			    goto RETURN;
			}
			case '0': case '1': case '2': case '3': case '4':
			case '5': case '6': case '7': case '8': case '9':
			{
			    CxpDistMatrixAppendC(self, c);
			    break;
			}
			case '.':
			{
			    CxpDistMatrixAppendC(self, c);
			    CxDistMatrixTokenState
				= CxDistMatrixTokenStateDec;
			    break;
			}
			case 'e': case 'E':
			{
			    CxpDistMatrixAppendC(self, c);
			    CxDistMatrixTokenState
				= CxDistMatrixTokenStateExp;
			    break;
			}
			default:
			{
			    CxpDistMatrixAppendC(self, c);
			    CxDistMatrixTokenState
				= CxDistMatrixTokenStateLabel;
			}
		    }
		    break;
		}
		case CxDistMatrixTokenStateDec:
		{
		    switch (c)
		    {
			case ' ': case '\t': case '\n': case '\r':
			{
			    /* Accept. */
			    *rTokenType = CxtDistMatrixTokenDec;
			    goto RETURN;
			}
			case '0': case '1': case '2': case '3': case '4':
			case '5': case '6': case '7': case '8': case '9':
			{
			    CxpDistMatrixAppendC(self, c);
			    break;
			}
			case 'e': case 'E':
			{
			    CxpDistMatrixAppendC(self, c);
			    CxDistMatrixTokenState
				= CxDistMatrixTokenStateExp;
			    break;
			}
			default:
			{
			    CxpDistMatrixAppendC(self, c);
			    CxDistMatrixTokenState
				= CxDistMatrixTokenStateLabel;
			}
		    }
		    break;
		}
		case CxDistMatrixTokenStateExp:
		{
		    switch (c)
		    {
			case ' ': case '\t': case '\n': case '\r':
			{
			    /* Accept. */
			    *rTokenType = CxtDistMatrixTokenDec;
			    goto RETURN;
			}
			case '+': case '-':
			case '0': case '1': case '2': case '3': case '4':
			case '5': case '6': case '7': case '8': case '9':
			{
			    CxpDistMatrixAppendC(self, c);
			    break;
			}
			default:
			{
			    CxpDistMatrixAppendC(self, c);
			    CxDistMatrixTokenState
				= CxDistMatrixTokenStateLabel;
			}
		    }
		    break;
		}
		case CxDistMatrixTokenStateExpCont:
		{
		    switch (c)
		    {
			case ' ': case '\t': case '\n': case '\r':
			{
			    /* Accept. */
			    *rTokenType = CxtDistMatrixTokenDec;
			    goto RETURN;
			}
			case '0': case '1': case '2': case '3': case '4':
			case '5': case '6': case '7': case '8': case '9':
			{
			    CxpDistMatrixAppendC(self, c);
			    break;
			}
			default:
			{
			    CxpDistMatrixAppendC(self, c);
			    CxDistMatrixTokenState
				= CxDistMatrixTokenStateLabel;
			}
		    }
		    break;
		}
		case CxDistMatrixTokenStateLabel:
		{
		    switch (c)
		    {
			case ' ': case '\t': case '\n': case '\r':
			{
			    /* Accept. */
			    *rTokenType = CxtDistMatrixTokenLabel;
			    goto RETURN;
			}
			default:
			{
			    CxpDistMatrixAppendC(self, c);
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
    }
    CxmNotReached();

    RETURN:
    return retval;
}

static bool
CxpDistMatrixTokenToDistance(CxtDistMatrixObject *self, double *rDistance)
{
    bool retval;
    double distance;

    CxpDistMatrixAppendC(self, '\0');

    errno = 0;
    distance = strtod(self->buf, NULL);
    if (errno == ERANGE)
    {
	retval = true;
	goto RETURN;
    }

    *rDistance = distance;
    retval = false;
    RETURN:
    return retval;
}

static bool
CxpDistMatrixSetLabel(CxtDistMatrixObject *self, long index,
		      const char *matrixFormat)
{
    bool retval, eof;
    CxtDistMatrixTokenType tokenType;
    long line, column;

    eof = CxpDistMatrixGetToken(self, &tokenType, &line, &column);
    switch (tokenType)
    {
	case CxtDistMatrixTokenNone:
	{
	    CxpDistMatrixError(CxgDistMatrixSyntaxError,
			       "At %ld:%ld (%s matrix format): End of input"
			       " reached",
			       self->line, self->column, matrixFormat);
	    retval = true;
	    goto RETURN;
	}
	case CxtDistMatrixTokenInt:
	case CxtDistMatrixTokenDec:
	{
	    CxpDistMatrixError(CxgDistMatrixSyntaxError,
			       "At %ld:%ld (%s matrix format): Missing taxon"
			       " label",
			       line, column, matrixFormat);
	    retval = true;
	    goto RETURN;
	}
	case CxtDistMatrixTokenLabel:
	{
	    /* Insert label into map. */
	    PyEval_CallMethod(self->map, "map", "s#i", self->buf,
			      self->tokenLen, index);
	    break;
	}
	default:
	{
	    CxmNotReached();
	}
    }
    if (eof)
    {
	CxpDistMatrixError(CxgDistMatrixSyntaxError,
			   "At %ld:%ld (%s matrix format): End of input reached"
			   " while reading label %ld",
			   self->line, self->column, matrixFormat, index);
	retval = true;
	goto RETURN;
    }

    retval = false;
    RETURN:
    return retval;
}

static bool
CxpDistMatrixSetDistance(CxtDistMatrixObject *self, long x, long y,
			 const char *matrixFormat)
{
    bool retval, eof;
    CxtDistMatrixTokenType tokenType;
    long line, column;

    eof = CxpDistMatrixGetToken(self, &tokenType, &line, &column);
    switch (tokenType)
    {
	case CxtDistMatrixTokenNone:
	{
	    CxpDistMatrixError(CxgDistMatrixSyntaxError,
			       "At %ld:%ld (%s matrix format): End of input"
			       " reached (%ld, %ld)",
			       self->line, self->column, matrixFormat, x, y);
	    retval = true;
	    goto RETURN;
	}
	case CxtDistMatrixTokenInt:
	case CxtDistMatrixTokenDec:
	{
	    if (CxpDistMatrixTokenToDistance(self,
					     &self->matrix[x * self->ntaxa
							   + y]))
	    {
		CxpDistMatrixError(CxgDistMatrixSyntaxError,
				   "At %ld:%ld (%s matrix format): Distance"
				   " conversion error (%ld, %ld) '%s'",
				   line, column, matrixFormat,
				   x, y, self->buf);
		retval = true;
		goto RETURN;
	    }
	    break;
	}
	case CxtDistMatrixTokenLabel:
	{
	    CxpDistMatrixError(CxgDistMatrixSyntaxError,
			       "At %ld:%ld (%s matrix format):"
			       " Missing distance (%ld, %ld)",
			       line, column, matrixFormat, x, y);
	    retval = true;
	    goto RETURN;
	}
	default:
	{
	    CxmNotReached();
	}
    }
    if (eof)
    {
	CxpDistMatrixError(CxgDistMatrixSyntaxError,
			   "At %ld:%ld (%s matrix format): End of input reached"
			   " (%ld, %ld)",
			   self->line, self->column, x, y);
	retval = true;
	goto RETURN;
    }

    retval = false;
    RETURN:
    return retval;
}

static bool
CxpDistMatrixParse(CxtDistMatrixObject *self)
{
    bool retval, eof;
    enum
    {
	CxDistMatrixFormatUnknown,
	CxDistMatrixFormatFull,
	CxDistMatrixFormatUpper,
	CxDistMatrixFormatLower
    } matrixFormat;
    CxtDistMatrixTokenType tokenType;
    long line, column, x, y;

    self->buf = NULL;
    self->bufLen = 0;
    self->line = 1;
    self->column = 0;

    matrixFormat = CxDistMatrixFormatUnknown;
    
    /* Get the number of taxa. */
    eof = CxpDistMatrixGetToken(self, &tokenType, &line, &column);
    switch (tokenType)
    {
	case CxtDistMatrixTokenNone:
	{
	    CxpDistMatrixError(CxgDistMatrixSyntaxError,
			       "End of input reached before matrix size"
			       " specification");
	    retval = true;
	    goto RETURN;
	}
	case CxtDistMatrixTokenInt:
	{
	    CxpDistMatrixAppendC(self, '\0');
	    self->ntaxa = strtol(self->buf, NULL, 10);
	    if (self->ntaxa < 2)
	    {
		CxpDistMatrixError(CxgDistMatrixSyntaxError,
				   "At %ld:%ld: Too few taxa (%ld)",
				   line, column, self->ntaxa);
		retval = true;
		goto RETURN;
	    }
	    self->matrix = (double *) CxmMalloc(sizeof(double)
						* self->ntaxa * self->ntaxa);
	    break;
	}
	case CxtDistMatrixTokenDec:
	case CxtDistMatrixTokenLabel:
	{
	    CxpDistMatrixError(CxgDistMatrixValueError,
			       "At %ld:%ld: Unspecified number of taxa", line,
			       column);
	    retval = true;
	    goto RETURN;
	}
	default:
	{
	    CxmNotReached();
	}
    }
    if (eof)
    {
	CxpDistMatrixError(CxgDistMatrixSyntaxError,
			   "At %ld:%ld: End of input reached",
			   self->line, self->column);
	retval = true;
	goto RETURN;
    }

    /* Get the first taxon label. */
    if (CxpDistMatrixSetLabel(self, 0, "unknown"))
    {
	retval = true;
	goto RETURN;
    }

    /* Get the next token; if it is a taxon label, then this matrix is in lower
     * triangle format. */
    eof = CxpDistMatrixGetToken(self, &tokenType, &line, &column);
    switch (tokenType)
    {
	case CxtDistMatrixTokenNone:
	{
	    CxpDistMatrixError(CxgDistMatrixSyntaxError,
			       "At %ld:%ld: End of input reached",
			       self->line, self->column);
	    retval = true;
	    goto RETURN;
	}
	case CxtDistMatrixTokenInt:
	case CxtDistMatrixTokenDec:
	{
	    /* This matrix is either in full or upper triangle format. */
	    break;
	}
	case CxtDistMatrixTokenLabel:
	{
	    /* The matrix is in lower triangle format. */
	    matrixFormat = CxDistMatrixFormatLower;
	    break;
	}
	default:
	{
	    CxmNotReached();
	}
    }
    if (eof)
    {
	CxpDistMatrixError(CxgDistMatrixSyntaxError,
			   "At %ld:%ld: End of input reached",
			   self->line, self->column);
	retval = true;
	goto RETURN;
    }
    
    if (matrixFormat == CxDistMatrixFormatLower)
    {
	/* Insert label into map. */
	PyEval_CallMethod(self->map, "map", "s#i", self->buf, self->tokenLen,
			  1);

	/* Get second row of distances. */
	if (CxpDistMatrixSetDistance(self, 1, 0, "lower"))
	{
	    retval = true;
	    goto RETURN;
	}

	/* Get remaining rows. */
	for (x = 2; x < self->ntaxa; x++)
	{
	    /* Get taxon label. */
	    if (CxpDistMatrixSetLabel(self, x, "lower"))
	    {
		retval = true;
		goto RETURN;
	    }

	    /* Get distances. */
	    for (y = 0; y < x; y++)
	    {
		if (CxpDistMatrixSetDistance(self, x, y, "lower"))
		{
		    retval = true;
		    goto RETURN;
		}
	    }
	}

	/* Reflect matrix contents. */
	for (x = 0; x < self->ntaxa; x++)
	{
	    for (y = 0; y < self->ntaxa; y++)
	    {
		self->matrix[x * self->ntaxa + y]
		    = self->matrix[y * self->ntaxa + x];
	    }
	}

	/* Initialize diagonal. */
	for (x = 0; x < self->ntaxa; x++)
	{
	    self->matrix[x * self->ntaxa + x] = 0.0;
	}
    }
    else
    {
	/* Full or upper triangle format. */

	/* Insert the first distance as though parsing an upper-triangle
	 * matrix. */
	if (CxpDistMatrixTokenToDistance(self,
					 &self->matrix[0 * self->ntaxa
						       + 1]))
	{
	    CxpDistMatrixError(CxgDistMatrixSyntaxError,
			       "At %ld:%ld: Distance conversion error (%s)",
			       line, column, self->buf);
	    retval = true;
	    goto RETURN;
	}
	/* Insert the rest of the first row of distances into the matrix as
	 * though parsing an upper-triangle matrix. */
	for (y = 2; y < self->ntaxa; y++)
	{
	    if (CxpDistMatrixSetDistance(self, 0, y, "full/upper"))
	    {
		retval = true;
		goto RETURN;
	    }
	}

	/* Determine whether this is a full or upper-triangle matrix. */
	eof = CxpDistMatrixGetToken(self, &tokenType, &line, &column);
	switch (tokenType)
	{
	    case CxtDistMatrixTokenNone:
	    {
		CxpDistMatrixError(CxgDistMatrixSyntaxError,
				   "At %ld:%ld: End of input reached",
				   self->line, self->column);
		retval = true;
		goto RETURN;
	    }
	    case CxtDistMatrixTokenInt:
	    case CxtDistMatrixTokenDec:
	    {
		/* This matrix is in full-matrix format. */
		matrixFormat = CxDistMatrixFormatFull;
		break;
	    }
	    case CxtDistMatrixTokenLabel:
	    {
		/* The matrix is in upper-triangle format. */
		matrixFormat = CxDistMatrixFormatUpper;
		break;
	    }
	    default:
	    {
		CxmNotReached();
	    }
	}
	if (eof)
	{
	    CxpDistMatrixError(CxgDistMatrixSyntaxError,
			       "At %ld:%ld: End of input reached",
			       self->line, self->column);
	    retval = true;
	    goto RETURN;
	}

	if (matrixFormat == CxDistMatrixFormatUpper)
	{
	    /* This is an upper-triangle matrix. */
	    
	    /* Insert label into map. */
	    PyEval_CallMethod(self->map, "map", "s#i", self->buf,
			      self->tokenLen, 1);

	    /* Get second row of distances. */
	    for (y = 2; y < self->ntaxa; y++)
	    {
		if (CxpDistMatrixSetDistance(self, 1, y, "upper"))
		{
		    retval = true;
		    goto RETURN;
		}
	    }

	    /* Get remaining rows. */
	    for (x = 2; x < self->ntaxa; x++)
	    {
		/* Get taxon label. */
		if (CxpDistMatrixSetLabel(self, x, "upper"))
		{
		    retval = true;
		    goto RETURN;
		}

		/* Get distances. */
		for (y = x + 1; y < self->ntaxa; y++)
		{
		    if (CxpDistMatrixSetDistance(self, x, y, "upper"))
		    {
			retval = true;
			goto RETURN;
		    }
		}
	    }

	    /* Reflect matrix contents. */
	    for (x = 0; x < self->ntaxa; x++)
	    {
		for (y = 0; y < self->ntaxa; y++)
		{
		    self->matrix[y * self->ntaxa + x]
			= self->matrix[x * self->ntaxa + y];
		}
	    }

	    /* Initialize diagonal. */
	    for (x = 0; x < self->ntaxa; x++)
	    {
		self->matrix[x * self->ntaxa + x] = 0.0;
	    }
	}
	else
	{
	    /* This is a full matrix. */

	    /* Shift the contents of the first row back one position. */
	    for (y = 0; y < self->ntaxa - 1; y++)
	    {
		self->matrix[0 * self->ntaxa + y]
		    = self->matrix[0 * self->ntaxa + y + 1];
	    }

	    /* Set the last distance on the first row. */
	    if (CxpDistMatrixTokenToDistance(self,
					     &self->matrix[0 * self->ntaxa
							   + self->ntaxa - 1]))
	    {
		CxpDistMatrixError(CxgDistMatrixSyntaxError,
				   "At %ld:%ld: Distance conversion"
				   " error (%ld, %ld) '%s'",
				   line, column,
				   0, self->ntaxa - 1, self->buf);
		retval = true;
		goto RETURN;
	    }

	    /* Get remaining rows. */
	    for (x = 1; x < self->ntaxa; x++)
	    {
		/* Get taxon label. */
		if (CxpDistMatrixSetLabel(self, x, "full"))
		{
		    retval = true;
		    goto RETURN;
		}

		/* Get Distances. */
		for (y = 0; y < self->ntaxa; y++)
		{
		    if (CxpDistMatrixSetDistance(self, x, y, "full"))
		    {
			retval = true;
			goto RETURN;
		    }
		}
	    }
	}
    }

    retval = false;
    RETURN:
    return retval;
}

static PyObject *
CxpDistMatrixNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtDistMatrixObject *self;

    /* Allocate object. */
    self = (CxtDistMatrixObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    self->buf = NULL;

    self->map = Py_None;
    Py_INCREF(self->map);

    retval = (PyObject *) self;
    RETURN:
    return retval;
}

static int
CxpDistMatrixTraverse(CxtDistMatrixObject *self, visitproc visit, void *arg)
{
    int retval;

    if (visit((PyObject *) self->map, arg) < 0)
    {
	retval = -1;
	goto RETURN;
    }

    retval = 0;
    RETURN:
    return retval;
}

static int
CxpDistMatrixClear(CxtDistMatrixObject *self)
{
    Py_XDECREF(self->map);
    self->map = NULL;

    return 0;
}

static void
CxpDistMatrixDelete(CxtDistMatrixObject *self)
{
    CxpDistMatrixClear(self);
    if (self->buf != NULL)
    {
	free(self->buf);
	self->buf = NULL;
    }
    self->ob_type->tp_free((PyObject*) self);
}

PyObject *
CxDistMatrixParse(CxtDistMatrixObject *self, PyObject *args)
{
    PyObject *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    PyObject *input, *map;

    /* Parse arguments. */
    if (PyArg_ParseTuple(args, "OO", &input, &map) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    /* Determine input type. */
    if (PyFile_Check(input))
    {
	self->inputType = CxDistMatrixInputFile;
	self->i.f.file = PyFile_AsFile(input);
    }
    else if (PyString_Check(input))
    {
	self->inputType = CxDistMatrixInputString;
	self->i.s.string = PyString_AsString(input);
	self->i.s.offset = 0;
    }
    else /* TaxonMap. */
    {
	self->inputType = CxDistMatrixInputTaxonMap;
    }

    Py_DECREF(self->map);
    self->map = map;
    Py_INCREF(self->map);

    /* Handle input. */
    retval = Py_None;
    Py_INCREF(retval);
    CxmXepBegin();
    CxmXepTry
    {
	switch (self->inputType)
	{
	    case CxDistMatrixInputFile:
	    case CxDistMatrixInputString:
	    {
		/* Parse. */
		if (CxpDistMatrixParse(self))
		{
		    Py_DECREF(retval);
		    retval = NULL;
		}
		break;
	    }
	    case CxDistMatrixInputTaxonMap:
	    {
		PyObject *result;
		int ntaxa;
		long i;

		/* Create an uninitialized (zero-filled) distance matrix of the
		 * appropriate size, given the number of taxa in the
		 * TaxonMap. */

		result = PyEval_CallMethod(self->map, "ntaxaGet", "()");
		if (PyArg_ParseTuple(result, "i", &ntaxa) == 0)
		{
		    Py_DECREF(retval);
		    retval = NULL;
		    break;
		}

		self->matrix = (double *) CxmMalloc(sizeof(double)
						    * ntaxa * ntaxa);

		for (i = 0; i < ntaxa * ntaxa; i++)
		{
		    self->matrix[i] = 0.0;
		}

		break;
	    }
	    default:
	    {
		CxmNotReached();
	    }
	}
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	PyErr_NoMemory();
	self = NULL;
    }
    CxmXepEnd();

    RETURN:
    return retval;
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
    int fr, to;

    /* Parse arguments. */
    if (PyArg_ParseTuple(args, "ii", &fr, &to) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (fr < 0 || fr >= self->ntaxa || to < 0 || to >= self->ntaxa)
    {
	CxpDistMatrixError(CxgDistMatrixValueError,
			   "Out of bounds matrix access: (%d, %d)", fr, to);
	retval = NULL;
	goto RETURN;
    }

    /* Get distance. */
    retval = Py_BuildValue("d", self->matrix[fr * self->ntaxa + to]);
    RETURN:
    return retval;
}

PyObject *
CxDistMatrixDistanceSet(CxtDistMatrixObject *self, PyObject *args)
{
    PyObject *retval;
    int fr, to;
    double distance;

    /* Parse arguments. */
    if (PyArg_ParseTuple(args, "ii", &fr, &to, &distance) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (fr < 0 || fr >= self->ntaxa || to < 0 || to >= self->ntaxa)
    {
	CxpDistMatrixError(CxgDistMatrixValueError,
			   "Out of bounds matrix access: (%d, %d)", fr, to);
	retval = NULL;
	goto RETURN;
    }

    /* Set distance. */
    self->matrix[fr * self->ntaxa + to] = distance;
    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

static PyMethodDef CxpDistMatrixMethods[] =
{
    {
	"_parse",
	(PyCFunction) CxDistMatrixParse,
	METH_VARARGS,
	"_parse"
    },
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
