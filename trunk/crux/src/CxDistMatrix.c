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

static PyTypeObject CxtDistMatrix;

// Convert from row/column matrix coordinates to array offsets for a symmetric
// matrix.  See CxTreeNj.c for more information about the underlying algorithms.
CxmpInline unsigned long
CxpDistMatrixXy2i(CxtDistMatrixObject *self, unsigned long aX, unsigned long aY)
{
    unsigned long retval;

    CxmAssert(aX < self->ntaxa);
    CxmAssert(aY < self->ntaxa);

    if (self->symmetric)
    {
	if (aX > aY)
	{
	    unsigned long t;

	    t = aX;
	    aX = aY;
	    aY = t;
	}

	retval = self->ntaxa * aX + aY - (((aX + 3) * aX) >> 1) - 1;
    }
    else
    {
	retval = self->ntaxa * aX + aY;
    }

    return retval;
}

CxmInline void
CxpDistMatrixAppendC(CxtDistMatrixObject *self, char c)
{
    // Expand the token buffer, if necessary.
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

    // Append the character.
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

    // Update line and column info.
    if (c == '\n')
    {
	self->line++;
	self->column = 0;
    }
    else
    {
	self->column++;
    }

    // Set returns.
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
	    // End of input.  Accept final token.
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
			    // Stay in the start state.
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
			    // Accept.
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
			    // Accept.
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
			    // Accept.
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
			    // Accept.
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
			    // Accept.
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

static size_t
CxpDistMatrixNtaxaAccept(CxtDistMatrixObject *self)
{
    size_t retval;

    if (self->symmetric)
    {
	retval = sizeof(float) * (CxpDistMatrixXy2i(self,
						    self->ntaxa - 2,
						    self->ntaxa - 1)
				  + 1);
	self->matrix = (float *) CxmMalloc(retval);
    }
    else
    {
	retval = sizeof(float) * self->ntaxa * self->ntaxa;
	self->matrix = (float *) CxmMalloc(retval);
    }

    return retval;
}

static bool
CxpDistMatrixLabelAccept(CxtDistMatrixObject *self, long index)
{
    bool retval;
    PyObject *result;

    result = PyEval_CallMethod(self->map, "map",
			       "s#i", self->buf, self->tokenLen, index);
    if (result == NULL)
    {
	retval = true;
	goto RETURN;
    }
    Py_DECREF(result);

    retval = false;
    RETURN:
    return retval;
}

static bool
CxpDistMatrixTokenToDistance(CxtDistMatrixObject *self, float *rDistance)
{
    bool retval;
    float distance;

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
CxpDistMatrixTokenSetDistance(CxtDistMatrixObject *self, long x, long y)
{
    bool retval;
    float distance;

    if (CxpDistMatrixTokenToDistance(self, &distance))
    {
	retval = true;
	goto RETURN;
    }

    CxDistMatrixDistanceSet(self, x, y, distance);

    retval = false;
    RETURN:
    return retval;
}

static bool
CxpDistMatrixProcessLabel(CxtDistMatrixObject *self, long index,
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
	    CxError(CxgDistMatrixSyntaxError,
		    "At %ld:%ld (%s matrix format): End of input reached",
		    self->line, self->column, matrixFormat);
	    retval = true;
	    goto RETURN;
	}
	case CxtDistMatrixTokenInt:
	case CxtDistMatrixTokenDec:
	{
	    CxError(CxgDistMatrixSyntaxError,
		    "At %ld:%ld (%s matrix format): Missing taxon label",
		    line, column, matrixFormat);
	    retval = true;
	    goto RETURN;
	}
	case CxtDistMatrixTokenLabel:
	{
	    // Insert label into map.
	    if (CxpDistMatrixLabelAccept(self, index))
	    {
		retval = true;
		goto RETURN;
	    }
	    break;
	}
	default:
	{
	    CxmNotReached();
	}
    }
    if (eof)
    {
	CxError(CxgDistMatrixSyntaxError,
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
CxpDistMatrixProcessDistance(CxtDistMatrixObject *self, long x, long y,
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
	    CxError(CxgDistMatrixSyntaxError,
		    "At %ld:%ld (%s matrix format): End of input"
		    " reached (%ld, %ld)",
		    self->line, self->column, matrixFormat, x, y);
	    retval = true;
	    goto RETURN;
	}
	case CxtDistMatrixTokenInt:
	case CxtDistMatrixTokenDec:
	{
	    if (CxpDistMatrixTokenSetDistance(self, x, y))
	    {
		CxError(CxgDistMatrixSyntaxError,
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
	    CxError(CxgDistMatrixSyntaxError,
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
	CxError(CxgDistMatrixSyntaxError,
		"At %ld:%ld (%s matrix format): End of input reached"
		" (%ld, %ld)",
		self->line, self->column, matrixFormat, x, y);
	retval = true;
	goto RETURN;
    }

    retval = false;
    RETURN:
    return retval;
}

static bool
CxpDistMatrixStashDistance(CxtDistMatrixObject *self, float *tdists, long y)
{
    bool retval, eof;
    CxtDistMatrixTokenType tokenType;
    long line, column;

    eof = CxpDistMatrixGetToken(self, &tokenType, &line, &column);
    switch (tokenType)
    {
	case CxtDistMatrixTokenNone:
	{
	    CxError(CxgDistMatrixSyntaxError,
		    "At %ld:%ld (full/upper matrix format): End of input"
		    " reached (0, %ld)",
		    self->line, self->column, y);
	    retval = true;
	    goto RETURN;
	}
	case CxtDistMatrixTokenInt:
	case CxtDistMatrixTokenDec:
	{
	    if (CxpDistMatrixTokenToDistance(self, &tdists[y]))
	    {
		CxError(CxgDistMatrixSyntaxError,
			"At %ld:%ld (full/upper matrix format): Distance"
			" conversion error (0, %ld) '%s'",
			line, column, y, self->buf);
		retval = true;
		goto RETURN;
	    }
	    break;
	}
	case CxtDistMatrixTokenLabel:
	{
	    CxError(CxgDistMatrixSyntaxError,
		    "At %ld:%ld (full/upper matrix format):"
		    " Missing distance (0, %ld)",
		    line, column, y);
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
	CxError(CxgDistMatrixSyntaxError,
		"At %ld:%ld (full/upper matrix format): End of input reached"
		" (0, %ld)",
		self->line, self->column, y);
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

    // Get the number of taxa.
    eof = CxpDistMatrixGetToken(self, &tokenType, &line, &column);
    switch (tokenType)
    {
	case CxtDistMatrixTokenNone:
	{
	    CxError(CxgDistMatrixSyntaxError,
		    "End of input reached before matrix size specification");
	    retval = true;
	    goto RETURN;
	}
	case CxtDistMatrixTokenInt:
	{
	    CxpDistMatrixAppendC(self, '\0');
	    self->ntaxa = strtol(self->buf, NULL, 10);
	    if (self->ntaxa < 2)
	    {
		CxError(CxgDistMatrixSyntaxError,
			"At %ld:%ld: Too few taxa (%ld)",
			line, column, self->ntaxa);
		retval = true;
		goto RETURN;
	    }
	    break;
	}
	case CxtDistMatrixTokenDec:
	case CxtDistMatrixTokenLabel:
	{
	    CxError(CxgDistMatrixValueError,
		    "At %ld:%ld: Unspecified number of taxa",
		    line, column);
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
	CxError(CxgDistMatrixSyntaxError,
		"At %ld:%ld: End of input reached",
		self->line, self->column);
	retval = true;
	goto RETURN;
    }

    // Get the first taxon label.
    if (CxpDistMatrixProcessLabel(self, 0, "unknown"))
    {
	retval = true;
	goto RETURN;
    }

    // Get the next token; if it is a taxon label, then this matrix is in lower
    // triangle format.
    eof = CxpDistMatrixGetToken(self, &tokenType, &line, &column);
    switch (tokenType)
    {
	case CxtDistMatrixTokenNone:
	{
	    CxError(CxgDistMatrixSyntaxError,
		    "At %ld:%ld: End of input reached",
		    self->line, self->column);
	    retval = true;
	    goto RETURN;
	}
	case CxtDistMatrixTokenInt:
	case CxtDistMatrixTokenDec:
	{
	    // This matrix is either in full or upper triangle format.
	    break;
	}
	case CxtDistMatrixTokenLabel:
	{
	    // The matrix is in lower triangle format.
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
	CxError(CxgDistMatrixSyntaxError,
		"At %ld:%ld: End of input reached",
		self->line, self->column);
	retval = true;
	goto RETURN;
    }
    
    if (matrixFormat == CxDistMatrixFormatLower)
    {
	// Allocate a symmetric matrix.
	self->symmetric = true;
	CxpDistMatrixNtaxaAccept(self);

	// Insert label into map.
	if (CxpDistMatrixLabelAccept(self, 1))
	{
	    retval = true;
	    goto RETURN;
	}

	// Get second row of distances.
	if (CxpDistMatrixProcessDistance(self, 1, 0, "lower"))
	{
	    retval = true;
	    goto RETURN;
	}

	// Get remaining rows.
	for (x = 2; x < self->ntaxa; x++)
	{
	    // Get taxon label.
	    if (CxpDistMatrixProcessLabel(self, x, "lower"))
	    {
		retval = true;
		goto RETURN;
	    }

	    // Get distances.
	    for (y = 0; y < x; y++)
	    {
		if (CxpDistMatrixProcessDistance(self, x, y, "lower"))
		{
		    retval = true;
		    goto RETURN;
		}
	    }
	}
    }
    else
    {
	float *tdists;

	// Full or upper triangle format.

	// Allocate a temporary array that is large enough to store the
	// distances until it's possible to tell whether this is a symmetric
	// matrix.
	tdists = (float *) CxmMalloc(sizeof(float) * (self->ntaxa - 1));

	// Insert the first distance as though parsing an upper-triangle
	// matrix.
	if (CxpDistMatrixTokenToDistance(self, &tdists[0]))
	{
	    CxError(CxgDistMatrixSyntaxError,
		    "At %ld:%ld: Distance conversion error (%s)",
		    line, column, self->buf);
	    CxmFree(tdists);
	    retval = true;
	    goto RETURN;
	}
	// Insert the rest of the first row of distances into the matrix as
	// though parsing an upper-triangle matrix.
	for (y = 2; y < self->ntaxa; y++)
	{
	    if (CxpDistMatrixStashDistance(self, tdists, y - 1))
	    {
		CxmFree(tdists);
		retval = true;
		goto RETURN;
	    }
	}

	// Determine whether this is a full or upper-triangle matrix.
	eof = CxpDistMatrixGetToken(self, &tokenType, &line, &column);
	switch (tokenType)
	{
	    case CxtDistMatrixTokenNone:
	    {
		CxError(CxgDistMatrixSyntaxError,
			"At %ld:%ld: End of input reached",
			self->line, self->column);
		CxmFree(tdists);
		retval = true;
		goto RETURN;
	    }
	    case CxtDistMatrixTokenInt:
	    case CxtDistMatrixTokenDec:
	    {
		// This matrix is in full-matrix format.
		matrixFormat = CxDistMatrixFormatFull;
		break;
	    }
	    case CxtDistMatrixTokenLabel:
	    {
		// The matrix is in upper-triangle format.
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
	    CxError(CxgDistMatrixSyntaxError,
		    "At %ld:%ld: End of input reached",
		    self->line, self->column);
	    CxmFree(tdists);
	    retval = true;
	    goto RETURN;
	}

	if (matrixFormat == CxDistMatrixFormatUpper)
	{
	    // This is an upper-triangle matrix.

	    // Allocate a symmetric matrix.
	    self->symmetric = true;
	    CxpDistMatrixNtaxaAccept(self);

	    // Move distances from tdists to matrix.
	    for (y = 0; y < self->ntaxa - 1; y++)
	    {
		CxDistMatrixDistanceSet(self, 0, y + 1, tdists[y]);
	    }

	    CxmFree(tdists);

	    // Insert label into map.
	    if (CxpDistMatrixLabelAccept(self, 1))
	    {
		retval = true;
		goto RETURN;
	    }

	    // Get second row of distances.
	    for (y = 2; y < self->ntaxa; y++)
	    {
		if (CxpDistMatrixProcessDistance(self, 1, y, "upper"))
		{
		    retval = true;
		    goto RETURN;
		}
	    }

	    // Get remaining rows.
	    for (x = 2; x < self->ntaxa; x++)
	    {
		// Get taxon label.
		if (CxpDistMatrixProcessLabel(self, x, "upper"))
		{
		    retval = true;
		    goto RETURN;
		}

		// Get distances.
		for (y = x + 1; y < self->ntaxa; y++)
		{
		    if (CxpDistMatrixProcessDistance(self, x, y, "upper"))
		    {
			retval = true;
			goto RETURN;
		    }
		}
	    }
	}
	else
	{
	    // This is a full matrix.

	    // Allocate a full matrix.
	    self->symmetric = false;
	    CxpDistMatrixNtaxaAccept(self);

	    // Move distances from tdists to matrix.
	    for (y = 0; y < self->ntaxa - 1; y++)
	    {
		CxDistMatrixDistanceSet(self, 0, y, tdists[y]);
	    }
	    CxmFree(tdists);

	    // Set the last distance on the first row.
	    if (CxpDistMatrixTokenSetDistance(self, 0, self->ntaxa - 1))
	    {
		CxError(CxgDistMatrixSyntaxError,
			"At %ld:%ld: Distance conversion error (%ld, %ld) '%s'",
			line, column, 0, self->ntaxa - 1, self->buf);
		retval = true;
		goto RETURN;
	    }

	    // Get remaining rows.
	    for (x = 1; x < self->ntaxa; x++)
	    {
		// Get taxon label.
		if (CxpDistMatrixProcessLabel(self, x, "full"))
		{
		    retval = true;
		    goto RETURN;
		}

		// Get Distances.
		for (y = 0; y < self->ntaxa; y++)
		{
		    if (CxpDistMatrixProcessDistance(self, x, y, "full"))
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

    // Allocate object.
    self = (CxtDistMatrixObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	retval = NULL;
	goto RETURN;
    }

    self->buf = NULL;

    self->map = Py_None;
    Py_INCREF(self->map);

    self->matrix = NULL;

    retval = (PyObject *) self;
    RETURN:
    return retval;
}

static int
CxpDistMatrixTraverse(CxtDistMatrixObject *self, visitproc visit, void *arg)
{
    int retval;

    if (self->map != NULL)
    {
	if (visit((PyObject *) self->map, arg) < 0)
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
    if (self->matrix != NULL)
    {
	free(self->matrix);
	self->matrix = NULL;
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
    int symmetric;

    // Parse arguments.
    if (PyArg_ParseTuple(args, "OOi", &input, &map, &symmetric) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (symmetric != 0 && symmetric != 1)
    {
	CxError(CxgDistMatrixValueError,
		"symmetric: False or True expected");
	retval = NULL;
	goto RETURN;
    }

    // Determine input type.
    if (PyFile_Check(input))
    {
	if (symmetric)
	{
	    CxError(CxgDistMatrixValueError,
		    "symmetric: Automatic for file input");
	    retval = NULL;
	    goto RETURN;
	}
	self->inputType = CxDistMatrixInputFile;
	self->i.f.file = PyFile_AsFile(input);
    }
    else if (PyString_Check(input))
    {
	if (symmetric)
	{
	    CxError(CxgDistMatrixValueError,
		    "symmetric: Automatic for string input");
	    retval = NULL;
	    goto RETURN;
	}
	self->inputType = CxDistMatrixInputString;
	self->i.s.string = PyString_AsString(input);
	self->i.s.offset = 0;
    }
    else // TaxonMap.
    {
	self->inputType = CxDistMatrixInputTaxonMap;
    }

    Py_DECREF(self->map);
    self->map = map;
    Py_INCREF(self->map);

    // Handle input.
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
		// Parse.
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
		long i, nelms;

		// Create an uninitialized (zero-filled) distance matrix of the
		// appropriate size, given the number of taxa in the TaxonMap.

		result = PyEval_CallMethod(self->map, "ntaxaGet", "()");
		if (result == NULL)
		{
		    retval = NULL;
		    break;
		}
		if (PyInt_Check(result) == false)
		{
		    CxError(CxgDistMatrixTypeError,
			    "Integer expected from distMatrix.ntaxaGet()");
		    Py_DECREF(retval);
		    retval = NULL;
		    break;
		}
		self->ntaxa = PyInt_AsLong(result);
		Py_DECREF(result);

		if (symmetric)
		{
		    nelms = CxpDistMatrixXy2i(self,
					      self->ntaxa - 2,
					      self->ntaxa - 1) + 1;
		    self->symmetric = true;
		}
		else
		{
		    nelms = self->ntaxa * self->ntaxa;
		    self->symmetric = false;
		}

		self->matrix = (float *) CxmMalloc(sizeof(float) * nelms);
		for (i = 0; i < nelms; i++)
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

static PyObject *
CxpDistMatrixDup(CxtDistMatrixObject *self, PyObject *args)
{
    PyObject *retval, *map;
    CxtDistMatrixObject *orig;

    if (PyArg_ParseTuple(args, "OO", &orig, &map) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    Py_DECREF(self->map);
    self->map = map;
    Py_INCREF(self->map);

    self->ntaxa = orig->ntaxa;
    self->symmetric = orig->symmetric;

    Py_INCREF(Py_None);
    retval = Py_None;
    CxmXepBegin();
    CxmXepTry
    {
	if (self->ntaxa != 0)
	{
	    size_t nbytes;

	    nbytes = CxpDistMatrixNtaxaAccept(self);
	    memcpy(self->matrix, orig->matrix, nbytes);
	}
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	Py_DECREF(retval);
	retval = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
    return retval;
}

static int
CxpDistMatrixSampleCompare(const void *aA, const void *aB)
{
    long *a = (long *) aA;
    long *b = (long *) aB;

    if (*a < *b)
    {
	return -1;
    }
    else if (*a > *b)
    {
	return 1;
    }
    else
    {
	return 0;
    }
}

static PyObject *
CxpDistMatrixSample(CxtDistMatrixObject *self, PyObject *args)
{
    PyObject *retval, *map, *rows;
    CxtDistMatrixObject *orig;
    long i, j, *rowTab;
    float dist;

    if (PyArg_ParseTuple(args, "OOO!", &orig, &map, &PyList_Type, &rows) == 0)
    {
	retval = NULL;
	goto RETURN;
    }

    Py_DECREF(self->map);
    self->map = map;
    Py_INCREF(self->map);

    self->ntaxa = PyList_Size(rows);
    self->symmetric = orig->symmetric;
    
    Py_INCREF(Py_None);
    retval = Py_None;
    CxmXepBegin();
    CxmXepTry
    {
	CxpDistMatrixNtaxaAccept(self);

	// Create a row table, then sort it so that orig matrix access will be
	// linear.
	rowTab = (long *) CxmMalloc((sizeof(long) << 1) * self->ntaxa);

	for (i = 0; i < self->ntaxa; i++)
	{
	    rowTab[i << 1] = PyInt_AsLong(PyList_GetItem(rows, i));
	    rowTab[(i << 1) + 1] = i;
	}
	qsort(rowTab, self->ntaxa, sizeof(long) << 1,
	      CxpDistMatrixSampleCompare);

	if (self->symmetric)
	{
	    for (i = 0; i < self->ntaxa - 1; i++)
	    {
		for (j = i + 1; j < self->ntaxa; j++)
		{
		    dist = CxDistMatrixDistanceGet(orig,
						   rowTab[i << 1],
						   rowTab[j << 1]);

		    CxDistMatrixDistanceSet(self,
					    rowTab[(i << 1) + 1],
					    rowTab[(j << 1) + 1],
					    dist);
		}
	    }
	}
	else
	{
	    for (i = 0; i < self->ntaxa; i++)
	    {
		for (j = 0; j < self->ntaxa; j++)
		{
		    dist = CxDistMatrixDistanceGet(orig,
						   rowTab[i << 1],
						   rowTab[j << 1]);

		    CxDistMatrixDistanceSet(self,
					    rowTab[(i << 1) + 1],
					    rowTab[(j << 1) + 1],
					    dist);
		}
	    }
	}
	CxmFree(rowTab);
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	Py_DECREF(retval);
	retval = PyErr_NoMemory();
    }
    CxmXepEnd();

    RETURN:
    return retval;
}

PyObject *
CxDistMatrixNtaxaGet(CxtDistMatrixObject *self)
{
    return Py_BuildValue("l", self->ntaxa);
}

PyObject *
CxDistMatrixIsSymmetric(CxtDistMatrixObject *self)
{
    return Py_BuildValue("i", self->symmetric);
}

PyObject *
CxDistMatrixTaxonMapGet(CxtDistMatrixObject *self)
{
    Py_INCREF(self->map);
    return self->map;
}

float
CxDistMatrixDistanceGet(CxtDistMatrixObject *self, long x, long y)
{
    return self->matrix[CxpDistMatrixXy2i(self, x, y)];
}

PyObject *
CxDistMatrixDistanceGetPargs(CxtDistMatrixObject *self, PyObject *args)
{
    PyObject *retval;
    long fr, to;
    float distance;

    // Parse arguments.
    if (PyArg_ParseTuple(args, "ll", &fr, &to) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (fr < 0 || fr >= self->ntaxa || to < 0 || to >= self->ntaxa)
    {
	CxError(CxgDistMatrixValueError,
		"Out of bounds matrix access: (%ld, %ld)",
		fr, to);
	retval = NULL;
	goto RETURN;
    }
    if (self->symmetric && fr == to)
    {
	distance = 0.0;
    }
    else
    {
	distance = CxDistMatrixDistanceGet(self, fr, to);
    }

    // Get distance.
    retval = Py_BuildValue("f", distance);
    RETURN:
    return retval;
}

void
CxDistMatrixDistanceSet(CxtDistMatrixObject *self, long x, long y, float dist)
{
    self->matrix[CxpDistMatrixXy2i(self, x, y)] = dist;
}

PyObject *
CxDistMatrixDistanceSetPargs(CxtDistMatrixObject *self, PyObject *args)
{
    PyObject *retval;
    long fr, to;
    float distance;

    // Parse arguments.
    if (PyArg_ParseTuple(args, "llf", &fr, &to, &distance) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    if (fr < 0 || fr >= self->ntaxa || to < 0 || to >= self->ntaxa)
    {
	CxError(CxgDistMatrixValueError,
		"Out of bounds matrix access: (%ld, %ld)",
		fr, to);
	retval = NULL;
	goto RETURN;
    }
    if (self->symmetric && fr == to)
    {
	CxError(CxgDistMatrixValueError,
		"Out of bounds matrix access for symmetric matrix: (%ld, %ld)",
		fr, to);
	retval = NULL;
	goto RETURN;
    }

    // Set distance.
    CxDistMatrixDistanceSet(self, fr, to, distance);
    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

CxmpInline void
CxpDistMatrixRowsSwap(CxtDistMatrixObject *self, long aA, long aB)
{
    long i;
    float distA, distB;

    if (self->symmetric)
    {
	// Swap rows.  Avoid the diagonal.
	for (i = 0; i < self->ntaxa; i++)
	{
	    if (i != aA && i != aB)
	    {
		distA = CxDistMatrixDistanceGet(self, i, aA);
		distB = CxDistMatrixDistanceGet(self, i, aB);

		CxDistMatrixDistanceSet(self, i, aA, distB);
		CxDistMatrixDistanceSet(self, i, aB, distA);
	    }
	}
    }
    else
    {
	// Swap rows.
	for (i = 0; i < self->ntaxa; i++)
	{
	    distA = CxDistMatrixDistanceGet(self, i, aA);
	    distB = CxDistMatrixDistanceGet(self, i, aB);

	    CxDistMatrixDistanceSet(self, i, aA, distB);
	    CxDistMatrixDistanceSet(self, i, aB, distA);
	}

	// Swap columns.
	for (i = 0; i < self->ntaxa; i++)
	{
	    distA = CxDistMatrixDistanceGet(self, aA, i);
	    distB = CxDistMatrixDistanceGet(self, aB, i);

	    CxDistMatrixDistanceSet(self, aA, i, distB);
	    CxDistMatrixDistanceSet(self, aB, i, distA);
	}
    }
}

static PyObject *
CxpDistMatrixShuffle(CxtDistMatrixObject *self, PyObject *args)
{
    PyObject *retval, *order, *tobj;
    long i, a, b, t, fromRow, *curOrder, *rowTab;

    // Parse arguments.
    if (PyArg_ParseTuple(args, "O!", &PyList_Type, &order) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    CxmAssert(PyList_Size(order) == self->ntaxa);

    // Create a lookup table that maps original row to the current row of the
    // matrix (rowTab), as well as a table that represents the current row order
    // of the matrix (curOrder).  This is needed to keep track of where rows end
    // up as repeated row swaps are done.
    rowTab = (long *) CxmMalloc(sizeof(long) * self->ntaxa);
    curOrder = (long *) CxmMalloc(sizeof(long) * self->ntaxa);
    for (i = 0; i < self->ntaxa; i++)
    {
	rowTab[i] = i;
	curOrder[i] = i;
    }

    // Iteratively swap the correct row into row i.  The last row need not be
    // swapped with itself.
    for (i = 0; i < self->ntaxa - 1; i++)
    {
	tobj = PyList_GetItem(order, i);
	CxmAssert(PyInt_Check(tobj));
	fromRow = PyInt_AsLong(tobj);

	a = i;
	b = rowTab[fromRow];

	CxpDistMatrixRowsSwap(self, a, b);

	// Update curOrder.
	t = curOrder[a];
	curOrder[a] = curOrder[b];
	curOrder[b] = t;

	// Update rowTab.
	t = rowTab[curOrder[a]];
	rowTab[curOrder[a]] = rowTab[curOrder[b]];
	rowTab[curOrder[b]] = t;
    }

    CxmFree(rowTab);
    CxmFree(curOrder);

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

CxmpInline bool
CxpDistMatrixFileLabelRender(CxtDistMatrixObject *self, FILE *aF, long aX)
{
    bool retval;
    PyObject *result;

    result = PyEval_CallMethod(self->map, "labelGet",
			       "(i)", (int) aX);
    if (result == NULL)
    {
	retval = true;
	goto RETURN;
    }
    if (PyString_Check(result) == 0)
    {
	Py_DECREF(result);
	CxError(CxgDistMatrixValueError,
		"map.labelGet(): String expected");
	retval = true;
	goto RETURN;
    }
    if (fprintf(aF, "%-10s", PyString_AsString(result)) < 0)
    {
	Py_DECREF(result);
	CxError(CxgDistMatrixIOError, "Error in fprintf()");
	retval = true;
	goto RETURN;
    }
    Py_DECREF(result);

    retval = false;
    RETURN:
    return retval;
}

static PyObject *
CxpDistMatrixFileRender(CxtDistMatrixObject *self, PyObject *args)
{
    PyObject *retval, *outFile;
    char *format, *distFormat;
    FILE *f;
    long x, y;

    // Parse arguments.
    if (PyArg_ParseTuple(args, "ssO!",
			 &format,
			 &distFormat,
			 &PyFile_Type, &outFile) == 0)
    {
	retval = NULL;
	goto RETURN;
    }
    f = PyFile_AsFile(outFile);

    if (fprintf(f, "%ld\n", self->ntaxa) < 0)
    {
	CxError(CxgDistMatrixIOError, "Error in fprintf()");
	retval = NULL;
	goto RETURN;
    }

    if (strcmp(format, "full") == 0)
    {
	for (x = 0; x < self->ntaxa; x++)
	{
	    if (CxpDistMatrixFileLabelRender(self, f, x))
	    {
		retval = NULL;
		goto RETURN;
	    }

	    for (y = 0; y < self->ntaxa; y++)
	    {
		if (fprintf(f, distFormat, CxDistMatrixDistanceGet(self, x, y))
		    < 0)
		{
		    CxError(CxgDistMatrixIOError, "Error in fprintf()");
		    retval = NULL;
		    goto RETURN;
		}
	    }
	    if (fprintf(f, "\n") < 0)
	    {
		CxError(CxgDistMatrixIOError, "Error in fprintf()");
		retval = NULL;
		goto RETURN;
	    }
	}
    }
    else if (strcmp(format, "upper") == 0)
    {
	for (x = 0; x < self->ntaxa; x++)
	{
	    if (CxpDistMatrixFileLabelRender(self, f, x))
	    {
		retval = NULL;
		goto RETURN;
	    }

	    for (y = 0; y < x + 1; y++)
	    {
		if (fprintf(f, "%8s", "") < 0)
		{
		    CxError(CxgDistMatrixIOError, "Error in fprintf()");
		    retval = NULL;
		    goto RETURN;
		}
	    }
	    for (; y < x; y++)
	    {
		if (fprintf(f, distFormat, CxDistMatrixDistanceGet(self, x, y))
		    < 0)
		{
		    CxError(CxgDistMatrixIOError, "Error in fprintf()");
		    retval = NULL;
		    goto RETURN;
		}
	    }
	    if (fprintf(f, "\n") < 0)
	    {
		CxError(CxgDistMatrixIOError, "Error in fprintf()");
		retval = NULL;
		goto RETURN;
	    }
	}
    }
    else if (strcmp(format, "lower") == 0)
    {
	for (x = 0; x < self->ntaxa; x++)
	{
	    if (CxpDistMatrixFileLabelRender(self, f, x))
	    {
		retval = NULL;
		goto RETURN;
	    }

	    for (y = 0; y < x; y++)
	    {
		if (fprintf(f, distFormat, CxDistMatrixDistanceGet(self, x, y))
		    < 0)
		{
		    CxError(CxgDistMatrixIOError, "Error in fprintf()");
		    retval = NULL;
		    goto RETURN;
		}
	    }
	    if (fprintf(f, "\n") < 0)
	    {
		CxError(CxgDistMatrixIOError, "Error in fprintf()");
		retval = NULL;
		goto RETURN;
	    }
	}
    }

    Py_INCREF(Py_None);
    retval = Py_None;
    RETURN:
    return retval;
}

// Hand off an array of floats that represent an upper-triangle distance matrix,
// and clean up such that this DistMatrix no longer refers to the data (though
// the TaxonMap continues to be referred to by the DistMatrix).
//
// This is not a terribly clean interface; ideally DistMatrix would be
// subclassed and initialization methods overridden, but the Python overhead for
// such a solution is unacceptably high.
void
CxDistMatrixUpperHandoff(CxtDistMatrixObject *self, float **rMatrix,
			 long *rNtaxa)
{
    if (self->symmetric == false)
    {
	unsigned long i;

	// Discard the lower triangle of the matrix by doing a series of
	// memmove() calls, followed by reallocating the matrix.
	self->symmetric = true;
	for (i = 0; i < self->ntaxa - 1; i++)
	{
	    memmove(&self->matrix[CxpDistMatrixXy2i(self, i, i + 1)],
		    &self->matrix[self->ntaxa * i + (i + 1)],
		    sizeof(float) * (self->ntaxa - (i + 1)));
	}

	self->matrix
	    = (float *) CxmRealloc(self->matrix, sizeof(float)
				   * (CxpDistMatrixXy2i(self,
							self->ntaxa - 2,
							self->ntaxa - 1)
				      + 1));
    }

    *rMatrix = self->matrix;
    *rNtaxa = self->ntaxa;

    self->matrix = NULL;
    self->ntaxa = 0;
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
	"isSymmetric",
	(PyCFunction) CxDistMatrixIsSymmetric,
	METH_NOARGS,
	"isSymmetric"
    },
    {
	"_dup",
	(PyCFunction) CxpDistMatrixDup,
	METH_VARARGS,
	"_dup"
    },
    {
	"_sample",
	(PyCFunction) CxpDistMatrixSample,
	METH_VARARGS,
	"_sample"
    },
    {
	"taxonMapGet",
	(PyCFunction) CxDistMatrixTaxonMapGet,
	METH_NOARGS,
	"taxonMapGet"
    },
    {
	"distanceGet",
	(PyCFunction) CxDistMatrixDistanceGetPargs,
	METH_VARARGS,
	"distanceGet"
    },
    {
	"distanceSet",
	(PyCFunction) CxDistMatrixDistanceSetPargs,
	METH_VARARGS,
	"distanceSet"
    },
    {
	"_matrixShuffle",
	(PyCFunction) CxpDistMatrixShuffle,
	METH_VARARGS,
	"_matrixShuffle"
    },
    {
	"_fileRender",
	(PyCFunction) CxpDistMatrixFileRender,
	METH_VARARGS,
	"_fileRender"
    },
    {NULL, NULL}
};

static PyTypeObject CxtDistMatrix =
{
    PyObject_HEAD_INIT(NULL)
    0,			// int ob_size
    "C_DistMatrix.C_DistMatrix",	// char *tp_name
    sizeof(CxtDistMatrixObject),	// int tp_basicsize
    0,			// int tp_itemsize
    (destructor) CxpDistMatrixDelete,	// destructor tp_dealloc
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
    "DistMatrix(): Create the C portion of a tree.",	// char *tp_doc
    (traverseproc) CxpDistMatrixTraverse,	// traverseproc tp_traverse
    (inquiry) CxpDistMatrixClear,	// inquiry tp_clear
    0,			// richcmpfunc tp_richcompare
    0,			// long tp_weaklistoffset
    0,			// getiterfunc tp_iter
    0,			// iternextfunc tp_iternext
    CxpDistMatrixMethods,	// struct PyMethodDef *tp_methods
    0,			// struct PyMemberDef *tp_members
    0,			// struct PyGetSetDef *tp_getset
    0,			// struct _typeobject *tp_base
    0,			// PyObject *tp_dict
    0,			// descrgetfunc tp_descr_get
    0,			// descrsetfunc tp_descr_set
    0,			// long tp_dictoffset
    0,			// initproc tp_init
    0,			// allocfunc tp_alloc
    CxpDistMatrixNew,		// newfunc tp_new
    _PyObject_Del,	// freefunc tp_free
    0			// inquiry tp_is_gc
};

static PyMethodDef CxpDistMatrixFuncs[] =
{
    {NULL}
};

PyObject *CxgDistMatrixException;
PyObject *CxgDistMatrixValueError;
PyObject *CxgDistMatrixTypeError;
PyObject *CxgDistMatrixSyntaxError;
PyObject *CxgDistMatrixIOError;

void
CxDistMatrixInit(void)
{
    PyObject *m;

    // Create new type.
    if (PyType_Ready(&CxtDistMatrix) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_DistMatrix", CxpDistMatrixFuncs,
		       "DistMatrix extensions");
    Py_INCREF(&CxtDistMatrix);
    PyModule_AddObject(m, "C_DistMatrix", (PyObject *) &CxtDistMatrix);

    // Create exception objects.
    // Exception.
    CxgDistMatrixException = PyErr_NewException("C_DistMatrix.Exception",
						 CxgException,
						 NULL);
    Py_INCREF(CxgDistMatrixException);
    PyModule_AddObject(m, "Exception", CxgDistMatrixException);

    // ValueError.
    CxgDistMatrixValueError = PyErr_NewException("C_DistMatrix.ValueError",
						  CxgDistMatrixException,
						  NULL);
    Py_INCREF(CxgDistMatrixValueError);
    PyModule_AddObject(m, "ValueError", CxgDistMatrixValueError);

    // TypeError.
    CxgDistMatrixTypeError = PyErr_NewException("C_DistMatrix.TypeError",
						  CxgDistMatrixException,
						  NULL);
    Py_INCREF(CxgDistMatrixTypeError);
    PyModule_AddObject(m, "TypeError", CxgDistMatrixTypeError);

    // SyntaxError.
    CxgDistMatrixSyntaxError = PyErr_NewException("C_DistMatrix.SyntaxError",
						   CxgDistMatrixException,
						   NULL);
    Py_INCREF(CxgDistMatrixSyntaxError);
    PyModule_AddObject(m, "SyntaxError", CxgDistMatrixSyntaxError);

    // IOError.
    CxgDistMatrixIOError = PyErr_NewException("C_DistMatrix.IOError",
					      CxgDistMatrixException,
					      NULL);
    Py_INCREF(CxgDistMatrixIOError);
    PyModule_AddObject(m, "IOError", CxgDistMatrixIOError);
}
