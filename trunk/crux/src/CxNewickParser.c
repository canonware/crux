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
//
// Recursive descent parser that implements Newick Tree Format parsing.  This
// implementation expands from Gary Olsen's interpretation of the Newick Tree
// Format Standard, available at:
//
//   http://evolution.genetics.washington.edu/phylip/newick_doc.html
//
// Joe Felsenstein's description of the Newick Tree Format is at:
//
//   http://evolution.genetics.washington.edu/phylip/newicktree.html
//
// The two descriptions directly conflict in their definitions of the top level
// production.  Felsenstein claims that the following is a valid tree:
//
//   A;
//
// However, Olsen claims that the descendant_list component of the tree
// production is mandatory, and that the descendant_list production must not be
// empty.  This implementation follows Felsenstein's interpretation, since it is
// a superset of Olsen's interpretation.
//
// Similarly, Felsenstein claims that zero length labels are allowable, whereas
// Olsen claims that a label is a "string of printing characters".  This
// implementation allows "zero length strings of printing characters".  This has
// the implication that the label productions always accept.
//
//==============================================================================
//
// Notation:
//
//   Production: production_name ::= value
//
//   Production reference: <production>
//
//   Optional (0-1 occurrences): [optional]
//
//   Alternation: <this> | <that>
//
//   Regular expression: /regex/
//
//   Special character: "[", "]", etc.
//
//   Whitespace between components of a production alternative indicates
//   optional whitespace (<ws>), which consists of whitespace and/or comments.
//
// Notes:
//
//   "_" in <unquoted_label> is converted to " ".
//
// tree ::= [ <descendant_list> ][ <root_label> ][ : <branch_length> ] ;
//
// descendant_list ::= ( <subtree> , <subtree> [ , <subtree> ] )
//
// subtree ::= <descendant_list> [ <internal_label> ][ : <branch_length> ]
//           | <leaf_label> [ : <branch_length> ]
//
// root_label ::= <label>
//
// internal_label ::= <label>
//
// leaf_label ::= <label>
//
// label ::= <unquoted_label>
//         | <quoted_label>
//         | <e>
//
// unquoted_label ::= <ulabel_char>[<unquoted_label>]
//
// ulabel_char ::= /[^ ()[\]':;,]/
//
// quoted_label ::= '<qlabel_char>[<quoted_label>]'
//
// qlabel_char ::= /[^']|('')/
//
// branch_length ::= <int>[.<int>][<exp>]
//                 | /[-+]/<int>[.<int>][<exp>]
//
// exp ::= /[eE]/[/[-+]/]<int>
//
// int ::= <digit>[<int>]
//
// digit ::= /[0-9]/
//
// ws ::= [<comment>][<whitespace>][<ws>]
//      | <e>
//
// whitespace ::= / \t\r\n/
//
// comment ::= \[ [<comment_chars>] [<comment>] [<comment_chars>] \]
//
// comment_chars ::= /[^[]/
//                 | <e>
//
// e ::= <epsilon (empty production)>
//
//==============================================================================

#include "../include/_cruxmodule.h"

static PyTypeObject CxtNewickParser;

static PyObject *
CxpNewickParserNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtNewickParserObject *self;

    self = (CxtNewickParserObject *) type->tp_alloc(type, 0);
    if (self == NULL)
    {
	rVal = NULL;
	goto RETURN;
    }

    self->buf = NULL;
    self->bufLen = 0;
    self->tokenLen = 0;

    rVal = (PyObject *) self;
    RETURN:
    return rVal;
}

static int
CxpNewickParserTraverse(CxtNewickParserObject *self, visitproc visit, void *arg)
{
    return 0;
}

static int
CxpNewickParserClear(CxtNewickParserObject *self)
{
    return 0;
}

static void
CxpNewickParserDelete(CxtNewickParserObject *self)
{
    CxpNewickParserClear(self);
    if (self->buf != NULL)
    {
	free(self->buf);
    }
    self->ob_type->tp_free((PyObject*) self);
}

CxmpInline void
CxpNewickParserAppendC(CxtNewickParserObject *self)
{
    // Expand the token buffer, if necessary.
    if (self->tokenLen + 1 > self->bufLen)
    {
	if (self->buf == NULL)
	{
#define CxmpNewickParserBufLenStart 1024
	    self->buf = (char *) CxmMalloc(CxmpNewickParserBufLenStart);
	    self->bufLen = CxmpNewickParserBufLenStart;
#undef CxmpNewickParserBufLenStart
	}
	else
	{
	    self->buf = (char *) CxmRealloc(self->buf, self->bufLen * 2);
	    self->bufLen *= 2;
	}
    }

    // Append the character.
    self->buf[self->tokenLen] = self->c;
    self->tokenLen++;
}

CxmpInline bool
CxpNewickParserGetC(CxtNewickParserObject *self)
{
    bool rVal;

    if (self->lookaheadValid)
    {
	self->c = self->lookaheadC;
	self->lookaheadValid = false;
    }
    else
    {
	// Update offset.
	self->offset++;

	if (self->fileInput)
	{
	    if (fread(&self->c, 1, 1, self->i.f.file) == 0)
	    {
		CxError(CxgNewickParserSyntaxError,
			"At offset %d (token '%.*s'):"
			" End of input reached",
			self->offset - 1, self->tokenLen, self->buf);
		rVal = true;
		goto RETURN;
	    }
	}
	else
	{
	    self->c = self->i.s.string[self->i.s.offset];
	    self->i.s.offset++;
	    if (self->c == '\0')
	    {
		CxError(CxgNewickParserSyntaxError,
			"At offset %d (token '%.*s'):"
			" End of input reached",
			self->offset - 1, self->tokenLen, self->buf);
		rVal = true;
		goto RETURN;
	    }
	}
    }

    // Append character to token string.
    CxpNewickParserAppendC(self);

    rVal = false;
    RETURN:
    return rVal;
}

CxmpInline void
CxpNewickParserUngetC(CxtNewickParserObject *self)
{
    CxmAssert(self->tokenLen > 0);
    CxmAssert(self->lookaheadValid == false);

    self->tokenLen--;
    self->lookaheadC = self->c;
    self->lookaheadValid = true;
}

CxmpInline bool
CxpNewickParserTokenAccept(CxtNewickParserObject *self,
			   char *aMethodName)
{
    bool rVal;
    PyObject *result;

    result = PyEval_CallMethod((PyObject *) self, aMethodName, "()");
    if (result == NULL)
    {
	rVal = true;
	goto RETURN;
    }
    Py_DECREF(result);

    // Clear token.
    self->tokenLen = 0;

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpNewickParserProdWhitespace(CxtNewickParserObject *self, bool *rAccepted)
{
    bool rVal;

    if (CxpNewickParserGetC(self))
    {
	rVal = true;
	goto RETURN;
    }

    switch (self->c)
    {
	case ' ': case '\r': case '\n': case '\t':
	{
	    while (true)
	    {
		if (CxpNewickParserGetC(self))
		{
		    rVal = true;
		    goto RETURN;
		}

		switch (self->c)
		{
		    case ' ': case '\r': case '\n': case '\t':
		    {
			break;
		    }
		    default:
		    {
			CxpNewickParserUngetC(self);
			goto OUT;
		    }
		}
	    }
	    OUT:

	    *rAccepted = true;
	    break;
	}
	default:
	{
	    CxpNewickParserUngetC(self);

	    *rAccepted = false;
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpNewickParserProdComment(CxtNewickParserObject *self, bool *rAccepted)
{
    bool rVal, accepted;

    if (CxpNewickParserGetC(self))
    {
	rVal = true;
	goto RETURN;
    }

    if (self->c == '[')
    {
	while (true)
	{
	    if (CxpNewickParserGetC(self))
	    {
		rVal = true;
		goto RETURN;
	    }

	    if (self->c == ']')
	    {
		break;
	    }
	    else
	    {
		CxpNewickParserUngetC(self);

		if (CxpNewickParserProdComment(self, &accepted))
		{
		    rVal = true;
		    goto RETURN;
		}

		if (accepted == false)
		{
		    if (CxpNewickParserGetC(self))
		    {
			rVal = true;
			goto RETURN;
		    }
		}
	    }
	}

	*rAccepted = true;
    }
    else
    {
	CxpNewickParserUngetC(self);

	*rAccepted = false;
    }

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpNewickParserProdWs(CxtNewickParserObject *self)
{
    bool rVal, accepted, again;

    for (again = true; again;)
    {
	again = false;

	// Match a comment.
	if (CxpNewickParserProdComment(self, &accepted))
	{
	    rVal = true;
	    goto RETURN;
	}

	if (accepted)
	{
	    if (CxpNewickParserTokenAccept(self, "commentAccept"))
	    {
		rVal = true;
		goto RETURN;
	    }

	    again = true;
	}

	// Match whitespace.
	if (CxpNewickParserProdWhitespace(self, &accepted))
	{
	    rVal = true;
	    goto RETURN;
	}

	if (accepted)
	{
	    if (CxpNewickParserTokenAccept(self, "whitespaceAccept"))
	    {
		rVal = true;
		goto RETURN;
	    }

	    again = true;
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpNewickParserProdBranchLength(CxtNewickParserObject *self)
{
    bool rVal;
    enum
    {
	CxNewickParserLengthStateIntSign,
	CxNewickParserLengthStateInt,
	CxNewickParserLengthStateIntCont,
	CxNewickParserLengthStateDec,
	CxNewickParserLengthStateDecCont,
	CxNewickParserLengthStateExpSign,
	CxNewickParserLengthStateExp,
	CxNewickParserLengthStateExpCont
    } state;

    state = CxNewickParserLengthStateIntSign;
    while (true)
    {
	if (CxpNewickParserGetC(self))
	{
	    rVal = true;
	    goto RETURN;
	}

	switch (state)
	{
	    case CxNewickParserLengthStateIntSign:
	    {
		switch (self->c)
		{
		    case '+': case '-':
		    {
			state = CxNewickParserLengthStateInt;
			break;
		    }
		    case '0': case '1': case '2': case '3': case '4':
		    case '5': case '6': case '7': case '8': case '9':
		    {
			state = CxNewickParserLengthStateIntCont;
			break;
		    }
		    default:
		    {
			CxError(CxgNewickParserSyntaxError,
				"At offset %d (token '%.*s', char '%c'):"
				" Expected digit",
				self->offset - 1, self->tokenLen, self->buf,
				self->c);
			rVal = true;
			goto RETURN;
		    }
		}
		break;
	    }
	    case CxNewickParserLengthStateInt:
	    {
		switch (self->c)
		{
		    case '0': case '1': case '2': case '3': case '4':
		    case '5': case '6': case '7': case '8': case '9':
		    {
			state = CxNewickParserLengthStateIntCont;
			break;
		    }
		    default:
		    {
			CxError(CxgNewickParserSyntaxError,
				"At offset %d (token '%.*s', char '%c'):"
				" Expected digit",
				self->offset - 1, self->tokenLen, self->buf,
				self->c);
			rVal = true;
			goto RETURN;
		    }
		}
		break;
	    }
	    case CxNewickParserLengthStateIntCont:
	    {
		switch (self->c)
		{
		    case '0': case '1': case '2': case '3': case '4':
		    case '5': case '6': case '7': case '8': case '9':
		    {
			break;
		    }
		    case '.':
		    {
			state = CxNewickParserLengthStateDec;
			break;
		    }
		    case 'e': case 'E':
		    {
			state = CxNewickParserLengthStateExpSign;
			break;
		    }
		    default:
		    {
			CxpNewickParserUngetC(self);

			if (CxpNewickParserTokenAccept(self, "lengthAccept"))
			{
			    rVal = true;
			    goto RETURN;
			}

			goto OUT;
		    }
		}
		break;
	    }
	    case CxNewickParserLengthStateDec:
	    {
		switch (self->c)
		{
		    case '0': case '1': case '2': case '3': case '4':
		    case '5': case '6': case '7': case '8': case '9':
		    {
			state = CxNewickParserLengthStateDecCont;
			break;
		    }
		    default:
		    {
			CxError(CxgNewickParserSyntaxError,
				"At offset %d (token '%.*s', char '%c'):"
				" Expected digit",
				self->offset - 1, self->tokenLen, self->buf,
				self->c);
			rVal = true;
			goto RETURN;
		    }
		}
		break;
	    }
	    case CxNewickParserLengthStateDecCont:
	    {
		switch (self->c)
		{
		    case '0': case '1': case '2': case '3': case '4':
		    case '5': case '6': case '7': case '8': case '9':
		    {
			break;
		    }
		    case 'e': case 'E':
		    {
			state = CxNewickParserLengthStateExpSign;
			break;
		    }
		    default:
		    {
			CxpNewickParserUngetC(self);

			if (CxpNewickParserTokenAccept(self, "lengthAccept"))
			{
			    rVal = true;
			    goto RETURN;
			}

			goto OUT;
		    }
		}
		break;
	    }
	    case CxNewickParserLengthStateExpSign:
	    {
		switch (self->c)
		{
		    case '+': case '-':
		    {
			state = CxNewickParserLengthStateExp;
			break;
		    }
		    case '0': case '1': case '2': case '3': case '4':
		    case '5': case '6': case '7': case '8': case '9':
		    {
			state = CxNewickParserLengthStateExpCont;
			break;
		    }
		    default:
		    {
			CxError(CxgNewickParserSyntaxError,
				"At offset %d (token '%.*s', char '%c'):"
				" Expected digit",
				self->offset - 1, self->tokenLen, self->buf,
				self->c);
			rVal = true;
			goto RETURN;
		    }
		}
		break;
	    }
	    case CxNewickParserLengthStateExp:
	    {
		switch (self->c)
		{
		    case '0': case '1': case '2': case '3': case '4':
		    case '5': case '6': case '7': case '8': case '9':
		    {
			state = CxNewickParserLengthStateIntCont;
			break;
		    }
		    default:
		    {
			CxError(CxgNewickParserSyntaxError,
				"At offset %d (token '%.*s', char '%c'):"
				" Expected digit",
				self->offset - 1, self->tokenLen, self->buf,
				self->c);
			rVal = true;
			goto RETURN;
		    }
		}
		break;
	    }
	    case CxNewickParserLengthStateExpCont:
	    {
		switch (self->c)
		{
		    case '0': case '1': case '2': case '3': case '4':
		    case '5': case '6': case '7': case '8': case '9':
		    {
			break;
		    }
		    default:
		    {
			CxpNewickParserUngetC(self);

			if (CxpNewickParserTokenAccept(self, "lengthAccept"))
			{
			    rVal = true;
			    goto RETURN;
			}

			goto OUT;
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
    OUT:

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpNewickParserProdQuotedLabelChar(CxtNewickParserObject *self, bool *rAccepted)
{
    bool rVal;

    if (CxpNewickParserGetC(self))
    {
	rVal = true;
	goto RETURN;
    }

    switch (self->c)
    {
	case '\'':
	{
	    if (CxpNewickParserGetC(self))
	    {
		rVal = true;
		goto RETURN;
	    }

	    if (self->c == '\'')
	    {
		*rAccepted = true;
	    }
	    else
	    {
		CxpNewickParserUngetC(self);

		*rAccepted = false;
	    }

	    break;
	}
	default:
	{
	    *rAccepted = true;
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpNewickParserProdQuotedLabel(CxtNewickParserObject *self, char *aMethodName,
			       bool *rAccepted)
{
    bool rVal, accepted;

    if (CxpNewickParserGetC(self))
    {
	rVal = true;
	goto RETURN;
    }

    if (self->c == '\'')
    {
	long i;

	while (true)
	{
	    if (CxpNewickParserProdQuotedLabelChar(self, &accepted))
	    {
		rVal = true;
		goto RETURN;
	    }

	    if (accepted == false)
	    {
		break;
	    }
	}

	// Remove quotes before accepting the token.
	if (self->tokenLen > 2)
	{
	    memmove(self->buf, &self->buf[1], self->tokenLen - 2);
	}
	self->tokenLen -= 2;

	// Collapse '' to ' before accepting the token.
	for (i = 0; i < self->tokenLen; i++)
	{
	    if (self->buf[i] == '\'')
	    {
		if (i < self->tokenLen - 1)
		{
		    memmove(&self->buf[i], &self->buf[i + 1],
			    self->tokenLen - 1 - i);
		}
		self->tokenLen--;
	    }
	}

	if (CxpNewickParserTokenAccept(self, aMethodName))
	{
	    rVal = true;
	    goto RETURN;
	}

	*rAccepted = true;
    }
    else
    {
	CxpNewickParserUngetC(self);

	*rAccepted = false;
    }

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpNewickParserProdUnquotedLabelChar(CxtNewickParserObject *self,
				     bool *rAccepted)
{
    bool rVal;

    if (CxpNewickParserGetC(self))
    {
	rVal = true;
	goto RETURN;
    }

    switch (self->c)
    {
	case ' ': case '(': case ')': case '[': case ']': case '\'': case ':':
	case ';': case ',':
	{
	    CxpNewickParserUngetC(self);

	    *rAccepted = false;
	    break;
	}
	default:
	{
	    *rAccepted = true;
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpNewickParserProdUnquotedLabel(CxtNewickParserObject *self, char *aMethodName,
				 bool *rAccepted)
{
    bool rVal, accepted;

    if (CxpNewickParserProdUnquotedLabelChar(self, &accepted))
    {
	rVal = true;
	goto RETURN;
    }

    if (accepted)
    {
	long i;

	// Consume all ulabel_char characters.
	while (true)
	{
	    if (CxpNewickParserProdUnquotedLabelChar(self, &accepted))
	    {
		rVal = true;
		goto RETURN;
	    }

	    if (accepted == false)
	    {
		break;
	    }
	}

	// Convert '_' to ' ' before accepting the token.
	for (i = 0; i < self->tokenLen; i++)
	{
	    if (self->buf[i] == '_')
	    {
		self->buf[i] = ' ';
	    }
	}

	if (CxpNewickParserTokenAccept(self, aMethodName))
	{
	    rVal = true;
	    goto RETURN;
	}

	*rAccepted = true;
    }
    else
    {
	*rAccepted = false;
    }

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpNewickParserProdLabel(CxtNewickParserObject *self, char *aMethodName)
{
    bool rVal, accepted;

    if (CxpNewickParserProdQuotedLabel(self, aMethodName, &accepted))
    {
	rVal = true;
	goto RETURN;
    }

    if (accepted == false)
    {
	if (CxpNewickParserProdUnquotedLabel(self, aMethodName, &accepted))
	{
	    rVal = true;
	    goto RETURN;
	}

	if (accepted == false)
	{
	    // Accept a zero length label.
	    if (CxpNewickParserTokenAccept(self, aMethodName))
	    {
		rVal = true;
		goto RETURN;
	    }
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpNewickParserProdDescendantList(CxtNewickParserObject *self, bool *rAccepted);

static bool
CxpNewickParserProdSubtree(CxtNewickParserObject *self)
{
    bool rVal, accepted;

    if (CxpNewickParserProdWs(self)
	|| CxpNewickParserProdDescendantList(self, &accepted)
	|| CxpNewickParserProdWs(self))
    {
	rVal = true;
	goto RETURN;
    }
    if (accepted)
    {
	if (CxpNewickParserProdLabel(self, "internalLabelAccept"))
	{
	    rVal = true;
	    goto RETURN;
	}
    }
    else
    {
	if (CxpNewickParserProdLabel(self, "leafLabelAccept"))
	{
	    rVal = true;
	    goto RETURN;
	}
    }

    if (CxpNewickParserProdWs(self)
	|| (CxpNewickParserGetC(self)))
    {
	rVal = true;
	goto RETURN;
    }

    if (self->c == ':')
    {
	if (CxpNewickParserTokenAccept(self, "colonAccept")
	    || CxpNewickParserProdWs(self)
	    || CxpNewickParserProdBranchLength(self))
	{
	    rVal = true;
	    goto RETURN;
	}
    }
    else
    {
	CxpNewickParserUngetC(self);
    }

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpNewickParserProdDescendantList(CxtNewickParserObject *self, bool *rAccepted)
{
    bool rVal;

    if (CxpNewickParserGetC(self))
    {
	rVal = true;
	goto RETURN;
    }

    if (self->c == '(')
    {
	if (CxpNewickParserTokenAccept(self, "openParenAccept")
	    || CxpNewickParserProdWs(self)
	    || CxpNewickParserProdSubtree(self)
	    || CxpNewickParserProdWs(self))
	{
	    rVal = true;
	    goto RETURN;
	}

	if (CxpNewickParserGetC(self))
	{
	    rVal = true;
	    goto RETURN;
	}

	if (self->c != ',')
	{
	    CxError(CxgNewickParserSyntaxError,
		    "At offset %d (token '%.*s', char '%c'):"
		    " ',' expected",
		    self->offset - 1, self->tokenLen, self->buf, self->c);
	    rVal = true;
	    goto RETURN;
	}

	if (CxpNewickParserTokenAccept(self, "commaAccept")
	    || CxpNewickParserProdWs(self)
	    || CxpNewickParserProdSubtree(self)
	    || CxpNewickParserProdWs(self))
	{
	    rVal = true;
	    goto RETURN;
	}

	while (true)
	{
	    if (CxpNewickParserGetC(self))
	    {
		rVal = true;
		goto RETURN;
	    }

	    if (self->c == ',')
	    {
		if (CxpNewickParserTokenAccept(self, "commaAccept")
		    || CxpNewickParserProdWs(self)
		    || CxpNewickParserProdSubtree(self)
		    || CxpNewickParserProdWs(self))
		{
		    rVal = true;
		    goto RETURN;
		}
	    }
	    else
	    {
		CxpNewickParserUngetC(self);
		break;
	    }
	}

	if (CxpNewickParserGetC(self))
	{
	    rVal = true;
	    goto RETURN;
	}

	if (self->c != ')')
	{
	    CxError(CxgNewickParserSyntaxError,
		    "At offset %d (token '%.*s', char '%c'):"
		    " ',' or ')' expected",
		    self->offset - 1, self->tokenLen, self->buf, self->c);
	    rVal = true;
	    goto RETURN;
	}

	if (CxpNewickParserTokenAccept(self, "closeParenAccept"))
	{
	    rVal = true;
	    goto RETURN;
	}

	if (rAccepted != NULL)
	{
	    *rAccepted = true;
	}
    }
    else
    {
	CxpNewickParserUngetC(self);

	if (rAccepted != NULL)
	{
	    *rAccepted = false;
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpNewickParserProdTree(CxtNewickParserObject *self)
{
    bool rVal;

    if (CxpNewickParserProdWs(self)
	|| CxpNewickParserProdDescendantList(self, NULL)
	|| CxpNewickParserProdWs(self)
	|| CxpNewickParserProdLabel(self, "rootLabelAccept")
	|| CxpNewickParserProdWs(self)
	|| CxpNewickParserGetC(self))
    {
	rVal = true;
	goto RETURN;
    }

    if (self->c == ':')
    {
	if (CxpNewickParserTokenAccept(self, "colonAccept")
	    || CxpNewickParserProdWs(self)
	    || CxpNewickParserProdBranchLength(self)
	    || CxpNewickParserProdWs(self))
	{
	    rVal = true;
	    goto RETURN;
	}
    }
    else
    {
	CxpNewickParserUngetC(self);
    }

    if (CxpNewickParserGetC(self))
    {
	rVal = true;
	goto RETURN;
    }

    if (self->c == ';')
    {
	if (CxpNewickParserTokenAccept(self, "semicolonAccept"))
	{
	    rVal = true;
	    goto RETURN;
	}
    }
    else
    {
	CxError(CxgNewickParserSyntaxError,
		"At offset %d (token '%.*s', char '%c'): ';' expected",
		self->offset - 1, self->tokenLen, self->buf, self->c);
	rVal = true;
	goto RETURN;
    }

    rVal = false;
    RETURN:
    return rVal;
}

PyObject *
CxNewickParserParse(CxtNewickParserObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    PyObject *input;

    if (PyArg_ParseTuple(args, "O", &input) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    // Determine input type.
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
	CxError(CxgNewickParserTypeError, "input: file or string expected");
	rVal = NULL;
	goto RETURN;
    }

    // Parse.
    CxmXepBegin();
    CxmXepTry
    {
	enum
	{
	    CxpStateStart,
	    CxpStateLabel,
	    CxpStateComment,
	    CxpStateChars
	} state;

	state = CxpStateStart;
	self->offset = 0;
	self->tokenLen = 0;
	self->c = '\0';
	self->lookaheadValid = false;
	self->lookaheadC = '\0';

	if (CxpNewickParserProdTree(self))
	{
	    goto ERROR;
	}

	Py_INCREF(Py_None);
	rVal = Py_None;
	break;

	ERROR:
	rVal = NULL;
    }
    CxmXepCatch(CxmXepOOM)
    {
	CxmXepHandled();
	PyErr_NoMemory();
	rVal = NULL;
    }
    CxmXepEnd();

    RETURN:
    return rVal;
}

PyObject *
CxNewickParserToken(CxtNewickParserObject *self)
{
    PyObject *rVal;

    if (self->buf != NULL)
    {
	rVal = PyString_FromStringAndSize(self->buf, self->tokenLen);
    }
    else
    {
	rVal = PyString_FromString("");
    }

    return rVal;
}

PyObject *
CxNewickParserOffset(CxtNewickParserObject *self)
{
    return Py_BuildValue("i", self->offset);
}

static PyMethodDef CxpNewickParserMethods[] =
{
    {
	"_parse",
	(PyCFunction) CxNewickParserParse,
	METH_VARARGS,
	"_parse"
    },
    {
	"token",
	(PyCFunction) CxNewickParserToken,
	METH_NOARGS,
	"token"
    },
    {
	"offset",
	(PyCFunction) CxNewickParserOffset,
	METH_NOARGS,
	"offset"
    },
    {NULL, NULL}
};

static PyTypeObject CxtNewickParser =
{
    PyObject_HEAD_INIT(NULL)
    0,			// int ob_size
    "C_NewickParser.C_NewickParser",	// char *tp_name
    sizeof(CxtNewickParserObject),	// int tp_basicsize
    0,			// int tp_itemsize
    (destructor) CxpNewickParserDelete,	// destructor tp_dealloc
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
    "NewickParser(): Create the C portion of a tree.",	// char *tp_doc
    (traverseproc) CxpNewickParserTraverse,	// traverseproc tp_traverse
    (inquiry) CxpNewickParserClear,	// inquiry tp_clear
    0,			// richcmpfunc tp_richcompare
    0,			// long tp_weaklistoffset
    0,			// getiterfunc tp_iter
    0,			// iternextfunc tp_iternext
    CxpNewickParserMethods,	// struct PyMethodDef *tp_methods
    0,			// struct PyMemberDef *tp_members
    0,			// struct PyGetSetDef *tp_getset
    0,			// struct _typeobject *tp_base
    0,			// PyObject *tp_dict
    0,			// descrgetfunc tp_descr_get
    0,			// descrsetfunc tp_descr_set
    0,			// long tp_dictoffset
    0,			// initproc tp_init
    0,			// allocfunc tp_alloc
    CxpNewickParserNew,		// newfunc tp_new
    _PyObject_Del,	// freefunc tp_free
    0			// inquiry tp_is_gc
};

static PyMethodDef CxpNewickParserFuncs[] =
{
    {NULL}
};

PyObject *CxgNewickParserException;
PyObject *CxgNewickParserValueError;
PyObject *CxgNewickParserTypeError;
PyObject *CxgNewickParserSyntaxError;

void
CxNewickParserInit(void)
{
    PyObject *m;

    // Create new type.
    if (PyType_Ready(&CxtNewickParser) < 0)
    {
	return;
    }
    m = Py_InitModule3("C_NewickParser", CxpNewickParserFuncs,
		       "NewickParser extensions");
    Py_INCREF(&CxtNewickParser);
    PyModule_AddObject(m, "C_NewickParser", (PyObject *) &CxtNewickParser);

    // Create exception objects.
    // Exception.
    CxgNewickParserException = PyErr_NewException("C_NewickParser.Exception",
						  CxgException,
						  NULL);
    Py_INCREF(CxgNewickParserException);
    PyModule_AddObject(m, "Exception", CxgNewickParserException);

    // ValueError.
    CxgNewickParserValueError = PyErr_NewException("C_NewickParser.ValueError",
						   CxgNewickParserException,
						   NULL);
    Py_INCREF(CxgNewickParserValueError);
    PyModule_AddObject(m, "ValueError", CxgNewickParserValueError);

    // TypeError.
    CxgNewickParserTypeError = PyErr_NewException("C_NewickParser.TypeError",
						  CxgNewickParserException,
						  NULL);
    Py_INCREF(CxgNewickParserTypeError);
    PyModule_AddObject(m, "TypeError", CxgNewickParserTypeError);

    // SyntaxError.
    CxgNewickParserSyntaxError
	= PyErr_NewException("C_NewickParser.SyntaxError",
			     CxgNewickParserException,
			     NULL);
    Py_INCREF(CxgNewickParserSyntaxError);
    PyModule_AddObject(m, "SyntaxError", CxgNewickParserSyntaxError);
}
