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

typedef struct
{
    PyObject_HEAD
    char *buf;
    long bufLen;
    long tokenLen;
    long offset;

    char c;
    bool lookaheadValid;
    char lookaheadC;

    bool fileInput;
    union
    {
	struct
	{
	    char *string;
	    long offset;
	} s;
	struct
	{
	    FILE *file;
	} f;
    } i;
} CxtNewickParserObject;

extern PyObject *CxgNewickParserException;
extern PyObject *CxgNewickParserValueError;
extern PyObject *CxgNewickParserTypeError;
extern PyObject *CxgNewickParserSyntaxError;

void
CxNewickParserInit(void);

PyObject *
CxNewickParserParse(CxtNewickParserObject *self, PyObject *args);
PyObject *
CxNewickParserToken(CxtNewickParserObject *self);
PyObject *
CxNewickParserOffset(CxtNewickParserObject *self);
