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

typedef struct
{
    PyObject_HEAD
    char *buf;
    int bufLen;
    int tokenLen;
    int line;
    int column;
    int offset;

    bool fileInput;
    union
    {
	struct
	{
	    char *string;
	} s;
	struct
	{
	    FILE *file;
	} f;
    } i;
} CxtFastaParserObject;

extern PyObject *CxgFastaParserException;
extern PyObject *CxgFastaParserValueError;
extern PyObject *CxgFastaParserTypeError;
extern PyObject *CxgFastaParserSyntaxError;

void
CxFastaParserInit(void);

CxtTreeObject *
CxFastaParserNew(void);
PyObject *
CxFastaParserParse(CxtFastaParserObject *self, PyObject *args);
PyObject *
CxFastaParserToken(CxtFastaParserObject *self);
PyObject *
CxFastaParserLine(CxtFastaParserObject *self);
