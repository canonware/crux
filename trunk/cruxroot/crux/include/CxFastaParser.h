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
    long bufLen;
    long tokenLen;
    long line;
    long column;

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
} CxtFastaParserObject;

extern PyObject *CxgFastaParserException;
extern PyObject *CxgFastaParserValueError;
extern PyObject *CxgFastaParserTypeError;
extern PyObject *CxgFastaParserSyntaxError;

void
CxFastaParserInit(void);

PyObject *
CxFastaParserParse(CxtFastaParserObject *self, PyObject *args);
PyObject *
CxFastaParserToken(CxtFastaParserObject *self);
PyObject *
CxFastaParserLine(CxtFastaParserObject *self);
