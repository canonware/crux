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

typedef enum
{
    CxDistMatrixInputFile,
    CxDistMatrixInputString,
    CxDistMatrixInputTaxonMap
} CxDistMatrixInputType;

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

    int ntaxa;
    PyObject *map;
    double *matrix;
} CxtDistMatrixObject;

extern PyObject *CxgDistMatrixException;
extern PyObject *CxgDistMatrixValueError;
extern PyObject *CxgDistMatrixTypeError;
extern PyObject *CxgDistMatrixSyntaxError;

void
CxDistMatrixInit(void);

CxtTreeObject *
CxDistMatrixNew(void);
PyObject *
CxDistMatrixParse(CxtDistMatrixObject *self, PyObject *args);
PyObject *
CxDistMatrixToken(CxtDistMatrixObject *self);
PyObject *
CxDistMatrixLine(CxtDistMatrixObject *self);
