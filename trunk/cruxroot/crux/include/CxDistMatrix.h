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
} CxtDistMatrixInputType;

typedef enum
{
    CxtDistMatrixTokenNone,
    CxtDistMatrixTokenInt,
    CxtDistMatrixTokenDec,
    CxtDistMatrixTokenLabel
} CxtDistMatrixTokenType;

typedef struct
{
    PyObject_HEAD
    char *buf;
    long bufLen;
    long tokenLen;
    long line;
    long column;

    CxtDistMatrixInputType inputType;
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

    PyObject *map;
    long ntaxa;
    double *matrix;
} CxtDistMatrixObject;

extern PyObject *CxgDistMatrixException;
extern PyObject *CxgDistMatrixValueError;
extern PyObject *CxgDistMatrixTypeError;
extern PyObject *CxgDistMatrixSyntaxError;

void
CxDistMatrixInit(void);

PyObject *
CxDistMatrixParse(CxtDistMatrixObject *self, PyObject *args);
PyObject *
CxDistMatrixNtaxaGet(CxtDistMatrixObject *self);
PyObject *
CxDistMatrixTaxonMapGet(CxtDistMatrixObject *self);
PyObject *
CxDistMatrixDistanceGet(CxtDistMatrixObject *self, PyObject *args);
PyObject *
CxDistMatrixDistanceSet(CxtDistMatrixObject *self, PyObject *args);
