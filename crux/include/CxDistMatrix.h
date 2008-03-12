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
    bool symmetric;
} CxtDistMatrixObject;

extern PyObject *CxgDistMatrixException;
extern PyObject *CxgDistMatrixValueError;
extern PyObject *CxgDistMatrixTypeError;
extern PyObject *CxgDistMatrixSyntaxError;
extern PyObject *CxgDistMatrixIOError;

void
CxDistMatrixInit(void);

PyObject *
CxDistMatrixParse(CxtDistMatrixObject *self, PyObject *args);
PyObject *
CxDistMatrixNtaxaGet(CxtDistMatrixObject *self);
PyObject *
CxDistMatrixIsSymmetric(CxtDistMatrixObject *self);
PyObject *
CxDistMatrixTaxonMapGet(CxtDistMatrixObject *self);
double
CxDistMatrixDistanceGet(CxtDistMatrixObject *self, long x, long y);
PyObject *
CxDistMatrixDistanceGetPargs(CxtDistMatrixObject *self, PyObject *args);
void
CxDistMatrixDistanceSet(CxtDistMatrixObject *self, long x, long y, double dist);
PyObject *
CxDistMatrixDistanceSetPargs(CxtDistMatrixObject *self, PyObject *args);

void
CxDistMatrixUpperHandoff(CxtDistMatrixObject *self, double **rMatrix,
			 long *rNtaxa);
