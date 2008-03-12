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
} CxtNexusParserObject;

typedef struct
{
    PyObject_HEAD
} CxtNexusParserBaseObject;

void
CxNexusParserInit(void);
void
CxNexusParserBaseInit(void);

// CxNexusParser.
PyObject *
CxNexusParserInputName(CxtNexusParserObject *self);
PyObject *
CxNexusParserTokenGet(CxtNexusParserObject *self);
PyObject *
CxNexusParserTokenUnget(CxtNexusParserObject *self);
PyObject *
CxNexusParserTextRange(CxtNexusParserObject *self, PyObject *args);
PyObject *
CxNexusParserHandlerSearch(CxtNexusParserObject *self, PyObject *args);
PyObject *
CxNexusParserRoot(CxtNexusParserObject *self);
PyObject *
CxNexusParserPrint(CxtNexusParserObject *self);

// CxNexusParserBase.
PyObject *
CxNexusParserBaseInput(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseText(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseOffset(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseLine(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseColumn(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseLengthGet(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseLengthSet(CxtNexusParserBaseObject *self, PyObject *args);
PyObject *
CxNexusParserBaseCText(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseFinish(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseFinished(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBasePrint(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseChildAppend(CxtNexusParserBaseObject *self, PyObject *args);
PyObject *
CxNexusParserBaseNext(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBasePrev(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseParent(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseIndex(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseNChildren(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseChild(CxtNexusParserBaseObject *self, PyObject *args);
PyObject *
CxNexusParserBaseLeft(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseRight(CxtNexusParserBaseObject *self);
PyObject *
CxNexusParserBaseAccessCallback(CxtNexusParserBaseObject *self);
