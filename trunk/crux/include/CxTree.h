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

typedef struct CxsTreeObject CxtTreeObject;
typedef struct CxsNodeObject CxtNodeObject;
typedef struct CxsEdgeObject CxtEdgeObject;
typedef struct CxsRingObject CxtRingObject;

struct CxsTreeObject
{
    PyObject_HEAD
    CxtTr *tr;
};

struct CxsNodeObject
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtTrNode node;
};

struct CxsEdgeObject
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtRingObject *ringA;
    CxtRingObject *ringB;
    CxtTrEdge edge;
    bool GcDetached;
};

struct CxsRingObject
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtEdgeObject *edge;
    CxtTrRing ring;
    bool GcDetached;
};

extern PyObject *CxgTreeException;
extern PyObject *CxgTreeValueError;
extern PyObject *CxgTreeTypeError;

extern PyObject *CxgEdgeException;
extern PyObject *CxgEdgeValueError;

extern PyObject *CxgRingException;
extern PyObject *CxgRingValueError;

void
CxTreeInit(void);
void
CxNodeInit(void);
void
CxEdgeInit(void);
void
CxRingInit(void);

/* Tree. */
CxtTreeObject *
CxTreeNew(void);
PyObject *
CxTreeNtaxaGet(CxtTreeObject *self);
PyObject *
CxTreeNedgesCget(CxtTreeObject *self);
PyObject *
CxTreeBaseGet(CxtTreeObject *self);
void
CxTreeBaseSet(CxtTreeObject *self, CxtNodeObject *aNode);
PyObject *
CxTreeBaseSetPargs(CxtTreeObject *self, PyObject *args);

/* Node. */
CxtNodeObject *
CxNodeNew(CxtTreeObject *a_tree);
PyObject *
CxNodeTree(CxtNodeObject *self);
PyObject *
CxNodeTaxonNumGet(CxtNodeObject *self);
void
CxNodeTaxonNumSet(CxtNodeObject *self, uint32_t a_taxon_num);
PyObject *
CxNodeTaxonNumSetPargs(CxtNodeObject *self, PyObject *args);
PyObject *
CxNodeRing(CxtNodeObject *self);
PyObject *
CxNodeDegree(CxtNodeObject *self);

/* Edge. */
CxtEdgeObject *
CxEdgeNew(CxtTreeObject *a_tree);
PyObject *
CxEdgeTree(CxtEdgeObject *self);
PyObject *
CxEdgeRingsGet(CxtEdgeObject *self);
PyObject *
CxEdgeLengthGet(CxtEdgeObject *self);
void
CxEdgeLengthSet(CxtEdgeObject *self, double aLength);
PyObject *
CxEdgeLengthSetPargs(CxtEdgeObject *self, PyObject *args);
void
CxEdgeAttach(CxtEdgeObject *self, CxtNodeObject *aNodeA,
	     CxtNodeObject *aNodeB);
PyObject *
CxEdgeAttachPargs(CxtEdgeObject *self, PyObject *args);
PyObject *
CxEdgeDetach(CxtEdgeObject *self);

/* Ring. */
CxtRingObject *
CxRingNew(CxtEdgeObject *aEdge, uint32_t aEnd);
PyObject *
CxRingTree(CxtRingObject *self);
PyObject *
CxRingNode(CxtRingObject *self);
PyObject *
CxRingEdge(CxtRingObject *self);
PyObject *
CxRingOther(CxtRingObject *self);
PyObject *
CxRingNext(CxtRingObject *self);
PyObject *
CxRingPrev(CxtRingObject *self);
