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
    CxtTr *tr;
} CxtTreeObject;

typedef struct
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtTrNode node;
} CxtNodeObject;

typedef struct
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtTrEdge edge;
} CxtEdgeObject;

void
CxTreeInit(void);
void
CxNodeInit(void);
void
CxEdgeInit(void);

CxtTreeObject *
CxTreeNew(void);
PyObject *
CxTreeNtaxaGet(CxtTreeObject *self);
PyObject *
CxTreeNedgesCget(CxtTreeObject *self);
PyObject *
CxTreeBaseGet(CxtTreeObject *self);
void
CxTreeBaseSetCargs(CxtTreeObject *self, CxtNodeObject *aNode);
PyObject *
CxTreeBaseSet(CxtTreeObject *self, PyObject *args);

CxtNodeObject *
CxNodeNew(CxtTreeObject *a_tree);
PyObject *
CxNodeTree(CxtNodeObject *self);
PyObject *
CxNodeTaxonNumGet(CxtNodeObject *self);
void
CxNodeTaxonNumSetCargs(CxtNodeObject *self, uint32_t a_taxon_num);
PyObject *
CxNodeTaxonNumSet(CxtNodeObject *self, PyObject *args);
PyObject *
CxNodeEdge(CxtNodeObject *self);
PyObject *
CxNodeDegree(CxtNodeObject *self);

CxtEdgeObject *
CxEdgeNew(CxtTreeObject *a_tree);
PyObject *
CxEdgeTree(CxtEdgeObject *self);
PyObject *
CxEdgeNodeCargs(CxtEdgeObject *self, uint32_t aInd);
PyObject *
CxEdgeNode(CxtEdgeObject *self, PyObject *args);
void
CxEdgeNextCargs(CxtEdgeObject *self, uint32_t aInd, CxtEdgeObject **rEdge,
		uint32_t *rNextEnd);
PyObject *
CxEdgeNext(CxtEdgeObject *self, PyObject *args);
void
CxEdgePrevCargs(CxtEdgeObject *self, uint32_t aInd, CxtEdgeObject **rEdge,
		uint32_t *rPrevEnd);
PyObject *
CxEdgePrev(CxtEdgeObject *self, PyObject *args);
PyObject *
CxEdgeLengthGet(CxtEdgeObject *self);
void
CxEdgeLengthSetCargs(CxtEdgeObject *self, double aLength);
PyObject *
CxEdgeLengthSet(CxtEdgeObject *self, PyObject *args);
void
CxEdgeAttachCargs(CxtEdgeObject *self, CxtNodeObject *aNodeA,
		  CxtNodeObject *aNodeB);
PyObject *
CxEdgeAttach(CxtEdgeObject *self, PyObject *args);
PyObject *
CxEdgeDetach(CxtEdgeObject *self);
