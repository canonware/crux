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

typedef struct CxsTreeAux CxtTreeAux;
typedef struct CxsNodeAux CxtNodeAux;
typedef struct CxsEdgeAux CxtEdgeAux;
typedef struct CxsRingAux CxtRingAux;

typedef struct CxsTreeObject CxtTreeObject;
typedef struct CxsNodeObject CxtNodeObject;
typedef struct CxsEdgeObject CxtEdgeObject;
typedef struct CxsRingObject CxtRingObject;

extern PyTypeObject CxtTree;
extern PyTypeObject CxtNode;
extern PyTypeObject CxtEdge;
extern PyTypeObject CxtRing;

// The (void *) argument to these types of functions is the "data" pointer
// associated with each aux registration.
//
// The (unsigned) argument is the aux index 
typedef bool
CxtNodeAuxInit(CxtNodeObject *, unsigned);
typedef bool
CxtEdgeAuxInit(CxtEdgeObject *, unsigned);
typedef bool
CxtRingAuxInit(CxtRingObject *, unsigned);

typedef void
CxtTreeAuxCleanup(CxtTreeObject *, void *, unsigned);
typedef void
CxtNodeAuxCleanup(CxtNodeObject *, void *, unsigned);
typedef void
CxtEdgeAuxCleanup(CxtEdgeObject *, void *, unsigned);
typedef void
CxtRingAuxCleanup(CxtRingObject *, void *, unsigned);

struct CxsTreeAux
{
    // Used to associate types of aux data with particular indices.
    const char *key;

    // Opaque data associated with this registration.
    void *data;

    // Called during tree destruction.
    CxtTreeAuxCleanup *cleanupFinal;

    // Called during cleanup of tree aux pointers.
    CxtTreeAuxCleanup *cleanupTree;
};

struct CxsNodeAux
{
    // Used to associate types of aux data with particular indices.
    const char *key;

    // Opaque data associated with this registration.
    void *data;

    // Called during node construction.
    CxtNodeAuxInit *initNode;

    // Called during tree destruction.
    CxtTreeAuxCleanup *cleanupFinal;

    // Called during node destruction.
    CxtNodeAuxCleanup *cleanupNode;
};

struct CxsEdgeAux
{
    // Used to associate types of aux data with particular indices.
    const char *key;

    // Opaque data associated with this registration.
    void *data;

    // Called during node construction.
    CxtEdgeAuxInit *initEdge;

    // Called during tree destruction.
    CxtTreeAuxCleanup *cleanupFinal;

    // Called during edge destruction.
    CxtEdgeAuxCleanup *cleanupEdge;
};

struct CxsRingAux
{
    // Used to associate types of aux data with particular indices.
    const char *key;

    // Opaque data associated with this registration.
    void *data;

    // Called during node construction.
    CxtRingAuxInit *initRing;

    // Called during tree destruction.
    CxtTreeAuxCleanup *cleanupFinal;

    // Called during ring destruction.
    CxtRingAuxCleanup *cleanupRing;
};

struct CxsTreeObject
{
    PyObject_HEAD
    CxtTr *tr;
    uint64_t seq;
    bool GcCleared;

    // Aux control data structures for trees, nodes, edges, and rings.  Actual
    // aux pointers reside in the respective objects, but all control
    // information is centrally located here.
    CxtTreeAux *treeAux;
    unsigned nTreeAux;

    CxtNodeAux *nodeAux;
    unsigned nNodeAux;

    CxtEdgeAux *edgeAux;
    unsigned nEdgeAux;

    CxtRingAux *ringAux;
    unsigned nRingAux;

    // Aux data.
    void **aux;
    unsigned nAux;
};

struct CxsNodeObject
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtTrNode node;
    bool GcCleared;

    // Aux data.
    void **aux;
    unsigned nAux;
};

struct CxsEdgeObject
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtRingObject *ringA;
    CxtRingObject *ringB;
    CxtTrEdge edge;
    bool GcDetached:1;
    bool GcCleared:1;

    // Aux data.
    void **aux;
    unsigned nAux;
};

struct CxsRingObject
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtEdgeObject *edge;
    CxtTrRing ring;
    bool GcDetached:1;
    bool GcCleared:1;

    // Aux data.
    void **aux;
    unsigned nAux;
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

// Tree.
CxtTreeObject *
CxTreeNew(void);

unsigned
CxTreeNtaxaGet(CxtTreeObject *self);
PyObject *
CxTreeNtaxaGetPargs(CxtTreeObject *self);

PyObject *
CxTreeNedgesCget(CxtTreeObject *self);

CxtNodeObject *
CxTreeBaseGet(CxtTreeObject *self);
PyObject *
CxTreeBaseGetPargs(CxtTreeObject *self);

void
CxTreeBaseSet(CxtTreeObject *self, CxtNodeObject *aNode);
PyObject *
CxTreeBaseSetPargs(CxtTreeObject *self, PyObject *args);

typedef enum
{
    CxTreeIteratorStagePre,
    CxTreeIteratorStageIn,
    CxTreeIteratorStagePost
} CxtTreeIteratorStage;

typedef bool
CxtTreeIterateNodeCallback(CxtNodeObject *, CxtTreeIteratorStage, void *);
typedef bool
CxtTreeIterateEdgeCallback(CxtEdgeObject *, CxtTreeIteratorStage, void *);
typedef bool
CxtTreeIterateRingCallback(CxtRingObject *, CxtTreeIteratorStage, void *);

bool
CxTreeSubtreeIterate(CxtRingObject *aRing,
		     CxtTreeIterateNodeCallback *aNodeCallback,
		     CxtTreeIterateEdgeCallback *aEdgeCallback,
		     CxtTreeIterateRingCallback *aRingCallback,
		     void *aContext);

bool
CxTreeIterate(CxtTreeObject *aTree,
	      CxtTreeIterateNodeCallback *aNodeCallback,
	      CxtTreeIterateEdgeCallback *aEdgeCallback,
	      CxtTreeIterateRingCallback *aRingCallback,
	      void *aContext);

#if (!defined(CxmUseInlines))
uint64_t
CxTreeSeq(CxtTreeObject *self);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline uint64_t
CxTreeSeq(CxtTreeObject *self)
{
    return self->seq;
}
#endif

bool
CxTreeAuxRegister(CxtTreeObject *self, const char *aKey, void *aData,
		  CxtTreeAuxCleanup *aCleanupFinal,
		  CxtTreeAuxCleanup *aCleanupTree, unsigned *rInd);
bool
CxTreeAuxSearch(CxtTreeObject *self, const char *aKey, unsigned *rInd);

#if (!defined(CxmUseInlines))
void *
CxTreeAuxData(CxtTreeObject *self, unsigned aInd);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxTreeAuxData(CxtTreeObject *self, unsigned aInd)
{
    CxmAssert(aInd < self->nTreeAux);

    return self->treeAux[aInd].data;
}
#endif

bool
CxTreeNodeAuxRegister(CxtTreeObject *self, const char *aKey, void *aData,
		      CxtNodeAuxInit *aInitNode,
		      CxtTreeAuxCleanup *aCleanupFinal,
		      CxtNodeAuxCleanup *aCleanupNode, unsigned *rInd);
bool
CxTreeNodeAuxSearch(CxtTreeObject *self, const char *aKey, unsigned *rInd);

#if (!defined(CxmUseInlines))
void *
CxNodeAuxData(CxtTreeObject *self, unsigned aInd);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxNodeAuxData(CxtTreeObject *self, unsigned aInd)
{
    CxmAssert(aInd < self->nNodeAux);

    return self->nodeAux[aInd].data;
}
#endif

bool
CxTreeEdgeAuxRegister(CxtTreeObject *self, const char *aKey, void *aData,
		      CxtEdgeAuxInit *aInitEdge,
		      CxtTreeAuxCleanup *aCleanupFinal,
		      CxtEdgeAuxCleanup *aCleanupEdge, unsigned *rInd);
bool
CxTreeEdgeAuxSearch(CxtTreeObject *self, const char *aKey, unsigned *rInd);

#if (!defined(CxmUseInlines))
void *
CxEdgeAuxData(CxtTreeObject *self, unsigned aInd);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxEdgeAuxData(CxtTreeObject *self, unsigned aInd)
{
    CxmAssert(aInd < self->nEdgeAux);

    return self->edgeAux[aInd].data;
}
#endif

bool
CxTreeRingAuxRegister(CxtTreeObject *self, const char *aKey, void *aData,
		      CxtRingAuxInit *aInitRing,
		      CxtTreeAuxCleanup *aCleanupFinal,
		      CxtRingAuxCleanup *aCleanupRing, unsigned *rInd);
bool
CxTreeRingAuxSearch(CxtTreeObject *self, const char *aKey, unsigned *rInd);

#if (!defined(CxmUseInlines))
void *
CxRingAuxData(CxtTreeObject *self, unsigned aInd);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxRingAuxData(CxtTreeObject *self, unsigned aInd)
{
    CxmAssert(aInd < self->nRingAux);

    return self->ringAux[aInd].data;
}
#endif

#if (!defined(CxmUseInlines))
void *
CxTreeAuxGet(CxtTreeObject *self, unsigned aInd);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxTreeAuxGet(CxtTreeObject *self, unsigned aInd)
{
    void *rVal;

    CxmAssert(aInd < self->nTreeAux);

    if (self->nAux <= aInd)
    {
	rVal = NULL;
	goto RETURN;
    }

    rVal = self->aux[aInd];
    RETURN:
    return rVal;
}
#endif

bool
CxTreeAuxSet(CxtTreeObject *self, unsigned aInd, void *aAux);

// Node.
CxtNodeObject *
CxNodeNew(CxtTreeObject *aTree);

PyObject *
CxNodeTree(CxtNodeObject *self);

uint32_t
CxNodeTaxonNumGet(CxtNodeObject *self);
PyObject *
CxNodeTaxonNumGetPargs(CxtNodeObject *self);

void
CxNodeTaxonNumSet(CxtNodeObject *self, uint32_t aTaxonNum);
PyObject *
CxNodeTaxonNumSetPargs(CxtNodeObject *self, PyObject *args);

CxtRingObject *
CxNodeRing(CxtNodeObject *self);
PyObject *
CxNodeRingPargs(CxtNodeObject *self);

void
CxNodeRingSet(CxtNodeObject *self, CxtRingObject *aRing);

unsigned
CxNodeDegree(CxtNodeObject *self);
PyObject *
CxNodeDegreePargs(CxtNodeObject *self);

#if (!defined(CxmUseInlines))
void *
CxNodeAuxGet(CxtNodeObject *self, unsigned aInd);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxNodeAuxGet(CxtNodeObject *self, unsigned aInd)
{
    void *rVal;

    CxmAssert(aInd < self->tree->nNodeAux);

    if (self->nAux <= aInd)
    {
	rVal = NULL;
	goto RETURN;
    }

    rVal = self->aux[aInd];
    RETURN:
    return rVal;
}
#endif

bool
CxNodeAuxSet(CxtNodeObject *self, unsigned aInd, void *aAux);

// Edge.
CxtEdgeObject *
CxEdgeNew(CxtTreeObject *aTree);
PyObject *
CxEdgeTree(CxtEdgeObject *self);

void
CxEdgeRingsGet(CxtEdgeObject *self, CxtRingObject **rRingA,
	       CxtRingObject **rRingB);
PyObject *
CxEdgeRingsGetPargs(CxtEdgeObject *self);

double
CxEdgeLengthGet(CxtEdgeObject *self);
PyObject *
CxEdgeLengthGetPargs(CxtEdgeObject *self);

void
CxEdgeLengthSet(CxtEdgeObject *self, double aLength);
PyObject *
CxEdgeLengthSetPargs(CxtEdgeObject *self, PyObject *args);

void
CxEdgeAttach(CxtEdgeObject *self, CxtNodeObject *aNodeA,
	     CxtNodeObject *aNodeB);
PyObject *
CxEdgeAttachPargs(CxtEdgeObject *self, PyObject *args);

bool
CxEdgeDetach(CxtEdgeObject *self);
PyObject *
CxEdgeDetachPargs(CxtEdgeObject *self);

#if (!defined(CxmUseInlines))
void *
CxEdgeAuxGet(CxtEdgeObject *self, unsigned aInd);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxEdgeAuxGet(CxtEdgeObject *self, unsigned aInd)
{
    void *rVal;

    CxmAssert(aInd < self->tree->nEdgeAux);

    if (self->nAux <= aInd)
    {
	rVal = NULL;
	goto RETURN;
    }

    rVal = self->aux[aInd];
    RETURN:
    return rVal;
}
#endif

bool
CxEdgeAuxSet(CxtEdgeObject *self, unsigned aInd, void *aAux);

// Ring.
CxtRingObject *
CxRingNew(CxtEdgeObject *aEdge, uint32_t aEnd);

CxtTreeObject *
CxRingTree(CxtRingObject *self);
PyObject *
CxRingTreePargs(CxtRingObject *self);

CxtNodeObject *
CxRingNode(CxtRingObject *self);
PyObject *
CxRingNodePargs(CxtRingObject *self);

CxtEdgeObject *
CxRingEdge(CxtRingObject *self);
PyObject *
CxRingEdgePargs(CxtRingObject *self);

CxtRingObject *
CxRingOther(CxtRingObject *self);
PyObject *
CxRingOtherPargs(CxtRingObject *self);

CxtRingObject *
CxRingNext(CxtRingObject *self);
PyObject *
CxRingNextPargs(CxtRingObject *self);

CxtRingObject *
CxRingPrev(CxtRingObject *self);
PyObject *
CxRingPrevPargs(CxtRingObject *self);

#if (!defined(CxmUseInlines))
void *
CxRingAuxGet(CxtRingObject *self, unsigned aInd);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxRingAuxGet(CxtRingObject *self, unsigned aInd)
{
    void *rVal;

    CxmAssert(aInd < self->tree->nRingAux);

    if (self->nAux <= aInd)
    {
	rVal = NULL;
	goto RETURN;
    }

    rVal = self->aux[aInd];
    RETURN:
    return rVal;
}
#endif

bool
CxRingAuxSet(CxtRingObject *self, unsigned aInd, void *aAux);
