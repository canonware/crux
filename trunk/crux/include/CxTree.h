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

typedef struct CxsTreeObject CxtTreeObject;
typedef struct CxsNodeObject CxtNodeObject;
typedef struct CxsEdgeObject CxtEdgeObject;
typedef struct CxsRingObject CxtRingObject;

extern PyTypeObject CxtTree;
extern PyTypeObject CxtNode;
extern PyTypeObject CxtEdge;
extern PyTypeObject CxtRing;

struct CxsTreeObject
{
    PyObject_HEAD
    CxtTr *tr;
    uint64_t seq;
    bool GcCleared;

#define CxmTreeObjectAuxCount 1
    void *aux[CxmTreeObjectAuxCount];
};

struct CxsNodeObject
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtTrNode node;
    bool GcCleared;

#define CxmNodeObjectAuxCount 1
    void *aux[CxmNodeObjectAuxCount];
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

#define CxmEdgeObjectAuxCount 1
#define CxmEdgeObjectAuxMp 0
    void *aux[CxmEdgeObjectAuxCount];
};

struct CxsRingObject
{
    PyObject_HEAD
    CxtTreeObject *tree;
    CxtEdgeObject *edge;
    CxtTrRing ring;
    bool GcDetached:1;
    bool GcCleared:1;

#define CxmRingObjectAuxCount 1
#define CxmRingObjectAuxMp 0
    void *aux[CxmRingObjectAuxCount];
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
CxTreeIterate(CxtTreeObject *aTree,
	      CxtTreeIterateNodeCallback *aNodeCallback,
	      CxtTreeIterateEdgeCallback *aEdgeCallback,
	      CxtTreeIterateRingCallback *aRingCallback,
	      void *aContext);

#if (!defined(CxmUseInlines))
uint64_t
CxTreeSeq(CxtTreeObject *self);
void *
CxTreeAuxGet(CxtTreeObject *self, unsigned aInd);
void
CxTreeAuxSet(CxtTreeObject *self, unsigned aInd, void *aAux);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline uint64_t
CxTreeSeq(CxtTreeObject *self)
{
    return self->seq;
}

CxmInline void *
CxTreeAuxGet(CxtTreeObject *self, unsigned aInd)
{
    CxmAssert(aInd < CxmTreeObjectAuxCount);

    return self->aux[aInd];
}

CxmInline void
CxTreeAuxSet(CxtTreeObject *self, unsigned aInd, void *aAux)
{
    CxmAssert(aInd < CxmTreeObjectAuxCount);

    self->aux[aInd] = aAux;
}
#endif

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

unsigned
CxNodeDegree(CxtNodeObject *self);
PyObject *
CxNodeDegreePargs(CxtNodeObject *self);

void
CxNodeRingSet(CxtNodeObject *self, CxtRingObject *aRing);

#if (!defined(CxmUseInlines))
void *
CxNodeAuxGet(CxtNodeObject *self, unsigned aInd);
void
CxNodeAuxSet(CxtNodeObject *self, unsigned aInd, void *aAux);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxNodeAuxGet(CxtNodeObject *self, unsigned aInd)
{
    CxmAssert(aInd < CxmNodeObjectAuxCount);

    return self->aux[aInd];
}

CxmInline void
CxNodeAuxSet(CxtNodeObject *self, unsigned aInd, void *aAux)
{
    CxmAssert(aInd < CxmNodeObjectAuxCount);

    self->aux[aInd] = aAux;
}
#endif

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
void
CxEdgeAuxSet(CxtEdgeObject *self, unsigned aInd, void *aAux);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxEdgeAuxGet(CxtEdgeObject *self, unsigned aInd)
{
    CxmAssert(aInd < CxmEdgeObjectAuxCount);

    return self->aux[aInd];
}

CxmInline void
CxEdgeAuxSet(CxtEdgeObject *self, unsigned aInd, void *aAux)
{
    CxmAssert(aInd < CxmEdgeObjectAuxCount);

    self->aux[aInd] = aAux;
}
#endif

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
void
CxRingAuxSet(CxtRingObject *self, unsigned aInd, void *aAux);
#endif
#if (defined(CxmUseInlines) || defined(CxmTree_c))
CxmInline void *
CxRingAuxGet(CxtRingObject *self, unsigned aInd)
{
    CxmAssert(aInd < CxmRingObjectAuxCount);

    return self->aux[aInd];
}

CxmInline void
CxRingAuxSet(CxtRingObject *self, unsigned aInd, void *aAux)
{
    CxmAssert(aInd < CxmRingObjectAuxCount);

    self->aux[aInd] = aAux;
}
#endif
