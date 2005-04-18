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

#include "../include/_cruxmodule.h"

typedef struct CxsTreeMpHeld CxtTreeMpHeld;
typedef struct CxsTreeMpData CxtTreeMpData;
typedef struct CxsTreeMpPs CxtTreeMpPs;

// Neighbor tree.
struct CxsTreeMpHeld
{
    // Neighbor index for the tree.  This can be passed to CxTrTbrNeighborGet()
    // to get the associated TBR parameters.
    unsigned neighbor;

    // Fitch parsimony score for the neighboring tree.
    unsigned score;
};

struct CxsTreeMpData
{
    // XXX Use these.
    //PyObject *cTMatrix;
    //uint64_t cTMatrixSeq;

    unsigned edgeAuxInd;
    unsigned ringAuxInd;

    char **taxa;
    unsigned nTaxa;
    unsigned nCharsTotal;

    // Settings/data for masking of uninformative characters.
    bool elimUninform;
    bool *charsMask;
    unsigned nInformative;

    // Number of characters actually used (equal to nCharsTotal or nInformative,
    // depending on the value of elimUninform).
    unsigned nChars;

    // Tree holding data.
    CxtTreeMpHeld *held;
    unsigned heldLen;
    unsigned nHeld;
};

// Character (in the systematics sense of the word).
typedef char CxtTreeMpC;

// Partial parsimony score information.
struct CxsTreeMpPs
{
    // Parent which most recently used this node's partial score when caching
    // its results.  Both children must still point to the parent in order for
    // the cached results to be valid.
    CxtTreeMpPs *parent;

    // Sum of the subtree scores, and this node's score, given particular
    // children.  In order for this to be useful, both childrens' parent
    // pointers must still point to this node.
    unsigned subtreesScore;
    unsigned nodeScore;

    // chars points to an array of Fitch parsimony state sets.  Each element in
    // the array contains a bitmap representation of a subset of {ACGT} in the 4
    // least significant bits.  T is the least significant bit.  1 means that a
    // nucleotide is in the set.
    //
    // There are nchars character state sets.
    //
    // achars is the actual allocation, which is padded in order to
    // be able to guarantee that chars is 16 byte-aligned.
    CxtTreeMpC *chars;
    unsigned nChars;
    CxtTreeMpC *aChars;
};

// Specifies different tree holding strategies.
typedef enum
{
    CxeTreeHoldBest,
    CxeTreeHoldBetter,
    CxeTreeHoldAll
} CxtTreeHoldHow;

#define CxmTreeMpMaxScoreNone 0xffffffffU
#define CxmTreeMpHoldAll 0xffffffffU

static bool
CxpTreeMpDataGet(CxtTreeObject *self, CxtTreeMpData **rData);

//==============================================================================
//
// CxTreeMpPs.
//

CxmpInline CxtTreeMpPs *
CxpTreeMpPsNew(void)
{
    CxtTreeMpPs *rVal;

    rVal = (CxtTreeMpPs *) malloc(sizeof(CxtTreeMpPs));
    if (rVal == NULL)
    {
	PyErr_NoMemory();
	goto RETURN;
    }

    rVal->parent = NULL;
    rVal->chars = NULL;

    RETURN:
    return rVal;
}

CxmpInline void
CxpTreeMpPsDelete(CxtTreeMpPs *aPs)
{
    if (aPs->chars != NULL)
    {
	free(aPs->aChars);
    }

    free(aPs);
}

#if (0) // Unused (so far).
CxmpInline CxtTreeMpC
CxpTreePsCharGet(CxtTreeMpPs *aPs, unsigned aOffset)
{
    CxtTreeMpC rVal;

    CxmCheckPtr(aPs->chars);
    CxmAssert(aOffset < aPs->nchars);

    rVal = aPs->chars[aOffset >> 1];
    rVall >>= ((aOffset & 1) * 4);

    return rVal;
}
#endif

CxmpInline void
CxpTreePsCharSet(CxtTreeMpPs *aPs, CxtTreeMpC aChar, unsigned aOffset)
{
    CxmCheckPtr(aPs->chars);
    CxmAssert((aChar & 0xfU) == aChar);
    CxmAssert(aOffset < aPs->nChars + ((32 - (aPs->nChars & 0x1fU)) & 0x1fU));

    if ((aOffset & 1) == 0)
    {
	aPs->chars[aOffset >> 1]
	    = (aChar << 4)
	    | (aPs->chars[aOffset >> 1] & 0xfU);
    }
    else
    {
	aPs->chars[aOffset >> 1]
	    = (aPs->chars[aOffset >> 1] & 0xf0U)
	    | aChar;
    }
}

static bool
CxpTreePsPrepare(CxtTreeMpData *aData, CxtTreeMpPs *aPs)
{
    bool rVal;

    // Clean up old character vector if it isn't the right size for
    // aData->nChars characters.
    if (aPs->chars != NULL && aPs->nChars != aData->nChars)
    {
	free(aPs->aChars);
	aPs->chars = NULL;
    }

    // Allocate character vector if necessary.
    if (aPs->chars == NULL)
    {
	if (aData->nChars != 0)
	{
	    unsigned nPad, i;

	    // Calculate the number of pad bytes to append, such that the total
	    // number of bytes is a multiple of 16 (total number of taxonomical
	    // characters is a multiple of 32).
	    nPad = (32 - (aData->nChars & 0x1fU)) & 0x1fU;
	    CxmAssert(((aData->nChars + nPad) & 0x1fU) == 0);

	    // Tack on 8 bytes; all modern systems provide at least 8 byte
	    // alignment.
	    aPs->aChars = (CxtTreeMpC *) malloc(sizeof(CxtTreeMpC)
						* (((aData->nChars + nPad)
						    >> 1))
						+ 8);
	    if (aPs->aChars == NULL)
	    {
		PyErr_NoMemory();
		rVal = true;
		goto RETURN;
	    }

	    // Make sure that aPs->chars is 16 byte-aligned.  Assume that
	    // aPs->aChars is at least 8 byte-aligned.
	    aPs->chars = &aPs->aChars[((unsigned) aPs->aChars) & 0xfU];

	    aPs->nChars = aData->nChars + nPad;

	    // Set all pad characters to {ACGT}.  This allows the pad characters
	    // to be calculated along with the actual characters, without
	    // affecting the score.
	    for (i = aData->nChars; i < aData->nChars + nPad; i++)
	    {
		CxpTreePsCharSet(aPs, 0xfU, i);
	    }
	}
    }

    // Make sure that cached values are not used.
    aPs->parent = NULL;

    rVal = false;
    RETURN:
    return rVal;
}

//
// End CxTreeMpPs.
//
//==============================================================================
//
// CxTreeMp.
//

static void
CxpTreeMpCleanupFinal(CxtTreeObject *aTree, void *aData, unsigned aInd)
{
    CxtTreeMpData *data = (CxtTreeMpData *) aData;

    if (data->held != NULL)
    {
	free(data->held);
    }

    if (data->charsMask != NULL)
    {
	free(data->charsMask);
    }

    if (data->taxa != NULL)
    {
	free(data->taxa);
    }

    free(data);
}

static bool
CxpTreeMpInitEdge(CxtEdgeObject *aEdge, unsigned aInd)
{
    bool rVal;
    CxtTreeMpData *data;
    CxtTreeMpPs *ps;

    if (CxpTreeMpDataGet(CxEdgeTree(aEdge), &data))
    {
	rVal = true;
	goto RETURN;
    }

    ps = CxEdgeAuxGet(aEdge, aInd);
    if (ps == NULL)
    {
	ps = CxpTreeMpPsNew();
	if (ps == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}

	if (CxEdgeAuxSet(aEdge, aInd, ps))
	{
	    rVal = true;
	    goto RETURN;
	}
    }

    if (CxpTreePsPrepare(data, ps))
    {
	rVal = true;
	goto RETURN;
    }

    rVal = false;
    RETURN:
    return rVal;
}

static void
CxpTreeMpCleanupEdge(CxtEdgeObject *aEdge, void *aData, unsigned aInd)
{
    CxtTreeMpPs *ps = (CxtTreeMpPs *) aData;

    if (ps != NULL)
    {
	CxpTreeMpPsDelete(ps);
	CxEdgeAuxSet(aEdge, aInd, NULL);
    }
}

static bool
CxpTreeMpInitRing(CxtRingObject *aRing, unsigned aInd)
{
    bool rVal;
    CxtTreeMpData *data;
    CxtTreeMpPs *ps;
    uint32_t taxonNum;

    if (CxpTreeMpDataGet(CxRingTree(aRing), &data))
    {
	rVal = true;
	goto RETURN;
    }

    ps = CxRingAuxGet(aRing, aInd);
    if (ps == NULL)
    {
	ps = CxpTreeMpPsNew();
	if (ps == NULL)
	{
	    rVal = true;
	    goto RETURN;
	}

	if (CxRingAuxSet(aRing, aInd, ps))
	{
	    rVal = true;
	    goto RETURN;
	}
    }

    if (CxpTreePsPrepare(data, ps))
    {
	rVal = true;
	goto RETURN;
    }

    // If this is a leaf node, initialize the character state sets and scores.
    taxonNum = CxNodeTaxonNumGet(CxRingNode(aRing));
    if (taxonNum != CxmTrNodeNone)
    {
	unsigned i, j;
	char *chars;

	ps->subtreesScore = 0;
	ps->nodeScore = 0;

	CxmAssert(taxonNum < data->nTaxa);
	chars = data->taxa[taxonNum];
	for (i = j = 0; i < data->nCharsTotal; i++)
	{
	    // Ignore uninformative characters.
	    if (data->elimUninform && data->charsMask[i] == false)
	    {
		continue;
	    }

	    switch (chars[i])
	    {
		case 'N':
		case 'n':
		case 'X':
		case 'x':
		    // Treat gaps as uncertainty.  This isn't the only way to do
		    // things, and may need to be made configurable.
		case '-':
		{
		    CxpTreePsCharSet(ps, 0xf, j);
		    break;
		}
		case 'V':
		case 'v':
		{
		    CxpTreePsCharSet(ps, 0xe, j);
		    break;
		}
		case 'H':
		case 'h':
		{
		    CxpTreePsCharSet(ps, 0xd, j);
		    break;
		}
		case 'M':
		case 'm':
		{
		    CxpTreePsCharSet(ps, 0xc, j);
		    break;
		}
		case 'D':
		case 'd':
		{
		    CxpTreePsCharSet(ps, 0xb, j);
		    break;
		}
		case 'R':
		case 'r':
		{
		    CxpTreePsCharSet(ps, 0xa, j);
		    break;
		}
		case 'W':
		case 'w':
		{
		    CxpTreePsCharSet(ps, 0x9, j);
		    break;
		}
		case 'A':
		case 'a':
		{
		    CxpTreePsCharSet(ps, 0x8, j);
		    break;
		}
		case 'B':
		case 'b':
		{
		    CxpTreePsCharSet(ps, 0x7, j);
		    break;
		}
		case 'S':
		case 's':
		{
		    CxpTreePsCharSet(ps, 0x6, j);
		    break;
		}
		case 'Y':
		case 'y':
		{
		    CxpTreePsCharSet(ps, 0x5, j);
		    break;
		}
		case 'C':
		case 'c':
		{
		    CxpTreePsCharSet(ps, 0x4, j);
		    break;
		}
		case 'K':
		case 'k':
		{
		    CxpTreePsCharSet(ps, 0x3, j);
		    break;
		}
		case 'G':
		case 'g':
		{
		    CxpTreePsCharSet(ps, 0x2, j);
		    break;
		}
		case 'T':
		case 't':
		{
		    CxpTreePsCharSet(ps, 0x1, j);
		    break;
		}
		default:
		{
		    CxmNotReached();
		}
	    }
	    j++;
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

static void
CxpTreeMpCleanupRing(CxtRingObject *aRing, void *aData, unsigned aInd)
{
    CxtTreeMpPs *ps = (CxtTreeMpPs *) aData;

    if (ps != NULL)
    {
	CxpTreeMpPsDelete(ps);
	CxRingAuxSet(aRing, aInd, NULL);
    }
}

static bool
CxpTreeMpDataGet(CxtTreeObject *self, CxtTreeMpData **rData)
{
    bool rVal;
    unsigned treeAuxInd;
    CxtTreeMpData *data;

    // Get aux indices for MP data.
    if (CxTreeAuxSearch(self, "MP", &treeAuxInd))
    {
	// No aux registration.
	data = (CxtTreeMpData *) malloc(sizeof(CxtTreeMpData));
	if (data == NULL)
	{
	    PyErr_NoMemory();
	    rVal = true;
	    goto RETURN;
	}

	// Create MP-specific aux mappings for edges and rings.
	if (CxTreeEdgeAuxRegister(self, "MP", NULL,
				  CxpTreeMpInitEdge, NULL,
				  CxpTreeMpCleanupEdge,
				  &data->edgeAuxInd)
	    || CxTreeRingAuxRegister(self, "MP", NULL,
				     CxpTreeMpInitRing, NULL,
				     CxpTreeMpCleanupRing,
				     &data->ringAuxInd)
	    || CxTreeAuxRegister(self, "MP", (void *) data,
				 CxpTreeMpCleanupFinal, NULL, &treeAuxInd))
	{
	    free(data);
	    rVal = true;
	    goto RETURN;
	}

	// Initialize data.
	data->taxa = NULL;
	data->nTaxa = 0;
	data->nCharsTotal = 0;
	data->elimUninform = false;
	data->charsMask = NULL;
	data->nInformative = 0;
	data->nChars = 0;
	data->held = NULL;
	data->heldLen = 0;
	data->nHeld = 0;
    }
    else
    {
	data = (CxtTreeMpData *) CxTreeAuxData(self, treeAuxInd);
    }

    *rData = data;
    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpTreeMpPrepareRecurse(CxtTreeObject *self, CxtTreeMpData *aData,
			CxtRingObject *aRing)
{
    bool rVal;
    CxtEdgeObject *edge;
    CxtRingObject *ring;

    // Prepare aRing.
    if (CxpTreeMpInitRing(aRing, aData->ringAuxInd))
    {
	rVal = true;
	goto RETURN;
    }

    // Recurse into subtrees.
    for (ring = CxRingNext(aRing); ring != aRing; ring = CxRingNext(ring))
    {
	// Prepare edge before recursing.
	edge = CxRingEdge(ring);
	if (CxpTreeMpInitEdge(edge, aData->edgeAuxInd))
	{
	    rVal = true;
	    goto RETURN;
	}

	// Prepare ring.
	if (CxpTreeMpInitRing(ring, aData->ringAuxInd))
	{
	    rVal = true;
	    goto RETURN;
	}

	// Recurse.
	if (CxpTreeMpPrepareRecurse(self, aData, CxRingOther(ring)))
	{
	    rVal = true;
	    goto RETURN;
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

static bool
CxpTreeMpPrepare(CxtTreeObject *self, bool aElimUninform,
		 char **aTaxa, uint32_t aNTaxa, uint32_t aNChars)
{
    bool rVal;
    CxtTreeMpData *data;
    CxtNodeObject *node;

    if (CxpTreeMpDataGet(self, &data))
    {
	rVal = true;
	goto RETURN;
    }

    if ((node = CxTreeBaseGet(self)) != NULL)
    {
	CxtRingObject *ringFirst, *ring;
	CxtEdgeObject *edge;
	uint32_t i, j, k, x, y;
	uint32_t codes[15];

	// Set fields in data.
	if (data->taxa != NULL)
	{
	    free(data->taxa);
	}
	data->taxa = aTaxa;

	data->nCharsTotal = aNChars;
	data->nTaxa = aNTaxa;

	data->elimUninform = aElimUninform;
	if (data->charsMask != NULL)
	{
	    free(data->charsMask);
	}
	data->charsMask = (bool *) malloc(sizeof(bool) * aNChars);
	if (data->charsMask == NULL)
	{
	    PyErr_NoMemory();
	    rVal = true;
	    goto RETURN;
	}

	// Create a mask of informative characters.
	data->nInformative = 0;
	for (i = 0; i < aNChars; i++)
	{
	    for (k = 0; k < 15; k++)
	    {
		codes[k] = 0;
	    }

	    for (j = 0; j < aNTaxa; j++)
	    {
		switch (aTaxa[j][i])
		{
		    case 'N':
		    case 'n':
		    case 'X':
		    case 'x':
		    case '-':
		    {
			// Treat gaps as uncertainty.  This isn't the only way
			// to do things, and may need to be made configurable.
			break;
		    }
		    case 'V':
		    case 'v':
		    {
			codes[14]++;
			break;
		    }
		    case 'H':
		    case 'h':
		    {
			codes[13]++;
			break;
		    }
		    case 'M':
		    case 'm':
		    {
			codes[12]++;
			break;
		    }
		    case 'D':
		    case 'd':
		    {
			codes[11]++;
			break;
		    }
		    case 'R':
		    case 'r':
		    {
			codes[10]++;
			break;
		    }
		    case 'W':
		    case 'w':
		    {
			codes[9]++;
			break;
		    }
		    case 'A':
		    case 'a':
		    {
			codes[8]++;
			break;
		    }
		    case 'B':
		    case 'b':
		    {
			codes[7]++;
			break;
		    }
		    case 'S':
		    case 's':
		    {
			codes[6]++;
			break;
		    }
		    case 'Y':
		    case 'y':
		    {
			codes[5]++;
			break;
		    }
		    case 'C':
		    case 'c':
		    {
			codes[4]++;
			break;
		    }
		    case 'K':
		    case 'k':
		    {
			codes[3]++;
			break;
		    }
		    case 'G':
		    case 'g':
		    {
			codes[2]++;
			break;
		    }
		    case 'T':
		    case 't':
		    {
			codes[1]++;
			break;
		    }
		    default:
		    {
			CxmNotReached();
		    }
		}
	    }

	    // Count the number of states in which two or more taxa exist.
	    data->charsMask[i] = false;
	    for (x = 1; x < 15; x++)
	    {
		for (y = 1; y < 15; y++)
		{
		    if ((x & y) == 0 && codes[x] >= 2 && codes[y] >= 2)
		    {
			if (data->charsMask[i] == false)
			{
			    data->charsMask[i] = true;
			    data->nInformative++;
			}
		    }
		}
	    }
	}

	if (aElimUninform)
	{
	    data->nChars = data->nInformative;
	}
	else
	{
	    data->nChars = data->nCharsTotal;
	}

	// Prepare the tree.
	ringFirst = CxNodeRing(node);
	if (ringFirst != NULL)
	{
	    ring = ringFirst;
	    do
	    {
		// Prepare edge before recursing.
		edge = CxRingEdge(ring);
		if (CxpTreeMpInitEdge(edge, data->edgeAuxInd))
		{
		    rVal = true;
		    goto RETURN;
		}

		// Prepare ring.
		if (CxpTreeMpInitRing(ring, data->ringAuxInd))
		{
		    rVal = true;
		    goto RETURN;
		}

		// Recurse.
		if (CxpTreeMpPrepareRecurse(self, data, CxRingOther(ring)))
		{
		    rVal = true;
		    goto RETURN;
		}

		ring = CxRingNext(ring);
	    } while (ring != ringFirst);
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

PyObject *
CxTreeMpPrepare(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    PyObject *taxa, *tobj;
    uint32_t elim, ntaxa, i, j;
    char **tarr
#ifdef CxmCcSilence
	= NULL
#endif
	;

    if (PyArg_ParseTuple(args, "O!i", &PyList_Type, &taxa, &elim) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }

    ntaxa = PyList_Size(taxa);

    // Make sure that all taxa have the same number of characters.
    if (ntaxa > 0)
    {
	uint32_t nchars
#ifdef CxmCcSilence
	    = 0
#endif
	    ;
	
	tobj = PyList_GetItem(taxa, 0);
	// Don't worry about raising ValueError error here, since the for
	// loop below will do so.
	if (PyString_Check(tobj))
	{
	    nchars = PyString_Size(tobj);
	}

	// Create an array of string pointers.
	tarr = (char **) malloc(sizeof(char *) * ntaxa);
	if (tarr == NULL)
	{
	    rVal = PyErr_NoMemory();
	    goto RETURN;
	}

	for (i = 0; i < ntaxa; i++)
	{
	    tobj = PyList_GetItem(taxa, i);
	    if (PyString_Check(tobj) == 0 || PyString_Size(tobj) != nchars)
	    {
		free(tarr);
		CxError(CxgTreeValueError,
			"Character string expected");
		rVal = NULL;
		goto RETURN;
	    }
	    tarr[i] = PyString_AsString(tobj);

	    // Make sure characters are valid codes.
	    for (j = 0; j < nchars; j++)
	    {
		switch (tarr[i][j])
		{
		    case 'N':
		    case 'n':
		    case 'X':
		    case 'x':
		    case 'V':
		    case 'v':
		    case 'H':
		    case 'h':
		    case 'M':
		    case 'm':
		    case 'D':
		    case 'd':
		    case 'R':
		    case 'r':
		    case 'W':
		    case 'w':
		    case 'A':
		    case 'a':
		    case 'B':
		    case 'b':
		    case 'S':
		    case 's':
		    case 'Y':
		    case 'y':
		    case 'C':
		    case 'c':
		    case 'K':
		    case 'k':
		    case 'G':
		    case 'g':
		    case 'T':
		    case 't':
		    case '-':
		    {
			break;
		    }
		    default:
		    {
			free(tarr);
			CxError(CxgTreeValueError,
				"Invalid character '%c'", tarr[i][j]);
			rVal = NULL;
			goto RETURN;
		    }
		}
	    }
	}

	// XXX Recurse through the tree and make sure that the taxa are
	// numbered correctly.

	// Do preparation.  tarr is managed internally from here on, so it
	// shouldn't be freed here unless there is an error.
	if (CxpTreeMpPrepare(self, elim, tarr, ntaxa, nchars))
	{
	    free(tarr);
	    rVal = NULL;
	    goto RETURN;
	}
    }

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
}

static void
CxpTreeMpFinishRecurse(CxtTreeMpData *aData, CxtRingObject *aRing)
{
    CxtRingObject *ring;
    CxtEdgeObject *edge;

    // Clean up aRing.
    CxpTreeMpCleanupRing(aRing, CxRingAuxGet(aRing, aData->ringAuxInd),
			 aData->ringAuxInd);

    // Recurse into subtrees.
    for (ring = CxRingNext(aRing); ring != aRing; ring = CxRingNext(ring))
    {
	// Clean up edge before recursing.
	edge = CxRingEdge(ring);
	CxpTreeMpCleanupEdge(edge, CxEdgeAuxGet(edge, aData->edgeAuxInd),
			     aData->edgeAuxInd);

	// Clean up ring.
	CxpTreeMpCleanupRing(ring, CxRingAuxGet(ring, aData->ringAuxInd),
			     aData->ringAuxInd);

	// Recurse.
	CxpTreeMpFinishRecurse(aData, CxRingOther(ring));
    }
}

PyObject *
CxTreeMpFinish(CxtTreeObject *self)
{
    PyObject *rVal;
    CxtNodeObject *base;
    CxtTreeMpData *data;

    if (CxpTreeMpDataGet(self, &data))
    {
	rVal = NULL;
	goto RETURN;
    }

    base = CxTreeBaseGet(self);
    if (base != NULL)
    {
	CxtRingObject *ringFirst, *ring;
	CxtEdgeObject *edge;

	// Clean up the tree.
	ringFirst = CxNodeRing(base);
	if (ringFirst != NULL)
	{
	    ring = ringFirst;
	    do
	    {
		// Clean up edge before recursing.
		edge = CxRingEdge(ring);
		CxpTreeMpCleanupEdge(edge, CxEdgeAuxGet(edge, data->edgeAuxInd),
				     data->edgeAuxInd);

		// Clean up ring.
		CxpTreeMpCleanupRing(ring, CxRingAuxGet(ring, data->ringAuxInd),
				     data->ringAuxInd);

		// Recurse.
		CxpTreeMpFinishRecurse(data, CxRingOther(ring));

		ring = CxRingNext(ring);
	    } while (ring != ringFirst);
	}
    }

    if (data->held != NULL)
    {
	free(data->held);
	data->held = NULL;
	data->heldLen = 0;
	data->nHeld = 0;
    }

    if (data->charsMask != NULL)
    {
	free(data->charsMask);
	data->charsMask = NULL;
	data->nInformative = 0;
    }

    if (data->taxa != NULL)
    {
	free(data->taxa);
	data->taxa = NULL;
	data->nTaxa = 0;
	data->nCharsTotal = 0;
	data->nChars = 0;
    }

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
}

#ifdef CxmCpuIa32
CxmpInline void
CxpTreeMpIa32PScore(CxtTreeMpPs *aP, CxtTreeMpPs *aA, CxtTreeMpPs *aB)
{
    unsigned curlimit, i, nbytes, ns;
    CxtTreeMpC *charsP, *charsA, *charsB;

    // Calculate node score.
    ns = 0;

    // Calculate partial Fitch parsimony scores for each character.
    charsP = aP->chars;
    charsA = aA->chars;
    charsB = aB->chars;

    nbytes = (aP->nChars >> 1);

    // Initialize SSE2 registers.
    {
	static const unsigned char low[] __attribute__ ((aligned (16))) =
	    "\x0f\x0f\x0f\x0f\x0f\x0f\x0f\x0f"
	    "\x0f\x0f\x0f\x0f\x0f\x0f\x0f\x0f";

	asm volatile (
	    // Clear pns.
	    "pxor %%xmm4, %%xmm4;"

	    // Fill xmm5 with masks for the least significant four bits of each
	    // byte.
	    "movdqa %[low], %%xmm5;"

	    // Fill xmm6 with masks for the most significant four bits of each
	    // byte.
	    "pcmpeqb %%xmm6, %%xmm6;"
	    "pxor %%xmm5, %%xmm6;"

	    // Fill xmm7 with 16 1's.
	    "pxor %%xmm7, %%xmm7;"
	    "pcmpeqb %%xmm0, %%xmm0;"
	    "psubb %%xmm0, %%xmm7;"
	    :
	    : [low] "m" (*low)
	    : "%xmm4", "%xmm5", "%xmm6", "%xmm7"
	    );
    }

    // The inner loop can be run a maximum of 127 times before the partial node
    // score results (stored in %xmm4) are added to ns (otherwise, overflow
    // could occur).  Therefore, the outer loop calculates the upper bound for
    // the inner loop, thereby avoiding extra computation in the inner loop.
    curlimit = 127 * 16;
    if (curlimit > nbytes)
    {
	curlimit = nbytes;
    }
    for (i = 0;;)
    {
	// Use SSE2 to evaluate the characters.  This loop handles 32 characters
	// per iteration.
	for (; i < curlimit; i += 16)
	{
	    asm volatile (
		// Read character data, and'ing and or'ing them together.
		//
		// a = *charsA;
		// b = *charsB;
		// p = a & b;
		// d = a | b;
		"movdqa %[a], %%xmm0;"
		"movdqa %%xmm0, %%xmm1;"
		"por %[b], %%xmm0;" // xmm0 contains d.
		"pand %[b], %%xmm1;" // xmm1 contains p.

		//==============================================================
		// Most significant bits.

		// Create bitmasks according to whether the character state sets
		// are empty.
		//
		// c = p ? 0x00 : 0xff;
		// e = (c & d);
		// s = c ? 0 : 1;
		// p = (p | e);
		// ns += s;
		"pxor %%xmm2, %%xmm2;"
		"movdqa %%xmm1, %%xmm3;"
		"pand %%xmm6, %%xmm3;" // Mask out unwanted bits of p.
		"pcmpeqb %%xmm3, %%xmm2;" // xmm2 contains c.
		"movdqa %%xmm0, %%xmm3;"
		"pand %%xmm6, %%xmm3;" // Mask out unwanted bits of d.
		"pand %%xmm2, %%xmm3;" // xmm3 contains e.
		"por %%xmm3, %%xmm1;" // Update p.
		"pand %%xmm7, %%xmm2;" // xmm2 contains s.
		"paddusb %%xmm2, %%xmm4;" // Update ns (add s).

		//==============================================================
		// Least significant bits.

		// Create bitmasks according to whether the character state sets
		// are empty.
		//
		// c = p ? 0x00 : 0xff;
		// e = (c & d);
		// s = c ? 0 : 1;
		// p = (p | e);
		// ns += s;
		"pxor %%xmm2, %%xmm2;"
		"movdqa %%xmm1, %%xmm3;"
		"pand %%xmm5, %%xmm3;" // Mask out unwanted bits of p.
		"pcmpeqb %%xmm3, %%xmm2;" // xmm2 contains c.
		"pand %%xmm5, %%xmm0;" // Mask out unwanted bits of d.
		"pand %%xmm2, %%xmm0;" // xmm0 contains e.
		"por %%xmm0, %%xmm1;" // Update p.
		"pand %%xmm7, %%xmm2;" // xmm2 contains s.
		"paddusb %%xmm2, %%xmm4;" // Update ns (add s).

		// Store results.
		//
		// *charsP = p;
		"movdqa %%xmm1, %[p];"
		: [p] "=m" (charsP[i])
		: [a] "m" (charsA[i]), [b] "m" (charsB[i])
		: "memory"
		);
	}

	// Update ns and reset pns.
	{
	    uint32_t pns;

	    asm volatile (
		// Sum the upper 8 bytes and lower 8 bytes separately (that's
		// what psadbw does).
		"pxor %%xmm0, %%xmm0;"
		"psadbw %%xmm0, %%xmm4;"

		// Combine the results of psadbw.
		"movdqa %%xmm4, %%xmm3;"
		"punpckhqdq %%xmm0, %%xmm4;"
		"paddq %%xmm3, %%xmm4;"

		// Store the result.
		"movd %%xmm4, %[pns];"
		"pxor %%xmm4, %%xmm4;"
		: [pns] "=r" (pns)
		:
		: "memory"
		);

	    ns += pns;
	}

	// Break out of the loop if the bound for the inner loop was the maximum
	// possible.
	if (curlimit == nbytes)
	{
	    break;
	}
	// Update the bound for the inner loop, taking care not to exceed the
	// maximum possible bound.
	curlimit += 127 * 16;
	if (curlimit > nbytes)
	{
	    curlimit = nbytes;
	}
    }

    aP->nodeScore = ns;
}
#endif

static void
CxpTreeMpCPScore(CxtTreeMpPs *aP, CxtTreeMpPs *aA, CxtTreeMpPs *aB)
{
    uint32_t i, nwords, ns, a, b, m, r, un;
    uint32_t *charsP, *charsA, *charsB;
    static const uint32_t bitsTable[] =
    {
	2, 1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1,
	1, 0
    };

    // Calculate node score.
    ns = 0;

#define CxmTreeMpCPScoreInner()						\
    a = charsA[i];							\
    b = charsB[i];							\
									\
    /* Get 1's in the least significant bits of state sets that are	\
     * non-empty after the intersection operation. */			\
    r = m = a & b;							\
    m |= (m >> 1);							\
    m |= (m >> 1);							\
    m |= (m >> 1);							\
									\
    /* Mask out garbage. */						\
    m &= 0x11111111;							\
									\
    /* Count up changes. */						\
    ns += bitsTable[m & 0xff]						\
	+ bitsTable[(m >> 8) & 0xff]					\
	+ bitsTable[(m >> 16) & 0xff]					\
	+ bitsTable[(m >> 24) & 0xff];					\
									\
    /* Propagate 1's to make a bit mask. */				\
    m |= (m << 1);							\
    m |= (m << 1);							\
    m |= (m << 1);							\
									\
    /* Combine results of intersection and union operations. */		\
    r &= m;								\
    un = a | b;								\
    m = (~m);								\
    un &= m;								\
    r |= un;								\
									\
    /* Store result. */							\
    charsP[i] = r;							\
									\
    i++;

    // Calculate preliminary Fitch parsimony scores for each character.
    charsP = (uint32_t *) aP->chars;
    charsA = (uint32_t *) aA->chars;
    charsB = (uint32_t *) aB->chars;
    for (i = 0, nwords = (aP->nChars >> 3); i < nwords;)
    {
	CxmTreeMpCPScoreInner();
	CxmTreeMpCPScoreInner();
	CxmTreeMpCPScoreInner();
	CxmTreeMpCPScoreInner();
    }
#undef CxmTreeMpCPScoreInner

    aP->nodeScore = ns;
}

// Unconditionally calculate the partial score for aP, using aA and aB as
// children.
CxmpInline void
CxpTreeMpPScore(CxtTreeMpPs *aP, CxtTreeMpPs *aA, CxtTreeMpPs *aB)
{
    // Reset this node's parent pointer, to keep the parent from using an
    // invalid cached value.
    aP->parent = NULL;

    // Calculate sum of subtree scores.
    aP->subtreesScore
	= aA->subtreesScore + aA->nodeScore
	+ aB->subtreesScore + aB->nodeScore;

#ifdef CxmCpuIa32
    if (CxgIa32UseSse2)
    {
	CxpTreeMpIa32PScore(aP, aA, aB);
    }
    else
#endif
    {
	CxpTreeMpCPScore(aP, aA, aB);
    }
}

// The sole purpose of this function is to assure that the contents of
// CxpTreeMpPScore() are not inlined in CxpTreeMpCachePScore().  Most of the
// time, the cache should be usable, so the actual scoring code doesn't usually
// get called.
static void
CxpTreeMpNoInlinePScore(CxtTreeMpPs *aP, CxtTreeMpPs *aA, CxtTreeMpPs *aB)
{
    CxpTreeMpPScore(aP, aA, aB);
}

// Calculate the partial score for aP, using aA and aB as children.  However, do
// some extra bookkeeping in order to be able to cache the results, and later
// recognize that precisely the same calculation was cached.
CxmpInline void
CxpTreeMpCachePScore(CxtTreeMpPs *aP, CxtTreeMpPs *aA, CxtTreeMpPs *aB)
{
//#define CxmTreeMpCachePScoreValidate
#ifdef CxmTreeMpCachePScoreValidate
    bool cached;
    unsigned cachedNodeScore;
#endif

    CxmCheckPtr(aP);
    CxmCheckPtr(aA);
    CxmCheckPtr(aB);

    // Only calculate the parent's node score if the cached value is invalid.
    if (aA->parent != aP || aB->parent != aP)
#ifdef CxmTreeMpCachePScoreValidate
    {
	cached = false;
    }
    else
    {
	cached = true;
	cachedNodeScore = aP->nodeScore;

	if (aP->subtreesScore
	    != (aA->subtreesScore + aA->nodeScore
		+ aB->subtreesScore + aB->nodeScore))
	{
	    fprintf(stderr,
		    "%s:%d:%s(): subtreesScore %u (should be %u)\n",
		    __FILE__, __LINE__, __func__,
		    aP->subtreesScore,
		    aA->subtreesScore + aA->nodeScore
		    + aB->subtreesScore + aB->nodeScore);
	    abort();
	}
    }
#endif
    {
	// Set parent pointers, so that cached values may be used in future
	// runs.
	aA->parent = aP;
	aB->parent = aP;

	// Calculate the partial score.
	CxpTreeMpNoInlinePScore(aP, aA, aB);
    }

#ifdef CxmTreeMpCachePScoreValidate
    if (cached)
    {
	if (cachedNodeScore != aP->nodeScore)
	{
	    fprintf(stderr, "%s:%d:%s(): nodeScore %u (should be %u)\n",
		    __FILE__, __LINE__, __func__,
		    cachedNodeScore, aP->nodeScore);
	    abort();
	}
    }
#endif
}

CxmpInline void
CxpTreeMpCacheInvalidate(CxtTreeMpPs *aPs)
{
    CxmCheckPtr(aPs);

    // Reset this node's parent pointer, to keep the old parent from using an
    // invalid cached value.
    aPs->parent = NULL;
}

#ifdef CxmCpuIa32
CxmpInline unsigned
CxpTreeMpIa32FScore(CxtTreeMpPs *aA, CxtTreeMpPs *aB, unsigned aMaxScore)
{
    unsigned rVal, i, nbytes, pns;
    CxtTreeMpC *charsA, *charsB;

    // Calculate sum of subtree scores.
    rVal
	= aA->subtreesScore + aA->nodeScore
	+ aB->subtreesScore + aB->nodeScore;

    // Calculate partial Fitch parsimony scores for each character.
    charsA = aA->chars;
    charsB = aB->chars;

    // Initialize SSE2 registers.
    {
	static const unsigned char low[] __attribute__ ((aligned (16))) =
	    "\x0f\x0f\x0f\x0f\x0f\x0f\x0f\x0f"
	    "\x0f\x0f\x0f\x0f\x0f\x0f\x0f\x0f";

	asm volatile (
	    // Clear pns.
	    "pxor %%xmm4, %%xmm4;"

	    // Fill xmm5 with masks for the least significant four bits of each
	    // byte.
	    "movdqa %[low], %%xmm5;"

	    // Fill xmm6 with masks for the most significant four bits of each
	    // byte.
	    "pcmpeqb %%xmm6, %%xmm6;"
	    "pxor %%xmm5, %%xmm6;"

	    // Fill xmm7 with 16 1's.
	    "pxor %%xmm7, %%xmm7;"
	    "pcmpeqb %%xmm0, %%xmm0;"
	    "psubb %%xmm0, %%xmm7;"
	    :
	    : [low] "m" (*low)
	    : "%xmm4", "%xmm5", "%xmm6", "%xmm7"
	    );
    }

    // Use SSE2 to evaluate the characters.  This loop handles 32 characters per
    // iteration.
    for (i = 0, nbytes = (aA->nChars >> 1); i < nbytes; i += 16)
    {
	asm volatile (
	    // Read character data, and'ing and or'ing them together.
	    //
	    // a = *charsA;
	    // b = *charsB;
	    // p = a & b;
	    "movdqa %[a], %%xmm1;"
	    "pand %[b], %%xmm1;" // xmm1 contains p.

	    //==================================================================
	    // Most significant bits.

	    // Create bitmasks according to whether the character state sets are
	    // empty.
	    //
	    // c = p ? 0x00 : 0xff;
	    // s = c ? 0 : 1;
	    // rVal += s;
	    "pxor %%xmm2, %%xmm2;"
	    "movdqa %%xmm1, %%xmm3;"
	    "pand %%xmm6, %%xmm3;" // Mask out unwanted bits of p.
	    "pcmpeqb %%xmm3, %%xmm2;" // xmm2 contains c.
	    "pand %%xmm7, %%xmm2;" // xmm2 contains s.
	    "paddusb %%xmm2, %%xmm4;" // Update rVal (add s).

	    //==================================================================
	    // Least significant bits.

	    // Create bitmasks according to whether the character state sets are
	    // empty.
	    //
	    // c = p ? 0x00 : 0xff;
	    // s = c ? 0 : 1;
	    // rVal += s;
	    "pxor %%xmm2, %%xmm2;"
	    "pand %%xmm5, %%xmm1;" // Mask out unwanted bits of p.
	    "pcmpeqb %%xmm1, %%xmm2;" // xmm2 contains c.
	    "pand %%xmm7, %%xmm2;" // xmm2 contains s.
	    "paddusb %%xmm2, %%xmm4;" // Update rVal (add s).

	    // Sum the upper 8 bytes and lower 8 bytes separately (there's no
	    // choice in the matter -- that's what psadbw does).
	    "pxor %%xmm0, %%xmm0;"
	    "psadbw %%xmm0, %%xmm4;"

	    // Combine the results of psadbw.
	    "movdqa %%xmm4, %%xmm3;"
	    "punpckhqdq %%xmm0, %%xmm4;"
	    "paddq %%xmm3, %%xmm4;"

	    // Store the result.
	    "movd %%xmm4, %[pns];"
	    "pxor %%xmm4, %%xmm4;"

	    : [pns] "=r" (pns)
	    : [a] "m" (charsA[i]), [b] "m" (charsB[i])
	    : "memory"
	    );

	// Update rVal and terminate if the max score was exceeded.
	rVal += pns;
	if (rVal > aMaxScore)
	{
	    rVal = UINT_MAX;
	    break;
	}
    }

    return rVal;
}
#endif

static unsigned
CxpTreeMpCFScore(CxtTreeMpPs *aA, CxtTreeMpPs *aB, unsigned aMaxScore)
{
    uint32_t rVal, i, nwords, a, b, m;
    uint32_t *charsA, *charsB;
    static const uint32_t bitsTable[] =
    {
	2, 1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1,
	1, 0
    };

    // Calculate sum of subtree scores.
    rVal
	= aA->subtreesScore + aA->nodeScore
	+ aB->subtreesScore + aB->nodeScore;

#define CxmTreeMpCFScoreInner()						\
    a = charsA[i];							\
    b = charsB[i];							\
									\
    /* Get 1's in the least significant bits of state sets that are	\
     * non-empty after the intersection operation. */			\
    m = a & b;								\
    m |= (m >> 1);							\
    m |= (m >> 1);							\
    m |= (m >> 1);							\
									\
    /* Mask out garbage. */						\
    m &= 0x11111111;							\
									\
    /* Count up changes. */						\
    rVal += bitsTable[m & 0xff]						\
	+ bitsTable[(m >> 8) & 0xff]					\
	+ bitsTable[(m >> 16) & 0xff]					\
	+ bitsTable[(m >> 24) & 0xff];					\
									\
    if (rVal > aMaxScore)						\
    {									\
	rVal = UINT_MAX;						\
	break;								\
    }									\
									\
    i++;

    // Calculate partial Fitch parsimony scores for each character.
    charsA = (uint32_t *) aA->chars;
    charsB = (uint32_t *) aB->chars;
    for (i = 0, nwords = (aA->nChars >> 3); i < nwords;)
    {
	CxmTreeMpCFScoreInner();
	CxmTreeMpCFScoreInner();
	CxmTreeMpCFScoreInner();
	CxmTreeMpCFScoreInner();
    }
#undef CxmTreeMpCFScoreInner

    return rVal;
}

// Unconditionally calculate the final score of a tree, using aA and aB as
// children.
CxmpInline unsigned
CxpTreeMpFScore(CxtTreeMpPs *aA, CxtTreeMpPs *aB, unsigned aMaxScore)
{
    unsigned rVal;

#ifdef CxmCpuIa32
    if (CxgIa32UseSse2)
    {
	rVal = CxpTreeMpIa32FScore(aA, aB, aMaxScore);
    }
    else
#endif
    {
	rVal = CxpTreeMpCFScore(aA, aB, aMaxScore);
    }

    return rVal;
}

static CxtTreeMpPs *
CxpTreeMpScoreRecurse(CxtTreeMpData *aData, CxtRingObject *aRing,
		      CxtEdgeObject *aBisect)
{
    CxtTreeMpPs *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    unsigned degree;
    bool adjacent;
    CxtRingObject *ring;

    // Get the degree of the node.  Don't count the bisection edge (only an
    // issue if this node is adjacent to the bisection).
    degree = 1;
    adjacent = false;
    for (ring = CxRingNext(aRing); ring != aRing; ring = CxRingNext(ring))
    {
	if (CxRingEdge(ring) != aBisect)
	{
	    degree++;
	}
	else
	{
	    adjacent = true;
	}
    }

    switch (degree)
    {
	case 1:
	{
	    // Leaf node.  Do nothing.
	    rVal = (CxtTreeMpPs *) CxRingAuxGet(aRing, aData->ringAuxInd);
	    break;
	}
	case 2:
	{
	    // This is a trifurcating node that is adjacent to the bisection.
	    // Return the child node's ps, since this node's ps is
	    // irrelevant.
	    CxmAssert(adjacent);

	    // Clear the cache for the view that is being bypassed.  This is
	    // critical to correctness of the caching machinery, since each view
	    // should never be claimed as the parent of more than two other
	    // views.
	    CxpTreeMpCacheInvalidate((CxtTreeMpPs *)
				     CxRingAuxGet(aRing, aData->ringAuxInd));

	    // Get the ring element that connects to the other portion of the
	    // subtree on this side of the bisection.
	    for (ring = CxRingNext(aRing);
		 ring != aRing;
		 ring = CxRingNext(ring))
	    {
		if (CxRingEdge(ring) != aBisect)
		{
		    rVal = CxpTreeMpScoreRecurse(aData, CxRingOther(ring),
						 aBisect);
		    break;
		}
	    }
	    break;
	}
	case 3:
	{
	    if (adjacent == false)
	    {
		CxtTreeMpPs *psA, *psB;

		// This is a normal trifurcating node.  This is the common case,
		// and is handled separately from the code below for performance
		// reasons.

		// Recursively calculate partial scores for the subtrees.
		ring = CxRingNext(aRing);
		psA = CxpTreeMpScoreRecurse(aData, CxRingOther(ring), aBisect);

		ring = CxRingNext(ring);
		psB = CxpTreeMpScoreRecurse(aData, CxRingOther(ring), aBisect);

		// Calculate the partial score for this node.
		rVal = (CxtTreeMpPs *) CxRingAuxGet(aRing, aData->ringAuxInd);
		CxpTreeMpCachePScore(rVal, psA, psB);

		break;
	    }
	    // Fall through if this node is adjacent to the bisection.
	}
	default:
	{
	    // This is a multifurcating node.
	    CxmError("XXX Not implemented");
	}
    }

    return rVal;
}

CxmpInline bool
CxpTreeMpScore(CxtTreeObject *self, unsigned *rScore)
{
    bool rVal;
    unsigned score;
    CxtNodeObject *base;
    CxtRingObject *ring;

    base = CxTreeBaseGet(self);
    if (base != NULL && (ring = CxNodeRing(base)) != NULL)
    {
	CxtTreeMpData *data;
	CxtEdgeObject *edge;
	CxtRingObject *ringA, *ringB;
	CxtTreeMpPs *psA, *psB;

	if (CxpTreeMpDataGet(self, &data))
	{
	    rVal = true;
	    goto RETURN;
	}

	edge = CxRingEdge(ring);
	CxEdgeRingsGet(edge, &ringA, &ringB);

	// Calculate partial scores for the subtrees on each end of edge.
	psA = CxpTreeMpScoreRecurse(data, CxRingOther(ringA), NULL);
	psB = CxpTreeMpScoreRecurse(data, CxRingOther(ringB), NULL);

	// Calculate the final score.
	score = CxpTreeMpFScore(psA, psB, UINT_MAX);
    }
    else
    {
	score = 0;
    }

    *rScore = score;
    rVal = false;
    RETURN:
    return rVal;
}

CxmpInline void
CxpTreeMpViewsRecurse(CxtTreeMpData *aData, CxtRingObject *aRing,
		      CxtTreeMpPs *aPs, CxtEdgeObject *aBisect)
{
    unsigned degree;
    bool adjacent;
    CxtRingObject *ring;

    // Get the degree of the node.  Don't count the bisection edge (only an
    // issue if this node is adjacent to the bisection).
    degree = 1;
    adjacent = false;
    for (ring = CxRingNext(aRing); ring != aRing; ring = CxRingNext(ring))
    {
	if (CxRingEdge(ring) != aBisect)
	{
	    degree++;
	}
	else
	{
	    adjacent = true;
	}
    }

    switch (degree)
    {
	case 1:
	{
	    // Leaf node.  Do nothing.
	    break;
	}
	case 2:
	{
	    // This is a trifurcating node that is adjacent to the bisection.
	    // Pass the parent's ps when recursing, since this node's ps is
	    // irrelevant.
	    CxmAssert(adjacent);

	    // Get the ring element that connects to the other portion of the
	    // subtree on this side of the bisection.
	    for (ring = CxRingNext(aRing);
		 ring != aRing;
		 ring = CxRingNext(ring))
	    {
		if (CxRingEdge(ring) != aBisect)
		{
		    // Clear the cache for the view that is being bypassed.
		    // This is critical to correctness of the caching machinery,
		    // since each view should never be claimed as the parent of
		    // more than two other views.
		    CxpTreeMpCacheInvalidate((CxtTreeMpPs *)
					     CxRingAuxGet(ring,
							  aData->ringAuxInd));

		    // Recurse.
		    CxpTreeMpViewsRecurse(aData, CxRingOther(ring), aPs,
					  aBisect);
		    break;
		}
	    }
	    break;
	}
	case 3:
	{
	    if (adjacent == false)
	    {
		CxtRingObject *ringA, *ringB;
		CxtRingObject *ringAOther, *ringBOther;
		CxtTreeMpPs *psA, *psB;
		CxtTreeMpPs *psAOther, *psBOther;
		CxtEdgeObject *edgeA, *edgeB;
		CxtTreeMpPs *psEdgeA, *psEdgeB;

		// This is a normal trifurcating node.  This is the common case,
		// and is handled separately from the code below for performance
		// reasons.

		// Get all variables that are necessary for view calculation and
		// recursion.
		ringA = CxRingNext(aRing);
		psA = (CxtTreeMpPs *) CxRingAuxGet(ringA, aData->ringAuxInd);
		ringAOther = CxRingOther(ringA);
		psAOther = (CxtTreeMpPs *) CxRingAuxGet(ringAOther,
							aData->ringAuxInd);
		edgeA = CxRingEdge(ringA);
		psEdgeA = (CxtTreeMpPs *) CxEdgeAuxGet(edgeA,
						       aData->edgeAuxInd);

		ringB = CxRingNext(ringA);
		psB = (CxtTreeMpPs *) CxRingAuxGet(ringB, aData->ringAuxInd);
		ringBOther = CxRingOther(ringB);
		psBOther = (CxtTreeMpPs *) CxRingAuxGet(ringBOther,
							aData->ringAuxInd);
		edgeB = CxRingEdge(ringB);
		psEdgeB = (CxtTreeMpPs *) CxEdgeAuxGet(edgeB,
						       aData->edgeAuxInd);

		// Calculate views and edges, and recurse.
		CxpTreeMpPScore(psA, aPs, psBOther);
		CxpTreeMpPScore(psEdgeA, psA, psAOther);
		CxpTreeMpViewsRecurse(aData, ringAOther, psA, aBisect);

		CxpTreeMpPScore(psB, aPs, psAOther);
		CxpTreeMpPScore(psEdgeB, psB, psBOther);
		CxpTreeMpViewsRecurse(aData, ringBOther, psB, aBisect);

		break;
	    }
	    // Fall through if this node is adjacent to the bisection.
	}
	default:
	{
	    // This is a multifurcating node.
	    CxmError("XXX Not implemented");
	}
    }
}

// Calculate the partial score for each edge in aEdges.  aEdges[0] must either
// be NULL, or the edge connected to the node that is in turn connected
// to the bisection edge.
CxmpInline bool
CxpTreeMpBisectionEdgeList(CxtTreeMpData *aData,
			   CxtEdgeObject **aEdges, unsigned aNEdges,
			   CxtEdgeObject *aBisect, unsigned aMaxScore)
{
    bool rVal;

    if (aEdges[0] != NULL)
    {
	CxtRingObject *ringA, *ringB;
	CxtTreeMpPs *ps, *psA, *psB;

	CxEdgeRingsGet(aEdges[0], &ringA, &ringB);

	// Recursively (post-order traversal) calculate the partial score at
	// each node, as viewed from the first edge in aEdges.  This leaves one
	// valid view at each node, which then makes it possible to calculate
	// the rest of the views during a pre-order traversal of the tree.
	psA = CxpTreeMpScoreRecurse(aData, ringA, aBisect);
	psB = CxpTreeMpScoreRecurse(aData, ringB, aBisect);

	// The first edge must be calculated using psA and psB as children,
	// rather than using the ps's at the ends of the edge.  This is because
	// one of the connected nodes is in turn connected to the bisection
	// edge, which means that the node does not have a useful ps.  The first
	// edge is the only one for which this is an issue, so it is handled
	// here.
	ps = (CxtTreeMpPs *) CxEdgeAuxGet(aEdges[0], aData->edgeAuxInd);
	CxpTreeMpPScore(ps, psA, psB);
	if (ps->subtreesScore + ps->nodeScore > aMaxScore)
	{
	    // Don't bother calculating other views or edge states, since this
	    // subtree exceeds the maximum score that's of interest.
	    rVal = true;
	    goto RETURN;
	}

	// Perform the pre-order traversal, calculating the remaining views that
	// were not calculated by the above post-order traversal, as well as
	// calculating the state sets for the edges along the way.  Take care to
	// pass the appropriate ps's.
	CxpTreeMpViewsRecurse(aData, ringA, psB, aBisect);
	CxpTreeMpViewsRecurse(aData, ringB, psA, aBisect);

#ifdef CxmDebug
	// Validate per-edge partial scores.
	{
	    unsigned i;

	    for (i = 1; i < aNEdges; i++)
	    {
		// All edge partial scores should have the same value, since the
		// location of the root is irrelevant to the score.
		psA = (CxtTreeMpPs *) CxEdgeAuxGet(aEdges[i],
						   aData->edgeAuxInd);
		if (psA->subtreesScore + psA->nodeScore
		    != ps->subtreesScore + ps->nodeScore)
		{
		    fprintf(stderr,
			    "%s:%d:%s(): Expected %u (%u + %u),"
			    " got %u (%u + %u)\n",
			    __FILE__, __LINE__, __func__,
			    ps->subtreesScore + ps->nodeScore,
			    ps->subtreesScore, ps->nodeScore,
			    psA->subtreesScore + psA->nodeScore,
			    psA->subtreesScore, psA->nodeScore);
		    abort();
		}
	    }
	}
#endif
    }

    rVal = false;
    RETURN:
    return rVal;
}

// Hold a tree.  If aMaxHeld is exceeded, the tree is not held.
CxmpInline bool
CxpTreeMpHold(CxtTreeMpData *aData, unsigned aMaxHold, unsigned aNeighbor,
	      unsigned aScore)
{
    bool rVal;

    if (aData->nHeld < aMaxHold)
    {
	CxtTreeMpHeld *held;

	// Make sure there is space to store another held tree.
	if (aData->held == NULL)
	{
	    // Allocate.
	    aData->held = (CxtTreeMpHeld *) malloc(sizeof(CxtTreeMpHeld));
	    if (aData->held == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }
	    aData->heldLen = 1;
	}
	else if (aData->nHeld == aData->heldLen)
	{
	    CxtTreeMpHeld *tHeld;

	    // Reallocate.
	    tHeld = (CxtTreeMpHeld *) realloc(aData->held,
					      sizeof(CxtTreeMpHeld)
					      * aData->heldLen * 2);
	    if (tHeld == NULL)
	    {
		rVal = true;
		goto RETURN;
	    }
	    aData->held = tHeld;
	    aData->heldLen *= 2;
	}

	// Hold this tree.
	held = &aData->held[aData->nHeld];
	held->neighbor = aNeighbor;
	held->score = aScore;

	aData->nHeld++;
    }

    rVal = false;
    RETURN:
    return rVal;
}

// Calculate the Fitch parsimony scores for all TBR neighbors of the tree, and
// hold results according to the function parameters.
static bool
CxpTreeMpTbrNeighbors(CxtTreeObject *self, unsigned aMaxHold,
		      unsigned aMaxScore, CxtTreeHoldHow aHow)
{
    bool rVal;
    CxtTreeMpData *data;
    unsigned neighbor, nEdges, nEdgesA, nEdgesB, i, j, k, curMax, score;
    CxtEdgeObject **edgesA, **edgesB, *bisect, *edgeA, *edgeB;
    CxtRingObject *ringA, *ringB;
    CxtTreeMpPs *psA, *psB;

    if (CxpTreeMpDataGet(self, &data))
    {
	rVal = true;
	goto RETURN;
    }

    curMax = aMaxScore;

    // Set up tree holding data structures.
    data->nHeld = 0;

    if (CxTreeTbrNEdgesGet(self, &nEdges))
    {
	rVal = true;
	goto RETURN;
    }
    for (i = neighbor = 0; i < nEdges; i++)
    {
	if (CxTreeTbrEdgeGet(self, i, &bisect))
	{
	    rVal = true;
	    goto RETURN;
	}

	// Determine which edges are in each subtree.
	if (CxTreeTbrBEdgeSetsGet(self, bisect,
				  &edgesA, &nEdgesA, &ringA,
				  &edgesB, &nEdgesB, &ringB))
	{
	    rVal = true;
	    goto RETURN;
	}

	// Calculate the partial score for each edge in the edge lists.  Don't
	// bother scoring the trees if either subtree exceeds the max score.
	if (CxpTreeMpBisectionEdgeList(data, edgesA, nEdgesA, bisect, curMax)
	    || CxpTreeMpBisectionEdgeList(data, edgesB, nEdgesB, bisect,
					  curMax))
	{
	    unsigned offsetA, offsetB;

	    if (CxTreeTbrEdgeOffset(self, i, &offsetA)
		|| CxTreeTbrEdgeOffset(self, i + 1, &offsetB))
	    {
		rVal = true;
		goto RETURN;
	    }

	    neighbor += offsetB - offsetA;
	    continue;
	}

	// Iteratively (logically) reconnect every legitimate pairing of edges
	// between the two subtrees and calculate final parsimony scores.
	for (j = 0; j < nEdgesA; j++)
	{
	    edgeA = edgesA[j];
	    if (edgeA != NULL)
	    {
		psA = (CxtTreeMpPs *) CxEdgeAuxGet(edgeA, data->edgeAuxInd);
	    }
	    else
	    {
		psA = (CxtTreeMpPs *) CxRingAuxGet(ringA, data->ringAuxInd);
	    }

	    for (k = 0; k < nEdgesB; k++)
	    {
		// Skip this iteration if the reconnection would result in
		// reversing the bisection.
		if (j == 0 && k == 0)
		{
		    continue;
		}

		edgeB = edgesB[k];
		if (edgeB != NULL)
		{
		    psB = (CxtTreeMpPs *) CxEdgeAuxGet(edgeB, data->edgeAuxInd);
		}
		else
		{
		    psB = (CxtTreeMpPs *) CxRingAuxGet(ringB, data->ringAuxInd);
		}

		// Calculate the final parsimony score for this reconnection.
		score = CxpTreeMpFScore(psA, psB, curMax);

		// Hold the tree, if appropriate.
		switch (aHow)
		{
		    case CxeTreeHoldBest:
		    {
			if (score == curMax)
			{
			    // This tree is as good as those currently held.
			    if (CxpTreeMpHold(data, aMaxHold, neighbor, score))
			    {
				rVal = true;
				goto RETURN;
			    }
			}
			else if (score < curMax)
			{
			    // No trees held, or this tree is better than those
			    // currently held.  Throw away previously held trees
			    // and hold this one.
			    data->nHeld = 0;
			    if (CxpTreeMpHold(data, aMaxHold, neighbor, score))
			    {
				rVal = true;
				goto RETURN;
			    }
			    curMax = score;
			}
			break;
		    }
		    case CxeTreeHoldBetter:
		    {
			if (score <= curMax)
			{
			    // No trees held, or this (neighboring) tree is
			    // better than the tree whose neighbors are being
			    // evaluated.
			    if (CxpTreeMpHold(data, aMaxHold, neighbor, score))
			    {
				rVal = true;
				goto RETURN;
			    }
			}
			break;
		    }
		    case CxeTreeHoldAll:
		    {
			// Hold all trees.
			if (CxpTreeMpHold(data, aMaxHold, neighbor, score))
			{
			    rVal = true;
			    goto RETURN;
			}
			break;
		    }
		    default:
		    {
			CxmNotReached();
		    }
		}

		// Increment the neighbor index here.  Due to the possibility of
		// loop short-circuting above, this must happen at the end of
		// the loop body, rather than in the 'for' statement.
		neighbor++;
	    }
	}
    }

    rVal = false;
    RETURN:
    return rVal;
}

PyObject *
CxTreeMp(CxtTreeObject *self)
{
    PyObject *rVal;
    unsigned score;

    if (CxpTreeMpScore(self, &score))
    {
	rVal = NULL;
	goto RETURN;
    }

    rVal = Py_BuildValue("i", score);
    RETURN:
    return rVal;
}

PyObject *
CxTreeTbrBestNeighborsMp(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    int maxHold;

    maxHold = CxmTrHoldAll;
    if (PyArg_ParseTuple(args, "|i", &maxHold) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }
    if (maxHold < 0 && maxHold != CxmTreeMpHoldAll)
    {
	CxError(CxgTreeValueError,
		"maxHold: non-negative integer expected");
	rVal = NULL;
	goto RETURN;
    }

    if (CxpTreeMpTbrNeighbors(self, maxHold, CxmTreeMpMaxScoreNone,
			      CxeTreeHoldBest))
    {
	rVal = NULL;
	goto RETURN;
    }

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
}

PyObject *
CxTreeTbrBetterNeighborsMp(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    int maxHold;
    unsigned score;

    maxHold = CxmTrHoldAll;
    if (PyArg_ParseTuple(args, "|i", &maxHold) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }
    if (maxHold < 0 && maxHold != CxmTrHoldAll)
    {
	CxError(CxgTreeValueError,
		"maxHold: non-negative integer expected");
	rVal = NULL;
	goto RETURN;
    }

    if (CxpTreeMpScore(self, &score))
    {
	rVal = NULL;
	goto RETURN;
    }

    if (CxpTreeMpTbrNeighbors(self, maxHold,
			      score > 0 ? score - 1 : 0,
			      CxeTreeHoldBetter))
    {
	rVal = NULL;
	goto RETURN;
    }

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
}

PyObject *
CxTreeTbrAllNeighborsMp(CxtTreeObject *self)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;

    if (CxpTreeMpTbrNeighbors(self, CxmTreeMpHoldAll, CxmTreeMpMaxScoreNone,
			      CxeTreeHoldAll))
    {
	rVal = NULL;
	goto RETURN;
    }

    Py_INCREF(Py_None);
    rVal = Py_None;
    RETURN:
    return rVal;
}

PyObject *
CxTreeNHeldGet(CxtTreeObject *self)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    CxtTreeMpData *data;

    if (CxpTreeMpDataGet(self, &data))
    {
	rVal = NULL;
	goto RETURN;
    }

    rVal = Py_BuildValue("i", data->nHeld);

    RETURN:
    return rVal;
}

PyObject *
CxTreeHeldGet(CxtTreeObject *self, PyObject *args)
{
    PyObject *rVal
#ifdef CxmCcSilence
	= NULL
#endif
	;
    int held;
    CxtEdgeObject *bisect, *reconnectA, *reconnectB;
    CxtTreeMpData *data;

    if (CxpTreeMpDataGet(self, &data))
    {
	rVal = NULL;
	goto RETURN;
    }

    if (PyArg_ParseTuple(args, "i", &held) == 0)
    {
	rVal = NULL;
	goto RETURN;
    }
    if (held >= data->nHeld)
    {
	Py_INCREF(PyExc_ValueError);
	rVal = PyExc_ValueError;
	goto RETURN;
    }

    if (CxTreeTbrNeighborGet(self, data->held[held].neighbor,
			     &bisect, &reconnectA, &reconnectB))
    {
	rVal = NULL;
	goto RETURN;
    }

    Py_INCREF(bisect);

    // reconnect[AB] may be NULL, which must be translated to None.
    if (reconnectA != NULL)
    {
	Py_INCREF(reconnectA);
    }
    else
    {
	Py_INCREF(Py_None);
	reconnectA = (CxtEdgeObject *) Py_None;
    }

    if (reconnectB != NULL)
    {
	Py_INCREF(reconnectB);
    }
    else
    {
	Py_INCREF(Py_None);
	reconnectB = (CxtEdgeObject *) Py_None;
    }

    rVal = Py_BuildValue("(i(OOO))",
			 data->held[held].score,
			 (PyObject *) bisect,
			 (PyObject *) reconnectA,
			 (PyObject *) reconnectB);

    RETURN:
    return rVal;
}
