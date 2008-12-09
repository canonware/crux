#ifndef CxDistMatrix_h
#define CxDistMatrix_h

#include "../Cx.h"

typedef float CxtDMDist;
typedef unsigned long CxtDMSize;

typedef enum {
    CxeDistMatrixLexerInputModeFd,
    CxeDistMatrixLexerInputModeStr
} CxtDistMatrixLexerInputMode;

typedef enum {
    CxeDistMatrixLexerTokTypeInt,
    CxeDistMatrixLexerTokTypeDist,
    CxeDistMatrixLexerTokTypeLabel
} CxtDistMatrixLexerTokType;

typedef struct {
    unsigned line;
    int ioerror;
    CxtDistMatrixLexerInputMode inputMode;
    union {
        int fd;
        struct {
            char *s;
            size_t len;
        } s;
    } input;

    CxtDistMatrixLexerTokType tokType;
    union {
        CxtDMDist dist;
        char *label;
    } tok;
} CxtDistMatrixLexerExtra;

#ifndef CxmUseInlines
CxtDMSize
CxDistMatrixNxy2i(CxtDMSize n, CxtDMSize x, CxtDMSize y);
#endif

#if (defined(CxmUseInlines) || defined(CxDistMatrix_c))
CxmInline CxtDMSize
CxDistMatrixNxy2i(CxtDMSize n, CxtDMSize x, CxtDMSize y) {
    CxmAssert(x < n);
    CxmAssert(y < n);
    CxmAssert(x != y);

    if (x > y) {
	CxtDMSize t;

	t = x;
	x = y;
	y = t;
    }

    return n*x + y - (((x+3)*x) >> 1) - 1;
}
#endif

#endif // CxDistMatrix_h
