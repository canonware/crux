%top{
#include <stdlib.h>

#include "../Cx.h"

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
        float dist;
        char *label;
    } tok;
} CxtDistMatrixLexerExtra;
}
%{
#include <errno.h>

#define YY_INPUT(buf, result, max_size) do {				\
    switch (yyextra->inputMode) {					\
	case CxeDistMatrixLexerInputModeFd: {				\
	    ssize_t r = read(yyextra->input.fd, buf, max_size);		\
	    if (r == 0) {						\
		result = YY_NULL;					\
		break;							\
	    } else if (r == -1) {					\
		yyextra->ioerror = 1;					\
		result = YY_NULL;					\
		break;							\
	    }								\
	    result = r;							\
	    break;							\
	} case CxeDistMatrixLexerInputModeStr: {			\
	    if (yyextra->input.s.len == 0) {				\
		result = YY_NULL;					\
		break;							\
	    }								\
	    result = (yyextra->input.s.len < max_size)			\
	      ? yyextra->input.s.len : max_size;			\
	    memcpy(buf, yyextra->input.s.s, result);			\
	    yyextra->input.s.s += result;				\
	    yyextra->input.s.len -= result;				\
	    break;							\
	} default: {							\
	    CxmError("Not reached");					\
	}								\
    }									\
} while (0)
%}

%option prefix="CxDistMatrixLexer_"
%option reentrant
%option extra-type="CxtDistMatrixLexerExtra *"
%option nointeractive
%option nounput

%%

[1-9][0-9]* {
    char *endptr;

    errno = 0;
    yyextra->tok.dist = strtof(yytext, &endptr);
    if (endptr == yytext || errno != 0) {
        return -1;
    }
    yyextra->tokType = CxeDistMatrixLexerTokTypeInt;
    return 1;
}

[-+]?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)? {
    char *endptr;

    errno = 0;
    yyextra->tok.dist = strtof(yytext, &endptr);
    if (endptr == yytext || errno != 0) {
        return -1;
    }
    yyextra->tokType = CxeDistMatrixLexerTokTypeDist;
    return 1;
}

[A-Za-z_][A-Za-z0-9_\-.]*([ \t\r\f\v]+[A-Za-z_][A-Za-z0-9_\-.]*)* {
    yyextra->tok.label = yytext;
    yyextra->tokType = CxeDistMatrixLexerTokTypeLabel;
    return 1;
}

[ \t\r\f\v]+ {}
[\n]+ {
    yyextra->line += yyleng;
}

. {
    return -1;
}

<<EOF>> {
    return 0;
}

%%

int
CxDistMatrixLexer_wrap(yyscan_t yyscanner) {
    return 1;
}
