%top{
#include <stdlib.h>

#include "Cx.h"

typedef enum {
    CxeNewickLexerInputModeFd,
    CxeNewickLexerInputModeStr
} CxtNewickLexerInputMode;

typedef enum {
    CxeNewickLexerTokTypeComment,
    CxeNewickLexerTokTypeLparen,
    CxeNewickLexerTokTypeRparen,
    CxeNewickLexerTokTypeComma,
    CxeNewickLexerTokTypeColon,
    CxeNewickLexerTokTypeSemicolon,
    CxeNewickLexerTokTypeBranchLength,
    CxeNewickLexerTokTypeUnquotedLabel,
    CxeNewickLexerTokTypeQuotedLabel,
    CxeNewickLexerTokTypeWhitespace
} CxtNewickLexerTokType;

typedef struct {
    unsigned line;
    int ioerror;
    CxtNewickLexerInputMode inputMode;
    union {
        int fd;
        struct {
            char *s;
            size_t len;
        } s;
    } input;

    CxtNewickLexerTokType tokType;
    struct {
	char *s;
	size_t len;
    } tok;
} CxtNewickLexerExtra;
}
%{
#include <errno.h>

#define YY_INPUT(buf, result, max_size) do {				\
    switch (yyextra->inputMode) {					\
	case CxeNewickLexerInputModeFd: {				\
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
	} case CxeNewickLexerInputModeStr: {			\
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

%option prefix="CxNewickLexer_"
%option reentrant
%option extra-type="CxtNewickLexerExtra *"
%option nointeractive
%option noinput
%option nounput

%option stack
%x comment

%%

\[[^[\]\n]*\] { // [...]
    yyextra->tok.s = yytext;
    yyextra->tok.len = yyleng;
    yyextra->tokType = CxeNewickLexerTokTypeComment;
    return 1;
}

\[[^[\]\n]*\[ { // [...[
    yy_push_state(comment, yyscanner);
    yy_push_state(comment, yyscanner);
    yymore();
}

<comment>\[[^[\]\n]*\n { // [...\n
    yyextra->line++;
    yy_push_state(comment, yyscanner);
    yymore();
}

<comment>\[ { // [
    yy_push_state(comment, yyscanner);
    yymore();
}

<comment>\] { // ]
    yy_pop_state(yyscanner);
    if (yy_top_state(yyscanner) == comment) {
	yymore();
    } else {
	yyextra->tok.s = yytext;
	yyextra->tok.len = yyleng;
	yyextra->tokType = CxeNewickLexerTokTypeComment;
	return 1;
    }
}

<comment>\n { // \n
    yyextra->line++;
    yymore();
}

<comment>[^\[\]\n]+ { // text
    yymore();
}

[(] {
    errno = 0;
    yyextra->tokType = CxeNewickLexerTokTypeLparen;
    return 1;
}

[)] {
    errno = 0;
    yyextra->tokType = CxeNewickLexerTokTypeRparen;
    return 1;
}

[,] {
    errno = 0;
    yyextra->tokType = CxeNewickLexerTokTypeComma;
    return 1;
}

[:] {
    errno = 0;
    yyextra->tokType = CxeNewickLexerTokTypeColon;
    return 1;
}

[;] {
    errno = 0;
    yyextra->tokType = CxeNewickLexerTokTypeSemicolon;
    return 1;
}

[-+]?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)? {
    errno = 0;
    yyextra->tok.s = yytext;
    yyextra->tok.len = yyleng;
    yyextra->tokType = CxeNewickLexerTokTypeBranchLength;
    return 1;
}

[^ \t\n\r\f\v()[\]':;,]+ {
    yyextra->tok.s = yytext;
    yyextra->tok.len = yyleng;
    yyextra->tokType = CxeNewickLexerTokTypeUnquotedLabel;
    return 1;
}

'(''|[^'])*' {
    yyextra->tok.s = yytext;
    yyextra->tok.len = yyleng;
    yyextra->tokType = CxeNewickLexerTokTypeQuotedLabel;
    return 1;
}

[ \t\r\f\v]+ {
    yyextra->tok.s = yytext;
    yyextra->tok.len = yyleng;
    yyextra->tokType = CxeNewickLexerTokTypeWhitespace;
    return 1;
}
[\n]+ {
    yyextra->line += yyleng;
    yyextra->tok.s = yytext;
    yyextra->tok.len = yyleng;
    yyextra->tokType = CxeNewickLexerTokTypeWhitespace;
    return 1;
}

. {
    return -1;
}

<<EOF>> {
    return 0;
}

%%

int
CxNewickLexer_wrap(yyscan_t yyscanner) {
    return 1;
}
