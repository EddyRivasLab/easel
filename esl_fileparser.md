# esl_fileparser - token-based data file input

The `fileparser` module parses simple input text data files that
consist of whitespace-delimited tokens.

Data files can contain blank lines and comments. Comments are defined
by a single character; for instance, a `#` character commonly means
that everything following the `#` on the line is a comment.

Two different styles of token input are supported. The simplest style
reads tokens one at a time, regardless of what line they occur on,
until the file ends. You can also read in a line-oriented way, in
which you get one data line at a time, then read all the tokens on that
line; this style lets you count how many tokens occur on a data line,
which allows better checking of your input.

The module implements one object, an `ESL_FILEPARSER`, that holds the
open input stream and the state of the parser.

## example

The `eslFILEPARSER_EXAMPLE` example in `esl_fileparser.c` opens a
file, reads all its tokens one at a time, and prints out token number,
token length, and the token itself.

A single character can be defined to serve as a comment character
(often `#`), using the `esl_fileparser_SetCommentChar()` call. The
parser will ignore the comment character, and the remainder of any line
following a comment character.

Each call to `esl_fileparser_GetToken()` retrieves one
whitespace-delimited token from the input stream; the call returns
`eslOK` if a token is parsed, and `eslEOF` when there are no more
tokens in the file. Whitespace is defined as space, tab, newline, or
carriage return (`" \t\n\r"`).

When the caller is done, the fileparser is closed with
`esl_fileparser_Close()`.

## a second example: line-oriented parsing

The `esl_fileparser_GetToken()` call provides a simple style of parsing
a file: read one token at a time until the file ends, regardless of
what line the tokens are on. However, you may want to know how many
tokens are on a given data line, either because you know how many there
should be (and you want to verify) or because you don't (and you need
to allocate some variable-size data structure appropriately). The
`eslFILEPARSER_EXAMPLE2` example in `esl_fileparser.c`
reads a file line by line.

The output from this example is, for each data line, the actual line
number (starting from 1), the data line number (a count that excludes
comments and blank lines), and the number of tokens on the line.

Note the use of `efp->linenumber` to obtain the current line in the
file. You can use this to produce informative error messages. If a
token is not what you expected, you probably want to provide some
diagnostic output to the user, and `efp->linenumber` lets you direct
the user to the line that the failure occurred at.
