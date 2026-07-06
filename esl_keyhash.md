# esl_keyhash - associative hashes

The `keyhash` module provides a semblance of associative arrays (for
example, Perl hashes), by associating keywords with an integer array
index, and storing the association in an internal hash table for rapid
access.

The module implements one object: the `ESL_KEYHASH`.

The application maintains data in normal C-style arrays that
are indexed by an integer index value, and it uses the keyhash to
associate a specific key with that integer index. To store info, you
first store the keyword and obtain a new index value (this simply
starts at 0 and counts up, as you store successive keys), then you
store the info your arrays at that index. To look up info, you look up
the keyword to obtain the index, then you access the info by indexing
into your arrays.

This is the moral equivalent of Perl's associative arrays, as in
`$foo{$key} = whatever; $bar{$key} = whatever`.

The `eslKEYHASH_EXAMPLE` example at the end of `esl_keyhash.c` shows a
contrived example of storing the keywords obtained from a list in one
file, then looking up keywords listed in a second file. It doesn't
demonstrate the idea of using the index to store and retrieve
additional info associated with the keyword, but it demonstrates the
essentials of the `keyhash` API.

