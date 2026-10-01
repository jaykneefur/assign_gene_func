def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    n, m = len(seq1), len(seq2)

    # Score table: (n+1) rows, (m+1) columns, all zeros to start
    D = [[0.0] * (m + 1) for _ in range(n + 1)]

    # First column: seq1 letters aligned to gaps
    for i in range(1, n + 1):
        D[i][0] = D[i-1][0] + scoring_function(seq1[i-1], "-")

    # First row: seq2 letters aligned to gaps
    for j in range(1, m + 1):
        D[0][j] = D[0][j-1] + scoring_function("-", seq2[j-1])

    # TODO: fill the rest of the table using the slide's formula
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            D[i][j] = max(
                D[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1]),
                D[i-1][j]   + scoring_function(seq1[i-1], "-"),
                D[i][j-1]   + scoring_function("-", seq2[j-1]),
            )

    return D  # temporary, just so we can look at it

def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    raise NotImplementedError()


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)
