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

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            D[i][j] = max(
                D[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1]),
                D[i-1][j]   + scoring_function(seq1[i-1], "-"),
                D[i][j-1]   + scoring_function("-", seq2[j-1]),
            )

    # Traceback from bottom-right
    i, j = n, m
    aligned1 = ""
    aligned2 = ""

    while i > 0 or j > 0:
        # Diagonal: letter matched/mismatched with letter
        if i > 0 and j > 0 and D[i][j] == D[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1]):
            aligned1 = seq1[i-1] + aligned1
            aligned2 = seq2[j-1] + aligned2
            i -= 1
            j -= 1

        # Up: seq1 letter aligned to a gap
        elif i > 0 and D[i][j] == D[i-1][j] + scoring_function(seq1[i-1], "-"):
            aligned1 = seq1[i-1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1

        # Left: gap aligned to seq2 letter
        else:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1

    return aligned1, aligned2, float(D[n][m])

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
    n, m = len(seq1), len(seq2)

    # Score table: (n+1) rows, (m+1) columns, all zeros to start
    D = [[0.0] * (m + 1) for _ in range(n + 1)]

    best_score, best_i, best_j = 0, 0, 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            D[i][j] = max(
                D[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1]),
                D[i-1][j]   + scoring_function(seq1[i-1], "-"),
                D[i][j-1]   + scoring_function("-", seq2[j-1]),
                0
            )
            if D[i][j] > best_score:
                best_score, best_i, best_j = D[i][j], i, j

    # # Traceback from the highest-scoring cell.
    i, j = best_i, best_j    
    aligned1 = ""
    aligned2 = ""

    while i > 0 and j > 0 and D[i][j] > 0:
        # Diagonal: letter matched/mismatched with letter
        if i > 0 and j > 0 and D[i][j] == D[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1]):
            aligned1 = seq1[i-1] + aligned1
            aligned2 = seq2[j-1] + aligned2
            i -= 1
            j -= 1

        # Up: seq1 letter aligned to a gap
        elif i > 0 and D[i][j] == D[i-1][j] + scoring_function(seq1[i-1], "-"):
            aligned1 = seq1[i-1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1

        # Left: gap aligned to seq2 letter
        else:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1

    return aligned1, aligned2, float(best_score)


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)

from Bio.Align import substitution_matrices

BLOSUM62 = substitution_matrices.load("BLOSUM62")
GAP_PENALTY = -4  # at least as bad as BLOSUM62's worst mismatch

# Convert to a plain dict once, for speed
_blosum_dict = {(a, b): BLOSUM62[a][b] for a in BLOSUM62.alphabet for b in BLOSUM62.alphabet}

def scoring_function_blosum62(aa_i, aa_j):
    if aa_i == "-" or aa_j == "-":
        return GAP_PENALTY
    return _blosum_dict[(aa_i, aa_j)]