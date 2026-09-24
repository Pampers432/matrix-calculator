from fractions import Fraction


def transpose(M):
    return [[M[i][j] for i in range(len(M))] for j in range(len(M[0]))]


def submat(M, i, j):
    return [row[:j] + row[j + 1:] for k, row in enumerate(M) if k != i]


def eye(n):
    return [[Fraction(1) if r == c else Fraction(0) for c in range(n)] for r in range(n)]


def is_square(M):
    return len(M) > 0 and len(M) == len(M[0])