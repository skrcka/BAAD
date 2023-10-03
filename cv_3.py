import random
import sys
import time


sys.setrecursionlimit(1500)


def generate_random_string(length: int) -> str:
    return (''.join([random.choice('acgt') for _ in range(length)]) for _ in range(2))


def get_hamming_distance(str1: str, str2: str) -> int:
    return sum([1 for i in range(len(str1)) if str1[i] != str2[i]])


def get_same_chars(str1: str, str2: str) -> int:
    return sum([1 for i in range(len(str1)) if str1[i] == str2[i]])


def get_edit_distance(str1: str, str2: str) -> int:
    if len(str1) == 0:
        return len(str2)
    if len(str2) == 0:
        return len(str1)

    if str1[-1] == str2[-1]:
        return get_edit_distance(str1[:-1], str2[:-1])
    else:
        return min(
            get_edit_distance(str1[:-1], str2[:-1]),
            get_edit_distance(str1[:-1], str2),
            get_edit_distance(str1, str2[:-1])
        ) + 1


def get_edit_distance_dp(str1: str, str2: str) -> int:
    matrix = [[0 for _ in range(len(str2) + 1)] for _ in range(len(str1) + 1)]

    for i in range(len(str1) + 1):
        matrix[i][0] = i
    for j in range(len(str2) + 1):
        matrix[0][j] = j

    for i in range(1, len(str1) + 1):
        for j in range(1, len(str2) + 1):

            if str1[i - 1] == str2[j - 1]:
                matrix[i][j] = matrix[i-1][j-1]

            else:
                matrix[i][j] = min(
                    matrix[i-1][j-1],
                    matrix[i-1][j],
                    matrix[i][j-1]
                ) + 1

    return matrix[len(str1)][len(str2)]


def global_alignment(str1: str, str2: str) -> (int, str):
    rows, cols = len(str1) + 1, len(str2) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    backtrack = [['' for _ in range(cols)] for _ in range(rows)]

    for i in range(1, rows):
        matrix[i][0] = matrix[i-1][0] - 2
        backtrack[i][0] = 'U'

    for j in range(1, cols):
        matrix[0][j] = matrix[0][j-1] - 2
        backtrack[0][j] = 'L'

    for i in range(1, rows):
        for j in range(1, cols):
            match = matrix[i-1][j-1] + (1 if str1[i-1] == str2[j-1] else -1)
            delete = matrix[i-1][j] - 2
            insert = matrix[i][j-1] - 2

            matrix[i][j] = max(match, delete, insert)
            if matrix[i][j] == match:
                backtrack[i][j] = 'D'
            elif matrix[i][j] == delete:
                backtrack[i][j] = 'U'
            else:
                backtrack[i][j] = 'L'

    i, j = len(str1), len(str2)
    cigar = []
    while i > 0 or j > 0:
        if backtrack[i][j] == 'D':
            if str1[i-1] == str2[j-1]:
                cigar.append('M')
            else:
                cigar.append('X')
            i -= 1
            j -= 1
        elif backtrack[i][j] == 'U':
            cigar.append('D')
            i -= 1
        else:
            cigar.append('I')
            j -= 1

    cigar = ''.join(cigar[::-1])
    return matrix[len(str1)][len(str2)], cigar


def local_alignment(str1: str, str2: str) -> (int, str):
    rows, cols = len(str1) + 1, len(str2) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    backtrack = [['' for _ in range(cols)] for _ in range(rows)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):
            match = matrix[i-1][j-1] + (1 if str1[i-1] == str2[j-1] else -1)
            delete = matrix[i-1][j] - 2
            insert = matrix[i][j-1] - 2
            matrix[i][j] = max(0, match, delete, insert)

            if matrix[i][j] == match:
                backtrack[i][j] = 'D'
            elif matrix[i][j] == delete:
                backtrack[i][j] = 'U'
            elif matrix[i][j] == insert:
                backtrack[i][j] = 'L'

            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)

    i, j = max_pos
    cigar = []
    while matrix[i][j] != 0:
        if backtrack[i][j] == 'D':
            if str1[i-1] == str2[j-1]:
                cigar.append('M')
            else:
                cigar.append('X')
            i -= 1
            j -= 1
        elif backtrack[i][j] == 'U':
            cigar.append('D')
            i -= 1
        else:
            cigar.append('I')
            j -= 1

    cigar = ''.join(cigar[::-1])
    return max_score, cigar


def semi_global_alignment(str1: str, str2: str) -> (int, str):
    rows, cols = len(str1) + 1, len(str2) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    backtrack = [['' for _ in range(cols)] for _ in range(rows)]

    for i in range(1, rows):
        matrix[i][0] = 0
        backtrack[i][0] = 'U'

    for j in range(1, cols):
        matrix[0][j] = matrix[0][j-1] - 2
        backtrack[0][j] = 'L'

    for i in range(1, rows):
        for j in range(1, cols):
            match = matrix[i-1][j-1] + (1 if str1[i-1] == str2[j-1] else -1)
            delete = matrix[i-1][j] - 2
            insert = matrix[i][j-1] - 2
            matrix[i][j] = max(match, delete, insert)

            if matrix[i][j] == match:
                backtrack[i][j] = 'D'
            elif matrix[i][j] == delete:
                backtrack[i][j] = 'U'
            else:
                backtrack[i][j] = 'L'

    max_j = max(range(cols), key=lambda j: matrix[rows-1][j])
    cigar = []
    i, j = rows - 1, max_j
    while i > 0 or j > 0:
        if backtrack[i][j] == 'D':
            if str1[i-1] == str2[j-1]:
                cigar.append('M')
            else:
                cigar.append('X')
            i -= 1
            j -= 1
        elif backtrack[i][j] == 'U':
            cigar.append('D')
            i -= 1
        else:
            cigar.append('I')
            j -= 1

    cigar = ''.join(cigar[::-1])
    return matrix[rows-1][max_j], cigar



def main():
    num_char = 2180
    str1, str2 = generate_random_string(num_char)
    print(f'str1: {str1}')
    print(f'str2: {str2}')

    # Call the Global Alignment function and display its results
    score_global, cigar_global = global_alignment(str1, str2)
    print(f"Global Alignment Score: {score_global}")
    print(f"Global Alignment CIGAR: {cigar_global}")

    # Call the Local Alignment function and display its results
    score_local, cigar_local = local_alignment(str1, str2)
    print(f"\nLocal Alignment Score: {score_local}")
    print(f"Local Alignment CIGAR: {cigar_local}")

    # Call the Semi-Global Alignment function and display its results
    score_semi_global, cigar_semi_global = semi_global_alignment(str1, str2)
    print(f"\nSemi-Global Alignment Score: {score_semi_global}")
    print(f"Semi-Global Alignment CIGAR: {cigar_semi_global}")



if __name__ == '__main__':
    main()
