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


def main():
    num_char = 2180
    str1, str2 = generate_random_string(num_char)
    #print(str1, str2)

    print('Hamming distance: ', get_hamming_distance(str1, str2))
    print(f'Same chars: {get_same_chars(str1, str2)}, in %: {get_same_chars(str1, str2) / num_char * 100}%')
    # cca 25 % stejnych znaku

    start = time.time()
    print('Edit distance (DP): ', get_edit_distance_dp(str1, str2))
    end = time.time()
    print('Time: ', end - start)
    # pro num_char >= 2180 je cas > 1 sekunda

    start = time.time()
    print('Edit distance: ', get_edit_distance(str1, str2))
    end = time.time()
    print('Time: ', end - start)
    # pro num_char > 12 je cas > 1 sekunda, rekurze je zabijak
    # pro simpanze a cloveka by tabulka mela 400 000 000x400 000 000


if __name__ == '__main__':
    main()
