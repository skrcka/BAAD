import time


def build_suffix_array(string: str) -> list[int]:
    suffixes = [(string[i:], i) for i in range(len(string))]
    suffixes.sort(key=lambda x: x[0])
    suffix_array = [suffix[1] for suffix in suffixes]

    return suffix_array


def construct_bwt(string: str, suffix_array: list[int]) -> str:
    return ''.join([string[i-1] for i in suffix_array])


def binary_search_on_suffix_array(string: str, suffix_array: str, pattern: str) -> int | None:
    l = 0
    r = len(suffix_array) - 1

    while l <= r:
        mid = (l + r) // 2
        suffix = string[suffix_array[mid]:]

        if pattern == suffix[:len(pattern)]:
            return suffix_array[mid]
        elif pattern < suffix:
            r = mid - 1
        else:
            l = mid + 1

    return None


def build_fl_columns(bwt: str) -> (str, str):
    return ''.join(sorted(bwt)), bwt


def compute_rank(bwt: str) -> dict[str, list[int]]:
    rank = {}
    char_count = {}

    for char in bwt:
        char_count[char] = char_count.get(char, 0) + 1
        if char not in rank:
            rank[char] = [0]
        rank[char].append(char_count[char])

    return rank


def fm_index_search(bwt: str, rank: dict[str, list[int]], pattern: str) -> int:
    char_to_column = {char: idx for idx, char in enumerate(sorted(set(bwt)))}
    top, bottom = 0, len(bwt) - 1

    for char in reversed(pattern):
        if char in rank:
            top = char_to_column[char] + rank[char][top]
            bottom = char_to_column[char] + rank[char][bottom + 1] - 1
            if top > bottom:
                return None
        else:
            return None

    return top


def main():
    #with open('chr17.fa', 'r', encoding='utf-8') as fna:
    #    genome_sequence = ''.join(fna.read().splitlines())
    genome_sequence = 'abaaba$'

    genome_sequence = genome_sequence.upper()

    suffix_array = build_suffix_array(genome_sequence)
    bwt = construct_bwt(genome_sequence, suffix_array)
    rank = compute_rank(bwt)

    #with open("fastq_R1.fastq.txt", "r") as file:
    #    next(file)
    #    sequence = file.readline().strip().upper()
    sequence = 'aba'

    start_time = time.time()
    position = fm_index_search(bwt, rank, sequence)
    end_time = time.time()
    print(f"Time to find sequence using FM-index: {end_time - start_time:.8f} seconds")



if __name__ == '__main__':
    main()
