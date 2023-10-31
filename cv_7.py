import itertools


def get_overlap(a: str, b: str, min_length: int = 3):
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1


def greedy_scs(reads: str, k: int, l: int = 3):
    reads: list[str] = [reads[i:i+k] for i in range(0, len(reads) - k + 1)]
    read_a, read_b = None, None
    while True:
        max_olen = 0
        for a, b in itertools.permutations(reads, 2):
            olen = get_overlap(a, b, l)
            if olen > max_olen:
                read_a, read_b = a, b
                max_olen = olen
        if max_olen == 0:
            break
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[max_olen:])
    return ''.join(reads)


def main():
    sequence = "ATGCGATGACGGT"
    print(greedy_scs(sequence, 4, 3))
    sequence = "ATGCGATGACGGTATGCGATGACGGTATGCGATGACGGTATGCGATGACGGT"
    print(greedy_scs(sequence, 4, 3))
    sequence = "ACGTACTGACCCTGACCTGACCTGACCT"
    print(greedy_scs(sequence, 4, 3))


if __name__ == '__main__':
    main()
