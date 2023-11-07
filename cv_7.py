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


def greedy_scs(sequences: list[str], l: int = 3):
    reads = sequences
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


def get_sequences(sequnce: str, k: int):
    return [sequnce[i:i+k] for i in range(0, len(sequnce) - k + 1)]


def main():
    sequences = ['GTACGT', 'TACGTA']
    print(greedy_scs(sequences, 4))

    sequences = get_sequences('ACCTGACTGACCTGACTGACCTGACTGACCTGACTGACCTGACTG', 6)
    print(greedy_scs(sequences, 4))


if __name__ == '__main__':
    main()
