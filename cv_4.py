import time


def build_suffix_array(string: str):
    suffixes = [(string[i:], i) for i in range(len(string))]
    suffixes.sort(key=lambda x: x[0])
    suffix_array = [suffix[1] for suffix in suffixes]

    return suffix_array


def binary_search_on_suffix_array(string: str, suffix_array: str, pattern: str):
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


def main():
    with open('chr17.fa', 'r', encoding='utf-8') as fna:
        genome_sequence = ''.join(fna.read().splitlines())

    genome_sequence.upper()

    start_time = time.time()
    suffix_array = build_suffix_array(genome_sequence)
    end_time = time.time()
    print(f"Time taken to build suffix array: {end_time - start_time:.2f} seconds")

    with open("fastq_R1.fastq.txt", "r") as file:
        next(file)
        sequence = file.readline().strip().upper()

    start_time = time.time()
    index = binary_search_on_suffix_array(genome_sequence, suffix_array, sequence)
    print(f"Time to find sequence in suffix array: {end_time - start_time:.2f} seconds")


if __name__ == '__main__':
    main()
