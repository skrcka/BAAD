import time


def build_suffix_array(string: str):
    """
    Build a suffix array for the given string.
    This function will return the suffix array which is essentially a list of starting indices of sorted suffixes.
    """
    # Creating a list of suffixes
    suffixes = [(string[i:], i) for i in range(len(string))]

    # Sorting the list of suffixes
    suffixes.sort(key=lambda x: x[0])

    # Extracting only the indices of the sorted suffixes
    suffix_array = [suffix[1] for suffix in suffixes]

    return suffix_array


def binary_search_on_suffix_array(string: str, suffix_array: str, pattern: str):
    """
    Perform a binary search on the suffix array to find the pattern.
    This function returns the starting index of the pattern if found, otherwise None.
    """
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
    # Load the sequence of chromosome 17 (assuming you have it in a plain text file, one continuous sequence)
    with open('covid.fna', 'r', encoding='utf-8') as fna:
        genome_sequence = ''.join(fna.read().splitlines())

    # Measure time taken to build the suffix array
    start_time = time.time()
    suffix_array = build_suffix_array(genome_sequence)
    end_time = time.time()
    print(f"Time taken to build suffix array: {end_time - start_time:.2f} seconds")

    # Assuming you have a FASTQ file (you'll need to parse it properly if it has multiple sequences)
    with open("sample.fastq", "r") as file:
        # Skipping the first line (sequence identifier)
        next(file)
        # Reading the sequence
        sequence = file.readline().strip().upper()

    # Measure time taken to find the sequence in the suffix array
    start_time = time.time()
    index = binary_search_on_suffix_array(genome_sequence, suffix_array, sequence)
    end_time


if __name__ == '__main__':
    main()
