from dataclasses import dataclass
import random


@dataclass
class SubsequenceOccurrence:
    subsequence: str
    location: tuple[int, int]
    neighbourhood_parent: str


@dataclass
class ExpandedSubsequence:
    best_extended_subseq: str
    score: int
    query_seq: str
    query_subsequence: str
    db_subsequence: str
    db_sequence: str
    location_db: tuple[int, int]
    location_query: tuple[int, int]


BLOSUM62 = {
    'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'B': -2, 'J': -1, 'Z': -1, 'X': -1, '*': -4},
    'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3, 'B': -1, 'J': -2, 'Z': 0, 'X': -1, '*': -4},
    'N': {'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3, 'B': 4, 'J': -3, 'Z': 0, 'X': -1, '*': -4},
    'D': {'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'J': -3, 'Z': 1, 'X': -1, '*': -4},
    'C': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1, 'B': -3, 'J': -1, 'Z': -3, 'X': -1, '*': -4},
    'Q': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2, 'B': 0, 'J': -2, 'Z': 4, 'X': -1, '*': -4},
    'E': {'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'J': -3, 'Z': 4, 'X': -1, '*': -4},
    'G': {'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3, 'B': -4, 'J': -4, 'Z': -2, 'X': -1, '*': -4},
    'H': {'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3, 'B': 0, 'J': -3, 'Z': 1, 'X': -1, '*': -4},
    'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3, 'B': -3, 'J': 3, 'Z': -3, 'X': -1, '*': -4},
    'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1, 'B': -4, 'J': 1, 'Z': -3, 'X': -1, '*': -4},
    'K': {'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'J': -1, 'Z': 1, 'X': -1, '*': -4},
    'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1, 'B': -3, 'J': 2, 'Z': -1, 'X': -1, '*': -4},
    'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1, 'B': -3, 'J': 0, 'Z': -3, 'X': -1, '*': -4},
    'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2, 'B': -3, 'J': -3, 'Z': -1, 'X': -2, '*': -4},
    'S': {'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2, 'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'J': -2, 'Z': 0, 'X': 0, '*': -4},
    'T': {'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0, 'B': -1, 'J': -1, 'Z': -1, 'X': 0, '*': -4},
    'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3, 'B': -3, 'J': -3, 'Z': -2, 'X': -2, '*': -4},
    'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1, 'B': -2, 'J': -1, 'Z': -2, 'X': -1, '*': -4},
    'V': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -1, 'V': -4, 'B': -1, 'J': -1, 'Z': -1, 'X': -4, '*': -1},
    'B': {'A': -2, 'R': 0, 'N': 0, 'D': 1, 'C': 4, 'Q': 4, 'E': 4, 'G': -1, 'H': 0, 'I': -3, 'L': -3, 'K': -3, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 0, 'J': -3, 'Z': -1, 'X': -4, '*': 1},
    'J': {'A': -3, 'R': -3, 'N': -3, 'D': -2, 'C': 0, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': 3, 'L': 3, 'K': 3, 'M': 3, 'F': 3, 'P': 2, 'S': 2, 'T': 2, 'W': 1, 'Y': -1, 'V': 0, 'B': -2, 'J': 2, 'Z': -2, 'X': -4, '*': -3},
    'Z': {'A': -3, 'R': -3, 'N': -3, 'D': -1, 'C': 0, 'Q': 1, 'E': 1, 'G': -2, 'H': 0, 'I': 3, 'L': 3, 'K': 3, 'M': 3, 'F': 3, 'P': 2, 'S': 2, 'T': 2, 'W': 1, 'Y': -1, 'V': 0, 'B': -2, 'J': 2, 'Z': 2, 'X': -4, '*': -3},
    'X': {'A': -1, 'R': -1, 'N': -1, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -1, 'H': -1, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -1, 'P': -1, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': -1, 'B': -1, 'J': -1, 'Z': -1, 'X': -1, '*': -1},
    '': {'A': -4, 'R': -4, 'N': -4, 'D': -4, 'C': -4, 'Q': -4, 'E': -4, 'G': -4, 'H': -4, 'I': -4, 'L': -4, 'K': -4, 'M': -4, 'F': -4, 'P': -4, 'S': -4, 'T': -4, 'W': -4, 'Y': -4, 'V': -4, 'B': -4, 'J': -4, 'Z': -4, 'X': -4, '': 1},
}

BLOSUM62_SYMBOLS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'J']


def generate_db_string(length: int) -> str:
    return ''.join(random.choice(BLOSUM62_SYMBOLS) for _ in range(length))


def generate_subsequences(seq: str, w: int) -> list[str]:
    return [seq[p:p+w] for p in range(len(seq) - w + 1)]


def test_generate_subsequences() -> None:
    seq = 'YANCLEHKMGS'
    w = 3
    subseqs = generate_subsequences(seq, w)
    print(subseqs)


def get_sequence_to_sequence_score(seq1: str, seq2: str) -> int:
    assert len(seq1) == len(seq2), "Sequences must be of equal length"
    return sum(get_blosum62_score(seq1[i], seq2[i]) for i in range(len(seq1)))


def get_blosum62_score(a: str, b: str) -> int:
    return BLOSUM62[a][b]


def generate_neighborhood_words_recurs(
    seq: str, neighb_seq: str, w: int, i: int, res_list: list[str], threshold: int
) -> None:
    if i < w:
        for sym in BLOSUM62_SYMBOLS:
            generate_neighborhood_words_recurs(seq, neighb_seq + sym, w, i + 1, res_list, threshold)
    else:
        if get_sequence_to_sequence_score(seq, neighb_seq) >= threshold:
            res_list.append(neighb_seq)


def generate_neighborhood_words(seq: str, window: int, threshold: int) -> list[str]:
    neighb = []
    generate_neighborhood_words_recurs(seq, '', window, 0, neighb, threshold)
    return neighb


def find_subsequence_in_database_str(
    subsequence_neighb_words: dict[str, list[str]], database_seq: str, window: int
) -> list[SubsequenceOccurrence]:
    occurrences = []
    for i in range(len(database_seq) - window + 1):
        subseq = database_seq[i:i+window]
        for orig_w, neigh_words in subsequence_neighb_words.items():
            if subseq in neigh_words:
                occurrences.append(SubsequenceOccurrence(subseq, (i, i + window - 1), orig_w))
    return occurrences


def expand_subsequences(
    occurrences: list[SubsequenceOccurrence], query_seq: str, database_seq: str, window: int, drop_off_score: int
) -> list[ExpandedSubsequence]:
    scoring_segment_pairs = []

    for occurrence in occurrences:
        db_start, db_end = occurrence.location
        q_start = query_seq.find(occurrence.neighbourhood_parent)
        q_end = q_start + window - 1

        # Init
        extended_subseq, best_extended_subseq = occurrence.subsequence, occurrence.subsequence
        current_score = best_score = get_sequence_to_sequence_score(extended_subseq, extended_subseq)

        # Right
        for i in range(1, min(len(query_seq)-q_end, len(database_seq)-db_end)):
            symbol_right_query = query_seq[q_end+i]
            symbol_right_database = database_seq[db_end+i]
            change_in_score = get_blosum62_score(symbol_right_query, symbol_right_database)
            if change_in_score < (drop_off_score)*-1:
                break
            extended_subseq = extended_subseq + symbol_right_database
            current_score += change_in_score
            if current_score > best_score:
                best_extended_subseq = extended_subseq
                best_score = current_score

        # Reset
        extended_subseq = best_extended_subseq
        current_score = best_score

        # Left
        for i in range(1, min(q_start, db_start)+1):
            symbol_left_query = query_seq[q_start-i]
            symbol_left_database = database_seq[db_start-i]
            change_in_score = get_blosum62_score(symbol_left_query, symbol_left_database)
            if change_in_score < (drop_off_score)*-1:
                break
            extended_subseq = symbol_left_database + extended_subseq
            current_score += change_in_score
            if current_score > best_score:
                best_extended_subseq = extended_subseq
                best_score = current_score

        scoring_segment_pairs.append(ExpandedSubsequence(
            best_extended_subseq, best_score, query_seq,
            occurrence.neighbourhood_parent, occurrence.subsequence,
            database_seq, (db_start, db_end), (q_start, q_end)
        ))

    return scoring_segment_pairs


def generate_db_with_sequence(sequence: str):
    length = len(sequence) * 20
    generated_db = ''.join(random.choice(BLOSUM62_SYMBOLS) for _ in range(length))
    random_insertion_pos = random.randint(0, len(generated_db) - len(sequence))
    generated_db = generated_db[:random_insertion_pos] + sequence + generated_db[random_insertion_pos:]
    print(f"Random insert location: {random_insertion_pos}")
    return generated_db


def main():
    # Priklad z prednasky
    query_sequence = "YANCLEHKMGS"
    subsequences = generate_subsequences(query_sequence, 3)

    subsequences_neigbhourhood_words = {}
    for s in subsequences:
        subsequences_neigbhourhood_words[s] = generate_neighborhood_words(s, 3, 11)
    print("\nSubsequences neighborhoods")
    for k, v in subsequences_neigbhourhood_words.items():
        print(k, v, sep='\n')
    database_sequence = "DAPCQEHKRGWPNDC"
    occurences = find_subsequence_in_database_str(subsequences_neigbhourhood_words, database_sequence, 3)
    print("\nOccurences")
    for o in occurences:
        print(o)
    best = expand_subsequences(occurences, query_sequence, database_sequence, 3, 2)
    print("\nResult")
    for b in best:
        print(b)

    print()

    # Random priklad
    query_sequence = "AHPWYJ"
    database_sequence = generate_db_with_sequence(query_sequence)
    subsequences = generate_subsequences(query_sequence, 3)
    subsequences_neigbhourhood_words = {}
    for s in subsequences:
        subsequences_neigbhourhood_words[s] = generate_neighborhood_words(s, 3, 11)
    print("\nSubsequences neighborhoods")
    for k, v in subsequences_neigbhourhood_words.items():
        print(k, v, sep='\n')
    occurences = find_subsequence_in_database_str(subsequences_neigbhourhood_words, database_sequence, 3)
    print("\nOccurences")
    for o in occurences:
        print(o)
    best = expand_subsequences(occurences, query_sequence, database_sequence, 3, 2)
    print("\nResult")
    for b in best:
        print(b)


if __name__ == "__main__":
    main()
