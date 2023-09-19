complement_mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A'}
map_genes_to_proteins = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': 'STOP', 'UAG': 'STOP',
    'UGU': 'C', 'UGC': 'C', 'UGA': 'STOP', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def check_if_rna(arr: str) -> bool:
    return 'U' in arr


def translate_gen_to_protein(arr: str) -> list[str]:
    return ''.join(
        [map_genes_to_proteins[arr[i:i+3]] for i in range(0, len(arr), 3)] if len(arr) % 3 == 0 else ''
    ).split('STOP')


def find_gens(arr: str) -> list[str]:
    return [x.split('UAA')[0].split('UAG')[0].split('UGA')[0] for x in arr.split('AUG') if len(x) > 0]


def convert_seq_to_complement(arr: str) -> str:
    return ''.join([complement_mapping[x] for x in arr])


def main():
    with open('covid.fna', 'r', encoding='utf-8') as fna:
        lines = ''.join(fna.read().splitlines())

    complemented_lines = convert_seq_to_complement(lines)
    print('Complemented lines: ', complemented_lines)

    print('Is RNA: ', check_if_rna(lines))

    if not check_if_rna(lines):
        lines = lines.replace('T', 'U')

    gens = find_gens(lines)
    gens = [x for x in gens if len(x) > 0]
    #print('Gens', gens)
    print('Gen count: ', len(gens))

    proteins = []
    for gen in gens:
        proteins += translate_gen_to_protein(gen)
    proteins = [x for x in proteins if len(x) > 0]
    #print('Proteins: ', proteins)
    print('Protein count: ', len(proteins))


if __name__ == '__main__':
    main()
