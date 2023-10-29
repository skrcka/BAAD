def calculate_bayesian_probability(ref_base, read_base, h=0.001, e_err=0.01):
    if ref_base == read_base:
        return (1 - h) * (1 - e_err)
    else:
        return h * (1 - e_err) + (1 - h) * e_err


def detect_snp(sam_file, ref_seq):
    variants = []
    with open(sam_file, 'r', encoding='utf-8') as file:
        for line in file:
            if not line.startswith('@'):
                fields = line.strip().split()
                read_seq = fields[9]
                pos = int(fields[3])

                for i, base in enumerate(read_seq):
                    ref_base = ref_seq[pos + i - 1]
                    prob = calculate_bayesian_probability(ref_base, base)
                    if prob < 0.5:
                        variants.append((pos+i, ref_base, base))
    return variants


def main():
    with open('chr17.fa', 'r', encoding='utf-8') as fna:
        ref_seq = ''.join(fna.read().splitlines())

    snps = detect_snp("bio.sam", ref_seq)

    with open("bio.vcf", "r", encoding='utf-8') as f:
        for line in f:
            if not line.startswith("#"):
                fields = line.strip().split()
                pos, ref, alt = int(fields[1]), fields[3], fields[4]
                if (pos, ref, alt) in snps:
                    print(f"Shoda varianty na pozici {pos} - {ref} -> {alt}")

    print("Jedná se o bodové záměny (SNP)")


if __name__ == "__main__":
    main()
