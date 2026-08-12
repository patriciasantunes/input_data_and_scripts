#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
import collections, random, csv
import argparse
import itertools, textwrap


def read_fasta(file: str):
    with open(file, "r") as handle:
        seqs: list[str] = []
        seq_ids: list[str] = []
        for line in handle.read().splitlines():
            if line.startswith(">"):
                seqs.append("")
                seq_ids.append(line[1:].strip())
            else:
                seqs[-1] += line
    return [*zip(seq_ids, seqs)]


def write_tabular(file: str, chosen: dict[str, list[str]], sep: str):
    with open(file, mode="w") as handle:
        handle.writelines(
            f"{seq}{sep}{cluster}\n" for cluster, seqs in chosen.items() for seq in seqs
        )


def get_clusters(file: str):
    clusters = collections.defaultdict(list)
    with open(file, newline="") as handle:
        data = [*csv.DictReader(handle, delimiter=",")]
        for line in data:
            clusters[line["clustername"]].append(line["leafname"])

    if "0" in clusters:
        print(f"Removing cluster 0 with {len(clusters['0'])} sequences")
        del clusters["0"]

    return clusters


def ratio_breakdown(total: int, ratios: list[int], bins: list[int]):
    """Splits total in a list proportional to ratios, while respecting maximum bin sizes.

    Args:
        total (int): Final sum.
        ratios (list[int]): Proportions to use (order irrelevant).
        bins (list[int]): Maximum size for proportions (order irrelevant).

    Returns:
        list[int]: Sorted list of integers that sum to `total`. 
    """
    ratios.sort(reverse=True)
    bins.sort(reverse=True)

    values = [max(round(r * total / sum(ratios)), 1) for r in ratios]

    # add to biggest
    pos = itertools.cycle(range(len(values) - 1, -1, -1))
    while sum(values) < total:
        idx = next(pos)
        values[idx] += 1

    # take from biggest
    pos = itertools.cycle(range(len(values) - 1, -1, -1))
    while sum(values) > total:
        idx = next(pos)
        if values[idx] > 1:
            values[idx] -= 1

    # give from biggest to idx
    pos = itertools.cycle(range(len(values) - 1, -1, -1))
    bins.sort()
    values.sort()
    for idx, (b, v) in enumerate(zip(bins, values)):
        if b < v:
            values[idx] = b
            for _ in range(v-b):
                while (p := next(pos)) <= idx:
                    continue
                values[p] += 1

    values.sort()
    return values


def select_sequences(clusters: dict[str, list[str]], criteria: int, n: int = 0):
    sizes = {cluster: len(seqs) for cluster, seqs in clusters.items()}
    order = sorted(sizes, key=sizes.__getitem__)

    if criteria != "exact" and n < (n_clusters := len(clusters)):
        raise ValueError(
            f"Number of clusters ({n_clusters}) is bigger than number of requested sequences ({n})"
        )
    if criteria != "exact" and n > (n_chosen := sum(sizes.values())):
        raise ValueError(
            f"Total number of sequences ({n_chosen}) is smaller than number of requested sequences ({n})"
        )
    if criteria == "exact" and n > (min_clusters := min(sizes.values())):
        raise ValueError(
            f"One cluster has less sequences ({min_clusters}) than requested ({n})"
        )

    if criteria == "exact":
        # n sequences per cluster
        sample_sizes = [n] * len(clusters)

    elif criteria == "equal":
        # n sequences total, same number per cluster, prioritize bigger
        sample_sizes = ratio_breakdown(n, [1] * len(clusters), [*sizes.values()])

    else:
        # n sequences total, proportional to cluster size, prioritize bigger
        sample_sizes = ratio_breakdown(n, [*sizes.values()], [*sizes.values()])

    chosen = {cl: random.sample(clusters[cl], sz) for cl, sz in zip(order, sample_sizes)}
    assert (n_chosen := sum(len(c) for c in chosen.values())) == n

    print(n_chosen, "sequences chosen")
    return chosen


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-if", "--input-fasta", required=True)
    parser.add_argument("-ic", "--input-csv", required=True)

    parser.add_argument(
        "-of",
        "--output-fasta",
        default="selected.fasta",
        help="Arquivo FASTA (default: selected.fasta)",
    )
    parser.add_argument("-oc", "--output-csv")
    parser.add_argument("-ot", "--output-tsv")
    parser.add_argument(
        "-c",
        "--criteria",
        default="exact", metavar="CRITERIA",
        choices=["exact", "equal", "proportional"],
        help=textwrap.dedent("""
                             Critério de escolha (default: exact)
                                exact           N sequences per cluster
                                equal           N sequences total, same number per cluster
                                proportional    N sequences total, proportional to cluster size"""
        ).lstrip(),
    )
    parser.add_argument("-n", default=1, type=int, help="Número de sequências")

   
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    fasta = read_fasta(args.input_fasta)
    clusters = get_clusters(args.input_csv)

    chosen = select_sequences(clusters, criteria=args.criteria, n=args.n)

    if args.output_tsv:
        write_tabular(args.output_tsv, chosen, sep="\t")
    if args.output_csv:
        write_tabular(args.output_csv, chosen, sep=",")
    with open(args.output_fasta, mode="w") as handle:
        chosen_ids = {seq for seqs in chosen.values() for seq in seqs}
        handle.writelines(
            f">{seqid}\n{seq}\n" for seqid, seq in fasta if seqid in chosen_ids
        )
