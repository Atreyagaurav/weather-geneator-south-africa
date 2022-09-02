import pandas as pd


df = pd.read_csv("csvs/weather-regimes-scott-6.csv")

marcov_len = 1
total_patterns = range(1, 9)
patterns = dict()
clusters = df.clusters.to_numpy(dtype=int)

for i in range(len(clusters) - marcov_len):
    prev = tuple(clusters[i: i+marcov_len+1])
    if prev in patterns:
        patterns[prev] += 1
    else:
        patterns[prev] = 1

sorted_patterns = sorted(patterns.items(), key=lambda x: x[1], reverse=True)

sorted_patterns[:10]

sequence = dict()

for pat, count in sorted_patterns:
    start = pat[:-1]
    if start in sequence:
        sequence[start][pat[-1]] = count
    else:
        sequence[start] = {pat[-1]: count}

list(sequence.items())[-1]

for k, v in sequence.items():
    total = sum(v.values())
    print(f"IF {k} THEN: ")
    for k2, v2 in v.items():
        print(f"p({k2})={v2/total:.3f}", end="; ")
    print()

prob_table = pd.DataFrame(sequence)
prob_table.sort_index(inplace=True)
prob_table.sort_index(axis=1, inplace=True)

# cols is previous, row is next
prob_table.columns = (f"F{c}" for c in prob_table.columns)

prob_table.transpose().to_csv(f"csvs/freq-table-scott-6-{marcov_len}.csv")

for i, total in prob_table.sum().iteritems():
    prob_table.loc[:, i] /= total

prob_table.transpose().to_csv(f"csvs/prob-table-scott-6-{marcov_len}.csv")
