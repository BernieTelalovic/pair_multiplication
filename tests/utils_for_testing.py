from itertools import combinations_with_replacement

def partitions(n):
    result = []
    # Generate partitions by trying all lengths from 1 up to n
    for length in range(1, n + 1):
        for combo in combinations_with_replacement(range(1, n + 1), length):
            if sum(combo) == n:
                result.append(list(combo))
    return result

# Filter for only distinct partitions (unique combinations of numbers)
def unique_partitions(n):
    all_parts = partitions(n)
    unique_parts = {tuple(sorted(part, reverse=True)) for part in all_parts}
    return [list(part) for part in sorted(unique_parts, reverse=True)]


