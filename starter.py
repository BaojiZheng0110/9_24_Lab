# starter.py — FFT integer multiplication lab (students fill TODOs)
import math
from typing import List

def next_power_of_two(m: int) -> int:
    return 1 << (m - 1).bit_length()

def split_digits(x: str, base: int) -> List[int]:
    if x == "0":
        return [0]
    k = len(str(base)) - 1
    out = []
    i = len(x)
    while i > 0:
        j = max(0, i - k)
        out.append(int(x[j:i]))
        i = j
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out

def join_digits(ds: List[int], base: int) -> str:
    i = len(ds) - 1
    while i > 0 and ds[i] == 0:
        i -= 1
    ds = ds[:i+1]
    k = len(str(base)) - 1
    head = str(ds[-1])
    body = "".join(f"{d:0{k}d}" for d in reversed(ds[:-1]))
    return head + body

# TODO: Implement fft (iterative Cooley–Tukey)
def fft(a: List[complex], invert: bool) -> None:
    n = len(a)
    j = 0
    # Bit-reversal permutation
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i < j:
            a[i], a[j] = a[j], a[i]
    length = 2
    while length <= n:
        ang = 2 * math.pi / length * (-1 if invert else 1)
        wlen = complex(math.cos(ang), math.sin(ang))
        for i in range(0, n, length):
            w = 1
            for jj in range(i, i + length // 2):
                u = a[jj]
                v = a[jj + length // 2] * w
                a[jj] = u + v
                a[jj + length // 2] = u - v
                w *= wlen
        length <<= 1
    if invert:
        inv_n = 1.0 / n
        for i in range(n):
            a[i] *= inv_n

# Convolution using FFT
def conv_fft(a: List[int], b: List[int]) -> List[int]:
    if not a or not b:
        return []
    n = len(a) + len(b) - 1
    size = next_power_of_two(n)
    fa = [complex(x, 0.0) for x in a] + [0j] * (size - len(a))
    fb = [complex(x, 0.0) for x in b] + [0j] * (size - len(b))
    fft(fa, invert=False)
    fft(fb, invert=False)
    for i in range(size):
        fa[i] *= fb[i]
    fft(fa, invert=True)
    res = [int(round(fa[i].real)) for i in range(n)]
    # Remove trailing zeros
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

# TODO: Implement carry propagation and bigmul_fft
def carry_base(c: List[int], base: int) -> List[int]:
    res = c[:]
    for i in range(len(res)):
        if res[i] >= base or res[i] < 0:
            if i + 1 == len(res):
                res.append(0)
            res[i + 1] += res[i] // base
            res[i] %= base
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res  

def bigmul_fft(a: str, b: str, k: int = 3) -> str:
    base = 10 ** k
    da = split_digits(a, base)
    db = split_digits(b, base)
    c = conv_fft(da, db)
    c = carry_base(c, base)
    return join_digits(c, base)

if __name__ == "__main__":

    print("Run pytest to check your implementation.")

# --- FFT Integer Multiplication Lab Report ---
#
# Runtime Analysis:
# FFT-based multiplication (bigmul_fft) uses the Cooley-Tukey FFT algorithm for fast polynomial convolution.
# For two n-digit numbers, the complexity is O(n log n) due to FFT and inverse FFT steps.
# Naive grade-school multiplication is O(n^2), which is much slower for large n.
#
# Empirical results (for large n, e.g., n=2000):
# - FFT-based: completes in milliseconds to seconds, depending on hardware.
# - Naive: can take seconds to minutes for very large numbers.
#
# Theoretical Comparison:
# - FFT-based: O(n log n)
# - Naive: O(n^2)
#
# FFT-based multiplication is significantly faster for large integers, making it suitable for cryptography and scientific computing.
