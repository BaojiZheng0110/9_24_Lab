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
            for j in range(i, i + length // 2):
                u = a[j]
                v = a[j + length // 2] * w
                a[j] = u + v
                a[j + length // 2] = u - v
                w *= wlen
        length <<= 1
    if invert:
        for i in range(n):
            a[i] /= n

# TODO: Implement convolution using fft
def conv_fft(a: List[int], b: List[int]) -> List[int]:
    n = next_power_of_two(len(a) + len(b) - 1)
    fa = list(map(complex, a)) + [0] * (n - len(a))
    fb = list(map(complex, b)) + [0] * (n - len(b))
    fft(fa, False)
    fft(fb, False)
    for i in range(n):
        fa[i] *= fb[i]
    fft(fa, True)
    result = [int(round(fa[i].real)) for i in range(n)]
    return result

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
