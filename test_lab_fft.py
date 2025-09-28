# test_lab_fft.py — Pytest tests for FFT lab
import solution


def test_small_cases():
    cases = [
        ("0","0"),
        ("0","123456"),
        ("1","999999"),
        ("123","456"),
        ("3141592653589793","2718281828459045")
    ]
    for x,y in cases:
        assert solution.bigmul_fft(x,y,k=3) == str(int(x)*int(y))

def test_large_case():
    import random
    n = 2000
    x = ''.join(str(random.randint(0,9)) for _ in range(n))
    y = ''.join(str(random.randint(0,9)) for _ in range(n))
    got = solution.bigmul_fft(x,y,k=3)
    exp = str(int(x)*int(y))
    assert got == exp
