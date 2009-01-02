"""
Various basic math functions.
"""

# Compute the number of 1 bits in x.
cdef inline unsigned pop(unsigned x):
    assert sizeof(unsigned) == 4

    x = x - ((x >> 1) & 0x55555555)
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
    x = (x + (x >> 4)) & 0x0f0f0f0f
    x = x + (x >> 8)
    x = x + (x >> 16)
    return x & 0x0000003f

# Compute the greatest common divisor of u and v.
#
# Use the binary GCD algorithm, attributed by Knuth to Stein (1961).  See
# Knuth's TAOCP Vol. 2, 3rd Ed., pg 338 for details.
cdef inline long gcd(long u, long v) except *: # XXX
    cdef long k, t

    if u <= 0 or v <= 0:
        print "XXX gcd(%d, %d)" % (u, v)
    assert u > 0
    assert v > 0

    # Find power of 2.
    k = 0
    while (u & 1) == 0 and (v & 1) == 0:
        k += 1
        u >>= 1
        v >>= 1

    # Initialize.
    if (u & 1):
        t = -v
    else:
        t = (u >> 1)

    # Reduce t.
    while (t & 1) == 0:
        t >>= 1

    # Reset max(u, v).
    if t > 0:
        u = t
    else:
        v = -t

    # Subtract.
    t = u - v
    while t != 0:
        # Repeat various earlier steps.
        t >>= 1
        while (t & 1) == 0:
            t >>= 1

        if t > 0:
            u = t
        else:
            v = -t

        t = u - v

    return (u << k)

# Compute n choose k, such that overflow is avoided if at all possible.  On
# overflow, return 0.
cdef inline unsigned long choose(unsigned long n, unsigned long k):
    cdef unsigned long ret, l, a, b, i, j, g, i2, b2, j2, a2
    cdef bint overflow

    # (n k) and (n (n-k)) are the same; assure that k is as small as possible.
    l = n - k
    if l < k:
        l = k
        k = n - k

    # We maintain accumulators for the numerator and denominator: a/b.
    a = 1
    b = 1

    i = n
    j = k
    while i > l or j > 1:
        # Prepare to detect overflow, which results in no progress for an
        # entire iteration.
        overflow = True

        # Accumulate into a (numerator).
        while i > l:
            g = gcd(i, b)
            i2 = i / g
            b2 = b / g
            if a > 0x7fffffffffffffff / i2:
                # Overflowed a; switch accumulation modes.
                break
            a *= i2
            b = b2
            overflow = False
            i -= 1

        # Accumulate into b (denominator).
        while j > 1:
            g = gcd(j, a)
            j2 = j / g
            a2 = a / g
            if b > 0x7fffffffffffffff / j2:
                # Overflowed b; switch accumulation modes.
                break
            b *= j2
            a = a2
            overflow = False
            j -= 1

        # Check for overflow.
        if overflow:
            return 0

    assert b == 1
    return a
