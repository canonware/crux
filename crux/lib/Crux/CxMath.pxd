"""
Various basic math functions.
"""

from libm cimport *

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
cdef inline long long gcd(long long u, long long v):
    cdef long long k, t

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
cdef inline unsigned long long choose(unsigned long long n,
  unsigned long long k):
    cdef unsigned long long ret, l, a, b, i, j, g, i2, b2, j2, a2
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

cdef inline double lnGamma(double x):
    """
        Compute the natural log of Gamma(x), accurate to 10 decimal places.

        This implementation is based on:

          Pike, M.C., I.D. Hill (1966) Algorithm 291: Logarithm of Gamma
          function [S14].  Communications of the ACM 9(9):684.
    """
    cdef double f, z

    assert x > 0.0

    if x < 7.0:
        f = 1.0
        z = x
        while z < 7.0:
            f *= z
            z += 1.0
        x = z
        f = -log(f)
    else:
        f = 0.0

    z = 1.0 / (x * x);

    return f + (x-0.5) * log(x) - x + 0.918938533204673 + \
      (((-0.000595238095238 * z + 0.000793650793651) * z - 0.002777777777778) \
      * z + 0.083333333333333) / x

cdef inline double gammaI(double x, double p, double lnGammaP):
    """
        Compute the incomplete Gamma ratio for [0..x], where p is the shape
        parameter, and lnGammaP is lnGamma(p).

        This implementation is based on:

          Bhattacharjee, G.P. (1970) Algorithm AS 32: The incomplete Gamma
          integral.  Applied Statistics 19:285-287.
    """
    cdef double acu, factor, oflo, gin, term, rn, a, b, an, dif
    cdef double pn[6]
    cdef unsigned i

    assert p > 0.0
    assert x >= 0.0

    if x == 0.0:
        return 0.0

    acu = 1.0e-10
    oflo = 1.0e30
    gin = 0.0
    factor = exp(p * log(x) - x - lnGammaP)

    if x <= 1.0 or x < p:
        # Calculation by series expansion.
        gin = 1.0
        term = 1.0
        rn = p

        while True:
            rn += 1.0
            term *= x / rn
            gin += term
            if term <= acu:
                gin *= factor / p
                return gin
    else:
        # Calculation by continued fraction.
        a = 1.0 - p
        b = a + x + 1.0
        term = 0.0
        pn[0] = 1.0
        pn[1] = x
        pn[2] = x + 1.0
        pn[3] = x * b
        gin = pn[2] / pn[3]

        while True:
            a += 1.0
            b += 2.0
            term += 1.0
            an = a * term
            for 0 <= i < 2:
                pn[i+4] = b * pn[i+2] - an * pn[i]
            if pn[5] != 0.0:
                rn = pn[4] / pn[5]
                dif = fabs(gin - rn)
                if dif <= acu and dif <= acu * rn:
                    gin = 1.0 - factor * gin
                    return gin
                gin = rn
            for 0 <= i < 4:
                pn[i] = pn[i+2]

            if fabs(pn[4]) >= oflo:
                for 0 <= i < 4:
                    pn[i] /= oflo

cdef inline double ptNorm(double p):
    """
        Given a value p in [0..1] of the lower tail area of the normal
        distribution, compute the limit on the definite integral from [-inf..z]
        that satisfies p, accurate to 16 decimal places.

        This implementation is based on:

          Wichura, M.J. (1988) Algorithm AS 241: The percentage points of the
          normal distribution.  Applied Statistics 37(3):477-484.
    """
    cdef double q, r, ret

    assert p > 0.0 and p < 1.0

    q = p - 0.5
    if fabs(q) <= 0.425:
        # p close to 1/2.
        r = 0.180625 - q * q
        return q * (((((((2.5090809287301226727e3 * r + \
          3.3430575583588128105e4) * r + 6.7265770927008700853e4) * r + \
          4.5921953931549871457e4) * r + 1.3731693765509461125e4) * r + \
          1.9715909503065514427e3) * r + 1.3314166789178437745e2) * r + \
          3.3871328727963666080e0) / (((((((5.2264952788528545610e3 * r + \
          2.8729085735721942674e4) * r + 3.9307895800092710610e4) * r + \
          2.1213794301586595867e4) * r + 5.3941960214247511077e3) * r + \
          6.8718700749205790830e2) * r + 4.2313330701600911252e1) * r + 1.0)
    else:
        if q < 0.0:
            r = p
        else:
            r = 1.0 - p
        assert r > 0.0

        r = sqrt(-log(r))
        if r <= 5.0:
            # p neither close to 1/2 nor 0 or 1.
            r -= 1.6
            ret = (((((((7.74545014278341407640e-4 * r + \
              2.27238449892691845833e-2) * r + 2.41780725177450611770e-1) * r \
              + 1.27045825245236838258e0) * r + 3.64784832476320460504e0) * r \
              + 5.76949722146069140550e0) * r + 4.63033784615654529590e0) * r \
              + 1.42343711074968357734e0) / (((((((1.05075007164441684324e-9 * \
              r + 5.47593808499534494600e-4) * r + 1.51986665636164571966e-2) \
              * r + 1.48103976427480074590e-1) * r + \
              6.89767334985100004550e-1) * r + 1.67638483018380384940e0) * r + \
              2.05319162663775882187e0) * r + 1.0)
        else:
            # p near 0 or 1.
            r -= 5.0
            ret = (((((((2.01033439929228813265e-7 * r + \
              2.71155556874348757815e-5) * r + 1.24266094738807843860e-3) * r \
              + 2.65321895265761230930e-2) * r + 2.96560571828504891230e-1) * \
              r + 1.78482653991729133580e0) * r + 5.46378491116411436990e0) * \
              r + 6.65790464350110377720e0) / \
              (((((((2.04426310338993978564e-15 * r + \
              1.42151175831644588870e-7) * r + 1.84631831751005468180e-5) * r \
              + 7.86869131145613259100e-4) * r + 1.48753612908506148525e-2) * \
              r + 1.36929880922735805310e-1) * r + 5.99832206555887937690e-1) \
              * r + 1.0)
        if q < 0.0:
            ret = -ret
        return ret

cdef inline double ptChi2(double p, double v, double lnGammaV2):
    """
        Given a value p in [0..1] of the lower tail area of the Chi^2
        distribution with v degrees of freedom, where lnGammaV2 is
        lnGamma(v/2.0), compute the upper limit on the definite integral from
        [0..z] that satisfies p, accurate to 12 decimal places.

        This implementation is based on:

          Best, D.J., D.E. Roberts (1975) Algorithm AS 91: The percentage
          points of the Chi^2 distribution.  Applied Statistics 24(3):385-388.

          Shea, B.L. (1991) Algorithm AS R85: A remark on AS 91: The percentage
          points of the Chi^2 distribution.  Applied Statistics 40(1):233-235.
    """
    cdef double e, aa, xx, c, ch, a, q, p1, p2, t, x, b, s1, s2, s3, s4, s5, s6
    cdef unsigned i

    assert p >= 0.0 and p < 1.0
    assert v > 0.0

    e = 5.0e-7
    aa = 0.6931471805

    xx = 0.5 * v
    c = xx - 1.0

    if v < -1.24 * log(p):
        # Starting approximation for small Chi^2.
        ch = pow(p * xx * exp(lnGammaV2 + xx * aa), 1.0 / xx)
        if ch - e < 0.0:
            return ch
    else:
        if v > 0.32:
            x = ptNorm(p)
            # Starting approximation using Wilson and Hilferty estimate.
            p1 = 0.222222 / v
            ch = v * pow(x * sqrt(p1) + 1.0 - p1, 3.0)
            # Starting approximation for p tending to 1.
            if ch > 2.2 * v + 6.0:
                ch = -2.0 * (log(1.0 - p) - c * log(0.5 * ch) + lnGammaV2)
        else:
            ch = 0.4
            a = log(1.0 - p)
            while True:
                q = ch
                p1 = 1.0 + ch * (4.67 + ch)
                p2 = ch * (6.73 + ch * (6.66 + ch))
                t = -0.5 + (4.67 + 2.0 * ch) / p1 - (6.73 + ch * (13.32 + 3.0 \
                  * ch)) / p2
                ch -= (1.0 - exp(a + lnGammaV2 + 0.5 * ch + c * aa) * p2 / p1) \
                  / t
                if fabs(q / ch - 1.0) - 0.01 <= 0.0:
                    break

    for 0 <= i < 20:
        # Calculation of seven-term Taylor series.
        q = ch
        p1 = 0.5 * ch
        if p1 < 0.0:
            return -1.0
        p2 = p - gammaI(p1, xx, lnGammaV2)
        t = p2 * exp(xx * aa + lnGammaV2 + p1 - c * log(ch))
        b = t / ch
        a = 0.5 * t - b * c
        s1 = (210.0 + a * (140.0 + a * (105.0 + a * (84.0 + a * (70.0 + 60.0 * \
          a))))) / 420.0
        s2 = (420.0 + a * (735.0 + a * (966.0 + a * (1141.0 + 1278.0 * a)))) / \
          2520.0
        s3 = (210.0 + a * (462.0 + a * (707.0 + 932.0 * a))) / 2520.0
        s4 = (252.0 + a * (672.0 + 1182.0 * a) + c * (294.0 + a * (889.0 + \
          1740.0 * a))) / 5040.0
        s5 = (84.0 + 264.0 * a + c * (175.0 + 606.0 * a)) / 2520.0
        s6 = (120.0 + c * (346.0 + 127.0 * c)) / 5040.0
        ch += t * (1.0 + 0.5 * t * s1 - b * c * (s1 - b * (s2 - b * (s3 - b * \
          (s4 - b * (s5 - b * s6))))))
        if fabs(q / ch - 1.0) <= e:
            break

    return ch
