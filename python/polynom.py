#!/usr/bin/python3

"""
This module contains a Python class, named Poly, which implements polynomials
with numerical coefficients and their usual ring operations, as well as
evaluation at numerical values, composition and symbolic differentiation.
"""

# For type annotations including Poly
from __future__ import annotations
# Abstract class for Numbers
from numbers import Number
# Arbitrary precision floating point numbers
from mpmath import mpf, mp, fabs


# Register mpmath floats as part of the abstract Number class
Number.register(mpf)

# Largest degree when polynomial multiplication is computed by Cauchy product
POLY_MUL_CUTOFF = 16


class Poly:
    """
    Python class implementing polynomials with Number coefficients and their
    ring operations, namely addition, subtraction, multiplication (implemented
    with Cauchy products or Karatsuba's algorithm, depending on the degree) and
    powers by non-negative integers, as well as evaluation at Number values,
    composition, and differentiation.
    """

    def __init__(self, *coeffs: list[Number], deg: int = None):
        """
        Initializes polynomial `coeffs[0] + ... + x**deg * coeffs[deg]`.
        """

        # If coefficients are not passed to `__init__`, set `self` to `0`
        if not coeffs:
            self.__coeffs__ = [0]
            self.__deg__ = 0
            return

        # Otherwise, coeffs must be a list of `Number`s
        elif not all(isinstance(coeff, Number) for coeff in coeffs):
            raise TypeError

        # If degree is not passed to `__init__`, set to default value
        if deg is None:
            deg = max(0, len(coeffs) - 1)

        # Degree of polynomial must be a non-negative integer
        elif not isinstance(deg, int):
            raise TypeError
        elif deg < 0:
            raise ValueError

        # Invoke setter for `self.coeffs` and `self.deg`
        self.coeffs = list(coeffs)
        self.deg = deg


    @property
    def deg(self) -> int:
        """Getter for degree. Returns internal degree variable."""
        return self.__deg__

    @deg.setter
    def deg(self, deg: int) -> None:
        """
        Setter for degree. Truncates coefficients if new degree is smaller,
        and does nothing if it is larger or equal. 
        """

        # `deg` must be a non-negative integer
        if not isinstance(deg, int):
            raise TypeError
        if deg < 0:
            raise ValueError

        # If `deg` is equal or greater than degree of `self`, do nothing
        if deg >= self.deg:
            return
        # Otherwise, truncate `self.__coeffs__` to match `deg`
        self.coeffs = self.__coeffs__[:deg + 1]

    @property
    def coeffs(self) -> list[Number]:
        """Getter for coefficients. Returns internal list of coefficients."""
        return self.__coeffs__.copy()

    @coeffs.setter
    def coeffs(self, coeffs: list[Number]) -> None:
        """
        Setter for list of coefficients. Prunes any leading zero coefficients
        and updates degree accordingly.
        """

        # `coeffs` must be a list of `Number`s
        if not isinstance(coeffs, list):
            raise TypeError
        if not all(isinstance(item, Number) for item in coeffs):
            raise TypeError

        # If `coeffs` is empty list, set `self` to `0`
        if not coeffs:
            self.__coeffs__ = [0]
            self.__deg__ = 0
            return

        # Otherwise, set `self.__coeffs__` to deep copy of `coeffs`
        self.__coeffs__ = coeffs.copy()
        # Prune leading zero coefficients from `self.__coeffs__`
        while len(self.__coeffs__) > 1 and self.__coeffs__[-1] == 0:
            self.__coeffs__.pop()
        # Update `deg` to match number of coefficients
        self.__deg__ = len(self.__coeffs__) - 1


    def __getitem__(self, key: int | slice) -> Number | Poly:
        """
        Getter for coefficients. Returns coefficient at index `key` if `key` is
        an `int`, or `Poly` of given coefficients if `key` is a `slice` object.
        """

        # Slice `self` if `key` is a `slice` object 
        if isinstance(key, slice):
            coeffs = self.__coeffs__[key]
            return Poly(*coeffs)

        # Otherwise, `key` must be a non-negative integer
        if not isinstance(key, int):
            raise TypeError
        if key < 0:
            raise KeyError

        # Always return `0` if `key` is larger than degree
        if key > self.deg:
            return 0
        
        # Otherwise, return coefficient at index `key`
        return self.__coeffs__[key]

    def __setitem__(self, key: int, value: Number) -> None:
        """
        Setter for coefficients. Sets coefficient at index `key` to `value`.
        """

        # `key` must be a non-negative integer less than or equal to degree
        if not isinstance(key, int):
            raise TypeError
        if key < 0 or key > self.deg:
            raise KeyError

        # `value` must be a Number
        if not isinstance(value, Number):
            raise ValueError
                
        # Set coefficient at index `key` to `value`
        self.__coeffs__[key] = value


    def __call__(self, other: Number | Poly) -> Number | Poly:
        """
        If `other` is a Number, then evaluate `self` at `other`. Otherwise,
        if `other` is a Poly, then compose polynomials `self` and `other`.
        """

        # If `other` is a Number, compute `self(other)` by Horner's method
        if isinstance(other, Number):
            # If `P(x) = coeffs[0] + ... + x**deg * coeffs[deg]` then we
            # have `P(x) = P(deg, x)`, where `P(i, x)` is defined by
            #
            #         P(0, x) == coeffs[deg]
            #     P(i + 1, x) == coeffs[deg - i] + x * P(i, x)
            #
            # for all 0 <= i <= deg.
            acc = 0
            for coeff in reversed(self.__coeffs__):
                acc = coeff + acc * other
            return acc

        # Otherwise, `other` must be a polynomial
        if not isinstance(other, Poly):
            return NotImplemented

        # Base case: if `self`` is constant, return its only coefficient
        if self.deg == 0:
            return Poly(self.__coeffs__[0])

        # Base case: if `other` is constant, evaluate `self` at `other[0]`
        if other.deg == 0:
            return self(other[0])

        # Otherwise, set d to the least power of 2 that is >= self.deg + 1
        d = 1 << (self.deg).bit_length()

        # Divide and conquer algorithm for computing P(Q(x)).
        #
        # To compute the composition of the polynomials P(x) and Q(x) of degree
        # d and n, where both are powers of 2, is equivalent to evaluating P(t)
        # at t = Q(x). First, we must find two polynomials P0(t) and P1(t) of
        # degree d/2 such that
        #
        #   P(t) = P0(t) + (t**(d/2)) * P1(t).
        #
        # Thus, we have reduced the problem of evaluating a degree d polynomial
        # to that of evaluating two degree d/2 polynomials. In addition, we
        # raise t = Q(x) to the power of d/2, which involves on the order of 
        # log(d/2) multiplications of degree (d/2) * n polynomials, and then
        # multiply the result by P1(t). By Karatsuba's algorithm for computing
        # products of polynomials, it takes on the order of (nd)**log2(3)
        # operationsto compute P(t) at t = Q(x), where log2(3) is approximately
        # equal to 1.585.

        # Define `P0`, `P1` such that `self(x) = P0(x) + (x ** (d/2)) * P1(x)`
        m = d >> 1
        P0 = self[:m]
        P1 = self[m:]

        # Evaluate `self(x) = P0(x) + (other ** (d/2)) * P1(x)` at `x = other`
        return P0(other) + (other ** m) * P1(other)


    def __str__(self) -> str:
        """
        Converts polynomial to string.
        """

        # If self is constant, return the str of its only coefficient
        if self.deg == 0:
            return str(self.__coeffs__[0])

        # List of formatted strings for each term
        terms = []
        # Terms are printed in order of decreasing degree
        for i in reversed(range(self.deg + 1)):
            coeff = self.__coeffs__[i]

            # Ignore zero coefficients when converting to string
            if coeff == 0:
                continue

            sign = ""
            if terms:
                # Add padded plus sign to non-leading positive terms
                if coeff > 0:
                    sign = " + "
                # Add padded minus sign to non-leading negative terms
                else:
                    sign = " - "
                    coeff = -coeff
            else:
                # Add minus sign to leading negative term
                if coeff < 0:
                    sign = "-"
                    coeff = -coeff

            # Omit non-leading unit coefficients
            if coeff == 1 and i > 0:
                coeff_str = ""
            else:
                coeff_str = str(coeff)

            # Format for powers of formal variable
            if i == 0:
                x_str = ""
            elif i == 1:
                x_str = "x"
            else:
                x_str = f"x**{i}"

            terms.append(f"{sign}{coeff_str}{x_str}")

        if terms:
            return "".join(terms)
        # Return "0" if all terms were null
        else:
            return "0"

    def __repr__(self) -> str:
        """
        Representation of Poly object (for debugging).
        """

        return f"Poly({self.__coeffs__}, deg={self.__deg__})"


    def __neg__(self) -> Poly:
        """
        Takes the negation of a polynomial.
        """

        # `-self(x)` is the same as `(-1) * self(x)`.
        return (-1) * self

    def __pos__(self) -> Poly:
        """
        Takes the positive of a polynomial.
        """
        # Returns itself
        return self

    def __add__(self, other: Poly) -> Poly:
        """
        Adds two polynomials coefficient by coefficient.
        """

        # `other` must be a polynomial
        if not isinstance(other, Poly):
            return NotImplemented

        deg = max(self.deg, other.deg)
        coeffs = [0] * (deg + 1)
        for i in range(deg + 1):
            if i <= self.deg:
                coeffs[i] += self[i]
            if i <= other.deg:
                coeffs[i] += other[i]
        return Poly(*coeffs)

    def __radd__(self, other: Poly) -> Poly:
        """
        Adds two polynomials coefficient by coefficient.
        """

        # `other` must be a polynomial
        if not isinstance(other, Poly):
            return NotImplemented

        return self.__add__(other)

    def __sub__(self, other: Poly) -> Poly:
        """
        Subtracts two polynomials coefficient by coefficient.
        """

        # `other` must be a polynomial
        if not isinstance(other, Poly):
            return NotImplemented

        deg = max(self.deg, other.deg)
        coeffs = [0] * (deg + 1)
        for i in range(deg + 1):
            if i <= self.deg:
                coeffs[i] += self[i]
            if i <= other.deg:
                coeffs[i] -= other[i]
        return Poly(*coeffs)

    def __rsub__(self, other: Poly) -> Poly:
        """
        Subtracts two polynomials coefficient by coefficient.
        """

        # `other` must be a polynomial
        if not isinstance(other, Poly):
            return NotImplemented
        
        
        return other.__sub__(self)


    def __mul__(self, other : Poly | Number) -> Poly:
        """
        Multiplities two polynomials or a polynomial by a scalar. Implements
        the Cauchy product formula for formal power series when the degree of
        `self` and `other` are low enough, and Karatsuba's algorithm when they
        are above a certain threshold.
        """

        # If other is a Number, return scalar product of `self` by `other`
        if isinstance(other, Number):
            coeffs = [other * self[i] for i in range(self.deg + 1)]
            return Poly(*coeffs)

        # Otherwise, `other` must be a polynomial
        if not isinstance(other, Poly):
            return NotImplemented

        # Set d to length of polynomial of highest degree
        d = max(self.deg, other.deg) + 1
        
        # Base case: if `self` and `other` are both constant, return product of
        # coefficients
        if d == 1:
            return Poly(self[0] * other[0])

        # Otherwise, compute `self * other` by the Cauchy product formula if d
        # is below a given threshold
        if d <= POLY_MUL_CUTOFF:
            deg = self.deg + other.deg
            coeffs = [0] * (deg + 1)

            # If `self` and `other` are given by
            # 
            #   self(x) = a[0] + ... + x**m * a[m],
            #   other(x) = b[0] + ... + x**n * b[n],
            #
            # then their product is given by
            #
            #   (self * other)(x) = c[0] + ... + x**(m + n) * c[m + n],
            #
            # where the coefficients `c[0]`, ..., `c[m + n]` obey the formula
            #
            #   c[k] = sum(a[i] b[k - i] for i from 0 to k)
            #
            # for every `0 <= k <= m + n`, setting `a[i] = 0`` and `b[j] = 0`
            # when `i > m` and `j > n`. In particular, we can prune the values
            # of `i` in the range from `0` to `k` if we know `a[i]` and `b[i]`
            # to always be 0 by this rule, i.e. when k - i > n and i > m.

            for k in range(deg + 1):
                i_min = k - other.deg if k > other.deg else 0
                i_max = min(k, self.deg)
                for i in range(i_min, i_max + 1):
                    coeffs[k] += self[i] * other[k - i]
            return Poly(*coeffs)

        # Otherwise, set d to least power of two g.t. or equal to d, and
        # compute `self * other` by Karatsuba's algorithm
        d = 1 << (d - 1).bit_length()
        
        # Define `P0`, `P1` such that `self(x) = P0(x) + (x ** (d/2)) * P1(x)`,
        # and `Q0`, `Q1` such that `other(x) = Q0(x) + (x ** (d/2)) * Q1(x)`
        m = d >> 1
        P0 = self[:m]
        P1 = self[m:]
        Q0 = other[:m]
        Q1 = other[m:]

        # Given `P0`, `P1`, `Q0` and `Q1` as above, first, note that
        #
        #     self * other
        #   = (P0 + (x**m) * P1) * (Q0 + (x**m) * Q1)
        #   = (P0 * Q0) + (x**(d/2)) * ((P0 * Q1) + (P1 * Q0)) + (x**d) * (P1 * Q1);
        #
        # and secondly, note that
        #
        #     (P0 * Q1) + (P1 * Q0)
        #   = (P0 + P1) * (Q0 + Q1) - (P0 * Q0) - (P1 * Q1);
        #
        # therefore, we can compute the product of two polynomials of
        # degree d with just 3 products of polynomials of degree d / 2.

        z0 = P0 * Q0
        z1 = (P0 + P1) * (Q0 + Q1)
        z2 = P1 * Q1
        middle = z1 - z0 - z2

        return z0 + (middle << m) + (z2 << d)

    def __rmul__(self, other : Poly | Number) -> Poly:
        """
        Multiplities two polynomials or a polynomial by a scalar. Implements
        the Cauchy product formula for formal power series when the degree of
        `self` and `other` are low enough, and Karatsuba's algorithm when they
        are above a certain threshold.
        """

        # If other is a Number, return scalar product of `self` by `other` 
        if isinstance(other, Number):
            return self.__mul__(other)

        # Otherwise, `other` is not a Number and not a polynomial
        return NotImplemented

    def __lshift__(self, other : int) -> Poly:
        """
        Multiplies polynomial by `x**other`.
        """

        # `other` must be a non-negative integer
        if not isinstance(other, int):
            raise TypeError
        if other < 0:
            raise ValueError

        # Degenerate case: if `self` is `0`, return itself
        if self.deg == 0 and self[0] == 0:
            return self

        # Otherwise, shift coefficients foward by `other` places
        coeffs = [0] * other + self.__coeffs__
        return Poly(*coeffs)
    
    def __pow__(self, other: int) -> Poly:
        """
        Takes the power of a polynomial by a non-negative integer. Implements
        exponentiation by squaring, using `Poly.__mul__` for multiplication.
        """

        # `other` must be a non-negative integer
        if not isinstance(other, int):
            raise TypeError
        if other < 0:
            raise ValueError

        acc = Poly(1)

        # Base cases: if `other` is zero, return unit polynomial
        if other == 0:
            return acc

        base = Poly(*self.__coeffs__)
        # ...and if `other` is one, return deep copy of `self``
        if other == 1:
            return base

        # Otherwise, compute `self(x) ** other` using exponentiation by
        # squaring. In essence, `self(x) ** other` reduces to either
        # 
        #   (self(x) ** 2) ** (other << 1)          if `other` is even,
        #   self(x) * (self(x) ** (other << 1))     if `other` is odd.
        while other > 0:
            if other % 2 == 1:
                acc = acc * base
            other >>= 1
            if other > 0:
                base = base * base
        return acc


    @property
    def d(self) -> Poly:
        """
        Takes the derivative of a polynomial.
        """

        # If self is constant, return 0
        if self.deg == 0:
            return Poly()

        # Otherwise, compute the derivative by symbolic differentiation
        coeffs = [i * self.__coeffs__[i] for i in range(1, self.deg + 1)]
        return Poly(*coeffs)

    def root(self, init: Number) -> Number:
        """
        Computes an approximation of the root of a polynomialis close to init.
        """

        # Initialise x to init
        x = mpf(init)
        # Initialise derivative of P
        d_self = self.d
        
        while True:
            step = self(x) / d_self(x)
            tol = mp.eps * fabs(x)

            # Break loop if step is less than the tolerance for precision errors
            if fabs(step) <= tol:
                break

            # New value of x, given by the Newton-Raphson algorithm
            x_new = x - step
            
            # Break loop if x will not update to a new and distinct value
            if x_new == x:
                break

            # Update x to x_new
            x = x_new

        return x

