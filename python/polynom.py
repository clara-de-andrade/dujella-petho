#!/usr/bin/python3

"""
This module contains a class Poly which implements polynomials with numerical
coefficients, as well as their ring operations, evaluation at numerical and
polynomial values, composition with other polynomials, and differentiation.
"""

from __future__ import annotations
from numbers import Number

from math import log2, floor
from mpmath import mpf

Number.register(mpf)

class Poly:
    """
    Python class implementing polynomials with Number coefficients and their
    ring operations, namely addition, subtraction, multiplication, powers by
    non-negative integers, and Euclidean division, as well as evaluation at
    Number values, composition with other polynomials, and differentiation.
    """

    def __init__(self, deg: int = 0, coeffs: list[Number] = None):
        if not isinstance(deg, int):
            raise TypeError
        if deg < 0:
            raise ValueError
        
        self.__deg__ = deg
        len_self = deg + 1

        # By default, set each coefficient to zero
        if coeffs is None:
            coeffs = [0] * len_self

        if not isinstance(coeffs, list):
            raise TypeError
        
        if not all(isinstance(coeff, Number) for coeff in coeffs):
            raise TypeError
        
        self.__coeffs__ = coeffs.copy()
        len_coeffs = len(coeffs)

        # If there are less coefficients than degree, pad with zeroes
        if len_coeffs < len_self:
            self.__coeffs__.extend([0] * (len_self - len_coeffs))
        # Otherwise, set degree to number of coefficients
        else:
            self.__deg__ = len_coeffs - 1


    @property
    def deg(self) -> int:
        return self.__deg__

    @deg.setter
    def deg(self, value: int) -> None:
        if not isinstance(value, int):
            raise TypeError
        if value < 0:
            raise ValueError

        # If there are less coefficients than value, pad with zeroes
        if self.__deg__ <= value:
            self.__coeffs__.extend([0] * (value - self.__deg__))
        # Otherwise, truncate coefficients to match value
        else:
            self.__coeffs__ = self.__coeffs__[:value + 1]

        # Set degree to value
        self.__deg__ = value

    @property
    def coeffs(self) -> list[Number]:
        return self.__coeffs__.copy()

    @coeffs.setter
    def coeffs(self, value: list[Number]) -> None:
        if not isinstance(value, list):
            raise TypeError
        if not all(isinstance(item, Number) for item in value):
            raise TypeError

        # If value is non-empty list, update coeffs and deg
        if value:
            self.__coeffs__ = value.copy()
            self.__deg__ = len(value) - 1
        # Otherwise, set self to null polynomial
        else:
            self.__coeffs__ = [0]
            self.__deg__ = 0


    def __len__(self) -> int:
        # Return number of coefficients
        return self.deg + 1

    def __setitem__(self, key: int, value: Number) -> None:
        if not isinstance(key, int):
            raise TypeError
        if key < 0 or key > self.deg:
            raise KeyError

        if not isinstance(value, Number):
            raise ValueError
                
        # Set coefficient at index key to value
        self.__coeffs__[key] = value

    def __getitem__(self, key: int | slice) -> Number | Poly:
        # Slice list of coefficients 
        if isinstance(key, slice):
            coeffs = self.__coeffs__[key]
            return Poly(max(0, len(coeffs) - 1), coeffs)

        if not isinstance(key, int):
            raise TypeError
        if key < 0:
            raise KeyError

        if key > self.deg:
            return 0
        
        # Return coefficient at index key
        return self.__coeffs__[key]


    def __call__(self, other: Poly | Number) -> Poly | Number:
        # Evaluate polynomial self(x) at x = other using Horner's method,
        # where self(x) = P_deg(x) and for all 0 <= n < deg we have that
        #
        #     P_0(x) = coeffs[deg],
        #     P_(n + 1)(x) = coeffs[deg - n] + x * P_n(x)
        if isinstance(other, Number):
            acc = 0
            for coeff in reversed(self.__coeffs__):
                acc = coeff + other * acc
            return acc

        if not isinstance(other, Poly):
            return NotImplemented

        # If self is a constant polynomial, return its only coefficient
        if self.deg == 0:
            return Poly(0, [self.__coeffs__[0]])

        # Set d to least power of two g.t. or equal to number of coefficients
        d = self.deg + 1
        if not log2(d).is_integer():
            d = 1 << (floor(log2(d)) + 1)

        # Set P0 and P1 such that self(x) = P0(x) + (x ** m) * P1(x)
        m = d >> 1
        P0 = self[:m]
        P1 = self[m:]

        # Evaluate self(x) = P0(x) + (other ** m) * P1(x) at x = other
        return P0(other) + (other ** m) * P1(other)

    def __str__(self) -> str:
        # If self is a constant polynomial, return the str of its coefficient
        if self.deg == 0:
            return str(self.__coeffs__[0])

        terms = []
        for i in reversed(range(self.deg + 1)):
            coeff = self.__coeffs__[i]

            # Ignore null coefficients
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

            # Format powers of formal variable
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
        return f"Poly(deg={self.deg}, coeffs={self.__coeffs__})"


    def __add__(self, other: Poly) -> Poly:
        # Term-wise addition of polynomials
        if isinstance(other, Poly):
            deg = max(self.deg, other.deg)
            coeffs = [0] * (deg + 1)
            for i in range(deg + 1):
                if i <= self.deg:
                    coeffs[i] += self[i]
                if i <= other.deg:
                    coeffs[i] += other[i]
            return Poly(deg, coeffs)

        else:
            return NotImplemented

    def __radd__(self, other: Poly) -> Poly:
        if isinstance(other, Poly):
            return self.__add__(other)

        else:
            return NotImplemented

    def __sub__(self, other: Poly) -> Poly:
        # Term-wise subtraction of polynomials
        if isinstance(other, Poly):
            deg = max(self.deg, other.deg)
            coeffs = [0] * (deg + 1)
            for i in range(deg + 1):
                if i <= self.deg:
                    coeffs[i] += self[i]
                if i <= other.deg:
                    coeffs[i] -= other[i]
            return Poly(deg, coeffs)

        else:
            return NotImplemented

    def __rsub__(self, other: Poly) -> Poly:
        if isinstance(other, Poly):
            return other.__sub__(self)

        else:
            return NotImplemented


    def __mul__(self, other : Poly | Number) -> Poly:
        # If other is a Number, return scalar product of self(x) by other
        if isinstance(other, Number):
            return Poly(self.deg, [other * self[i] for i in range(len(self))])

        # Otherwise, if other is a polynomial, compute self(x) * other(x)
        # using Karatsuba's algorithm

        if not isinstance(other, Poly):
            return NotImplemented

        # If self and other are constant polynomials, return the product of
        # their coefficients. (Trivial case.)
        d = max(self.deg, other.deg) + 1
        if d == 1:
            return Poly(0, [self[0] * other[0]])


        # Set d to least power of two g.t. or equal to max number of coefficients
        if not log2(d).is_integer():
            d = 1 << (floor(log2(d)) + 1)

        
        # Set P0 and P1 such that self(x) = P0(x) + (x ** m) * P1(x), and
        # set Q0 and Q1 such that other(x) = Q0(x) + (x ** m) * Q1(x)
        m = d >> 1
        P0 = self[:m]
        P1 = self[m:]
        Q0 = other[:m]
        Q1 = other[m:]

        # Given the above polynomials, first, note that
        #
        #       self * other
        #     = (P0 + (x**m) * P1) * (Q0 + (x**m) * Q1)
        #     = (P0 * Q0) + (x**m) * (P0 * Q1 + P1 * Q0) + (x**d) * (P1 * Q1);
        #
        # and secondly, note that
        #
        #       (P0 * Q1 + P1 * Q0)
        #     = (P0 + P1) * (Q0 + Q1) - P0 * Q0) - (P1 * Q1);
        #
        # therefore, we can compute a product self * other of polynomials of
        # degree d with just 3 products of polynomials of degree d / 2.

        z0 = P0 * Q0
        z1 = (P0 + P1) * (Q0 + Q1)
        z2 = P1 * Q1
        z3 = z1 - z0 - z2

        deg = self.deg + other.deg
        coeffs = [0] * (2 * d - 1)
        for i in range(len(z0)):
            coeffs[i] += z0[i]
        for i in range(len(z3)):
            coeffs[m + i] += z3[i]
        for i in range(len(z2)):
            coeffs[d + i] += z2[i]

        return Poly(deg, coeffs[:deg + 1])

    def __rmul__(self, other : Poly | Number) -> Poly:
        if isinstance(other, Number):
            return self.__mul__(other)
    
        return NotImplemented

    def __pow__(self, exp: int) -> Poly:
        if not isinstance(exp, int):
            raise TypeError
        if exp < 0:
            raise ValueError

        # Compute self(x) ** exp using exponentiation by squaring
        acc = Poly(0, [1])
        base = Poly(self.deg, self.__coeffs__)
        while exp > 0:
            if exp & 1 == 1:
                acc = acc * base
            base = base * base
            exp = exp >> 1
        return acc

    def __neg__(self) -> Poly:
        return (-1) * self

    def __pos__(self) -> Poly:
        return self


    @property
    def d(self) -> Poly:
        # If self is a constant polynomial, the derivative is a null polynomial
        if self.deg == 0:
            return Poly()

        # Otherwise, compute the derivative by the power rule
        deg = self.deg - 1
        coeffs = [i * self.__coeffs__[i] for i in range(1, self.deg + 1)]
        return Poly(deg, coeffs)

