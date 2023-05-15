package Code;

import java.util.*;

public class Polynomial {

    private final int[] coefficients;
    private final GaloisField F;

    /**
     * Class to represent Polynomials over Fq for basis q.
     *
     * @param coefficients List of coefficients
     * @param field the galois field Fq
     * @pre foreach i, coefficients[i] = coefficients[i] % basis;
     */
    public Polynomial(int[] coefficients, GaloisField field) {
        if (coefficients.length == 0)
            this.coefficients = new int[]{0};
        else
            this.coefficients = coefficients;
        this.F = field;
    }

    public Polynomial(Polynomial other) {
        this.coefficients = other.getCoefficients();
        this.F = other.getField();
    }

    public static Polynomial ZERO(GaloisField field) {
        return new Polynomial(new int[]{0}, field);
    }

    public static Polynomial ONE(GaloisField field) {
        return new Polynomial(new int[]{1}, field);
    }

    /**
     * Returns the degree of the polynomial.
     * @return the degree of the polynomial
     */
    public int degree() {
        return coefficients.length - 1;
    }

    /**
     * Returns the field Fq this polynomial is over.
     * @return the field Fq for which this is in Fq[x]
     */
    public GaloisField getField() {
        return this.F;
    }

    /**
     * Returns the coefficients of the polynomial.
     * @return an integer array of the polynomial's coefficients
     */
    public int[] getCoefficients() {
        return coefficients;
    }

    /**
     * Returns the coefficient for the $degree$th term modulo q (positive value).
     * @param degree the degree of the term to retrieve the coefficient for
     * @return the coefficient for the $degree$th term
     */
    public int getCoefficient(int degree) {
        if (degree < 0 || degree > this.degree())
            return 0;
        return F.mod(this.coefficients[degree]);
    }

    /**
     * Returns a new polynomial that is the difference of this polynomial and the given polynomial.
     * @param other the polynomial to subtract from this polynomial
     * @return a new polynomial that is the difference of this polynomial and the given polynomial

     * @pre this.basis == other.getBasis();
     * @post result.getBasis() == this.basis; => new polynomial is over basis q
     * @post result.degree() <= Math.max(this.degree(), other.degree())
     */
    public Polynomial subtract(Polynomial other) {
        int[] result = new int[Math.max(this.degree(), other.degree()) + 1];

        // calculate subtraction
        for (int i = 0; i <= this.degree(); i++) {
            result[i] = F.add(result[i], this.getCoefficient(i));
        }
        for (int i = 0; i <= other.degree(); i++) {
            result[i] = F.subtract(result[i], other.getCoefficient(i));
        }

        // ensure result is over basis q
        int degreeOfResult = result.length - 1;
        for (int i = result.length - 1; i > 0; i--) {
            if (result[i] == 0)
                degreeOfResult--;
            else break;
        }
        if (degreeOfResult == result.length - 1)
            return new Polynomial(result, this.F);

        int[] resultOverQ = Arrays.copyOf(result, degreeOfResult + 1);
        return new Polynomial(resultOverQ, this.F);
    }

    /**
     * Given a polynomial returns the result of adding this polynomial to the given polynomial.
     * @param other given polynomial to add to this
     * @return Polynomial result of this+other
     * @pre this.getBasis() == other.getBasis()
     * @post this.getBasis() == result.getBasis()
     * @post result.degree() <= Math.max(this.degree(), other.degree())
     */
    public Polynomial add(Polynomial other) {
        int[] newCoeffs = new int[Math.max(this.degree(), other.degree()) + 1];
        for(int i = 0; i < newCoeffs.length; i++) {
            if (i > this.degree()) {
                newCoeffs[i] = other.getCoefficient(i);
                continue;
            }
            if (i > other.degree()) {
                newCoeffs[i] = this.getCoefficient(i);
                continue;
            }
            newCoeffs[i] = F.add(this.getCoefficient(i), other.getCoefficient(i));
        }
        return new Polynomial(newCoeffs, this.F);
    }

    /**
     * Returns a new polynomial that is the multiplication of this polynomial and the given polynomial.
     * @param other the polynomial to multiply this polynomial by
     * @return a new polynomial that is the result of multiplication of this polynomial and the given polynomial

     * @pre this.basis == other.getBasis();
     * @post result.getBasis() == this.basis;
     * @post result.degree() <= this.degree() + other.degree() + 1;
     */
    public Polynomial multiply(Polynomial other) {
        int[] result = new int[this.degree() + other.degree() + 1];
        // calculate multiplication over basis q
        for (int i = 0; i <= this.degree(); i++) {
            for (int j = 0; j <= other.degree(); j++) {
                result[i + j] = F.add(result[i + j], F.multiply(this.getCoefficient(i), other.getCoefficient(j)));
            }
        }

        // ensure result is over basis q
        int degreeOfResult = result.length - 1;
        for (int i = result.length - 1; i > 0; i--) {
            if (result[i] == 0)
                degreeOfResult--;
            else break;
        }
        if (degreeOfResult == result.length - 1)
            return new Polynomial(result, this.F);

        int[] resultOverQ = Arrays.copyOf(result, degreeOfResult + 1);
        return new Polynomial(resultOverQ, this.F);
    }

    /**
     * Calculates the remainder of this polynomial modulo a given divisor polynomial.
     * @param divisor the polynomial to use as divisor in the modulo operation
     * @return the remainder of this polynomial modulo the given divisor polynomial

     * @pre this.getBasis() == divisor.getBasis();
     * @post result.getBasis() == this.getBasis();
     * @post result.degree() < divisor.degree();
     */
    public Polynomial mod(Polynomial divisor) throws IllegalArgumentException {
        if (this.degree() < divisor.degree()) {
            return this;
        }

        if (divisor.degree() == 0 && divisor.getCoefficient(0) == 0)
            throw new IllegalArgumentException("Divisor can't be zero!!!!");

        Polynomial dividend = new Polynomial(this.coefficients, this.F);
        int divisorsLeadingCoefficient = divisor.getCoefficient(divisor.degree());

        while(dividend.degree() >= divisor.degree()) {
            int dividendsLeadingCoefficient = dividend.getCoefficient(dividend.degree());
            int multiplier = this.F.div(dividendsLeadingCoefficient, divisorsLeadingCoefficient);
            int degreeDifference = dividend.degree() - divisor.degree();
            int[] coeffsForMultiplierPolynomial = new int[degreeDifference + 1];
            coeffsForMultiplierPolynomial[coeffsForMultiplierPolynomial.length - 1] = multiplier;
            Polynomial multiplierPolynomial = new Polynomial(coeffsForMultiplierPolynomial, this.F);
            Polynomial multipliedDivisor = divisor.multiply(multiplierPolynomial);
            dividend = dividend.subtract(multipliedDivisor);
            if (dividend.equals(ZERO(this.F)))
                break;
        }
        // dividend once completed is the remainder
        return dividend;

    }

//    /**
//     * Given a polynomial divisor, divide this polynomial by divisor using long polynomial division and return the result.
//     * @param divisor the polynomial to use as divisor in the modulo operation
//     * @return the division of this polynomial by the given divisor polynomial
//
//     * @pre this.getBasis() == divisor.getBasis();
//     * @post result.getBasis() == this.getBasis();
//     * @post result.degree() < divisor.degree();
//     */
    public Polynomial div(Polynomial divisor) {
        if (divisor.equals(ZERO(this.F))) {
            throw new ArithmeticException("Division by zero polynomial");
        }

        Polynomial dividend = new Polynomial(this.coefficients, this.F);
        Polynomial result = ZERO(this.F);
        int divisorsLeadingCoefficient = divisor.getCoefficient(divisor.degree());

        while(dividend.degree() >= divisor.degree()) {
            int dividendsLeadingCoefficient = dividend.getCoefficient(dividend.degree());
            int multiplier = this.F.div(dividendsLeadingCoefficient, divisorsLeadingCoefficient);
            int degreeDifference = dividend.degree() - divisor.degree();
            int[] coeffsForMultiplierPolynomial = new int[degreeDifference + 1];
            coeffsForMultiplierPolynomial[coeffsForMultiplierPolynomial.length - 1] = multiplier;
            Polynomial multiplierPolynomial = new Polynomial(coeffsForMultiplierPolynomial, this.F);
            Polynomial multipliedDivisor = divisor.multiply(multiplierPolynomial);
            dividend = dividend.subtract(multipliedDivisor);
            result = result.add(multiplierPolynomial);
            if (dividend.equals(ZERO(this.F)))
                break;
        }
        return result;
    }


    /**
     * Returns the result of evaluating this polynomial at point x, modulo q.
     * @param x the value of the point to evaluate this polynomial at.
     * @return result of the evaluation of this polynomial at point x.
     */
    public int evaluatePolynomial(int x) {
        x = this.F.mod(x);
        int res = 0;
        for(int i = 0; i <= this.degree(); i++) {
            res = res + this.getCoefficient(i) * (int)Math.pow(x,i);
        }
        return this.F.mod(res);
    }

    @Override
    public String toString() {
        StringBuilder res = new StringBuilder();
        for (int i = degree(); i >= 0; i--) {
            if (this.getCoefficient(i) != 0) {
                if (res.length() > 0) {
                    res.append(" + ");
                }
                if (this.getCoefficient(i) != 1 || i == 0) {
                    res.append(getCoefficient(i));
                }
                if (i > 0) {
                    res.append("x");
                    if (i > 1) {
                        res.append("^").append(i);
                    }
                }
            } else if (this.degree() == 0)
                res.append("0");
        }
        return res.toString();
    }

    /**
     * Returns true if an identical polynomial to this is within c.
     * @param c collection to check whether this polynomial is in
     * @return true if found, false otherwise.
     */
    public boolean in(Collection<Polynomial> c) {
        for (Polynomial p : c) {
            if (!p.getField().equals(this.F))
                continue;
            if (Arrays.equals(this.coefficients, p.getCoefficients()))
                return true;
        }
        return false;
    }

    /**
     * Returns true if an identical polynomial to this is within arr.
     * @param arr array to check whether this polynomial is in
     * @return true if found, false otherwise.
     */
    public boolean in(Polynomial[] arr) {
        for (Polynomial p : arr) {
            if (!(p.getField().equals(this.F)))
                continue;
            if (Arrays.equals(this.coefficients, p.getCoefficients()))
                return true;
        }
        return false;
    }

    @Override
    public boolean equals(Object other) {
        if (this == other) return true;
        if (other == null || getClass() != other.getClass()) return false;
        Polynomial that = (Polynomial) other;
        if (!(this.F.equals(that.getField())) || this.degree() != that.degree()) return false;
        int k = this.degree() + 1;
        for(int i = 0; i < k; i++) {
            if (that.getCoefficient(i) != this.getCoefficient(i))
                return false;
        }
        return true;
    }

}