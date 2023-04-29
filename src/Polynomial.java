import java.util.Arrays;
import java.util.Random;

public class Polynomial {

    private final int[] coefficients;
    private final int basis;

    /**
     * Class to represent Polynomials over Fq for basis q.
     *
     * @param coefficients List of coefficients
     * @param basis
     * @pre foreach i, coefficients[i] = coefficients[i] % basis;
     */
    public Polynomial(int[] coefficients, int basis) {
        this.coefficients = coefficients;
        this.basis = basis;
    }

    /**
     * Returns the degree of the polynomial.
     *
     * @return the degree of the polynomial
     */
    public int degree() {
        return coefficients.length - 1;
    }

    /**
     * Returns the basis q
     *
     * @return the basis of the polynomial
     */
    public int getBasis() {
        return this.basis;
    }

    /**
     * Returns the coefficients of the polynomial.
     *
     * @return an integer array of the polynomial's coefficients
     */
    public int[] getCoefficients() {
        return coefficients;
    }

    /**
     * Returns the coefficient for the $degree$th term.
     *
     * @param degree the degree of the term to retrieve the coefficient for
     * @return the coefficient for the $degree$th term
     */
    public int getCoefficient(int degree) throws IllegalArgumentException {
        if (degree < 0 || degree > this.degree())
            throw new IllegalArgumentException("Invalid degree: " + degree);
        return coefficients[degree];
    }

    /**
     * Returns a new polynomial that is the difference of this polynomial and the given polynomial.
     *
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
            result[i] += this.getCoefficient(i);
        }
        for (int i = 0; i <= other.degree(); i++) {
            result[i] = Math.abs(result[i] - other.getCoefficient(i)) % this.basis;
        }

        // ensure result is over basis q
        int degreeOfResult = result.length - 1;
        for (int i = result.length - 1; i > 0; i--) {
            if (result[i] == 0)
                degreeOfResult--;
            else break;
        }
        if (degreeOfResult == result.length - 1)
            return new Polynomial(result, this.basis);

        int[] resultOverQ = Arrays.copyOf(result, degreeOfResult + 1);
        return new Polynomial(resultOverQ, this.basis);
    }

    /**
     * Returns a new polynomial that is the multiplication of this polynomial and the given polynomial.
     *
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
                result[i + j] = (result[i + j] + this.getCoefficient(i) * other.getCoefficient(j)) % this.basis;
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
            return new Polynomial(result, this.basis);

        int[] resultOverQ = Arrays.copyOf(result, degreeOfResult + 1);
        return new Polynomial(resultOverQ, this.basis);
    }

    /**
     * Calculates the remainder of this polynomial modulo a given divisor polynomial.
     *
     * @param divisor the polynomial to use as divisor in the modulo operation
     * @return the remainder of this polynomial modulo the given divisor polynomial
     * //
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

        int q = this.basis;
        int degreeDiff = this.degree() - divisor.degree();
        int[] resultCoeffs = Arrays.copyOf(this.coefficients, this.coefficients.length);

        for (int i = degreeDiff; i >= 0; i--) {
            int ratio = resultCoeffs[i + divisor.degree()] / divisor.coefficients[divisor.degree()];
            for (int j = 0; j <= divisor.degree(); j++) {
                resultCoeffs[i + j] = (resultCoeffs[i + j] - ratio * divisor.coefficients[j] + q) % q;
            }
        }

        // ensure result is over basis q
        int resultDegree = this.degree();
        while (resultDegree > 0 && resultCoeffs[resultDegree] == 0) {
            resultDegree--;
        }
        return new Polynomial(Arrays.copyOf(resultCoeffs, resultDegree + 1), this.basis);
    }

    public int evaluatePolynomial(int x) {
        int res = 0;
        for(int i = 0; i <= this.degree(); i++) {
            res = res + this.getCoefficient(i) * (int)Math.pow(x,i);
        }
        return res % basis;
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

    public static Polynomial findIrreduciblePoly(int basis, int degree) {
        int[] coefficients = new int[degree + 1];
        coefficients[0] = 1;
        coefficients[degree] = 1;
        Polynomial poly = new Polynomial(coefficients, basis);

        for (int i = 2; i <= degree / 2; i++) {
            for (int j = i; j < degree; j += i) {
                poly.getCoefficients()[j] = (poly.getCoefficients()[j] + 1) % basis;
            }
        }

        for (int i = 2; i <= degree / 2; i++) {
            if (degree % i != 0) continue;
            Polynomial factor = findIrreduciblePoly(basis, i);
            Polynomial potentialIrreducible = factor == null ? null : poly.mod(factor);
            if (potentialIrreducible == null || potentialIrreducible.degree() == 0) return null;
        }

        return poly;
    }


}