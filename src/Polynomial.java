import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

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

    public static Polynomial ZERO(int basis) {
        return new Polynomial(new int[]{0}, basis);
    }

    public static Polynomial ONE(int basis) {
        return new Polynomial(new int[]{1}, basis);
    }

    public static Polynomial SYMBOL_POLYNOMIAL(int[] coeffs) {
        return new Polynomial(coeffs, -1);
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
        int coeff = this.coefficients[degree] >= 0 ? this.coefficients[degree] :
                this.basis - Math.abs(this.coefficients[degree] % basis);
        return coeff;
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
        x = x % this.basis;
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
    public boolean in(Collection<Polynomial> c) {
        for (Polynomial p : c) {
            if (p.getBasis() != this.basis)
                continue;
            if (Arrays.equals(this.coefficients, p.getCoefficients()))
                return true;
        }
        return false;
    }

    public boolean in(Polynomial[] arr) {
        for (Polynomial p : arr) {
            if (p.getBasis() != this.basis)
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
        if (this.basis != that.getBasis() || this.degree() != that.degree()) return false;
        int k = this.degree() + 1;
        int q = this.basis;
        for(int i = 0; i < k; i++) {
            if (that.getCoefficient(i) != this.getCoefficient(i))
                return false;
        }
        return true;
    }

    /**
     * Returns the generator polynomial of Field Fq over (n,k).
     *
     * @param q the basis of Fq
     * @param n the length of the encoded message
     * @param k the length of the original message
     * @return the generator polynomial of field Fq over (n,k)
     */
    public static Polynomial computeGeneratorPolynomial(int q, int n, int k) {
        int[] alphaPowers = new int[n - k + 1];
        int alpha = Polynomial.findPrimitiveElement(q);

        for (int i = 0; i < alphaPowers.length; i++) {
            alphaPowers[i] = Polynomial.power(alpha, i + 1, q);
        }
        int[] coeffsForGenerator = Polynomial.findPolynomialFromRoots(alphaPowers, q);
        int degree = coeffsForGenerator.length;
        for(int i = coeffsForGenerator.length - 1; i >= 0; i--) {
            if (coeffsForGenerator[i] == 0) degree--;
            else break;
        }
        int[] newCoeffs = Arrays.copyOf(coeffsForGenerator, degree);
        return new Polynomial(newCoeffs, q);
    }

    /**
     * Returns the coefficients array of a desired polynomial over Fq based on the given roots.
     *
     * @param roots the to-be roots of the returned polynomial
     * @param q the basis of the field Fq
     * @return an integer array which holds the coefficient for a polynomial from Fq[X] whose roots are given.
     */
    public static int[] findPolynomialFromRoots(int[] roots, int q) {
        int[] polynomial = new int[roots.length + 1];
        polynomial[0] = 1;

        for (int root : roots) {
            int[] factor = {-root, 1};

            polynomial = Polynomial.multiplyPolynomials(polynomial, factor, q);
        }

        return polynomial;
    }

    public static int[] multiplyPolynomials(int[] p, int[] q, int mod) {
        int[] result = new int[p.length + q.length - 1];

        for (int i = 0; i < p.length; i++) {
            for (int j = 0; j < q.length; j++) {
                result[i+j] = (result[i+j] + p[i] * q[j]) % mod;
            }
        }

        return result;
    }

    public static int power(int base, int exponent, int mod) {
        int result = 1;
        while (exponent > 0) {
            if (exponent % 2 == 1) {
                result = (result * base) % mod;
            }
            base = (base * base) % mod;
            exponent /= 2;
        }
        return result;
    }

    public static int findPrimitiveElement(int q) {
        int[] factors = Polynomial.factor(q - 1);
        int alpha = 2;

        while (true) {
            boolean isPrimitive = true;
            for (int factor : factors) {
                int power = (q - 1) / factor;
                int alphaPower = Polynomial.power(alpha, power, q);
                if (alphaPower == 1) {
                    isPrimitive = false;
                    break;
                }
            }
            if (isPrimitive) {
                return alpha;
            }
            alpha++;
        }
    }

    public static int[] factor(int n) {
        List<Integer> factors = new ArrayList<>();
        int i = 2;

        while (i <= n) {
            if (n % i == 0) {
                factors.add(i);
                n /= i;
            } else {
                i++;
            }
        }

        return factors.stream().mapToInt(Integer::intValue).toArray();
    }





}