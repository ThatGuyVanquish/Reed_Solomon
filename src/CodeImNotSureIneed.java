import java.util.Arrays;

public class CodeImNotSureIneed {

    /**
     * Returns the syndrome polynomial, generated by evaluating the Polynomial R(x) obtained by
     * treating the encoded symbols as coefficients, at n-k distinct powers of the basis q
     *
     * @param encodedSymbols the encoded symbols received by the encoder
     * @param k              length of the original message
     * @return the syndrome Polynomial of the encoded symbols.
     */
    public static Polynomial getSyndromePolynomial(Polynomial encodedSymbols, int k) {
        int q = encodedSymbols.getBasis();
        int n = encodedSymbols.degree() + 1;
        int[] syndromeCoefficients = new int[n - k];
        for (int i = 0; i < n - k; i++) {
            syndromeCoefficients[i] = encodedSymbols.evaluatePolynomial(ReedSolomon.powerModQ(q, i, q));
        }
        return new Polynomial(syndromeCoefficients, q);
    }

    /**
     * Returns the polynomial greatest common divisor of p1 and p2
     *
     * @param p1 polynomial to be the dividend
     * @param p2 polynomial to be the divisor
     * @return a Polynomial which is the GCD of p1 and p2.
     * @pre p1.getBasis() == p2.getBasis();
     * @post res.degree() <= n-k/2
     */
    public static Polynomial GCD(Polynomial p1, Polynomial p2) {
        Polynomial zero = Polynomial.ZERO(p1.getBasis());
        Polynomial dividend = new Polynomial(p1);
        Polynomial divisor = new Polynomial(p2);
        while (divisor.degree() > 0 || !divisor.equals(zero)) {
            Polynomial remainder = dividend.mod(divisor);
            dividend = divisor;
            divisor = remainder;
        }
        return dividend;
    }



    // I have another implementation of this in ReedSolomon using a different implementation of the next method
    public static Polynomial computeGeneratorPolynomial(int q, int n, int k) {
        int[] alphaPowers = new int[n - k + 1];
        int alpha = ReedSolomon.findPrimitiveElement(q);

        for (int i = 0; i < alphaPowers.length; i++) {
            alphaPowers[i] = ReedSolomon.powerModQ(alpha, i + 1, q);
        }
        int[] coeffsForGenerator = findPolynomialFromRoots(alphaPowers, q);
        int degree = coeffsForGenerator.length;
        for (int i = coeffsForGenerator.length - 1; i >= 0; i--) {
            if (coeffsForGenerator[i] == 0) degree--;
            else break;
        }
        int[] newCoeffs = Arrays.copyOf(coeffsForGenerator, degree);
        return new Polynomial(newCoeffs, q);
    }

    // Have an implementation in ReedSolomon that uses and returns a polynomial
    public static int[] findPolynomialFromRoots(int[] roots, int q) {
        int[] polynomial = new int[roots.length + 1];
        polynomial[0] = 1;

        for (int root : roots) {
            int[] factor = {-root, 1};

            polynomial = multiplyPolynomials(polynomial, factor, q);
        }

        return polynomial;
    }

    public static int[] multiplyPolynomials(int[] p, int[] q, int mod) {
        int[] result = new int[p.length + q.length - 1];

        for (int i = 0; i < p.length; i++) {
            for (int j = 0; j < q.length; j++) {
                result[i + j] = (result[i + j] + p[i] * q[j]) % mod;
            }
        }

        return result;
    }
}
