import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class CodeImNotSureIneed {

//    /**
//     * Returns the syndrome polynomial, generated by evaluating the Polynomial R(x) obtained by
//     * treating the encoded symbols as coefficients, at n-k distinct powers of the basis q
//     *
//     * @param encodedSymbols the encoded symbols received by the encoder
//     * @param k              length of the original message
//     * @return the syndrome Polynomial of the encoded symbols.
//     */
//    public static Polynomial getSyndromePolynomial(Polynomial encodedSymbols, int k) {
//        int q = encodedSymbols.getBasis();
//        int n = encodedSymbols.degree() + 1;
//        int[] syndromeCoefficients = new int[n - k];
//        for (int i = 0; i < n - k; i++) {
//            syndromeCoefficients[i] = encodedSymbols.evaluatePolynomial(ReedSolomon.powerModQ(q, i, q));
//        }
//        return new Polynomial(syndromeCoefficients, q);
//    }

//    /**
//     * Returns the polynomial greatest common divisor of p1 and p2
//     *
//     * @param p1 polynomial to be the dividend
//     * @param p2 polynomial to be the divisor
//     * @return a Polynomial which is the GCD of p1 and p2.
//     * @pre p1.getBasis() == p2.getBasis();
//     * @post res.degree() <= n-k/2
//     */
//    public static Polynomial GCD(Polynomial p1, Polynomial p2) {
//        Polynomial zero = Polynomial.ZERO(p1.getBasis());
//        Polynomial dividend = new Polynomial(p1);
//        Polynomial divisor = new Polynomial(p2);
//        while (divisor.degree() > 0 || !divisor.equals(zero)) {
//            Polynomial remainder = dividend.mod(divisor);
//            dividend = divisor;
//            divisor = remainder;
//        }
//        return dividend;
//    }

    /**
     * Given basis q, finds the first primitive element of Fq.
     * @param q integer basis of field Fq
     * @return the first primitive element of field Fq
     */
    public static int findPrimitiveElement(int q) {
        int[] factors = CodeImNotSureIneed.factor(q - 1);
        int alpha = 2;

        while (true) {
            boolean isPrimitive = true;
            for (int factor : factors) {
                int power = (q - 1) / factor;
                int alphaPower = ReedSolomon.powerModQ(alpha, power, q);
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

    /**
     * Given an integer n, returns an array of all of n's factors.
     * @param n integer to find factors for
     * @return an integer array consisting of all of n's factors from 2 to n
     */
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

    /**
     * Returns the interpolation polynomial based on Lagrange Interpolation on the coordinates coords over basis q
     * @param coords coords to use for interpolation
     * @param q basis of the field F
     * @return a new Polynomial
     */
    public static Polynomial interpolate(int[][] coords, int q) {
        List<Polynomial> terms = new LinkedList<>();
        for (int[] coord : coords) {
            terms.add(new Polynomial(new int[]{-coord[0], 1}, q));
        }
        List<Polynomial> lagrangePolynomials = new LinkedList<>();
        for(int i = 0; i < coords.length; i++) {
            Polynomial l = Polynomial.ONE(q);
            for(int j = 0; j < coords.length; j++) {
                if (j == i) continue;
                l = l.multiply(terms.get(j));
            }
            int denominator = l.evaluatePolynomial(coords[i][0]);
        }
        return null;
    }

    // need to check if i have a better implementation of the following code:
//    /**
//     * Calculates the remainder of this polynomial modulo a given divisor polynomial.
//     * @param divisor the polynomial to use as divisor in the modulo operation
//     * @return the remainder of this polynomial modulo the given divisor polynomial
//
//     * @pre this.getBasis() == divisor.getBasis();
//     * @post result.getBasis() == this.getBasis();
//     * @post result.degree() < divisor.degree();
//     */
//    public Polynomial mod(Polynomial divisor) throws IllegalArgumentException {
//        if (this.degree() < divisor.degree()) {
//            return this;
//        }
//
//        if (divisor.degree() == 0 && divisor.getCoefficient(0) == 0)
//            throw new IllegalArgumentException("Divisor can't be zero!!!!");
//
//        int q = this.basis;
//        int degreeDiff = this.degree() - divisor.degree();
//        int[] resultCoeffs = Arrays.copyOf(this.coefficients, this.coefficients.length);
//
//        for (int i = degreeDiff; i >= 0; i--) {
//            int ratio = resultCoeffs[i + divisor.degree()] / divisor.coefficients[divisor.degree()];
//            for (int j = 0; j <= divisor.degree(); j++) {
//                resultCoeffs[i + j] = (resultCoeffs[i + j] - ratio * divisor.coefficients[j] + q) % q;
//            }
//        }
//
//        // ensure result is over basis q
//        int resultDegree = this.degree();
//        while (resultDegree > 0 && resultCoeffs[resultDegree] == 0) {
//            resultDegree--;
//        }
//        return new Polynomial(Arrays.copyOf(resultCoeffs, resultDegree + 1), this.basis);
//    }
}