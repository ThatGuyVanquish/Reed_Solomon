import Code.GaloisField;
import Code.Polynomial;
import Code.ReedSolomon;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class CodeImNotSureIneed {

    /**
     * Given basis q, finds the first primitive element of Fq.
     * @param F the galois field to operate over
     * @return the first primitive element of field Fq
     */
    public static int findPrimitiveElement(GaloisField F) {
        int q = F.getPrime();
        int[] factors = CodeImNotSureIneed.factor(q - 1);
        int alpha = 2;

        while (true) {
            boolean isPrimitive = true;
            for (int factor : factors) {
                int power = (q - 1) / factor;
                int alphaPower = ReedSolomon.powerModQ(alpha, power, F);
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
     * Given the encoded symbols polynomial, the basis q, the lengths k of the original message and n of the encoded message
     * calculate and return the syndrome polynomial received by evaluating
     * @param encodedSymbols the encoded symbols received by the encoder
     * @param k length of the original message
     * @param F the Galois Field to calculate over
     * @param n length of the encoded message
     * @return the syndrome Polynomial of the encoded symbols.
     */
    public static Polynomial getSyndromePolynomial(Polynomial encodedSymbols, GaloisField F, int k, int n) {
        int alpha = ReedSolomon.findPrimitiveElement(F);
        int t = (n - k) / 2; // Maximum number of errors
        int[] syndromeCoefficients = new int[n - t];
        int index = 0;
        for (int i = n - t; i <= n - 1; i++) {
            syndromeCoefficients[index++] = encodedSymbols.evaluatePolynomial(ReedSolomon.powerModQ(alpha, i, F));
        }
        int accurateLength = n - t;
        for(int i = accurateLength - 1; i > 0; i--) {
            if (syndromeCoefficients[i] == 0)
                accurateLength--;
            else break;
        }
        if (accurateLength == n - t)
            return new Polynomial(syndromeCoefficients, F);
        return new Polynomial(Arrays.copyOf(syndromeCoefficients, accurateLength), F);
    }

    /**
     * Given two polynomials, returns their greatest common divisor
     *
     * @param p1 polynomial to be the dividend
     * @param p2 polynomial to be the divisor
     * @return a Polynomial which is the GCD of p1 and p2.
     * @pre p1.getBasis() == p2.getBasis();
     * @post res.degree() <= n-k/2
     */
    public static Polynomial GCD(Polynomial p1, Polynomial p2) {
        Polynomial zero = Polynomial.ZERO(p1.getField());
        Polynomial dividend = new Polynomial(p1);
        Polynomial divisor = new Polynomial(p2);
        while (divisor.degree() > 0 || !divisor.equals(zero)) {
            Polynomial remainder = dividend.mod(divisor);
            dividend = divisor;
            divisor = remainder;
        }
        return dividend;
    }

    /**
     * Given the list of the encoded symbols and the generator polynomial for (q,n,k),
     * checks for errors in the encoded symbols polynomial using the generator polynomial.
     * Evaluates the generator polynomial at every power of the basis q from 1 to n and compares to the
     * @param encodedSymbols Polynomial whose coefficients are the encoded symbols
     * @param generatorPolynomial Polynomial which is a generator polynomial for Fq
     * @return a list of indices in which there is an error
     */
    public static List<Integer> checkForErrorsInSymbols(Polynomial encodedSymbols, Polynomial generatorPolynomial) {
        List<Integer> errorIndices = new LinkedList<>();
        GaloisField F = encodedSymbols.getField();
        int primitive = findPrimitiveElement(F);
        for (int i = 0; i < encodedSymbols.degree() + 1; i++) {
            if (generatorPolynomial.evaluatePolynomial(ReedSolomon.powerModQ(primitive, i + 1, F)) != encodedSymbols.getCoefficient(i))
                errorIndices.add(i);
        }
        return errorIndices;
    }
}
