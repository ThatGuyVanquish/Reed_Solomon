import Code.GaloisField;
import Code.Interpolation;
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

    public static Polynomial uniqueDecoder_e(Polynomial symbols, int k, int e) {
        int n = symbols.degree() + 1;
        GaloisField F = symbols.getField();
        int q = F.getPrime();

        int[][] equations = new int[n][n];
        int[] values;
        // Assuming there are at most maxNumOfErrors errors, values[0]...values[maxNumOfErrors - 1] are the error
        // variables, and the rest are the "correct" values variables
        int[] result = new int[n];
        // First calculate the results
        for(int i = 0; i < n; i++) {
            result[i] = F.mod(((int)Math.pow(i, 2) * (-1) * symbols.getCoefficient(i)));
        }
        // Generate the linear equations coefficients based on the Berlekamp-Welch algorithm
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                switch (j) {
                    case 0 -> equations[i][j] = symbols.getCoefficient(i);
                    case 1 -> equations[i][j] = (symbols.getCoefficient(i) * i) % q;
                    case 2 -> equations[i][j] = q - 1;
                    default -> equations[i][j] = (q - ((int) Math.pow(i, j - 2) % q)) % q;
                }
            }
        }

        values = F.gaussianElimination(equations, result);
        int[] errorCoeffs = Arrays.copyOfRange(values, 0, e);
        errorCoeffs = Arrays.copyOf(errorCoeffs, e + 1);
        errorCoeffs[e] = 1;
        Polynomial E = new Polynomial(errorCoeffs, F);
        int[] Qcoeffs = Arrays.copyOfRange(values, e, values.length);
        Polynomial Q = new Polynomial(Qcoeffs, F);
        if (Q.mod(E).equals(Polynomial.ZERO(F))) {
            return Q.div(E);
        }
        // Can't correct errors, return null
        return null;
    }

    public static Polynomial uniqueDecoder_Le(Polynomial symbols, int k, int e) {
        int n = symbols.degree() + 1;
        GaloisField F = symbols.getField();
        int q = F.getPrime();

        int[][] equations = new int[n][n];
        int[] values;
        // Assuming there are at most maxNumOfErrors errors, values[0]...values[maxNumOfErrors - 1] are the error
        // variables, and the rest are the "correct" values variables
        int[] result = new int[n];
        // First calculate the results
        for(int i = 0; i < n; i++) {
            result[i] = F.mod(((int)Math.pow(i, 2) * (-1) * symbols.getCoefficient(i)));
        }
        // Generate the linear equations coefficients based on the Berlekamp-Welch algorithm
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                switch (j) {
                    case 0 -> equations[i][j] = symbols.getCoefficient(i);
                    case 1 -> equations[i][j] = (symbols.getCoefficient(i) * i) % q;
                    case 2 -> equations[i][j] = q - 1;
                    default -> equations[i][j] = (q - ((int) Math.pow(i, j - 2) % q)) % q;
                }
            }
        }

        values = F.gaussianElimination(equations, result);
        Polynomial Q, E, M = null;
        int[] errorCoeffs = Arrays.copyOfRange(values, 0, e);
        errorCoeffs = Arrays.copyOf(errorCoeffs, e + 1);
        errorCoeffs[e] = 1;
        E = new Polynomial(errorCoeffs, F);

        int[] Qcoeffs = Arrays.copyOfRange(values, e, values.length);
        Q = new Polynomial(Qcoeffs, F);

        if (Q.mod(E).equals(Polynomial.ZERO(F))) {
            M = Q.div(E);
        }
        if (M == null)
            return null;

        Polynomial correctedSymbols;
        List<Integer> errorIndices = new LinkedList<>();
        if (e > 0) {
            for(int i = 0; i < n; i++) {
                if (E.evaluatePolynomial(i) == 0)
                    errorIndices.add(i);
            }
            int[] symbolsArr = symbols.getCoefficients();
            for(int index : errorIndices) {
                symbolsArr[index] = M.evaluatePolynomial(index);
            }
            correctedSymbols = new Polynomial(symbolsArr, F);
        }
        else correctedSymbols = symbols;
        int[][] coordsOfSymbols = Interpolation.getInterpolationCoordinates(symbols, errorIndices);
//        System.out.println("Error indices: " + errorIndices);
//        System.out.println("coords of symbols:\n" + Arrays.deepToString(coordsOfSymbols));
        int[] lagrangeCoeffsOfSymbols = Interpolation.lagrangeInterpolation(coordsOfSymbols, F);
        Polynomial lagrangeOfSymbols = new Polynomial(lagrangeCoeffsOfSymbols, F);
//        System.out.println("LAGRANGE OF SYMBOLS: " + lagrangeOfSymbols);
        int[] originalMessageCoeffs = new int[k];
        for(int i = 0; i < k; i++) {
            originalMessageCoeffs[i] = lagrangeOfSymbols.evaluatePolynomial(i);
        }
        return new Polynomial(originalMessageCoeffs, F);
    }

}
