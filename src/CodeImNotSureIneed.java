import Code.*;

import java.util.*;

public class CodeImNotSureIneed {
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
        Polynomial lagrangeOfSymbols = Interpolation.lagrangeInterpolation(coordsOfSymbols, F);
//        System.out.println("LAGRANGE OF SYMBOLS: " + lagrangeOfSymbols);
        int[] originalMessageCoeffs = new int[k];
        for(int i = 0; i < k; i++) {
            originalMessageCoeffs[i] = lagrangeOfSymbols.evaluatePolynomial(i);
        }
        return new Polynomial(originalMessageCoeffs, F);
    }

    /**
     * Given the basis q, desired length of encryption n and original message length k, computes the generator
     * polynomial for Fq over (n,k).
     * @param F the Galois Field to calculate over
     * @param n the length of the encoded message
     * @param k the length of the original message
     * @return the generator polynomial of field Fq over (n,k)
     */
    public static Polynomial computeGeneratorPolynomial(GaloisField F, int n, int k) {
        int[] alphaPowers = new int[n - k];
//        int alpha = ReedSolomon.findPrimitiveElement(F);
        int alpha = findPrimitiveElement(F);

        for (int i = 0; i < alphaPowers.length; i++) {
//            alphaPowers[i] = ReedSolomon.powerModQ(alpha, i + 1, F);
            alphaPowers[i] = powerModQ(alpha, i + 1, F);
        }

        return findPolynomialFromRoots(alphaPowers, F);
    }

    /**
     * Given the to-be roots of the returned polynomial and the basis of the field Fq,
     * computes the polynomial generated by multiplying all the terms (x-root).
     * @param roots the to-be roots of the returned polynomial
     * @param F the Galois Field to calculate over
     * @return an integer array which holds the coefficient for a polynomial from Fq[X] whose roots are given.
     */
    public static Polynomial findPolynomialFromRoots(int[] roots, GaloisField F) {
        Polynomial polynomial = Polynomial.ONE(F);

        for (int root : roots) {
            Polynomial termByRoot = new Polynomial(new int[]{-root, 1}, F);
            polynomial = polynomial.multiply(termByRoot);
        }

        return polynomial;
    }

    /**
     * Given the base, the exponent to power base by and the basis q of Fq, computes (base^exponent) modulo q.
     * @param base base to be powered
     * @param exponent exponent to power base by
     * @param F GaloisField to operate over
     * @return value of base^exponent % mod
     */
    public static int powerModQ(int base, int exponent, GaloisField F) {
        int mod = F.getPrime();
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

    public static int findPrimitiveElement(GaloisField F) {
        int q = F.getPrime();
        for(int i = 2; i < q; i++) {
            boolean[] nonZeroElements = new boolean[q];
            nonZeroElements[0] = true;
            nonZeroElements[i] = true;
            boolean found = true;
            for(int j = 2; j < q; j++) {
                int pow = (int)Math.pow(i, j) % q;
                if (nonZeroElements[pow]) {
                    found = false;
                    break;
                }
                nonZeroElements[pow] = true;
            }
            if (found)
                return i;
        }
        return -1;
    }

}
