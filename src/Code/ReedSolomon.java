package Code;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class ReedSolomon {

    /**
     * Given a message polynomial and desired length of encryption, uses Reed-Solomon to encrypt the message and
     * generate the encoded message, a list of encoded symbols and a generator polynomial for Fq.
     * @param msg   message to be encoded, a polynomial over Fp for a prime p
     * @param n desired length of the encoded message
     * @return a list of Polynomials of length 4:
     * #0 => Encoded message
     * #1 => Encoded symbols
     * #2 => k, Length of original message
     * #3 => Generator polynomial
     * @throws IllegalArgumentException if basis q of message has no irreducible polynomials
     */
    public static List<Polynomial> RSEncoder(Polynomial msg, int n) throws IllegalArgumentException {
        int k = msg.degree() + 1;
        GaloisField F = msg.getField();

        int q = F.getPrime();
        Polynomial generatorPolynomial = ReedSolomon.computeGeneratorPolynomial(F, n, k);
        int alpha = ReedSolomon.findPrimitiveElement(F);

        int[] symbolsArr = new int[n];
        for (int i = 0; i < n; i++) {
            symbolsArr[i] = msg.evaluatePolynomial(ReedSolomon.powerModQ(alpha, i + 1, F));
        }

        Polynomial encodedMsg = msg.multiply(generatorPolynomial);
        Polynomial encodedSymbols = new Polynomial(symbolsArr, F);
        Polynomial constantK = new Polynomial(new int[]{k}, F);

        List<Polynomial> res = new LinkedList<>();
        res.add(encodedMsg);
        res.add(encodedSymbols);
        res.add(constantK);
        res.add(generatorPolynomial);

        return res;
    }

    public static List<Polynomial> RSEncoder_L(Polynomial msg, int n) throws IllegalArgumentException {
        int k = msg.degree() + 1;
        GaloisField F = msg.getField();

        Polynomial generatorPolynomial = ReedSolomon.computeGeneratorPolynomial(F, n, k);

        int[] lagrangeCoeffs = Interpolation.lagrangeInterpolation(Interpolation.getInterpolationCoordinates(msg, new LinkedList<>()), F);
        Polynomial L = new Polynomial(lagrangeCoeffs, F);

        int[] symbolsArr = new int[n];
        for (int i = 0; i < n; i++) {
            symbolsArr[i] = L.evaluatePolynomial(i);
        }

        Polynomial encodedMsg = msg.multiply(generatorPolynomial);
        Polynomial encodedSymbols = new Polynomial(symbolsArr, F);
        Polynomial constantK = new Polynomial(new int[]{k}, F);

        List<Polynomial> res = new LinkedList<>();
        res.add(encodedMsg);
        res.add(encodedSymbols);
        res.add(constantK);
        res.add(generatorPolynomial);

        return res;
    }

    /**
     * Given a polynomial of encoded symbols and the original message length k
     * decodes the original message of length k using Reed-Solomon unique decoding algorithm of Berlekamp-Welch.
     * @param symbols Encoded symbols polynomial
     * @return the original message polynomial if it can be decoded, null otherwise
     */
    public static Polynomial uniqueDecoder(Polynomial symbols, int k) {
        int n = symbols.degree() + 1;
        GaloisField F = symbols.getField();
        int q = F.getPrime();

        int maxNumOfErrors = (n - k) / 2;

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
        int currentNumOfErrors = maxNumOfErrors;
        for(;currentNumOfErrors > 0; currentNumOfErrors--) {
            int[] errorCoeffs = Arrays.copyOfRange(values, 0, currentNumOfErrors);
            errorCoeffs = Arrays.copyOf(errorCoeffs, currentNumOfErrors + 1);
            errorCoeffs[currentNumOfErrors] = 1;
            Polynomial E = new Polynomial(errorCoeffs, F);

            int[] Qcoeffs = Arrays.copyOfRange(values, currentNumOfErrors, values.length);
            Polynomial Q = new Polynomial(Qcoeffs, F);

            System.out.println("Q: " + Q + "\nE: " + E);

            if (Q.mod(E).equals(Polynomial.ZERO(F))) {
                System.out.println("num of errors is " + currentNumOfErrors);
                return Q.div(E);
            }
        }
        // Can't correct errors, return null
        return null;
    }

    public static Polynomial uniqueDecoder_L(Polynomial symbols, int k) {
        int n = symbols.degree() + 1;
        GaloisField F = symbols.getField();
        int q = F.getPrime();

        int maxNumOfErrors = (n - k) / 2;

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
        int currentNumOfErrors = maxNumOfErrors;
        Polynomial Q, E = null, M = null;

        for(;currentNumOfErrors > 0; currentNumOfErrors--) {
            int[] errorCoeffs = Arrays.copyOfRange(values, 0, currentNumOfErrors);
            errorCoeffs = Arrays.copyOf(errorCoeffs, currentNumOfErrors + 1);
            errorCoeffs[currentNumOfErrors] = 1;
            E = new Polynomial(errorCoeffs, F);

            int[] Qcoeffs = Arrays.copyOfRange(values, currentNumOfErrors, values.length);
            Q = new Polynomial(Qcoeffs, F);

            if (Q.mod(E).equals(Polynomial.ZERO(F))) {
                M = Q.div(E);
                break;
            }
        }
        Polynomial correctedSymbols;
        List<Integer> errorIndices = new LinkedList<>();
        if (currentNumOfErrors > 0) {
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
        int[][] coordsOfSymbols = Interpolation.getInterpolationCoordinates(correctedSymbols, new LinkedList<>());
        int[] lagrangeCoeffsOfSymbols = Interpolation.lagrangeInterpolation(coordsOfSymbols, F);
        Polynomial lagrangeOfSymbols = new Polynomial(lagrangeCoeffsOfSymbols, F);
        int[] originalMessageCoeffs = new int[k];
        for(int i = 0; i < k; i++) {
            originalMessageCoeffs[i] = lagrangeOfSymbols.evaluatePolynomial(i);
        }
        return new Polynomial(originalMessageCoeffs, F);
    }


    public static String printMatrix(int[][] mat) {
        StringBuilder res = new StringBuilder();
        for (int[] ints : mat) {
            res.append(Arrays.toString(ints)).append("\n");
        }
        return res.toString();
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
        int alpha = ReedSolomon.findPrimitiveElement(F);

        for (int i = 0; i < alphaPowers.length; i++) {
            alphaPowers[i] = ReedSolomon.powerModQ(alpha, i + 1, F);
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