package Code;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class ReedSolomon {

    /**
     * Given a message polynomial and desired length of encryption, uses Reed-Solomon to encrypt the message and
     * generate the encoded message, a list of encoded symbols and a generator polynomial for Fq.
     * @param msg message to be encoded, a polynomial over Fp for a prime p
     * @param n desired length of the encoded message
     * @return a list of Polynomials of length 2 [4]:
     * #0 => Encoded message /DEPRECATED
     * #1 => Encoded symbols
     * #2 => k, Length of original message
     * #3 => Generator polynomial /DEPRECATED
     * @throws IllegalArgumentException if basis q of message has no irreducible polynomials
     */
    public static List<Polynomial> RSEncoder(Polynomial msg, int n) throws IllegalArgumentException {
        int k = msg.degree() + 1;
        GaloisField F = msg.getField();

//        int q = F.getPrime();
//        Polynomial generatorPolynomial = ReedSolomon.computeGeneratorPolynomial(F, n, k);
//        int alpha = ReedSolomon.findPrimitiveElement(F);

        int[] symbolsArr = new int[n];
        for (int i = 0; i < n; i++) {
//            symbolsArr[i] = msg.evaluatePolynomial(ReedSolomon.powerModQ(alpha, i + 1, F));
            symbolsArr[i] = msg.evaluatePolynomial(i);
        }

//        Polynomial encodedMsg = msg.multiply(generatorPolynomial);
        Polynomial encodedSymbols = new Polynomial(symbolsArr, F);
        Polynomial constantK = new Polynomial(new int[]{k}, F);

        List<Polynomial> res = new LinkedList<>();
//        res.add(encodedMsg);
        res.add(encodedSymbols);
        res.add(constantK);
//        res.add(generatorPolynomial);

        return res;
    }

    /**
     * Given a message polynomial and desired length of encryption, uses Reed-Solomon to encrypt the message and
     * generate the encoded message, a list of encoded symbols and a generator polynomial for Fq.
     * Encodes the original message by interpolating the given polynomial's coefficients over Fq and generating
     * the encoded symbols by evaluating the interpolated polynomial at every i in [0, n).
     * @param msg message to be encoded, a polynomial over Fp for a prime p
     * @param n desired length of the encoded message
     * @return a list of Polynomials of length 2 [4]:
     * #0 => Encoded message /DEPRECATED
     * #1 => Encoded symbols
     * #2 => k, Length of original message
     * #3 => Generator polynomial /DEPRECATED
     * @throws IllegalArgumentException if basis q of message has no irreducible polynomials
     */
    public static List<Polynomial> RSEncoder_L(Polynomial msg, int n) throws IllegalArgumentException {
        int k = msg.degree() + 1;
        GaloisField F = msg.getField();

//        Polynomial generatorPolynomial = ReedSolomon.computeGeneratorPolynomial(F, n, k);

        int[][] coords = Interpolation.getInterpolationCoordinates(msg, new LinkedList<>());
        Polynomial L = Interpolation.lagrangeInterpolation(coords, F);

        int[] symbolsArr = new int[n];
        for (int i = 0; i < n; i++) {
            symbolsArr[i] = L.evaluatePolynomial(i);
        }

//        Polynomial encodedMsg = msg.multiply(generatorPolynomial);
        Polynomial encodedSymbols = new Polynomial(symbolsArr, F);
        Polynomial constantK = new Polynomial(new int[]{k}, F);

        List<Polynomial> res = new LinkedList<>();
//        res.add(encodedMsg);
        res.add(encodedSymbols);
        res.add(constantK);
//        res.add(generatorPolynomial);

        return res;
    }

    /**
     * Given a polynomial of encoded symbols and the original message length k
     * decodes the original message of length k using Reed-Solomon unique decoding algorithm of Berlekamp-Welch.
     * @param symbols Encoded symbols polynomial
     * @param k Length of original message
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

//            System.out.println("Q: " + Q + "\nE: " + E);

            if (Q.mod(E).equals(Polynomial.ZERO(F))) {
                return Q.div(E);
            }
        }
        // Can't correct errors, return null
        return null;
    }

    /**
     * Given a polynomial of encoded symbols and the original message length k
     * decodes the original message of length k using Reed-Solomon unique decoding algorithm of Berlekamp-Welch
     * by interpolation.
     * @param symbols Encoded symbols polynomial
     * @param k Length of original message
     * @return the original message polynomial if it can be decoded, null otherwise
     */
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
     * Given a polynomial of encoded symbols and the original message length k
     * returns a list consisting of close matches to the original message decoded by the Guruswami-Sudan [40] algorithm.
     * @param symbols Encoded code word polynomial (y_i's)
     * @param alphas Alphas at which the original message is evaluated to generate the code word sent
     * @param k Length of original message
     * @return a list of codewords that should contain the original message polynomial
     */
    public static List<Polynomial> listDecoder(Polynomial symbols, Polynomial alphas, int k) {
        int D = (int)Math.sqrt(2 * k * symbols.degree());
        //Polynomial fValues = Polynomial.ONE(symbols.getField());
        //BivariatePolynomial Q = Interpolation.interpolateBivariatePolynomial(symbols, alphas, fValues, D);
        return null;
    }

    public static String printMatrix(int[][] mat) {
        StringBuilder res = new StringBuilder();
        for (int[] ints : mat) {
            res.append(Arrays.toString(ints)).append("\n");
        }
        return res.toString();
    }
}