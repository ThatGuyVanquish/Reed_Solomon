import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class ReedSolomon {

    public enum EncodedLength {
        NOT_ENOUGH,
        SHORT,
        MEDIUM,
        LONG
    }

    private static final int ENCODE_MULTIPLIER = 4;
    private static final double NOT_ENOUGH_MULTIPLIER = 1.5;

    public static int calculateEncodedLength(int k, EncodedLength nType) {
        switch (nType) {
            case NOT_ENOUGH -> {
                return (int) (NOT_ENOUGH_MULTIPLIER * k);
            }
            case SHORT -> {
                return 2 * k;
            }
            case MEDIUM -> {
                return ENCODE_MULTIPLIER * k;
            }
            case LONG -> {
                return 2 * ENCODE_MULTIPLIER * k;
            }
            default -> {
                return k;
            }
        }
    }

    /**
     * Given a message polynomial and desired length of encryption, uses Reed-Solomon to encrypt the message and
     * generate the encoded message, a list of encoded symbols and a generator polynomial for Fq.
     * @param msg   message to be encoded, a polynomial over Fp for a prime p
     * @param nType desired length type of the encoded message, used to get encoded message length
     * @return a list of Polynomials of length 4:
     * #0 => Encoded message
     * #1 => Encoded symbols
     * #2 => k, Length of original message
     * #3 => Generator polynomial
     * @throws IllegalArgumentException if basis q of message has no irreducible polynomials
     */
    public static List<Polynomial> RSEncoder(Polynomial msg, EncodedLength nType) throws IllegalArgumentException {
        int k = msg.degree() + 1;
        GaloisField F = msg.getField();
        if (nType == null) nType = EncodedLength.SHORT;
        int n = calculateEncodedLength(k, nType);

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

    /**
     * Given the encoded message (unnecessary), a polynomial of encoded symbols, the original message length k and the
     * generator polynomial for Fq, decodes the original message of length k using Reed-Solomon unique decoding algorithm.
     * @param encodedMsg Message encoded by RSEncoder
     * @param symbols Encoded symbols polynomial
     * @param ogMsgLength Polynomial of basis q whose coefficient is the value k - the length of the original message
     * @param generator generator polynomial for Fq used by the RSEncoder
     * @return the original message polynomial if it can be decoded, null otherwise
     */
    public static Polynomial uniqueDecoder(Polynomial encodedMsg, Polynomial symbols, Polynomial ogMsgLength, Polynomial generator) {
        int n = symbols.degree() + 1;
        int k = ogMsgLength.getCoefficient(0);
        GaloisField F = ogMsgLength.getField();
        int q = F.getPrime();

        int numErrors = (n - k) / 2;

        List<Integer> errorIndices = checkForErrorsInSymbols(symbols, generator);

        if (errorIndices.size() > numErrors) return null; // There are too many errors, can't decode

        int[][] coordsForInterpolation = Interpolation.getInterpolationCoordinates(symbols, errorIndices);

        Polynomial interpolation = new Polynomial(Interpolation.lagrangeInterpolation(coordsForInterpolation, F), F);

        int[] originalMessageCoeffs = new int[k];
        for(int i = 0; i < k; i++)
            originalMessageCoeffs[i] = interpolation.evaluatePolynomial(i);

        return new Polynomial(originalMessageCoeffs, F);
    }

    public static Polynomial uniqueDecoder2(Polynomial symbols, int k) {
        int n = symbols.degree() + 1;
        GaloisField F = symbols.getField();
        int q = F.getPrime();

        int maxNumOfErrors = (n - k) / 2;

        int[][] equations = new int[n][n];
        int[] values = new int[n];
        // Assuming there are maxNumOfErrors errors, values[0]...values[maxNumOfErrors - 1] are the error
        // variables, and the rest are the "correct" values variables
        int[] result = new int[n];
        // First calculate the results
        for(int i = 0; i < n; i++) {
            result[i] = ((int)Math.pow(i, 2) * (-1) * symbols.getCoefficient(i)) % q;
            if (result[i] < 0) result[i] = q + result[i];
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
        F.gaussianElimination(equations, result);


        return null;
    }

//    private static int[] GaussianElimination(int[][] equations, int[] result, int q) {
//
//        int n = result.length;
//
//
//        return null;
//    }

    public static String printMatrix(int[][] mat) {
        StringBuilder res = new StringBuilder();
        for (int[] ints : mat) {
            res.append(Arrays.toString(ints)).append("\n");
        }
        return res.toString();
    }

    public static int[] findUnknows(int[][] mat, int[] sol, GaloisField F) {
        F.gaussianElimination(mat, sol);
        int[] unknowns = new int[sol.length];
        for(int i = sol.length - 1; i >= 0; i--) {
            for(int j = i; j < sol.length - 1; j++) {
                unknowns[i] = F.add(unknowns[i], F.div(sol[j], mat[i][j]));
            }
        }
        return unknowns;
    }

    public static void main(String[] args) {
        Polynomial symbols = new Polynomial(new int[]{1, 5, 3, 6, 3, 2, 2}, new GaloisField(7));
        int q = 7;
        int[] result = new int[7];
        // First calculate the results
        for(int i = 0; i < 7; i++) {
            result[i] = ((int)Math.pow(i, 2) * (-1) * symbols.getCoefficient(i)) % q;
            if (result[i] < 0) result[i] = q + result[i];
        }
        System.out.println(Arrays.toString(result));

        int[][] equations = new int[7][7];
        for(int i = 0; i < 7; i++) {
            for(int j = 0; j < 7; j++) {
                switch (j) {
                    case 0 -> equations[i][j] = symbols.getCoefficient(i);
                    case 1 -> equations[i][j] = (symbols.getCoefficient(i) * i) % q;
                    case 2 -> equations[i][j] = q - 1;
                    default -> equations[i][j] = (q - ((int) Math.pow(i, j - 2) % q)) % q;
                }
            }
        }
//        System.out.println("MATRIX IS\n" + printMatrix(equations));
//        int[] values = GaussianElimination(equations, result, 7);
//        int[] coeffsForError = {values[0], values[1], 1};
//        int[] coeffsForQ = {values[2], values[3], values[4], values[5], values[6]};
//        Polynomial E = new Polynomial(coeffsForError, 7);
//        Polynomial Q = new Polynomial(coeffsForQ, 7);
//        System.out.println(Arrays.toString(values));
//        System.out.println("E: " + E);
//        System.out.println("Q: " + Q);
//        System.out.println(Q.div(E));

        System.out.println(new GaloisField(7).div(84, 5));

        //System.out.println(Arrays.toString(GaussianElimination(equations, result, 929)));
        //System.out.println(Arrays.toString(new GaloisField(929).gaussianElimination(equations, result)));
        //System.out.println(Arrays.toString(findUnknows(equations, result, 929)));
        //System.out.println(Arrays.deepToString(equations));
//        System.out.println(Interpolation.multiplicativeInverse(167, 7));
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