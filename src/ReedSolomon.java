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
        int q = msg.getBasis();
        if (nType == null) nType = EncodedLength.SHORT;
        int n = calculateEncodedLength(k, nType);

        Polynomial generatorPolynomial = ReedSolomon.computeGeneratorPolynomial(q, n, k);
        int alpha = ReedSolomon.findPrimitiveElement(q);

        int[] symbolsArr = new int[n];
        for (int i = 0; i < n; i++) {
            symbolsArr[i] = msg.evaluatePolynomial(ReedSolomon.powerModQ(alpha, i + 1, q));
        }

        Polynomial encodedMsg = msg.multiply(generatorPolynomial);
        Polynomial encodedSymbols = new Polynomial(symbolsArr, q);
        Polynomial constantK = new Polynomial(new int[]{k}, q);

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
        int q = ogMsgLength.getBasis();

        int numErrors = (n - k) / 2;

        List<Integer> errorIndices = checkForErrorsInSymbols(symbols, generator);

        if (errorIndices.size() > numErrors) return null; // There are too many errors, can't decode

        int[][] coordsForInterpolation = Interpolation.getInterpolationCoordinates(symbols, errorIndices);

        Polynomial interpolation = new Polynomial(Interpolation.lagrangeInterpolation(coordsForInterpolation, q), q);

        int[] originalMessageCoeffs = new int[k];
        for(int i = 0; i < k; i++)
            originalMessageCoeffs[i] = interpolation.evaluatePolynomial(i);

        return new Polynomial(originalMessageCoeffs, q);
    }

    public static Polynomial uniqueDecoder2(Polynomial symbols, int k) {
        int n = symbols.degree() + 1;
        int q = symbols.getBasis();

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
        GaussianElimination(equations, result, q);


        return null;
    }

    private static int[] GaussianElimination(int[][] equations, int[] result, int q) {
        int N = result.length;
        for (int k = 0; k < N; k++)
        {
            // find pivot row
            int max = k;
            for (int i = k + 1; i < N; i++)
                if (Math.abs(equations[i][k]) > Math.abs(equations[max][k]))
                    max = i;

            // swap row in A matrix
            int[] temp = equations[k];
            equations[k] = equations[max];
            equations[max] = temp;

            // swap corresponding values in constants matrix
            int t = result[k];
            result[k] = result[max];
            result[max] = t;

            // pivot within A and B
            for (int i = k + 1; i < N; i++)
            {
                int factor = Interpolation.divideModQ(equations[i][k], equations[k][k], q);
                result[i] = Interpolation.subtractModQ(result[i], Interpolation.multiplyModQ(factor, result[k], q), q);
                for (int j = k; j < N; j++)
                    equations[i][j] = Interpolation.subtractModQ(equations[i][j], Interpolation.multiplyModQ(factor, equations[k][j], q), q);
            }
        }

        // equations matrix is in row echelon form, reduce it
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                equations[i][j] = Interpolation.divideModQ(equations[i][j], equations[i][i], q);
            }
            result[i] = Interpolation.divideModQ(result[i], equations[i][i], q);
        }
        System.out.println(Arrays.deepToString(equations));

        // now that we have leading 1s, subtract row i from all rows <i
        for(int i = N - 1; i > 0; i--) {
            for(int j = i - 1; j >= 0; j--) {
                int multiplier = equations[j][i];
                result[j] = Interpolation.subtractModQ(result[j], Interpolation.multiplyModQ(result[i], multiplier, q), q);
                equations[j][i] = 0;
            }
        }
        System.out.println(Arrays.deepToString(equations));
        return result;
    }

    public static int[] findUnknows(int[][] mat, int[] sol, int q) {
        GaussianElimination(mat, sol, q);
        int[] unknowns = new int[sol.length];
        for(int i = sol.length - 1; i >= 0; i--) {
            for(int j = i; j < sol.length - 1; j++) {
                unknowns[i] = Interpolation.additionModQ(unknowns[i], Interpolation.divideModQ(sol[j], mat[i][j], q), q);
            }
        }
        return unknowns;
    }

    public static void main(String[] args) {
        Polynomial symbols = new Polynomial(new int[]{1, 6, 123, 456, 57, 86, 121}, 929);
        int q = 929;
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

        System.out.println(Arrays.toString(GaussianElimination(equations, result, 929)));
        //System.out.println(Arrays.toString(findUnknows(equations, result, 929)));
        //System.out.println(Arrays.deepToString(equations));
    }

    /**
     * Given the encoded symbols polynomial, the basis q, the lengths k of the original message and n of the encoded message
     * calculate and return the syndrome polynomial received by evaluating
     * @param encodedSymbols the encoded symbols received by the encoder
     * @param k length of the original message
     * @param q basis for Fq
     * @param n length of the encoded message
     * @return the syndrome Polynomial of the encoded symbols.
     */
    public static Polynomial getSyndromePolynomial(Polynomial encodedSymbols, int q, int k, int n) {
        int alpha = ReedSolomon.findPrimitiveElement(q);
        int t = (n - k) / 2; // Maximum number of errors
        int[] syndromeCoefficients = new int[n - t];
        int index = 0;
        for (int i = n - t; i <= n - 1; i++) {
            syndromeCoefficients[index++] = encodedSymbols.evaluatePolynomial(ReedSolomon.powerModQ(alpha, i, q));
        }
        int accurateLength = n - t;
        for(int i = accurateLength - 1; i > 0; i--) {
            if (syndromeCoefficients[i] == 0)
                accurateLength--;
            else break;
        }
        if (accurateLength == n - t)
            return new Polynomial(syndromeCoefficients, q);
        return new Polynomial(Arrays.copyOf(syndromeCoefficients, accurateLength), q);
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


    /**
     * Given the basis q, desired length of encryption n and original message length k, computes the generator
     * polynomial for Fq over (n,k).
     * @param q the basis of Fq
     * @param n the length of the encoded message
     * @param k the length of the original message
     * @return the generator polynomial of field Fq over (n,k)
     */
    public static Polynomial computeGeneratorPolynomial(int q, int n, int k) {
        int[] alphaPowers = new int[n - k];
        int alpha = ReedSolomon.findPrimitiveElement(q);

        for (int i = 0; i < alphaPowers.length; i++) {
            alphaPowers[i] = ReedSolomon.powerModQ(alpha, i + 1, q);
        }

        return findPolynomialFromRoots(alphaPowers, q);
    }

    /**
     * Given the to-be roots of the returned polynomial and the basis of the field Fq,
     * computes the polynomial generated by multiplying all the terms (x-root).
     * @param roots the to-be roots of the returned polynomial
     * @param q     the basis of the field Fq
     * @return an integer array which holds the coefficient for a polynomial from Fq[X] whose roots are given.
     */
    public static Polynomial findPolynomialFromRoots(int[] roots, int q) {
        Polynomial polynomial = Polynomial.ONE(q);

        for (int root : roots) {
            Polynomial termByRoot = new Polynomial(new int[]{-root, 1}, q);
            polynomial = polynomial.multiply(termByRoot);
        }

        return polynomial;
    }

    /**
     * Given the base, the exponent to power base by and the basis q of Fq, computes (base^exponent) modulo q.
     * @param base base to be powered
     * @param exponent exponent to power base by
     * @param mod basis of F to modulo result
     * @return value of base^exponent % mod
     */
    public static int powerModQ(int base, int exponent, int mod) {
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
        int q = encodedSymbols.getBasis();
        int primitive = findPrimitiveElement(q);
        for (int i = 0; i < encodedSymbols.degree() + 1; i++) {
            if (generatorPolynomial.evaluatePolynomial(ReedSolomon.powerModQ(primitive, i + 1, q)) != encodedSymbols.getCoefficient(i))
                errorIndices.add(i);
        }
        return errorIndices;
    }
}