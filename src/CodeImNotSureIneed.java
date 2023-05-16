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
        int[] lagrangeCoeffsOfSymbols = Interpolation.lagrangeInterpolation(coordsOfSymbols, F);
        Polynomial lagrangeOfSymbols = new Polynomial(lagrangeCoeffsOfSymbols, F);
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

// METHODS FROM BIVARIATE POLYNOMIAL CLASS, NOT SURE THAT ARE NECESSARY BUT I WROTE THEM BEFORE REALIZING THAT :FACEPALM
//    /**
//     * Given a Bivariate Polynomial other, returns the result of this - other.
//     * @param other Bivariate Polynomial to be subtracted from this.
//     * @return a new Bivariate Polynomial with the result of this - other.
//     */
//    public BivariatePolynomial subtract(BivariatePolynomial other) {
//        Map<Integer, Map<Integer, Integer>> newTerms = new HashMap<>();
//        // Subtract the terms of other from this if they have the same exponents.
//        for(Integer exponentX : this.terms.keySet()) {
//            Map<Integer, Integer> newTerm = new HashMap<>();
//            for(Integer exponentY : this.terms.get(exponentX).keySet()) {
//                int coefficient = F.mod(this.getCoefficient(exponentX, exponentY) - other.getCoefficient(exponentX, exponentY));
//                if (coefficient != 0)
//                    newTerm.put(exponentY, coefficient);
//            }
//            newTerms.put(exponentX, newTerm);
//        }
//
//        /* Some terms of other exist whose exponents are not in this.
//           For example, this = 2x^2*y^3 + 4xy, other = 6x^2*y^3 + 3x^2*y + 4
//           The first for loop above will return -4x^2*y^3, we're missing the rest of the terms of other.*/
//        for(Integer exponentX : other.terms.keySet()) {
//            Map<Integer, Integer> newTerm = new HashMap<>();
//            if (!this.terms.containsKey(exponentX)) { // Other contains an x exponent value that this does not.
//                for(Integer exponentY : other.terms.get(exponentX).keySet()) {
//                    int coefficient = F.mod(-other.getCoefficient(exponentX, exponentY));
//                    newTerm.put(exponentY,coefficient);
//                }
//            }
//            else { // Found an x exponent value that this and other share, need to subtract missing y exponent coefficients.
//                for(Integer exponentY : other.terms.get(exponentX).keySet()) {
//                    if (this.terms.get(exponentX).containsKey(exponentY)) // Already subtracted in the first for loop.
//                        continue;
//                    int coefficient = F.mod(-other.getCoefficient(exponentX, exponentY));
//                    newTerm.put(exponentY,coefficient);
//                }
//            }
//            newTerms.put(exponentX, newTerm);
//        }
//
//        return new BivariatePolynomial(newTerms, this.F);
//    }
//
//    /**
//     * Given a bivariate polynomial, returns the result of adding this polynomial to the given polynomial.
//     * @param other given polynomial to add to this
//     * @return Bivariate Polynomial result of this+other
//     * @pre this.getBasis() == other.getBasis()
//     * @post this.getBasis() == result.getBasis()
//     * @post result.weightedDegree() <= Math.max(this.weightedDegree(), other.weightedDegree())
//     */
//    public BivariatePolynomial add(BivariatePolynomial other) {
//        Map<Integer, Map<Integer, Integer>> newTerms = new HashMap<>();
//        // Add the terms of other from this if they have th-e same exponents.
//        for(Integer exponentX : this.terms.keySet()) {
//            Map<Integer, Integer> newTerm = new HashMap<>();
//            for(Integer exponentY : this.terms.get(exponentX).keySet()) {
//                int coefficient = F.mod(this.getCoefficient(exponentX, exponentY) + other.getCoefficient(exponentX, exponentY));
//                if (coefficient != 0)
//                    newTerm.put(exponentY, coefficient);
//            }
//            newTerms.put(exponentX, newTerm);
//        }
//
//        /* Some terms of other exist whose exponents are not in this.
//           For example, this = 2x^2*y^3 + 4xy, other = 6x^2*y^3 + 3x^2*y + 4
//           The first for loop above will return 8x^2*y^3, we're missing the rest of the terms of other.*/
//        for(Integer exponentX : other.terms.keySet()) {
//            Map<Integer, Integer> newTerm = new HashMap<>();
//            if (!this.terms.containsKey(exponentX)) { // Other contains an x exponent value that this does not.
//                for(Integer exponentY : other.terms.get(exponentX).keySet()) {
//                    int coefficient = F.mod(other.getCoefficient(exponentX, exponentY));
//                    newTerm.put(exponentY,coefficient);
//                }
//            }
//            else { // Found an x exponent value that this and other share, need to add missing y exponent coefficients.
//                for(Integer exponentY : other.terms.get(exponentX).keySet()) {
//                    if (this.terms.get(exponentX).containsKey(exponentY)) // Already added in the first for loop.
//                        continue;
//                    int coefficient = F.mod(other.getCoefficient(exponentX, exponentY));
//                    newTerm.put(exponentY,coefficient);
//                }
//            }
//            newTerms.put(exponentX, newTerm);
//        }
//
//        return new BivariatePolynomial(newTerms, this.F);
//    }
//
//    /**
//     * Returns a new polynomial that is the multiplication of this polynomial and the given polynomial.
//     * @param other the polynomial to multiply this polynomial by
//     * @return a new polynomial that is the result of multiplication of this polynomial and the given polynomial
//
//     * @pre this.F == other.getField();
//     * @post result.getField() == this.F;
//     */
//    public BivariatePolynomial multiply(BivariatePolynomial other) {
//        Map<Integer, Map<Integer, Integer>> newTerms = new HashMap<>();
//        for(Integer exponentX : this.terms.keySet()) {
//            for(Integer exponentY : this.terms.get(exponentX).keySet()) {
//                for(Integer exponentXOther : other.terms.keySet()) {
//                    for(Integer exponentYOther : other.terms.get(exponentXOther).keySet()) {
//                        int coefficient = F.mod(this.getCoefficient(exponentX, exponentY) * other.getCoefficient(exponentXOther, exponentYOther));
//                        if (coefficient != 0) {
//                            newTerms.computeIfAbsent(exponentX + exponentXOther, k -> new HashMap<>()).put(exponentY + exponentYOther, coefficient);
//                        }
//                    }
//                }
//            }
//        }
//        return new BivariatePolynomial(newTerms, this.F);
//    }



}
