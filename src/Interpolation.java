import java.util.Arrays;
import java.util.List;

public class Interpolation {
    /**
     * Given the Polynomial of encoded symbols and a list of the error indices, returns an array of arrays of length 2
     * such that each sub-array is the x,y coordinates that'll be used to interpolate the original message
     * polynomial
     * @param symbols Polynomial whose coefficients are the encoded symbols
     * @param errorIndices List of indices such that g(a^i) != symbol[i]
     * @return the coordinates to use to interpolate the original message
     */
    public static int[][] getInterpolationCoordinates(Polynomial symbols, List<Integer> errorIndices) {
        int[] symbolCoeffs = symbols.getCoefficients();
        int q = symbols.getBasis();
        int primitive = ReedSolomon.findPrimitiveElement(q);
        int numOfCoords = symbolCoeffs.length - errorIndices.size();
        int[][] coords = new int[numOfCoords][];

        int indexOfError = 0, indexOfSymbol = 0;
        for (int i = 0; i < symbolCoeffs.length; i++) {
            if (i == errorIndices.get(indexOfError)) {
                indexOfError++;
                continue;
            }
            coords[indexOfSymbol] = new int[]{ReedSolomon.powerModQ(primitive, i, q), symbolCoeffs[i]};
            indexOfSymbol++;
        }
        return coords;
    }

    public static int[] lagrangeInterpolation(int[][] coordinates, int q) {
        int n = coordinates.length;
        int[] lagrangeFactors = new int[n];
        Arrays.fill(lagrangeFactors, 1);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    int numerator = q - coordinates[j][0];
                    int denominator = coordinates[i][0] - coordinates[j][0];
                    int inverseDenominator = multiplicativeInverse(denominator, q);
                    lagrangeFactors[i] = multiplyModQ(lagrangeFactors[i], multiplyModQ(numerator, inverseDenominator, q), q);
                }
            }
            lagrangeFactors[i] = multiplyModQ(coordinates[i][1], lagrangeFactors[i], q);
        }
        return lagrangeFactors;
    }

    /**
     * Given two numbers x,y and a basis q, returns the positive value result of (x*y)%q
     * @param x factor
     * @param y multiplier
     * @param q basis of the field Fq
     * @return result of (x*y)%q
     */
    private static int multiplyModQ(int x, int y, int q) {
        return (int) (((long) x * y) % q);
    }

    /**
     * Given two numbers x,y and a basis q, returns the positive value result of (x-y)%q
     * @param x subtrahend
     * @param y subtracted
     * @param q basis of the field Fq
     * @return result of (x-y)%q
     */
    public static int subtractModQ(int x, int y, int q) {
        return ((x - y) % q + q) % q;
    }

    public static int multiplicativeInverse(int a, int q) {
        int t = 0, newt = 1;
        int r = q, newr = a;
        while (newr != 0) {
            int quotient = r / newr;
            int temp = newt;
            newt = subtractModQ(t, multiplyModQ(quotient, newt, q), q);
            t = temp;
            temp = newr;
            newr = subtractModQ(r, multiplyModQ(quotient, newr, q), q);
            r = temp;
        }
        if (r > 1) {
            throw new ArithmeticException("a is not invertible");
        }
        if (t < 0) {
            t += q;
        }
        return t;
    }

    public static int divide(int dividend, int divisor, int q) {
        int inverseOfDivisor = multiplicativeInverse(divisor, q);
        return (dividend * inverseOfDivisor) % q;
    }
}
