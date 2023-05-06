package Code;

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
        GaloisField q = symbols.getField();
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

    public static int[] lagrangeInterpolation(int[][] coordinates, GaloisField F) {
        int n = coordinates.length;
        int[] lagrangeFactors = new int[n];
        Arrays.fill(lagrangeFactors, 1);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    int numerator = F.getPrime() - coordinates[j][0];
                    int denominator = coordinates[i][0] - coordinates[j][0];
                    int inverseDenominator = F.modInverse(denominator);
                    lagrangeFactors[i] = F.multiply(lagrangeFactors[i], F.multiply(numerator, inverseDenominator));
                }
            }
            lagrangeFactors[i] = F.multiply(coordinates[i][1], lagrangeFactors[i]);
        }
        return lagrangeFactors;
    }
}
