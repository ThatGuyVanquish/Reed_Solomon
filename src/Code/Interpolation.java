package Code;

import java.util.LinkedList;
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
            if (errorIndices.size() != 0 && i == errorIndices.get(indexOfError)) {
                indexOfError++;
                continue;
            }
            coords[indexOfSymbol] = new int[]{i, symbolCoeffs[i]};
            indexOfSymbol++;
        }
        return coords;
    }

    public static int[] lagrangeInterpolation(int[][] coordinates, GaloisField F) {
        List<Polynomial> lagrangePolynomials = new LinkedList<>();
        for(int i = 0; i < coordinates.length; i++) {
            Polynomial l_i = Polynomial.ONE(F);
            for(int j = 0; j < coordinates.length; j++) {
                if (j == i)
                    continue;
                int negativeXValue = -coordinates[j][0];
                int diff = coordinates[i][0] + negativeXValue;
                Polynomial mul = new Polynomial(new int[]{negativeXValue, 1}, F);
                Polynomial div = new Polynomial(new int[]{diff}, F);
                l_i = l_i.multiply(mul).div(div);
            }
            lagrangePolynomials.add(l_i);
        }
        Polynomial lagrange = Polynomial.ZERO(F);
        for(int i = 0; i < coordinates.length; i++) {
            Polynomial l_i = lagrangePolynomials.get(i);
            Polynomial coeff = new Polynomial(new int[]{coordinates[i][1]}, F);
            lagrange = lagrange.add(l_i.multiply(coeff));
        }
        return lagrange.getCoefficients();
    }
}
