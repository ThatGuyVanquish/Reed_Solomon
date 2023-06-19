package Code;

import java.util.ArrayList;
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
//        GaloisField q = symbols.getField();
//        int primitive = ReedSolomon.findPrimitiveElement(q);
        int numOfCoords = symbolCoeffs.length - errorIndices.size();
        int[][] coords = new int[numOfCoords][];

        int indexOfError = 0, indexOfSymbol = 0;
        for (int i = 0; i < symbolCoeffs.length; i++) {
            if (indexOfError < errorIndices.size() && i == errorIndices.get(indexOfError)) {
                indexOfError++;
                continue;
            }
            coords[indexOfSymbol] = new int[]{i, symbolCoeffs[i]};
            indexOfSymbol++;
        }
        return coords;
    }

    public static Polynomial lagrangeInterpolation(int[][] coordinates, GaloisField F) {
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
        return lagrange;
    }

    /**
     * Given a polynomial of values, returns a list of polynomials obtained by multiplying
     * (x - values_i) for every i
     * @param values
     * @return list of Polynomials of the shape (X - values_0) * (X - values_1) * ... * (X - values_n)
     */
    public static List<Polynomial> getLagrangeTerms(Polynomial values)
    {
        int[] vArr = values.getCoefficients();
        GaloisField GF = values.getField();
        List<Polynomial> terms = new ArrayList<>();
        for(int i = 0; i < vArr.length; i++) {
            Polynomial term = Polynomial.ONE(GF);
            int denominator = 1;
            for(int j = 0; j < vArr.length; j++) {
                if (j == i) continue;
                Polynomial currentTerm = new Polynomial(new int[]{-vArr[j], 1}, GF);
                term = term.multiply(currentTerm);
                denominator *= (vArr[i] - vArr[j]);
            }
            Polynomial denomPoly = new Polynomial(new int[]{denominator}, GF);
            terms.add(term.div(denomPoly));
        }
        return terms;
    }

    /**
     * Given two vectors of alphas and f(alpha_i),
     * returns the bivariate Polynomial Q of (1, k)-weighted degree at most D
     * for which Q(alpha_i, y_i) = 0 for each i in [1, n]
     * @param xValues x values for which to calculate the interpolated polynomial
     * @param yValues y values for which to calculate the interpolated polynomial
     * @param f a bivariate polynomial at which to evaluate the coefficient for the lagrange terms
     * @param D maximum (1,k) weighted degree of result polynomial
     * @param k length of original message
     * @return the interpolated bivariate Polynomial Q
     */
    public static BivariatePolynomial interpolateBivariatePolynomial(Polynomial xValues, Polynomial yValues,
                                                                     BivariatePolynomial f, int D, int k) {
        GaloisField GF = xValues.getField();
        // Generate the lagrange terms
        List<Polynomial> xTerms = getLagrangeTerms(xValues), yTerms = getLagrangeTerms(yValues);
        List<BivariatePolynomial> xTermsBPs = new ArrayList<>(),
                yTermsBPs = new ArrayList<>();
        for(Polynomial p : xTerms) {
            xTermsBPs.add(BivariatePolynomial.univariateToBivariate(p, 'x'));
        }
        for(Polynomial p : yTerms) {
            yTermsBPs.add(BivariatePolynomial.univariateToBivariate(p, 'y'));
        }
        List<BivariatePolynomial> coeffs = new ArrayList<>();
        int[] alphasArr = xValues.getCoefficients();
        int[] yArr = yValues.getCoefficients();
        // Evaluate f at (alpha_i, received_word_i) for every i
        for(int i = 0; i < alphasArr.length; i++) {
            int coeff = f.evaluatePolynomial(alphasArr[i], yArr[i]);
            BivariatePolynomial coeffBP = BivariatePolynomial.univariateToBivariate(new Polynomial(new int[]{coeff}, GF), 'x');
            coeffs.add(coeffBP);
        }
        BivariatePolynomial Q = BivariatePolynomial.ZERO(GF);
        int index = 0;
        for(BivariatePolynomial xTerm : xTermsBPs) {
            BivariatePolynomial bp = BivariatePolynomial.ZERO(GF);
            for(BivariatePolynomial yTerm : yTermsBPs) {
                bp.add(yTerm);
            }
            BivariatePolynomial coeff = coeffs.get(index);
            index++;
            Q.add((xTerm.multiply(bp)).multiply(coeff));
        }
        if (Q.weightedDegree(k) <= D)
            return Q;

        return null;
    }
}
