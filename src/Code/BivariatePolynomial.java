package Code;

import java.util.*;

public class BivariatePolynomial {

    private Map<Integer, Map<Integer, Integer>> terms;
    private final GaloisField F;

    /**
     * Class to represent Bivariate over Fq for basis q.
     * Note that the terms are represented as:
     * <Key: exponent of X, Value: <Key: exponent of Y, Value: coefficient >>
     * @param terms the terms of the polynomial
     * @param field the field Fq for which this is in Fq[x]
     */
    public BivariatePolynomial(Map<Integer, Map<Integer, Integer>> terms, GaloisField field) {
        this.terms = terms;
        this.F = field;
    }

    public BivariatePolynomial(BivariatePolynomial other) {
        this.terms = other.getTerms();
        this.F = other.getField();
    }

    public static BivariatePolynomial ZERO(GaloisField field) {
        Map<Integer, Map<Integer, Integer>> terms = new HashMap<Integer, Map<Integer, Integer>>();
        Map<Integer, Integer> termZero = new HashMap<Integer, Integer>();
        termZero.put(0, 0);
        terms.put(0, termZero);
        return new BivariatePolynomial(terms, field);
    }

    public static BivariatePolynomial ONE(GaloisField field) {
        Map<Integer, Map<Integer, Integer>> terms = new HashMap<Integer, Map<Integer, Integer>>();
        Map<Integer, Integer> termOne = new HashMap<Integer, Integer>();
        termOne.put(0, 1);
        terms.put(0, termOne);
        return new BivariatePolynomial(terms, field);
    }

    /**
     * Returns the (1,k)-weighted degree of the polynomial.
     * @return the degree of the polynomial
     */
    public int weightedDegree(int k) {
        if (k < 0) return 0;
        int l, j;
        int weightedDegree = 0;
        for(Integer exponentX : this.terms.keySet()) {
            int currentMaxDegree = 0;
            for (Integer exponentY : this.terms.get(exponentX).keySet()) {
                currentMaxDegree = Math.max(currentMaxDegree, exponentX + k * exponentY);
            }
            weightedDegree = Math.max(weightedDegree, currentMaxDegree);
        }
        return weightedDegree;
    }

    /**
     * Returns the field Fq this polynomial is over.
     * @return the field Fq for which this is in Fq[x]
     */
    public GaloisField getField() {
        return this.F;
    }

    /**
     * Returns the coefficients of the polynomial.
     * @return a map of the polynomial's terms.
     */
    public Map<Integer, Map<Integer, Integer>> getTerms() {
        return this.terms;
    }

    /**
     * Returns the coefficient for the requested term modulo q (positive value).
     * @param exponentX the exponent of X
     * @param exponentY the exponent of Y
     * @return the coefficient for the X^exponentX * Y^exponentY term or 0 if such term does not exist.
     */
    public int getCoefficient(int exponentX, int exponentY) {
        if (exponentX < 0 || exponentY < 0)
            return 0;
        if (this.terms.containsKey(exponentX) && this.terms.get(exponentX).containsKey(exponentY))
            return F.mod(this.terms.get(exponentX).get(exponentY));
        return 0;
    }

    /**
     * Returns the result of evaluating this polynomial at point (x,y), modulo q.
     * @param x the value of the point to evaluate this polynomial at.
     * @param y the value of the point to evaluate this polynomial at.
     * @return result of the evaluation of this polynomial at point (x,y).
     */
    public int evaluatePolynomial(int x, int y) {
        int res = 0;
        for(Integer exponentX : this.terms.keySet()) {
            for(Integer exponentY : this.terms.get(exponentX).keySet()) {
                res += this.F.mod(this.getCoefficient(exponentX, exponentY) * (int)Math.pow(x, exponentX) * (int)Math.pow(y, exponentY));
            }
        }
        return this.F.mod(res);
    }

    /**
     * Given a polynomial, returns the single variable polynomial calculated by evaluating this bivariate
     * polynomial at (X, P(X)).
     * @param p polynomial to act as y value.
     * @return new polynomial whose coefficients are obtained by evaluating this bivariate polynomial at (X, P(X)).
     */
    public Polynomial evaluatePolynomial(Polynomial p) {
        int coefficientsLength = this.weightedDegree(p.degree()) + 1;
        int[] newCoefficients = new int[coefficientsLength];
        // Use memoization to avoid recalculation of exponents
        Map<Integer, Polynomial> exponentsOfP = new HashMap<>();
        exponentsOfP.put(0, Polynomial.ONE(this.F));
        exponentsOfP.put(1, p);

        for(Integer exponentX : this.terms.keySet()) {
            for(Integer exponentY : this.terms.get(exponentX).keySet()) {
                // Calculate p^exponentY
                Polynomial pExpY = null;
                if (exponentsOfP.containsKey(exponentY))
                    pExpY = exponentsOfP.get(exponentY);
                else {
                    for(int i = 2; i <= exponentY; i++) {
                        if (exponentsOfP.containsKey(i))
                            continue;
                        if (pExpY == null)
                            pExpY = exponentsOfP.get(i - 1);
                        pExpY = pExpY.multiply(p);
                        exponentsOfP.put(i, pExpY);
                    }
                }
                // Now that we have p^exponentY, calculate coeff(x,y) * x^exponentX * p^exponentY
                int[] xCoeffs = new int[exponentX + 1];
                xCoeffs[exponentX] = this.getCoefficient(exponentX, exponentY);
                Polynomial xExpX = new Polynomial(xCoeffs, this.F);
                Polynomial result = xExpX.multiply(pExpY);
                for(int i = 0; i < result.getCoefficients().length; i++) {
                    newCoefficients[i] = F.mod(newCoefficients[i] + result.getCoefficient(i));
                }
            }
        }
        return new Polynomial(newCoefficients, this.F);
    }

    @Override
    public String toString() {
        StringBuilder res = new StringBuilder();
        for(Integer exponentX : this.terms.keySet()) {
            for(Integer exponentY : this.terms.get(exponentX).keySet()) {
                int coefficient = this.getCoefficient(exponentX, exponentY);
                if (coefficient != 0) {
                    res.append(coefficient);
                    if (exponentX != 0 || exponentY != 0) {
                        res.append("x");
                        if (exponentX > 1) {
                            res.append("^");
                            res.append(exponentX);
                        }
                        if (exponentY > 1) {
                            res.append("y");
                            res.append("^");
                            res.append(exponentY);
                        }
                    }
                    res.append("+");
                }
            }
        }
        if (res.length() > 0) {
            res.deleteCharAt(res.length() - 1);
        }
        return res.toString();
    }


/*
    @Override
    public boolean equals(Object other) {
        if (this == other) return true;
        if (other == null || getClass() != other.getClass()) return false;
        Polynomial that = (Polynomial) other;
        if (!(this.F.equals(that.getField())) || this.degree() != that.degree()) return false;
        int k = this.degree() + 1;
        for(int i = 0; i < k; i++) {
            if (that.getCoefficient(i) != this.getCoefficient(i))
                return false;
        }
        return true;
    }
*/
}
