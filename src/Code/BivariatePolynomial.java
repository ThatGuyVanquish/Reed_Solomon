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
        x = F.mod(x);
        y = F.mod(y);
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
        // Reverse the order of the keys in the `terms` map.
        List<Integer> exponentXList = new ArrayList<>(this.terms.keySet());
        Collections.reverse(exponentXList);

        for(Integer exponentX : exponentXList) {
            List<Integer> exponentYList = new ArrayList<>(this.terms.get(exponentX).keySet());
            Collections.reverse(exponentYList);
            for(Integer exponentY : exponentYList) {
                int coefficient = this.getCoefficient(exponentX, exponentY);
                if (coefficient != 0) {
                    if (coefficient == 1 && exponentX == 0 && exponentY == 0) {
                        res.append("1");
                    }
                    if (coefficient > 1)
                        res.append(coefficient);
                    if (exponentX != 0) {
                        if (exponentX == 1)
                            res.append("x");
                        else {
                            res.append("(x^");
                            res.append(exponentX);
                            res.append(")");
                        }
                    }
                    if (exponentY != 0) {
                        if (exponentY == 1)
                            res.append("y");
                        else {
                            res.append("(y^");
                            res.append(exponentY);
                            res.append(")");
                        }
                    }
                    res.append(" + ");
                }
                else res.append("0");
            }
        }
        // Remove the last " + " character.
        if (res.length() > 3) {
            res.delete(res.length() - 3, res.length());
        }
        return res.toString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        BivariatePolynomial that = (BivariatePolynomial) o;

        if (this.terms.size() != that.terms.size()) return false;
        for (Integer exponentX : this.terms.keySet()) {
            for (Integer exponentY : this.terms.get(exponentX).keySet()) {
                if (this.getCoefficient(exponentX, exponentY) != that.getCoefficient(exponentX, exponentY)) {
                    return false;
                }
            }
        }

        return true;
    }


    /**
     * Given a Bivariate Polynomial other, returns the result of this - other.
     * @param other Bivariate Polynomial to be subtracted from this.
     * @return a new Bivariate Polynomial with the result of this - other.
     */
    public BivariatePolynomial subtract(BivariatePolynomial other) {
        Map<Integer, Map<Integer, Integer>> newTerms = new HashMap<>();
        // Subtract the terms of other from this if they have the same exponents.
        for(Integer exponentX : this.terms.keySet()) {
            Map<Integer, Integer> newTerm = new HashMap<>();
            for(Integer exponentY : this.terms.get(exponentX).keySet()) {
                int coefficient = F.mod(this.getCoefficient(exponentX, exponentY) - other.getCoefficient(exponentX, exponentY));
                if (coefficient != 0)
                    newTerm.put(exponentY, coefficient);
            }
            newTerms.put(exponentX, newTerm);
        }

        /* Some terms of other exist whose exponents are not in this.
           For example, this = 2x^2*y^3 + 4xy, other = 6x^2*y^3 + 3x^2*y + 4
           The first for loop above will return -4x^2*y^3, we're missing the rest of the terms of other.*/
        for(Integer exponentX : other.terms.keySet()) {
            Map<Integer, Integer> newTerm = new HashMap<>();
            if (!this.terms.containsKey(exponentX)) { // Other contains an x exponent value that this does not.
                for(Integer exponentY : other.terms.get(exponentX).keySet()) {
                    int coefficient = F.mod(-other.getCoefficient(exponentX, exponentY));
                    newTerm.put(exponentY,coefficient);
                }
            }
            else { // Found an x exponent value that this and other share, need to subtract missing y exponent coefficients.
                for(Integer exponentY : other.terms.get(exponentX).keySet()) {
                    if (this.terms.get(exponentX).containsKey(exponentY)) // Already subtracted in the first for loop.
                        continue;
                    int coefficient = F.mod(-other.getCoefficient(exponentX, exponentY));
                    newTerm.put(exponentY,coefficient);
                }
            }
            newTerms.put(exponentX, newTerm);
        }

        return new BivariatePolynomial(newTerms, this.F);
    }

    /**
     * Given a bivariate polynomial, returns the result of adding this polynomial to the given polynomial.
     * @param other given polynomial to add to this
     * @return Bivariate Polynomial result of this+other
     * @pre this.getBasis() == other.getBasis()
     * @post this.getBasis() == result.getBasis()
     * @post result.weightedDegree() <= Math.max(this.weightedDegree(), other.weightedDegree())
     */
    public BivariatePolynomial add(BivariatePolynomial other) {
        Map<Integer, Map<Integer, Integer>> newTerms = new HashMap<>();
        // Add the terms of other from this if they have th-e same exponents.
        for(Integer exponentX : this.terms.keySet()) {
            Map<Integer, Integer> newTerm = new HashMap<>();
            for(Integer exponentY : this.terms.get(exponentX).keySet()) {
                int coefficient = F.mod(this.getCoefficient(exponentX, exponentY) + other.getCoefficient(exponentX, exponentY));
                if (coefficient != 0)
                    newTerm.put(exponentY, coefficient);
            }
            newTerms.put(exponentX, newTerm);
        }

        /* Some terms of other exist whose exponents are not in this.
           For example, this = 2x^2*y^3 + 4xy, other = 6x^2*y^3 + 3x^2*y + 4
           The first for loop above will return 8x^2*y^3, we're missing the rest of the terms of other.*/
        for(Integer exponentX : other.terms.keySet()) {
            Map<Integer, Integer> newTerm = new HashMap<>();
            if (!this.terms.containsKey(exponentX)) { // Other contains an x exponent value that this does not.
                for(Integer exponentY : other.terms.get(exponentX).keySet()) {
                    int coefficient = F.mod(other.getCoefficient(exponentX, exponentY));
                    newTerm.put(exponentY,coefficient);
                }
            }
            else { // Found an x exponent value that this and other share, need to add missing y exponent coefficients.
                for(Integer exponentY : other.terms.get(exponentX).keySet()) {
                    if (this.terms.get(exponentX).containsKey(exponentY)) // Already added in the first for loop.
                        continue;
                    int coefficient = F.mod(other.getCoefficient(exponentX, exponentY));
                    newTerm.put(exponentY,coefficient);
                }
            }
            newTerms.put(exponentX, newTerm);
        }

        return new BivariatePolynomial(newTerms, this.F);
    }

    /**
     * Returns a new polynomial that is the multiplication of this polynomial and the given polynomial.
     * @param other the polynomial to multiply this polynomial by
     * @return a new polynomial that is the result of multiplication of this polynomial and the given polynomial

     * @pre this.F == other.getField();
     * @post result.getField() == this.F;
     */
    public BivariatePolynomial multiply(BivariatePolynomial other) {
        Map<Integer, Map<Integer, Integer>> newTerms = new HashMap<>();
        for(Integer exponentX : this.terms.keySet()) {
            for(Integer exponentY : this.terms.get(exponentX).keySet()) {
                for(Integer exponentXOther : other.terms.keySet()) {
                    for(Integer exponentYOther : other.terms.get(exponentXOther).keySet()) {
                        int coefficient = F.mod(this.getCoefficient(exponentX, exponentY) * other.getCoefficient(exponentXOther, exponentYOther));
                        if (coefficient != 0) {
                            newTerms.computeIfAbsent(exponentX + exponentXOther, k -> new HashMap<>()).put(exponentY + exponentYOther, coefficient);
                        }
                    }
                }
            }
        }
        return new BivariatePolynomial(newTerms, this.F);
    }

    public static BivariatePolynomial univariateToBivariate(Polynomial p, char c) {
        return c == 'x' ? univariateToBivariateX(p) :
                c == 'y' ? univariateToBivariateY(p) :
                        null;
    }
    public static BivariatePolynomial univariateToBivariateX(Polynomial p) {
        Map<Integer, Map<Integer, Integer>> terms = new HashMap<>();
        int[] coeffArr = p.getCoefficients();
        GaloisField GF = p.getField();
        for(int i = 0; i < coeffArr.length; i++) {
            int coeff = coeffArr[i];
            if (coeff == 0) continue;
            Map<Integer, Integer> coeffMap = new HashMap<>();
            coeffMap.put(0, coeff);
            terms.put(i, coeffMap);
        }
        return new BivariatePolynomial(terms, GF);
    }

    public static BivariatePolynomial univariateToBivariateY(Polynomial p) {
        Map<Integer, Map<Integer, Integer>> terms = new HashMap<>();
        int[] coeffArr = p.getCoefficients();
        GaloisField GF = p.getField();
        for(int i = 0; i < coeffArr.length; i++) {
            int coeff = coeffArr[i];
            if (coeff == 0) continue;
            Map<Integer, Integer> coeffMap = new HashMap<>();
            coeffMap.put(i, coeff);
            terms.put(0, coeffMap);
        }
        return new BivariatePolynomial(terms, GF);
    }
}
