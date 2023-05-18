package Tests;

import Code.BivariatePolynomial;
import Code.GaloisField;
import Code.Polynomial;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

public class BivariatePolynomialTests {

    BivariatePolynomial bp;
    GaloisField F = new GaloisField(929);
    @BeforeEach
    void setUp() {
        Map<Integer, Map<Integer, Integer>> terms = new HashMap<>();
        Map<Integer, Integer> t1 = new HashMap<>();
        t1.put(0,6);
        t1.put(2,4);
        terms.put(0, t1);
        Map<Integer, Integer> t2 = new HashMap<>();
        t2.put(0, 7);
        t2.put(1, 4);
        t2.put(3, 1);
        terms.put(2, t2);
        Map<Integer, Integer> t3 = new HashMap<>();
        t3.put(1,1);
        terms.put(1, t3);
        bp = new BivariatePolynomial(terms, F);
    }

    @Test
    void getTerms() {
        Map<Integer, Map<Integer, Integer>> terms = bp.getTerms();
        assert(terms.size() == 3);
        assertTrue(terms.containsKey(0));
        assertFalse(terms.containsKey(6));
        assertEquals(terms.get(0).get(2), 4);
        assertEquals(terms.get(2).get(0), 7);
    }

    @Test
    void getCoefficient() {
        assertEquals(bp.getCoefficient(0, 2), 4);
        assertEquals(bp.getCoefficient(2, 3), 1);
        assertEquals(bp.getCoefficient(2, 6), 0);
        assertEquals(bp.getCoefficient(6, 6), 0);
        assertEquals(bp.getCoefficient(-1, 6), 0);
    }

    @Test
    void ZERO() {
        BivariatePolynomial zero = BivariatePolynomial.ZERO(F);
        assertEquals(zero.getCoefficient(0, 2), 0);
        assertEquals(zero.getCoefficient(3, 2), 0);
        assertEquals(zero.getCoefficient(0, 0), 0);
        assertEquals(1, zero.getTerms().size());
        assertEquals(1, zero.getTerms().get(0).size());
    }

    @Test
    void ONE() {
        BivariatePolynomial one = BivariatePolynomial.ONE(F);
        assertEquals(one.getCoefficient(0, 2), 0);
        assertEquals(one.getCoefficient(3, 2), 0);
        assertEquals(one.getCoefficient(0, 0), 1);
        assertEquals(1, one.getTerms().size());
        assertEquals(1, one.getTerms().get(0).size());
    }

    @Test
    void weightedDegree() {

    }

    @Test
    void evaluatePolynomial() {
        assertEquals(bp.evaluatePolynomial(0, 0), 6);
        assertEquals(bp.evaluatePolynomial(1, 0), 13);
        assertEquals(bp.evaluatePolynomial(0, 1), 10);
        assertEquals(bp.evaluatePolynomial(1, 1), 23);
        assertEquals(bp.evaluatePolynomial(-928, 1), 23);
    }

    @Test
    void testEvaluatePolynomial() {
        // Asserting a single variable polynomial represented as a bivariate polynomial evaluates to itself
        Polynomial p1 = new Polynomial(new int[]{0, 1}, F); // x
        Map<Integer, Integer> coeffs = new HashMap<>();
        coeffs.put(0,1);
        Map<Integer, Map<Integer,Integer>> terms = new HashMap<>();
        terms.put(1, coeffs);
        BivariatePolynomial x = new BivariatePolynomial(terms, F);
        assertEquals(x.evaluatePolynomial(p1), p1);

        // Asserting the zero bivariate polynomial asserting to itself with any polynomial.
        assertEquals(BivariatePolynomial.ZERO(F).evaluatePolynomial(p1), Polynomial.ZERO(F));
        assertEquals(BivariatePolynomial.ZERO(F).evaluatePolynomial(Polynomial.ONE(F)), Polynomial.ZERO(F));

        // Asserting a bivariate polynomial evaluation with 0 resulting in only the terms without y.
        Polynomial p2 = new Polynomial(new int[]{6, 0, 7}, F); // 7x^2 + 6
        assertEquals(bp.evaluatePolynomial(Polynomial.ZERO(F)), p2);

        // Evaluating bp with the polynomial x.
        Polynomial bp1 = new Polynomial(new int[]{6,0,12,4,0,1}, F);
        assertEquals(bp.evaluatePolynomial(p1), bp1);
    }

    @Test
    void testToString() {
        assertEquals("0", BivariatePolynomial.ZERO(F).toString());
        assertEquals("1", BivariatePolynomial.ONE(F).toString());
        assertEquals("(x^2)(y^3) + 4(x^2)y + 7(x^2) + xy + 4(y^2) + 6", bp.toString());

        Map<Integer, Integer> coeffs = new HashMap<>();
        coeffs.put(0,1);
        Map<Integer, Map<Integer,Integer>> terms = new HashMap<>();
        terms.put(0, coeffs);
        terms.put(1, coeffs);
        BivariatePolynomial xPlusOne = new BivariatePolynomial(terms, F);
        assertEquals("x + 1", xPlusOne.toString());

        Map<Integer, Integer> coeffs2 = new HashMap<>();
        coeffs2.put(0,5);
        Map<Integer, Map<Integer,Integer>> terms2 = new HashMap<>();
        terms2.put(0, coeffs2);
        BivariatePolynomial five = new BivariatePolynomial(terms2, F);
        assertEquals("5", five.toString());
    }

    @Test
    void testEquals() {
        assertNotEquals(null, bp);
        assertNotEquals(new Object(), bp);

        Map<Integer, Map<Integer, Integer>> terms = new HashMap<>();
        Map<Integer, Integer> t1 = new HashMap<>();
        t1.put(0,6);
        t1.put(2,4);
        terms.put(0, t1);
        Map<Integer, Integer> t2 = new HashMap<>();
        t2.put(0, 7);
        t2.put(1, 4);
        t2.put(3, 1);
        terms.put(2, t2);
        BivariatePolynomial bp2 = new BivariatePolynomial(terms, F);
        assertNotEquals(bp, bp2);

        Map<Integer, Integer> t3 = new HashMap<>();
        t3.put(1,1);
        terms.put(1, t3);
        BivariatePolynomial bp3 = new BivariatePolynomial(terms, F);

        assertEquals(bp, bp3);

    }
}
