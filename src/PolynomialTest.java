import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.LinkedList;

import static org.junit.jupiter.api.Assertions.*;

class PolynomialTest {

    Polynomial a, b, zero, one;
    final int BASIS = 7;

    @BeforeEach
    void setUp() {
        a = new Polynomial(new int[]{0,0,0,0,6}, BASIS);
        b = new Polynomial(new int[]{5,2,0,3}, BASIS);
        zero = Polynomial.ZERO(BASIS);
        one = Polynomial.ONE(BASIS);
    }

    @Test
    void getCoefficient() {
        assertEquals(0, a.getCoefficient(3));
        assertEquals(6, a.getCoefficient(4));
        assertThrows(IllegalArgumentException.class, () -> a.getCoefficient(a.getCoefficients().length + 5));
    }

    @Test
    void subtract() {
        assertArrayEquals(a.subtract(a).getCoefficients(), zero.getCoefficients());
        assertArrayEquals(a.subtract(zero).getCoefficients(), a.getCoefficients());
    }

    @Test
    void multiply() {
        assertArrayEquals(a.multiply(one).getCoefficients(), a.getCoefficients());
        assertArrayEquals(one.multiply(a).getCoefficients(), a.getCoefficients());
        assertArrayEquals(a.multiply(zero).getCoefficients(), zero.getCoefficients());
    }

    @Test
    void mod() {
        assertArrayEquals(a.mod(one).getCoefficients(), zero.getCoefficients());
        assertArrayEquals(one.mod(a).getCoefficients(), one.getCoefficients());
        assertArrayEquals(zero.mod(one).getCoefficients(), zero.getCoefficients());
        assertThrows(IllegalArgumentException.class, () -> a.mod(zero));
    }


    @Test
    void evaluatePolynomial() {
        assertEquals(0, zero.evaluatePolynomial(1337));
        assertEquals(1, one.evaluatePolynomial(1337));
        assertEquals(2, b.evaluatePolynomial(4));
        /* TODO: If Polynomial is changed to have a Polynomial array of coefficients
                 then add an evaluation test for the polynomials which aren't constants
         */
    }

    @Test
    void findIrreduciblePoly() {
        LinkedList<Polynomial> irreducibleOverF2 = new LinkedList<>();
        int[][] coeffsForF2 = new int[][] {{1, 1}, {1, 1, 1}, {1, 1, 0, 1}, {1, 0, 1, 1}, {1, 0, 1, 0, 0, 1}, {1, 1, 1, 0, 1, 1}};
        for(int[] coeffs : coeffsForF2)
            irreducibleOverF2.add(new Polynomial(coeffs, 2));
        Polynomial deg1 = Polynomial.findIrreduciblePoly(2, 1);
//        System.out.println("f2 deg1: " + deg1);
        Polynomial deg2 = Polynomial.findIrreduciblePoly(2, 2);
//        System.out.println("f2 deg2: " + deg2);
        Polynomial deg3 = Polynomial.findIrreduciblePoly(2, 3);
//        System.out.println("f2 deg3: " + deg3);
        Polynomial deg5 = Polynomial.findIrreduciblePoly(2, 5);
//        System.out.println("f2 deg5: " + deg5);
        Polynomial deg6 = Polynomial.findIrreduciblePoly(2, 6);
//        System.out.println("f2 deg6: " + deg6);
        assertTrue(deg1.in(irreducibleOverF2));

        LinkedList<Polynomial> irreducibleOverF3 = new LinkedList<>();
        int[][] coeffsForF3 = new int[][] {{0, 1}, {1, 1}, {2, 1}, {1, 0, 1}, {2, 0, 1}, {2, 1, 1}, {2, 2, 1}, {1, 0, 0, 1}, {2, 2, 0, 1}, {1, 0, 1, 0, 1, 1}};
        for(int[] coeffs : coeffsForF3)
            irreducibleOverF3.add(new Polynomial(coeffs, 3));
        Polynomial f3Deg1 = Polynomial.findIrreduciblePoly(3, 1);
//        System.out.println("f3 deg1: " + f3Deg1);
        Polynomial f3Deg2 = Polynomial.findIrreduciblePoly(3, 2);
//        System.out.println("f3 deg2: " + f3Deg2);
        assertTrue(f3Deg2.in(irreducibleOverF3));
        Polynomial f3Deg3 = Polynomial.findIrreduciblePoly(3, 3);
//        System.out.println("f3 deg3: " + f3Deg3);
        assertTrue(f3Deg3.in(irreducibleOverF3));
        // no irreducible polynomial of degree 4
        Polynomial f3Deg5 = Polynomial.findIrreduciblePoly(3, 5);
//        System.out.println("f3 deg5: " + f3Deg5);
        assertTrue(f3Deg5.in(irreducibleOverF3));
    }

    @Test
    void in() {
        // Collections implementation
        LinkedList<Polynomial> lst = new LinkedList<>();
        lst.add(zero);
        lst.add(one);
        assertFalse(a.in(lst));
        assertTrue(one.in(lst));
        lst.add(b);
        lst.add(a);
        assertTrue(b.in(lst));
        lst.pop();
        assertFalse(zero.in(lst));

        // Arrays implementation
        Polynomial[] arr = new Polynomial[]{a, b, zero};
        assertTrue(a.in(arr));
        assertFalse(one.in(arr));
        arr[0] = one;
        assertFalse(a.in(arr));
        assertTrue(one.in(arr));
    }
}