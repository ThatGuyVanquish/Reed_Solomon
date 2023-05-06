package Tests;

import Code.Polynomial;
import Code.GaloisField;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertFalse;

import java.util.LinkedList;

class PolynomialTest {

    Polynomial a, b, zero, one;
    GaloisField F7 = new GaloisField(7);

    @BeforeEach
    void setUp() {
        a = new Polynomial(new int[]{0,0,0,0,6}, F7);
        b = new Polynomial(new int[]{5,2,0,3}, F7);
        zero = Polynomial.ZERO(F7);
        one = Polynomial.ONE(F7);
    }

    @Test
    void getCoefficient() {
        assertEquals(0, a.getCoefficient(3));
        assertEquals(6, a.getCoefficient(4));
        assertEquals(0, a.getCoefficient(a.getCoefficients().length + 5));
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