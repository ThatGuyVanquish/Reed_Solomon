import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

class PolynomialTest {

    Polynomial a, b, zero, one;
    final int BASIS = 7;

    @org.junit.jupiter.api.BeforeEach
    void setUp() {
        a = new Polynomial(new int[]{0,0,0,0,6}, BASIS);
        b = new Polynomial(new int[]{5,2,0,3}, BASIS);
        zero = new Polynomial(new int[]{0}, BASIS);
        one = new Polynomial(new int[]{1}, BASIS);
    }

    @org.junit.jupiter.api.Test
    void getCoefficient() {
        assertEquals(0, a.getCoefficient(3));
        assertEquals(6, a.getCoefficient(4));
        assertThrows(IllegalArgumentException.class, () -> a.getCoefficient(a.getCoefficients().length + 5));
    }

    @org.junit.jupiter.api.Test
    void subtract() {
        assertArrayEquals(a.subtract(a).getCoefficients(), zero.getCoefficients());
        assertArrayEquals(a.subtract(zero).getCoefficients(), a.getCoefficients());
    }

    @org.junit.jupiter.api.Test
    void multiply() {
        assertArrayEquals(a.multiply(one).getCoefficients(), a.getCoefficients());
        assertArrayEquals(one.multiply(a).getCoefficients(), a.getCoefficients());
        assertArrayEquals(a.multiply(zero).getCoefficients(), zero.getCoefficients());
    }

    @org.junit.jupiter.api.Test
    void mod() {
        assertArrayEquals(a.mod(one).getCoefficients(), zero.getCoefficients());
        assertArrayEquals(one.mod(a).getCoefficients(), one.getCoefficients());
        assertArrayEquals(zero.mod(one).getCoefficients(), zero.getCoefficients());
        assertThrows(IllegalArgumentException.class, () -> a.mod(zero));
    }
}