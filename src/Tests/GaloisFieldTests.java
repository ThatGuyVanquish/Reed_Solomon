package Tests;
import Code.GaloisField;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;


public class GaloisFieldTests {

    GaloisField GF7 = new GaloisField(7);

    @BeforeEach
    void setUp() {
    }

    @Test
    void getPrime() {
        assertEquals(GF7.getPrime(), 7);
    }

    @Test
    void add() {
        assertEquals(GF7.add(2, 4), 6);
        assertEquals(GF7.add(3,0), 3);
        assertEquals(GF7.add(3,6), 2);
        assertEquals(GF7.add(-3, 5), 2);
    }

    @Test
    void subtract() {
        assertEquals(GF7.subtract(6, 4), 2);
        assertEquals(GF7.subtract(3,0), 3);
        assertEquals(GF7.subtract(3,6), 4);
        assertEquals(GF7.subtract(1, -5), 6);
    }

    @Test
    void multiply() {
        assertEquals(GF7.multiply(2, 3), 6);
        assertEquals(GF7.multiply(3,0), 0);
        assertEquals(GF7.multiply(3,1), 3);
        assertEquals(GF7.multiply(3, 5), 1);
    }

    @Test
    void modInverse() {
        assertEquals(GF7.modInverse(3), 5);
        assertEquals(GF7.modInverse(-2), 3);
        assertThrows(ArithmeticException.class, () -> GF7.modInverse(0));
    }

    @Test
    void div() {
        assertEquals(GF7.div(4, 2), 2);
        assertEquals(GF7.div(3,1), 3);
        assertEquals(GF7.div(3, 5), 2);
        assertThrows(IllegalArgumentException.class, () -> GF7.div(3, 0));
    }

    @Test
    void mod() {
        assertEquals(GF7.mod(6), 6);
        assertEquals(GF7.mod(15), 1);
        assertEquals(GF7.mod(-4), 3);
    }

    @Test
    void gaussianElimination() {
        int[] sol = {1, 3, 5};
        int[][] mat = {{3, 4, 6}, {0, 1, 5}, {2, 1, 2}};
        int[] expected = {1, 3, 0};
        assertArrayEquals(expected, GF7.gaussianElimination(mat, sol));
    }

    @Test
    void testEquals() {
        assertEquals(GF7, new GaloisField(7));
    }
}
