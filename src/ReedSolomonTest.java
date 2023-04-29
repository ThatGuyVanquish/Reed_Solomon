import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class ReedSolomonTest {

    Polynomial polyOverF5;
    Polynomial polyOverF20;
    Polynomial polyOverF30;

    @BeforeEach
    void setUp() {
        polyOverF20 = new Polynomial(new int[]{2, 5, 0, 1, 7, 0, 0, 3}, 20);
        polyOverF30 = new Polynomial(new int[]{2, 5, 0, 1, 0, 9, 0, 3, 15, 21, 0, 0, 28}, 30);
    }

    @Test
    void calculateEncodedLength() {
        int k = polyOverF20.degree() + 1;
        assertEquals(12, ReedSolomon.calculateEncodedLength(k, ReedSolomon.EncodedLength.NOT_ENOUGH));
        assertEquals(16, ReedSolomon.calculateEncodedLength(k, ReedSolomon.EncodedLength.SHORT));
        assertEquals(32, ReedSolomon.calculateEncodedLength(k, ReedSolomon.EncodedLength.MEDIUM));
        assertEquals(128, ReedSolomon.calculateEncodedLength(k, ReedSolomon.EncodedLength.LONG));
    }

    @Test
    void RSEncoder() {
    }

    @Test
    void uniqueDecoder() {

    }

    @Test
    void listDecoder() {

    }
}