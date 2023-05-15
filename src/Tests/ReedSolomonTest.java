package Tests;

import Code.GaloisField;
import Code.Polynomial;
import Code.ReedSolomon;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.LinkedList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class ReedSolomonTest {

    GaloisField GF929 = new GaloisField(929);
    GaloisField GF7 = new GaloisField(7);
    Polynomial P929;
    Polynomial P7;

    @BeforeEach
    void setUp() {
        P929 = new Polynomial(new int[]{3,2,1}, GF929);
        P7 = new Polynomial(new int[]{3,2,1}, GF7);
    }

    @Test
    void RSEncoder() {
        List<Polynomial> encoded = ReedSolomon.RSEncoder(P7, 12);
        List<Polynomial> expected = new LinkedList<>();
        expected.add(new Polynomial(new int[]{3, 6, 4, 4, 2, 6, 4, 1, 3, 3, 5, 1}, GF7)); // encoded message
        expected.add(new Polynomial(new int[]{4, 4, 2, 6, 3, 6, 4, 4, 2, 6, 3, 6}, GF7));
        expected.add(new Polynomial(new int[]{3}, GF7));
        expected.add(new Polynomial(new int[]{1, 6, 4, 6, 0, 0, 6, 1, 3, 1}, GF7));

        assertEquals(expected.get(0), encoded.get(0));
        assertEquals(expected.get(1), encoded.get(1));
        assertEquals(expected.get(2), encoded.get(2));
        assertEquals(expected.get(3), encoded.get(3));
    }

    @Test
    void uniqueDecoder() {
        List<Polynomial> encodedMsg = ReedSolomon.RSEncoder(P7, 12);
        Polynomial symbols = encodedMsg.get(1);
        int k = encodedMsg.get(2).getCoefficient(0);
        // Should succeed decoding as there are 0 errors
        Polynomial decodedMsg = ReedSolomon.uniqueDecoder(symbols, k);
        assertEquals(decodedMsg, P7);

        int[] symbolsWithErrors = symbols.getCoefficients();

        // n = 12, therefore since max_errors = (n-k)/2 = 12-3/2 = 4 should fail
        symbolsWithErrors[1] = 1;
        symbolsWithErrors[3] = 3;
        symbolsWithErrors[5] = 1;
        symbolsWithErrors[6] = 3;
        symbolsWithErrors[9] = 1;
        symbolsWithErrors[11] = 3;
        Polynomial P7_with_errors = new Polynomial(symbolsWithErrors, GF7);
        Polynomial failedDecoding = ReedSolomon.uniqueDecoder(P7_with_errors, k);
        assertNull(failedDecoding);
    }

//    @Test
//    void listDecoder() {
//        List<Polynomial> encodedMsg = ReedSolomon.RSEncoder(P7, 12);
//        Polynomial symbols = encodedMsg.get(1);
//        int k = encodedMsg.get(2).getCoefficient(0);
//
//        List<Polynomial> possibleMsg = ReedSolomon.listDecoder(symbols, k);
//        assertTrue(P7.in(possibleMsg));
//
//        assertEquals(ReedSolomon.listDecode(symbols, k), P7);
//
//        int[] symbolsWithErrors = symbols.getCoefficients();
//
//        // n = 6, therefore since max_errors = (n-k)/2 = 6-3/2 = 1 should fail
//        symbolsWithErrors[1] = 1;
//        symbolsWithErrors[3] = 3;
//        Polynomial P7_with_errors = new Polynomial(symbolsWithErrors, GF7);
//        Polynomial failedDecoding = ReedSolomon.uniqueDecoder2(P7_with_errors, k);
//
//        possibleMsg = ReedSolomon.listDecode(symbolsWithErrors, k);
//        assertNull(possibleMsg);
//    }
}