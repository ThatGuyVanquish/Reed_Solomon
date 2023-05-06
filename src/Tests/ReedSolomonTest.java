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
    void calculateEncodedLength() {
        int k = P7.degree() + 1;
        assertEquals(4, ReedSolomon.calculateEncodedLength(k, ReedSolomon.EncodedLength.NOT_ENOUGH));
        assertEquals(6, ReedSolomon.calculateEncodedLength(k, ReedSolomon.EncodedLength.SHORT));
        assertEquals(12, ReedSolomon.calculateEncodedLength(k, ReedSolomon.EncodedLength.MEDIUM));
        assertEquals(24, ReedSolomon.calculateEncodedLength(k, ReedSolomon.EncodedLength.LONG));
    }

    @Test
    void RSEncoder() {
        testEncoderNOT_ENOUGH();
        testEncoderSHORT();
        testEncoderMEDIUM();
        testEncoderLONG();
    }

    private void testEncoderNOT_ENOUGH() {
        List<Polynomial> encodedNOT_ENOUGH = ReedSolomon.RSEncoder(P7, ReedSolomon.EncodedLength.NOT_ENOUGH);
        List<Polynomial> expectedNOT_ENOUGH = new LinkedList<>();
        expectedNOT_ENOUGH.add(new Polynomial(new int[]{5, 4, 6, 1}, GF7)); // encoded message
        expectedNOT_ENOUGH.add(new Polynomial(new int[]{4, 4, 2, 6}, GF7));
        expectedNOT_ENOUGH.add(new Polynomial(new int[]{3}, GF7));
        expectedNOT_ENOUGH.add(new Polynomial(new int[]{4, 1}, GF7));

        assertEquals(expectedNOT_ENOUGH.get(0), encodedNOT_ENOUGH.get(0));
        assertEquals(expectedNOT_ENOUGH.get(1), encodedNOT_ENOUGH.get(1));
        assertEquals(expectedNOT_ENOUGH.get(2), encodedNOT_ENOUGH.get(2));
        assertEquals(expectedNOT_ENOUGH.get(3), encodedNOT_ENOUGH.get(3));
    }

    private void testEncoderSHORT() {
        List<Polynomial> encodedSHORT = ReedSolomon.RSEncoder(P7, ReedSolomon.EncodedLength.SHORT);
        List<Polynomial> expectedSHORT = new LinkedList<>();
        expectedSHORT.add(new Polynomial(new int[]{4, 1, 3, 3, 5, 1}, GF7)); // encoded message
        expectedSHORT.add(new Polynomial(new int[]{4, 4, 2, 6, 3, 6}, GF7));
        expectedSHORT.add(new Polynomial(new int[]{3}, GF7));
        expectedSHORT.add(new Polynomial(new int[]{6, 1, 3, 1}, GF7));

        assertEquals(expectedSHORT.get(0), encodedSHORT.get(0));
        assertEquals(expectedSHORT.get(1), encodedSHORT.get(1));
        assertEquals(expectedSHORT.get(2), encodedSHORT.get(2));
        assertEquals(expectedSHORT.get(3), encodedSHORT.get(3));
    }

    private void testEncoderMEDIUM() {
        List<Polynomial> encodedMEDIUM = ReedSolomon.RSEncoder(P7, ReedSolomon.EncodedLength.MEDIUM);
        List<Polynomial> expectedMEDIUM = new LinkedList<>();
        expectedMEDIUM.add(new Polynomial(new int[]{3, 6, 4, 4, 2, 6, 4, 1, 3, 3, 5, 1}, GF7)); // encoded message
        expectedMEDIUM.add(new Polynomial(new int[]{4, 4, 2, 6, 3, 6, 4, 4, 2, 6, 3, 6}, GF7));
        expectedMEDIUM.add(new Polynomial(new int[]{3}, GF7));
        expectedMEDIUM.add(new Polynomial(new int[]{1, 6, 4, 6, 0, 0, 6, 1, 3, 1}, GF7));

        assertEquals(expectedMEDIUM.get(0), encodedMEDIUM.get(0));
        assertEquals(expectedMEDIUM.get(1), encodedMEDIUM.get(1));
        assertEquals(expectedMEDIUM.get(2), encodedMEDIUM.get(2));
        assertEquals(expectedMEDIUM.get(3), encodedMEDIUM.get(3));
    }

    private void testEncoderLONG() {
        List<Polynomial> encodedLONG = ReedSolomon.RSEncoder(P7, ReedSolomon.EncodedLength.LONG);
        List<Polynomial> expectedLONG = new LinkedList<>();
        expectedLONG.add(new Polynomial(new int[]{3, 6, 4, 4, 2, 6, 5, 3, 2, 2, 1, 3, 2, 4, 5, 5, 6, 4, 4, 1, 3, 3, 5, 1}, GF7)); // encoded message
        expectedLONG.add(new Polynomial(new int[]{4, 4, 2, 6, 3, 6, 4, 4, 2, 6, 3, 6, 4, 4, 2, 6, 3, 6, 4, 4, 2, 6, 3, 6}, GF7));
        expectedLONG.add(new Polynomial(new int[]{3}, GF7));
        expectedLONG.add(new Polynomial(new int[]{1, 6, 4, 6, 0, 0, 4, 3, 2, 3, 0, 0, 3, 4, 5, 4, 0, 0, 6, 1, 3, 1}, GF7));

        assertEquals(expectedLONG.get(0), encodedLONG.get(0));
        assertEquals(expectedLONG.get(1), encodedLONG.get(1));
        assertEquals(expectedLONG.get(2), encodedLONG.get(2));
        assertEquals(expectedLONG.get(3), encodedLONG.get(3));
    }

    @Test
    void uniqueDecoder() {
        uniqueDecoderNOT_ENOUGH();
        uniqueDecoderSHORT();
        uniqueDecoderMEDIUM();
        uniqueDecoderLONG();
    }

    void uniqueDecoderNOT_ENOUGH() {
        List<Polynomial> encodedMsgNOT_ENOUGH = ReedSolomon.RSEncoder(P7, ReedSolomon.EncodedLength.NOT_ENOUGH);
        Polynomial symbols = encodedMsgNOT_ENOUGH.get(1);
        int k = encodedMsgNOT_ENOUGH.get(2).getCoefficient(0);
        Polynomial decodedMsg = ReedSolomon.uniqueDecoder2(symbols, k);
        assertEquals(decodedMsg, P7);

        int[] symbolsWithErrors = symbols.getCoefficients();

        // n = 4, therefore since max_errors = (n-k)/2 = 4-3/2 = 0 should fail
        symbolsWithErrors[1] = 1;
        symbolsWithErrors[3] = 3;
        symbolsWithErrors[4] = 1;
        Polynomial P7_with_errors = new Polynomial(symbolsWithErrors, GF7);
        Polynomial failedDecoding = ReedSolomon.uniqueDecoder2(P7_with_errors, k);
        assertNull(failedDecoding);
    }

    void uniqueDecoderSHORT() {
        List<Polynomial> encodedMsgSHORT = ReedSolomon.RSEncoder(P7, ReedSolomon.EncodedLength.SHORT);
        Polynomial symbols = encodedMsgSHORT.get(1);
        int k = encodedMsgSHORT.get(2).getCoefficient(0);
        // Should succeed decoding as there are 0 errors
        Polynomial decodedMsg = ReedSolomon.uniqueDecoder2(symbols, k);
        assertEquals(decodedMsg, P7);

        int[] symbolsWithErrors = symbols.getCoefficients();

        // n = 6, therefore since max_errors = (n-k)/2 = 6-3/2 = 1 should fail
        symbolsWithErrors[1] = 1;
        symbolsWithErrors[3] = 3;
        Polynomial P7_with_errors = new Polynomial(symbolsWithErrors, GF7);
        Polynomial failedDecoding = ReedSolomon.uniqueDecoder2(P7_with_errors, k);
        assertNull(failedDecoding);
    }

    void uniqueDecoderMEDIUM() {
        List<Polynomial> encodedMsgMEDIUM = ReedSolomon.RSEncoder(P7, ReedSolomon.EncodedLength.MEDIUM);
        Polynomial symbols = encodedMsgMEDIUM.get(1);
        int k = encodedMsgMEDIUM.get(2).getCoefficient(0);
        // Should succeed decoding as there are 0 errors
        Polynomial decodedMsg = ReedSolomon.uniqueDecoder2(symbols, k);
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
        Polynomial failedDecoding = ReedSolomon.uniqueDecoder2(P7_with_errors, k);
        assertNull(failedDecoding);
    }
    void uniqueDecoderLONG() {
        List<Polynomial> encodedMsgLONG = ReedSolomon.RSEncoder(P7, ReedSolomon.EncodedLength.LONG);
        Polynomial symbols = encodedMsgLONG.get(1);
        int k = encodedMsgLONG.get(2).getCoefficient(0);
        // Should succeed decoding as there are 0 errors
        Polynomial decodedMsg = ReedSolomon.uniqueDecoder2(symbols, k);
        assertEquals(decodedMsg, P7);

        int[] symbolsWithErrors = symbols.getCoefficients();

        // n = 24, therefore since max_errors = (n-k)/2 = 24-3/2 = 10 should fail
        symbolsWithErrors[1] = 1;
        symbolsWithErrors[3] = 3;
        symbolsWithErrors[5] = 1;
        symbolsWithErrors[6] = 3;
        symbolsWithErrors[9] = 1;
        symbolsWithErrors[11] = 3;
        symbolsWithErrors[12] = 1;
        symbolsWithErrors[14] = 3;
        symbolsWithErrors[15] = 1;
        symbolsWithErrors[16] = 3;
        symbolsWithErrors[19] = 1;
        symbolsWithErrors[21] = 3;
        Polynomial P7_with_errors = new Polynomial(symbolsWithErrors, GF7);
        Polynomial failedDecoding = ReedSolomon.uniqueDecoder2(P7_with_errors, k);
        assertNull(failedDecoding);
    }

//    @Test
//    void listDecoder() {
//        List<Polynomial> encodedMsg = ReedSolomon.RSEncoder(P7, ReedSolomon.EncodedLength.SHORT);
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