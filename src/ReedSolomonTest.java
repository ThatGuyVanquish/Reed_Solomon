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
    Polynomial shortPolyOverF50;

    @BeforeEach
    void setUp() {
        P929 = new Polynomial(new int[]{3,2,1}, GF929);
        P7 = new Polynomial(new int[]{3,2,1}, GF7);
    }

    @Test
    void calculateEncodedLength() {
        int k = P7.degree() + 1;
        assertEquals(12, ReedSolomon.calculateEncodedLength(k, ReedSolomon.EncodedLength.NOT_ENOUGH));
        assertEquals(16, ReedSolomon.calculateEncodedLength(k, ReedSolomon.EncodedLength.SHORT));
        assertEquals(32, ReedSolomon.calculateEncodedLength(k, ReedSolomon.EncodedLength.MEDIUM));
        assertEquals(64, ReedSolomon.calculateEncodedLength(k, ReedSolomon.EncodedLength.LONG));
    }

    @Test
    void RSEncoder() {
        /* TODO:
            Calculate the generator polynomial for the 4 encode length
                1) Test for NOT_ENOUGH? Generally ill advised so IDK if I should test for this
                2) Test for SHORT
                3) Test for MEDIUM
                4) Test for LONG
         */
        System.out.println("F50 primitive element is 2");
        /*
            Basis is 50, therefore alpha (Primitive element) is 2
            The polynomial is 42x^2+6x+13
         */encodedSymbols = new Polynomial(new int[]{43, 9, 49, 11}, 50);
        testEncoderNOT_ENOUGH();
        testEncoderSHORT();
        testEncoderMEDIUM();
        testEncoderLONG();
    }

    private void testEncoderNOT_ENOUGH() {
        List<Polynomial> encodedNOT_ENOUGH = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.NOT_ENOUGH);
        Polynomial generatorNOT_ENOUGH =
                ReedSolomon.computeGeneratorPolynomial(50, ReedSolomon.calculateEncodedLength(3, ReedSolomon.EncodedLength.NOT_ENOUGH), 3);
        List<Polynomial> expectedNOT_ENOUGH = new LinkedList<>();
        expectedNOT_ENOUGH.add(new Polynomial(new int[]{4, 20, 13, 4, 42}, 50)); // encoded message
        expectedNOT_ENOUGH.add(encodedSymbols);
        expectedNOT_ENOUGH.add(new Polynomial(new int[]{3}, 50));
        expectedNOT_ENOUGH.add(generatorNOT_ENOUGH);

        assertEquals(expectedNOT_ENOUGH.get(0), encodedNOT_ENOUGH.get(0));
        assertEquals(expectedNOT_ENOUGH.get(1), encodedNOT_ENOUGH.get(1));
        assertEquals(expectedNOT_ENOUGH.get(2), encodedNOT_ENOUGH.get(2));
        assertEquals(expectedNOT_ENOUGH.get(3), encodedNOT_ENOUGH.get(3));
    }

    private void testEncoderSHORT() {
        List<Polynomial> encodedSHORT = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.SHORT);
        int n = ReedSolomon.calculateEncodedLength(3, ReedSolomon.EncodedLength.SHORT);
        Polynomial generatorSHORT =
                ReedSolomon.computeGeneratorPolynomial(50, n, 3);
        List<Polynomial> expectedSHORT = new LinkedList<>();
        expectedSHORT.add(new Polynomial(new int[]{12, 14, 38, 20, 43, 46, 42}, 50)); // encoded message
        encodedSymbols = new Polynomial(new int[]{43, 9, 49, 11, 13, 29}, 50);
        expectedSHORT.add(encodedSymbols);
        expectedSHORT.add(new Polynomial(new int[]{3}, 50));
        expectedSHORT.add(generatorSHORT);

        assertEquals(expectedSHORT.get(0), encodedSHORT.get(0));
        assertEquals(expectedSHORT.get(1), encodedSHORT.get(1));
        assertEquals(expectedSHORT.get(2), encodedSHORT.get(2));
        assertEquals(expectedSHORT.get(3), encodedSHORT.get(3));
    }

    private void testEncoderMEDIUM() {
        List<Polynomial> encodedMEDIUM = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.MEDIUM);
        int n = ReedSolomon.calculateEncodedLength(3, ReedSolomon.EncodedLength.MEDIUM);
        Polynomial generatorMEDIUM =
                ReedSolomon.computeGeneratorPolynomial(50, n, 3);
        List<Polynomial> expectedMEDIUM = new LinkedList<>();
        expectedMEDIUM.add(new Polynomial(new int[]{34, 40, 8, 34, 34, 30, 44, 42, 40, 30, 23, 24, 42}, 50)); // encoded message
        encodedSymbols = new Polynomial(new int[]{43, 9, 49, 11, 13, 29, 9, 11, 33, 49, 19, 11}, 50);
        expectedMEDIUM.add(encodedSymbols);
        expectedMEDIUM.add(new Polynomial(new int[]{3}, 50));
        expectedMEDIUM.add(generatorMEDIUM);

        assertEquals(expectedMEDIUM.get(0), encodedMEDIUM.get(0));
        assertEquals(expectedMEDIUM.get(1), encodedMEDIUM.get(1));
        assertEquals(expectedMEDIUM.get(2), encodedMEDIUM.get(2));
        assertEquals(expectedMEDIUM.get(3), encodedMEDIUM.get(3));
    }

    private void testEncoderLONG() {
        List<Polynomial> encodedLONG = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.LONG);
        int n = ReedSolomon.calculateEncodedLength(3, ReedSolomon.EncodedLength.LONG);
        Polynomial generatorLONG =
                ReedSolomon.computeGeneratorPolynomial(50, n, 3);
        List<Polynomial> expectedLONG = new LinkedList<>();
        expectedLONG.add(new Polynomial(new int[]{46, 30, 12, 46, 28, 0, 40, 20, 20, 0, 20, 10, 20, 0, 30, 40, 0, 0, 10, 30, 44, 20, 13, 4, 42}, 50)); // encoded message
        encodedSymbols = new Polynomial(new int[]{43, 9, 49, 11, 13, 29, 9, 11, 33, 49, 19, 11, 3, 19, 29, 11, 23, 39, 39, 11, 43, 9, 49, 11}, 50);
        expectedLONG.add(encodedSymbols);
        expectedLONG.add(new Polynomial(new int[]{3}, 50));
        expectedLONG.add(generatorLONG);

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

        List<Polynomial> encodedMsgSHORT = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.SHORT);
        List<Polynomial> encodedMsgMEDIUM = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.MEDIUM);
        List<Polynomial> encodedMsgLONG = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.LONG);

    }

    void uniqueDecoderNOT_ENOUGH() {
        List<Polynomial> encodedMsgNOT_ENOUGH = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.NOT_ENOUGH);
        /*
            change some values within the encoded symbols
         */

    }

    void uniqueDecoderSHORT() {
        List<Polynomial> encodedMsgSHORT = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.SHORT);
        /*
            change some values within the encoded symbols
         */

    }

    void uniqueDecoderMEDIUM() {
        List<Polynomial> encodedMsgMEDIUM = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.MEDIUM);
        /*
            change some values within the encoded symbols
         */

    }
    void uniqueDecoderLONG() {
        List<Polynomial> encodedMsgLONG = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.LONG);
        /*
            change some values within the encoded symbols
         */

    }


    @Test
    void listDecoder() {

    }
}