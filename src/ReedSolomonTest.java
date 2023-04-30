import com.sun.source.tree.ArrayAccessTree;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class ReedSolomonTest {

    Polynomial polyOverF5;
    Polynomial polyOverF20;
    Polynomial polyOverF30;
    Polynomial shortPolyOverF50;
    Polynomial encodedSymbols;

    @BeforeEach
    void setUp() {
        polyOverF20 = new Polynomial(new int[]{2, 5, 0, 1, 7, 0, 0, 3}, 20);
        polyOverF30 = new Polynomial(new int[]{2, 5, 0, 1, 0, 9, 0, 3, 15, 21, 0, 0, 28}, 30);
        shortPolyOverF50 = new Polynomial(new int[]{13, 6, 42}, 50);
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
         */
        encodedSymbols = Polynomial.SYMBOL_POLYNOMIAL(new int[]{43, 9, 49, 11});
        testEncoderNOT_ENOUGH();
        testEncoderSHORT();
        testEncoderMEDIUM();
        testEncoderLONG();

        List<Polynomial> encodedMEDIUM = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.MEDIUM);
        Polynomial generatorMEDIUM =
                Polynomial.computeGeneratorPolynomial(50, ReedSolomon.calculateEncodedLength(3, ReedSolomon.EncodedLength.MEDIUM), 3);
//        System.out.println("Generator polynomial for MEDIUM is " + generatorMEDIUM);
        List<Polynomial> encodedLONG = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.LONG);
        Polynomial generatorLONG =
                Polynomial.computeGeneratorPolynomial(50, ReedSolomon.calculateEncodedLength(3, ReedSolomon.EncodedLength.LONG), 3);
        //System.out.println("Generator polynomial for LONG is " + generatorLONG);

        
        
    }

    private void testEncoderNOT_ENOUGH() {
        List<Polynomial> encodedNOT_ENOUGH = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.NOT_ENOUGH);
        Polynomial generatorNOT_ENOUGH =
                Polynomial.computeGeneratorPolynomial(50, ReedSolomon.calculateEncodedLength(3, ReedSolomon.EncodedLength.NOT_ENOUGH), 3);
        List<Polynomial> expectedNOT_ENOUGH = new LinkedList<>();
        expectedNOT_ENOUGH.add(new Polynomial(new int[]{4, 20, 13, 4, 42}, 50)); // encoded message
        expectedNOT_ENOUGH.add(encodedSymbols);
        expectedNOT_ENOUGH.add(Polynomial.SYMBOL_POLYNOMIAL(new int[]{3}));
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
                Polynomial.computeGeneratorPolynomial(50, n, 3);
        List<Polynomial> expectedSHORT = new LinkedList<>();
        expectedSHORT.add(new Polynomial(new int[]{12, 14, 38, 20, 43, 46, 42}, 50)); // encoded message
        encodedSymbols = Polynomial.SYMBOL_POLYNOMIAL(new int[]{43, 9, 49, 11, 13, 29});
        expectedSHORT.add(encodedSymbols);
        expectedSHORT.add(Polynomial.SYMBOL_POLYNOMIAL(new int[]{3}));
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
                Polynomial.computeGeneratorPolynomial(50, n, 3);
        List<Polynomial> expectedMEDIUM = new LinkedList<>();
        expectedMEDIUM.add(new Polynomial(new int[]{12}, 50)); // encoded message
        encodedSymbols = Polynomial.SYMBOL_POLYNOMIAL(new int[]{43, 9, 49, 11, 13, 29, 9, 11, 33, 49, 19, 11});
        expectedMEDIUM.add(encodedSymbols);
        expectedMEDIUM.add(Polynomial.SYMBOL_POLYNOMIAL(new int[]{3}));
        expectedMEDIUM.add(generatorMEDIUM);

        //assertEquals(expectedMEDIUM.get(0), encodedMEDIUM.get(0));
        assertEquals(expectedMEDIUM.get(1), encodedMEDIUM.get(1));
        assertEquals(expectedMEDIUM.get(2), encodedMEDIUM.get(2));
        assertEquals(expectedMEDIUM.get(3), encodedMEDIUM.get(3));
    }

    private void testEncoderLONG() {
        List<Polynomial> encodedLONG = ReedSolomon.RSEncoder(shortPolyOverF50, ReedSolomon.EncodedLength.LONG);
        int n = ReedSolomon.calculateEncodedLength(3, ReedSolomon.EncodedLength.LONG);
        Polynomial generatorLONG =
                Polynomial.computeGeneratorPolynomial(50, n, 3);
        List<Polynomial> expectedLONG = new LinkedList<>();
        expectedLONG.add(new Polynomial(new int[]{12, 14, 38, 20, 43, 46, 42}, 50)); // encoded message
        encodedSymbols = Polynomial.SYMBOL_POLYNOMIAL(new int[]{43, 9, 49, 11, 13, 29, 9, 11, 33, 49, 19, 11, 3, 19, 29, 11, 23, 39, 39, 11, 43, 9, 49, 11});
        expectedLONG.add(encodedSymbols);
        expectedLONG.add(Polynomial.SYMBOL_POLYNOMIAL(new int[]{3}));
        expectedLONG.add(generatorLONG);

        //assertEquals(expectedLONG.get(0), encodedLONG.get(0));
        assertEquals(expectedLONG.get(1), encodedLONG.get(1));
        assertEquals(expectedLONG.get(2), encodedLONG.get(2));
        assertEquals(expectedLONG.get(3), encodedLONG.get(3));
    }
    
    @Test
    void uniqueDecoder() {

    }

    @Test
    void listDecoder() {

    }
}