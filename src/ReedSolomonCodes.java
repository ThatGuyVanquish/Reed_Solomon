import Code.GaloisField;
import Code.Polynomial;
import Code.ReedSolomon;
import Code.BivariatePolynomial;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

public class ReedSolomonCodes {

    public static void test() {
        System.out.println("Insert desired prime number as field base:");
        Scanner reader = new Scanner(System.in);
        int base = reader.nextInt();
        GaloisField field = new GaloisField(base);
        System.out.println("Insert an array of coefficients, comma separated:");
        reader.nextLine();
        String coeffsString = reader.nextLine();
        String[] coeffs = coeffsString.split(",");
        int[] coeffsList = new int[coeffs.length];
        for (int i = 0; i < coeffs.length; i++) {
            String coeff = coeffs[i];
            coeffsList[i] = field.mod(Integer.parseInt(coeff));
        }
        Polynomial pol = new Polynomial(coeffsList, field);
        System.out.println("The message to encode over GF" + base + " is: " + pol + "\n");
        int msgLength = coeffsList.length;

        System.out.println("Insert the length of the desired encoding:");
        int encodedLength = reader.nextInt();
        List<Polynomial> regEncoder = ReedSolomon.RSEncoder(pol, encodedLength);
        List<Polynomial> interpolEncoder = ReedSolomon.RSEncoder_L(pol, encodedLength);

        System.out.println("Regular encoded symbols: " + regEncoder.get(0) + "\n"); // encoded symbols
        System.out.println("Interpolation encoded symbols: " + interpolEncoder.get(0) + "\n"); // encoded symbols

        System.out.println("Insert the number of desired errors in the received symbols:");
        int errorCount = reader.nextInt();

        Polynomial RESymbols = regEncoder.get(0);
        int[] REScoeffs = RESymbols.getCoefficients();
        Polynomial IESymbols = interpolEncoder.get(0);
        int[] IEScoeffs = IESymbols.getCoefficients();

        int index = 0;
        // Randomize <=errorCount errors from index 0
        while (errorCount > 0) {
            int error = field.mod((int) (Math.random() * base));
            REScoeffs[index] = error;
            IEScoeffs[index] = error;
            index++;
            errorCount--;
        }

        Polynomial reReceivedSymbols = new Polynomial(REScoeffs, field);
        Polynomial ieReceivedSymbols = new Polynomial(IEScoeffs, field);

        System.out.println("Received symbols from regular encoder: " + reReceivedSymbols + "\n");
        System.out.println("Received symbols from interpolation encoder: " + ieReceivedSymbols + "\n");

        System.out.println("DECODING:\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        Polynomial re_rd = ReedSolomon.uniqueDecoder(reReceivedSymbols, msgLength);
        Polynomial re_id = ReedSolomon.uniqueDecoder_L(reReceivedSymbols, msgLength);
        Polynomial ie_rd = ReedSolomon.uniqueDecoder(ieReceivedSymbols, msgLength);
        Polynomial ie_id = ReedSolomon.uniqueDecoder_L(ieReceivedSymbols, msgLength);

        System.out.println("Regular Encoder => Regular Decoder: " + re_rd + "\n");
        System.out.println("Regular Encoder => Interpolation Decoder: " + re_id + "\n");
        System.out.println("Interpolation Encoder => Regular Decoder: " + ie_rd + "\n");
        System.out.println("Interpolation Encoder => Interpolation Decoder: " + ie_id + "\n");
    }


    public static void GF7Tests() {
        GaloisField GF7 = new GaloisField(7);
        Polynomial tPol = new Polynomial(new int[]{1,6,3,5,4}, GF7);
        System.out.println("Regular encoder:");
        List<Polynomial> regEnc = ReedSolomon.RSEncoder(tPol, 7);
        System.out.println("Encoded symbols: " + regEnc.get(0) + "\n"); // encoded symbols
        System.out.println("Interpolation encoder:");
        List<Polynomial> interpolationEnc = ReedSolomon.RSEncoder_L(tPol, 7);
        System.out.println("Encoded symbols: " + interpolationEnc.get(0) + "\n"); // encoded symbols
        System.out.println();

//      For some reason, the decoding the regular encoder produces terrible results
        System.out.println("Regular encoder => Regular decoder:");
        Polynomial regEncRegDec = ReedSolomon.uniqueDecoder(regEnc.get(0), 5);
        System.out.println(regEncRegDec);
        System.out.println();

        System.out.println("Regular encoder => Interpolation decoder:");
        Polynomial regEncInterpolDec = ReedSolomon.uniqueDecoder_L(regEnc.get(0), 5);
        System.out.println(regEncInterpolDec);
        System.out.println();

        /*
            On the other hand, any decoding of the interpolation encoder produces
            good results when there are no errors.
         */
        System.out.println("Interpolation encoder => Regular decoder:");
        Polynomial interpolationEncRegDec = ReedSolomon.uniqueDecoder_L(interpolationEnc.get(0), 5);
        System.out.println(interpolationEncRegDec);
        System.out.println();

        System.out.println("Interpolation encoder => Interpolation decoder:");
        Polynomial interpolationEncInterpolDec = ReedSolomon.uniqueDecoder_L(interpolationEnc.get(0), 5);
        System.out.println(interpolationEncInterpolDec);
    }

    public static void GF929Tests() {
        GaloisField GF929 = new GaloisField(929);
        Polynomial tPol = new Polynomial(new int[]{1,6,3,42, 13}, GF929);
        System.out.println("Regular encoder:");
        List<Polynomial> regEnc = ReedSolomon.RSEncoder(tPol, 7);
        System.out.println("Encoded symbols: " + regEnc.get(0) + "\n"); // encoded symbols
        System.out.println("Interpolation encoder:");
        List<Polynomial> interpolationEnc = ReedSolomon.RSEncoder_L(tPol, 7);
        System.out.println("Encoded symbols: " + interpolationEnc.get(0) + "\n"); // encoded symbols
        System.out.println();

//      For some reason, the decoding the regular encoder produces terrible results
//        System.out.println("Regular encoder => Regular decoder:");
//        Polynomial regEncRegDec = ReedSolomon.uniqueDecoder(regEnc.get(1), 3);
//        System.out.println(regEncRegDec);
//        System.out.println();
//
//        System.out.println("Regular encoder => Interpolation decoder:");
//        Polynomial regEncInterpolDec = ReedSolomon.uniqueDecoder_L(regEnc.get(1), 3);
//        System.out.println(regEncInterpolDec);
//        System.out.println();

        /*
            On the other hand, interpolation decoding of the interpolation encoder produces
            good results when there are no errors. regular decoding of the interpolation encoder seems
            to mostly produce null
         */
        System.out.println("Interpolation encoder => Regular decoder:");
        Polynomial interpolationEncRegDec = ReedSolomon.uniqueDecoder(interpolationEnc.get(0), 5);
        System.out.println(interpolationEncRegDec);
        System.out.println();

        System.out.println("Interpolation encoder => Interpolation decoder:");
        Polynomial interpolationEncInterpolDec = ReedSolomon.uniqueDecoder_L(interpolationEnc.get(0), 5);
        System.out.println(interpolationEncInterpolDec);
    }

    private static void BVPolyTest() {
        GaloisField GF7 = new GaloisField(7);
        Map<Integer, Map<Integer, Integer>> terms = new HashMap<>();
        terms.put(243, new HashMap<>());
        terms.get(243).put(1, 6);
        terms.put(1, new HashMap<>());
        terms.get(1).put(7, 4);
        // issue with coefficient of y = 1
        BivariatePolynomial poly = new BivariatePolynomial(terms, GF7);
        System.out.println(poly);
        System.out.println(poly.evaluatePolynomial(1, 2));
        System.out.println(poly.evaluatePolynomial(new Polynomial(new int[]{1, 1}, GF7)));
    }

    public static void main(String[] args) {
//        GF929Tests();
//        GF7Tests();
//        test();
        BVPolyTest();
    }
}
