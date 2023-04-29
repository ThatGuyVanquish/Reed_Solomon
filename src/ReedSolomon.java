import java.util.LinkedList;
import java.util.List;

public class ReedSolomon {

    public enum EncodedLength {
        NOT_ENOUGH,
        SHORT,
        MEDIUM,
        LONG
    }

    private static final int ENCODE_MULTIPLIER = 4;
    private static final double NOT_ENOUGH_MULTIPLIER = 1.5;

    public static int calculateEncodedLength(int k, EncodedLength nType) {
        switch (nType) {
            case NOT_ENOUGH -> {
                return (int)(NOT_ENOUGH_MULTIPLIER * k);
            }
            case SHORT -> {
                return 2 * k;
            }
            case MEDIUM -> {
                return ENCODE_MULTIPLIER * k;
            }
            case LONG -> {
                return (int)Math.pow(2, ENCODE_MULTIPLIER) * k;
            }
            default ->  {
                return k;
            }
        }
    }

    /**
     * Returns a pair of Polynomials which are {encoded message polynomial, encoded symbols}
     *
     * @param msg message to be endoded
     * @param nType desired length type of the encoded message, used to get encoded message length
     * @return a list of Polynomials of length 4:
     *  #0 => Encoded message
     *  #1 => Encoded symbols
     *  #2 => k, Length of original message
     *  #3 => Generator polynomial
     * @throws IllegalArgumentException if basis q of message has no irreducible polynomials
     */
    public static List<Polynomial> RSEncoder(Polynomial msg, EncodedLength nType) throws IllegalArgumentException{
        int k = msg.degree();
        int q = msg.getBasis();
        if (nType == null) nType = EncodedLength.SHORT;
        int n = calculateEncodedLength(k, nType);
//        if (n > basis)
//            return RSEncoder(msg, n, powerOfQ(basis, n));
        //int generatorDegree = n - k + 1;
        Polynomial generatorPolynomial = Polynomial.computeGeneratorPolynomial(q, n, k);
        int alpha = Polynomial.findPrimitiveElement(q);

        int[] symbolsArr = new int[n];
        for(int i = 0; i < n; i++) {
            symbolsArr[i] = msg.evaluatePolynomial((int)Math.pow(alpha, i + 1));
        }

        Polynomial encodedMsg = msg.multiply(generatorPolynomial);
        Polynomial encodedSymbols = Polynomial.SYMBOL_POLYNOMIAL(symbolsArr);
        Polynomial constantK = Polynomial.SYMBOL_POLYNOMIAL(new int[]{k});

        List<Polynomial> res = new LinkedList<>();
        res.add(encodedMsg);
        res.add(encodedSymbols);
        res.add(constantK);
        res.add(generatorPolynomial);

        return res;
    }

    private static int powerOfQ(int q, int minLength) {
        int res = q;
        int deg = 1;
        while(res < minLength) {
            res *= q;
            deg++;
        }
        return deg;
    }

}
