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
     * @return a list of Polynomials of length 2, with res.at(0) = encoded message and res.at(1) = encoded symbols
     * @throws IllegalArgumentException if basis q of message has no irreducible polynomials
     */
    public static List<Polynomial> RSEncoder(Polynomial msg, EncodedLength nType) throws IllegalArgumentException{
        int k = msg.degree();
        int basis = msg.getBasis();
        if (nType == null) nType = EncodedLength.SHORT;
        int n = calculateEncodedLength(k, nType);
//        if (n > basis)
//            return RSEncoder(msg, n, powerOfQ(basis, n));
        int generatorDegree = n - k + 1;
        Polynomial generatorPolynomial = Polynomial.findIrreduciblePoly(msg.getBasis(), generatorDegree);

        if (generatorPolynomial == null)
            throw new IllegalArgumentException("Message basis has no irreducible polynomials");

        // message should be multiplied by x^(n-k+1)
        int[] multiplierCoeffs = new int[n-k+1];
        multiplierCoeffs[n-k] = 1;
        Polynomial multiplier = new Polynomial(multiplierCoeffs, basis);
        Polynomial multipliedMsg = msg.multiply(multiplier);
        Polynomial encodedMsg = multipliedMsg.mod(generatorPolynomial);

        // obtain encoded symbols by evaluating encodedMsg in n distinct points over Fq

        return null;
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
