import Code.GaloisField;
import Code.Polynomial;
import Code.ReedSolomon;

import java.util.List;

public class ReedSolomonCodes {

    public static void main(String[] args) {
        GaloisField GF7 = new GaloisField(7);
        List<Polynomial> test = ReedSolomon.RSEncoder_L(new Polynomial(new int[]{1,6,3}, GF7), 7);
        Polynomial t = test.get(1);
        int[] t1 = t.getCoefficients();
        t1[1] = 5;
        t1[4] = 3;
        t = new Polynomial(t1, GF7);
        Polynomial res = ReedSolomon.uniqueDecoder_L(t, 3);
        System.out.println(res);
    }
}
