import Code.GaloisField;
import Code.Polynomial;
import Code.ReedSolomon;

public class ReedSolomonCodes {

    public static void main(String[] args) {
        Polynomial res = ReedSolomon.uniqueDecoder(new Polynomial(new int[]{1,4,3,6,1,2,2}, new GaloisField(7)), 3);
        System.out.println(res);
    }
}
