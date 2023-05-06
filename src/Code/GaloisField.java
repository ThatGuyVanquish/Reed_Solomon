package Code;

public class GaloisField {

    private final int prime;

    /**
     * Implementation of Galois Field class constructor
     * Used to centralize method over prime based fields
     * @param prime prime to be the basis of the field
     */
    public GaloisField(int prime) {
        this.prime = prime;
    }

    /**
     * Returns the prime this galois field is over
     * @return an integer which is the basis for this galois field
     */
    public int getPrime() {
        return this.prime;
    }

    /**
     * Given two integers, return their positive value sum over Fp
     * @param x
     * @param y
     * @return positive integer equal to (x + y) % p
     */
    public int add(int x, int y) {
        return mod(x + y);
    }

    /**
     * Given two integers, return their positive value subtraction over Fp
     * @param x
     * @param y
     * @return positive integer equal to (x - y) % p
     */
    public int subtract(int x, int y) {
        return mod(x - y);
    }

    /**
     * Given two integers, returns their positive value multiplication over Fp
     * @param x
     * @param y
     * @return positive integer equal to (x * y) % p
     */
    public int multiply(int x, int y) {
        return mod(x * y);
    }

    /**
     * Given an integer, returns its multiplicative inverse over Fp if one exists, otherwise throws an exception.
     * @param a
     * @return the multiplicative index of the given integer over Fp
     * @throws ArithmeticException if given integer has no multiplicative inverse over Fp
     */
    public int modInverse(int a) throws ArithmeticException{
        a = mod(a);
        for (int x = 1; x < prime; x++) {
            if (mod(a * x) == 1) {
                return x;
            }
        }
        throw new ArithmeticException("modular inverse does not exist");
    }

    /**
     * Given two integers, calculate the positive value division of x/y.
     * @param x
     * @param y
     * @return positive integer equal to (x/y) % q
     */
    public int div(int x, int y) {
        int inverseOfDivisor = modInverse(y);
        return (x * inverseOfDivisor) % prime;
    }

    /**
     * Given an integer, returns its positive representation over Fp
     * @param x
     * @return a positive integer equal to x % p
     */
    public int mod(int x) {
        return (x % prime + prime) % prime; // ensure positive result
    }

    public int[] gaussianElimination(int[][] mat, int[] sol) {
        int rows = mat.length;
        int cols = mat[0].length;

        // Append solution array to matrix
        int[][] augmentedMat = new int[rows][cols + 1];
        for (int i = 0; i < rows; i++) {
            System.arraycopy(mat[i], 0, augmentedMat[i], 0, cols);
            augmentedMat[i][cols] = sol[i];
        }

        // Perform row echelon form elimination
        for (int i = 0; i < rows; i++) {
            // Find pivot row
            int pivotRow = i;
            while (pivotRow < rows && augmentedMat[pivotRow][i] == 0) {
                pivotRow++;
            }

            if (pivotRow >= rows) {
                break; // No non-zero pivot found, so exit
            }

            if (pivotRow != i) {
                // Swap current row with pivot row
                int[] temp = augmentedMat[i];
                augmentedMat[i] = augmentedMat[pivotRow];
                augmentedMat[pivotRow] = temp;
            }

            // Perform row reduction
            int pivot = augmentedMat[i][i];
            int invPivot = modInverse(pivot);
            for (int j = i; j <= cols; j++) {
                augmentedMat[i][j] = mod(augmentedMat[i][j] * invPivot);
            }
            for (int k = 0; k < rows; k++) {
                if (k != i) {
                    int factor = augmentedMat[k][i];
                    for (int j = i; j <= cols; j++) {
                        augmentedMat[k][j] = mod(augmentedMat[k][j] - factor * augmentedMat[i][j]);
                    }
                }
            }
        }
        // Extract solution array from reduced row echelon form
        int[] solution = new int[rows];
        for (int i = 0; i < rows; i++) {
            solution[i] = augmentedMat[i][cols];
        }
        return solution;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GaloisField)) return false;
        GaloisField that = (GaloisField) o;
        return prime == that.prime;
    }
}
