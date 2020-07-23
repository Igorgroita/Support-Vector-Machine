package alg;

import java.math.BigDecimal;


public class Util{


public static double dotProduct(double[] a, double[] b){
	double s = 0;
	for(int i = 0; i < a.length; i++)
		s += a[i] * b[i];
	return s;
}


public static float dotProduct(float[] a, float[] b){
	float s = 0;
	for(int i = 0; i < a.length; i++)
		s += a[i] * b[i];
	return s;
}


public static double determinant(double[][] A){			// lucreaza bine cu dimensiuni mari
	return (new Determinant(A)).determinant().doubleValue();
}


public static double determinant(double[][] A, int N){  //metoda recursiva si lenta
	double det=0;
	if(N == 1)
		det = A[0][0];
	else if (N == 2)
		det = A[0][0]*A[1][1] - A[1][0]*A[0][1];
	else{
		det=0;
		for(int j1=0;j1<N;j1++){
			double[][] m = new double[N-1][];
			for(int k=0;k<(N-1);k++)
				m[k] = new double[N-1];
			for(int i=1;i<N;i++){
				int j2=0;
				for(int j=0;j<N;j++){
					if(j == j1)continue;
					m[i-1][j2] = A[i][j];
					j2++;
				}
			}
			det += Math.pow(-1.0,1.0+j1+1.0)* A[0][j1] * determinant(m,N-1);
		}
	}
	return det;
}


public static double[][] matrixInverse(double a[][]){
	int n = a.length;
	double x[][] = new double[n][n];
	double b[][] = new double[n][n];
	int index[] = new int[n];
	
	for (int i=0; i<n; ++i) 
		b[i][i] = 1;
			
	gaussian(a, index); 
			
	for (int i=0; i<n-1; ++i)
		for (int j=i+1; j<n; ++j)
			for (int k=0; k<n; ++k)
				b[index[j]][k] -= a[index[j]][i]*b[index[i]][k];

	for (int i=0; i<n; ++i){
		x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
		for (int j=n-2; j>=0; --j){
			x[j][i] = b[index[j]][i];
			for (int k=j+1; k<n; ++k) 
				x[j][i] -= a[index[j]][k]*x[k][i];
			x[j][i] /= a[index[j]][j];
		}
	}
	
	return x;
}


public static void gaussian(double a[][], int index[]){
	int n = index.length;
	double c[] = new double[n];
	for (int i=0; i<n; ++i) index[i] = i;
	for (int i=0; i<n; ++i){
		double c1 = 0;
		for (int j=0; j<n; ++j){
			double c0 = Math.abs(a[i][j]);
			if (c0 > c1) c1 = c0;
		}
		c[i] = c1;
	}
	int k = 0;
	for (int j=0; j<n-1; ++j){
		double pi1 = 0;
		for (int i=j; i<n; ++i){
			double pi0 = Math.abs(a[index[i]][j]);
			pi0 /= c[index[i]];
			if (pi0 > pi1){
				pi1 = pi0;
				k = i;
			}
		}
		int itmp = index[j];
		index[j] = index[k];
		index[k] = itmp;
		for (int i=j+1; i<n; ++i){
			double pj = a[index[i]][j]/a[index[j]][j];
			a[index[i]][j] = pj;
			for (int l=j+1; l<n; ++l)
				a[index[i]][l] -= pj*a[index[j]][l];
		}
	}
}


public static float[][] transposeMatrix(float[][] x){
	float[][] transpose = new float[x[0].length][x.length];
	for(int i = 0; i < x[0].length; i++)
		for(int j = 0; j < x.length; j++)
			transpose[i][j] = x[j][i];
	return transpose;
}

	
public static float findValue(float[][] G, float[][] x){
	float[][] tx = new float[x.length][1];
	for(int i = 0; i < x.length;i++){
		tx[i][0] = x[i][0];
	}
	float number = 0;
	for(int i = 0; i < G[0].length; i++){
		for(int j = 0; j < x.length; j++){
			number += (1/2)  * tx[j][0] * G[i][j] * x[j][0];
		}
	}
	return number;   	
}

	
public static float[][] computeX(float[][] x, float alpha, float[][] delta){
	float[][] compute = new float[x.length][1];
	for(int i = 0; i < x.length; i++)
		compute[i][0] = x[i][0] + alpha * delta[i][0];
	return compute;
}
	
	
public static float[][] oneMinusTau(float tau, float[][] y){
	float[][] compute = new float[y.length][1];
	for(int i = 0; i < y.length; i++)
		compute[i][0] = (1 - tau) * y[i][0];
	return compute;
}

	
public static float[][] compareIneq(float[][] first, float[][] second) {
	float[][] compute = new float[first.length][1];
	for(int i = 0; i < first.length; i++)
		if(first[i][0] <= second[i][0])
			compute[i][0] = 1;
		else compute[i][0] = 0;
	return compute;
}

	
public static boolean sumCompareIneq(float[][] first) {
	float sum = 0;
	for(int i = 0; i < first.length; i++) sum += first[i][0];		
	if(sum > 0)
		return true;
	else 
		return false;
}

	
public static float[][] addMatrix(float[][] A, float[][] B){
	float[][] compute = new float[A.length][A[0].length];
	for(int i = 0; i < A.length; i++) 
		for(int j = 0; j < A[0].length; j++)
			compute[i][j] = A[i][j] + B[i][j];
	return compute;
}

	
public static float[][] computeRD(float[][] G, float[][] x, float[][] A, float[][] lambda, float[][] c){
	float[][] compute = new float[G.length][1];
	float[][] transposeA = transposeMatrix(A);
	for(int i = 0; i < G.length; i++){
		float Gx = 0;
		float Alambda = 0;
		for(int j = 0; j < G[0].length; j++)
			Gx += G[i][j] * x[j][0];
		for(int k = 0; k < A.length; k++)
			Alambda += transposeA[i][k] * lambda[k][0];
		compute[i][0] = Gx - Alambda + c[i][0];
	}
	return compute;
}

	
public static float[][] computeRP(float[][] A, float[][] x, float[][] y, float[][] b){
	float[][] compute = new float[A.length][1];
	for(int i = 0; i < A.length; i++){
		float Ax = 0;
		for(int j = 0; j < A[0].length; j++)
			Ax += A[i][j] * x[j][0];
		compute[i][0] = Ax - y[i][0] - b[i][0];
	}
	return compute;
}
	
	
public static float[][] diag(float[][] a){
	float[][] dia = new float[a.length][a.length];
	for(int i = 0; i < a.length; i++)
		for(int j = 0; j < a.length; j++)
			if(i==j)
				dia[i][j] = a[i][0];
	return dia;
}

	
public static float[][] multiplyMatrixMinus(float[][] A, float[][] B){
	float[][] C = new float[A.length][B[0].length];
	for(int i = 0; i < A.length; i++)
		for(int j = 0; j < B[0].length; j++)
			for(int k = 0; k < B.length; k++)
			C[i][j] += (-1) * A[i][k] * B[k][j];
	return C;
}

	
public static float[][] multiplyMatrix(float[][] A, float[][] B){
	float[][] C = new float[A.length][B[0].length];
	for(int i = 0; i < A.length; i++)
		for(int j = 0; j < B[0].length; j++)
			for(int k = 0; k < B.length; k++)
			C[i][j] += A[i][k] * B[k][j];
	return C;
}

	
public static float[][] computeSum(float[][] y, float alpha, float[][] delta){
	float[][] compute = new float[y.length][1];
	for(int i = 0; i < y.length; i++) 
		compute[i][0] = y[i][0] + alpha * delta[i][0];
	return compute;
}

	
public static float[][] ScalarsMatrix(float a, float b, float[][] ones){
	float[][] compute = new float[ones.length][1];
	for(int i = 0; i < ones.length; i++)
		compute[i][0] = a * b * ones[i][0];
	return compute;
}

	
public static float maxSumLessThan0(float[][] y, float alpha, float[][] delta){
	float[][] compute = new float[y.length][1];
	for(int i = 0; i < y.length; i++) 
		compute[i][0] = y[i][0] + alpha * delta[i][0];
	for(int j = 0; j < y.length; j++)
		if(compute[j][0] < 0)
			compute[j][0] = 1;
		else compute[j][0] = 0;
	float maxVector = 0;
	for(int i = 0; i < y.length; i++)
		if(maxVector <= compute[i][0])
			maxVector = compute[i][0];
	return maxVector;
}

	
public static boolean compareTwoMax(float m1, float m2) {
	if(m1 == 0)
		if(m2 == 1)
			return true;
		else return false;
	if(m2 == 0)
		if(m1 == 1)
			return true;
		else return false;
	return true;
}





}


//=============================================


class Determinant {

private double[][] matrix;
private int sign = 1;


public Determinant(double[][] matrix) {
    this.matrix = matrix;
}

public int getSign() {
    return sign;
}

public BigDecimal determinant() {

    BigDecimal deter;
    if (isUpperTriangular() || isLowerTriangular())
        deter = multiplyDiameter().multiply(BigDecimal.valueOf(sign));

    else {
        makeTriangular();
        deter = multiplyDiameter().multiply(BigDecimal.valueOf(sign));

    }
    return deter;
}

/*  receives a matrix and makes it triangular using allowed operations
    on columns and rows
*/
public void makeTriangular() {

    for (int j = 0; j < matrix.length; j++) {
        sortCol(j);
        for (int i = matrix.length - 1; i > j; i--) {
            if (matrix[i][j] == 0)
                continue;

            double x = matrix[i][j];
            double y = matrix[i - 1][j];
            multiplyRow(i, (-y / x));
            addRow(i, i - 1);
            multiplyRow(i, (-x / y));
        }
    }
}

public boolean isUpperTriangular() {

    if (matrix.length < 2)
        return false;

    for (int i = 0; i < matrix.length; i++) {
        for (int j = 0; j < i; j++) {
            if (matrix[i][j] != 0)
                return false;

        }

    }
    return true;
}

public boolean isLowerTriangular() {

    if (matrix.length < 2)
        return false;

    for (int j = 0; j < matrix.length; j++) {
        for (int i = 0; j > i; i++) {
            if (matrix[i][j] != 0)
                return false;

        }

    }
    return true;
}

public BigDecimal multiplyDiameter() {

    BigDecimal result = BigDecimal.ONE;
    for (int i = 0; i < matrix.length; i++) {
        for (int j = 0; j < matrix.length; j++) {
            if (i == j)
                result = result.multiply(BigDecimal.valueOf(matrix[i][j]));

        }

    }
    return result;
}

// when matrix[i][j] = 0 it makes it's value non-zero
public void makeNonZero(int rowPos, int colPos) {

    int len = matrix.length;

    outer:
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            if (matrix[i][j] != 0) {
                if (i == rowPos) { // found "!= 0" in it's own row, so cols must be added
                    addCol(colPos, j);
                    break outer;

                }
                if (j == colPos) { // found "!= 0" in it's own col, so rows must be added
                    addRow(rowPos, i);
                    break outer;
                }
            }
        }
    }
}

//add row1 to row2 and store in row1
public void addRow(int row1, int row2) {
    for (int j = 0; j < matrix.length; j++)
        matrix[row1][j] += matrix[row2][j];
}

//add col1 to col2 and store in col1
public void addCol(int col1, int col2) {
    for (int i = 0; i < matrix.length; i++)
        matrix[i][col1] += matrix[i][col2];
}

//multiply the whole row by num
public void multiplyRow(int row, double num) {
    if (num < 0)
        sign *= -1;

    for (int j = 0; j < matrix.length; j++) {
        matrix[row][j] *= num;
    }
}

//multiply the whole column by num
public void multiplyCol(int col, double num) {
    if (num < 0)
        sign *= -1;

    for (int i = 0; i < matrix.length; i++)
        matrix[i][col] *= num;
}

// sort the cols from the biggest to the lowest value
public void sortCol(int col) {
    for (int i = matrix.length - 1; i >= col; i--) {
        for (int k = matrix.length - 1; k >= col; k--) {
            double tmp1 = matrix[i][col];
            double tmp2 = matrix[k][col];

            if (Math.abs(tmp1) < Math.abs(tmp2))
                replaceRow(i, k);
        }
    }
}

//replace row1 with row2
public void replaceRow(int row1, int row2) {
    if (row1 != row2)
        sign *= -1;

    double[] tempRow = new double[matrix.length];

    for (int j = 0; j < matrix.length; j++) {
        tempRow[j] = matrix[row1][j];
        matrix[row1][j] = matrix[row2][j];
        matrix[row2][j] = tempRow[j];
    }
}

//replace col1 with col2
public void replaceCol(int col1, int col2) {
    if (col1 != col2)
        sign *= -1;

    System.out.printf("replace col%d with col%d, sign = %d%n", col1, col2, sign);
    double[][] tempCol = new double[matrix.length][1];

    for (int i = 0; i < matrix.length; i++) {
        tempCol[i][0] = matrix[i][col1];
        matrix[i][col1] = matrix[i][col2];
        matrix[i][col2] = tempCol[i][0];
    }
} 

}