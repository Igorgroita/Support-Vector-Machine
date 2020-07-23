package alg;

import io.*;
import java.util.Random;

public class Optimization{


//====================== algoritmi de optimizare ======================

	/*
	programare patratica cu simplex(Wolfe 1959)
	Hang T. Lau, A Java Library of Graph Algorithms and Optimization, 2006

	parametri:
		n: nr de variabile 
		m: nr de constrangeri
		Q: matricea (simetrica) ce defineste partea patratica a functiei obiectiv, double[n+1][n+1], Q[i][j], i=1,..,n, j=1,..,n
		v: matricea care defineste partea liniara a functiei obiectiv, double[n+1], v[j], j=1,..,n	
		X: matricea constrangerilor, double[m+1][n+1], X[i][j], i=1,..,m, j=1,..,n
		c: matricea termen liber a constrangerilor, double[m+1], c[i], i=1,..,m
		sol: solutia optima gasita (vectorul w), double[n+1], sol[j], j=1,..,n, iar sol[0] este valoarea minima a functiei obiectiv
		
	metoda returneaza o valoare intreaga:
		0: solutia optima a fost gasita
		1: functia obiectiv este nemarginita
		2: a fost depasit numarul maxim de iteratii
	*/

	public static int Wolfe(int n, int m, double Q[][], double v[], double A[][], double c[], int maxiterations, double sol[]){
		int mplusn,m2n2,m2n3,m2n3plus1,big,unrestricted;
		int i,j,iterations,temp,col,row=0;
		int index[] = new int[n+m+1];
		double total,candidate,dividend;
		double table[][] = new double[n+m+1][2*m+3*n+2];
		double price[] = new double[2*m+3*n+2];
		double change[] = new double[2*m+3*n+2];
		double net[] = new double[2*m+3*n+2];
		double gain[] = new double[2*m+3*n+2];
		double fraction[] = new double[n+m+1];
		
		//initializari
		for (i=0; i<=n; i++) sol[i] = 0.;
		
		for (i=0; i<=n; i++)
			for (j=0; j<=n; j++)
				if (i != j) Q[i][j] /= 2.;
		
		big = Integer.MAX_VALUE;
		mplusn = m + n;
		m2n2 = m + m + n + n;
		m2n3 = m2n2 + n;
		m2n3plus1 = m2n3 + 1;
		
		for (i=1; i<=mplusn; i++)
			for (j=1; j<=m2n3plus1; j++)
				table[i][j] = 0.;
		for (i=1; i<=m; i++) table[i][1] = c[i];
		for (i=m+1; i<=mplusn; i++)	table[i][1] = -v[i-m];
		for (i=1; i<=m; i++)
			for (j=1; j<=n; j++)
				table[i][j+1] = A[i][j];
		for (i=1; i<=n; i++)
			for (j=1; j<=n; j++)
				table[i+m][j+1] = 2. * Q[i][j];
		for (i=m+1; i<=mplusn; i++)
			for (j=n+2; j<=mplusn+1; j++)
				table[i][j] = A[j-n-1][i-m];
		for (i=1; i<=mplusn; i++) {
			temp = i + mplusn + n + 1;
			for (j=m2n2+2; j<=m2n3plus1; j++)
				if (j == temp) table[i][j] = 1.;
		}
		for (i=m+1; i<=mplusn; i++) {
			temp = i - m + mplusn + 1;
			for (j=mplusn+2; j<=m2n3plus1; j++)
				if (j == temp) table[i][j] = -1.;
		}
		
		for (j=1; j<=m2n3; j++) price[j] = 0.;
		for (i=1; i<=m; i++) price[n+1+i] = table[i][1];
		for (j=m2n2+2; j<=m2n3plus1; j++) price[j] = big - 1;
		
		for (i=1; i<=mplusn; i++) index[i] = m2n3 - mplusn + i;
		
		iterations = 0;
		
		//calcul
		while (true) {
			iterations++;
			for (j=1; j<=m2n3plus1; j++) gain[j] = 0.;
			for (j=1; j<=m2n3plus1; j++) {
				total = 0.;
				for (i=1; i<=mplusn; i++)
					total += price[index[i]+1] * table[i][j];
				gain[j] = total;
				change[j] = price[j] - gain[j];
			}
			// cauta pivotul
			col = 0;
			candidate = 0.;
			// obtine variabila cu cel mai mare profit
			for (i=2; i<=m2n3plus1; i++)
				if (change[i] < candidate) {
					candidate = change[i];
					col = i;
				}
			if (col <= 0) break;
			unrestricted = 0;
			for (i=1; i<=mplusn; i++) {
				if (table[i][col] > 0)
					fraction[i] = table[i][1] / table[i][col];
				else {
					unrestricted++;
					if (unrestricted == mplusn)
						return 1; // functia obiectiv este nemarginita
					else
						fraction[i] = Double.MAX_VALUE;
					}
			}
			for (i=1; i<=mplusn; i++)
				if (fraction[i] >= 0) {
					if (fraction[i] > big) fraction[i] = big;
					candidate = fraction[i];
					row = i;
					break;
				}
			for (i=1; i<=mplusn; i++)
				if (candidate > fraction[i]) {
					candidate = fraction[i];
					row = i;
				}
			// pivotare si introducerea unei noi variabile
			dividend = table[row][col];
			for (j=1; j<=m2n3plus1; j++) table[row][j] /= dividend;
			for (i=1; i<=mplusn; i++)
				if (i != row) {
					for (j=1; j<=m2n3plus1; j++)
						net[j] = table[row][j] * table[i][col] / table[row][col];
					for (j=1; j<=m2n3plus1; j++)
						table[i][j] -= net[j];
				}
			price[row] = price[col];
			index[row] = col - 1;
			// recalcularea pretului
			for (j=1; j<=m2n2+1; j++) price[j] = 0.;
			for (i=1; i<=mplusn; i++) {
				if (index[i] <= mplusn)
					temp = index[i] + mplusn + 1;
				else {
					if (index[i] > m2n2) continue;
					temp = index[i] - (mplusn - 1);	
				}
				price[temp] = table[i][1];
			}
			if (iterations >= maxiterations) return 2; // depasirea nr maxim de iteratii
		}
		// returnarea solutiei optime
		total = 0.;
		for (i=1; i<=mplusn; i++)
			if (index[i] <= n) total += v[index[i]] * table[i][1];
		sol[0] = total;
		total =0.;
		for (i=1; i<=mplusn; i++)
			for (j=1; j<=mplusn; j++) {
				if (index[i] > n) continue;
				if (index[j] > n) continue;
				total += Q[index[i]][index[j]] * table[i][1] * table[j][1];
			}
		sol[0] += total;
		for (i=1; i<=mplusn; i++)
			if ((table[i][1] != 0) && (index[i] <= n))
				sol[index[i]] = table[i][1];
		return 0;
	}



	
	// Sequential Minimal Optimization Algorithm (SMO) 

	public static double[] SMO(int N, io.Vector[] V){
		int k = 0;
		int i, j, y_changed_i;
		double C = 10000; // this large number can provide a good model
		// Though, the model could overfit due to this C
		// C can be decreased. Shall try with 1e-5f, 1e-3f, 1e-1f, 1, 10, 100, 1000
		double tol = 0.01f;
		double bias = 0;
		int num_changed_alphas;
		int max_passes = 100;
		int passes = 0;
		int y_changed_j, y_changed_k;
		double old_alpha_i = 0, old_alpha_j = 0;
		double E_i = 0, E_j = 0;
		double L = 0;
		double H = 0;
		double eta;
		double b1, b2;
		double[] alpha = new double[N];	
		Random rand = new Random();

		for(i = 0; i < N; i++) alpha[i] = 0;
			
		while(passes < max_passes){
			passes++;
			num_changed_alphas = 0;
				
			for(i = 0; i < N; i++){
				E_i = 0;
				y_changed_i = V[i].cl.Y == 0? -1 : 1;
				for(k = 0; k < N; k++){
					y_changed_k = V[k].cl.Y == 0? -1 : 1;
					E_i += alpha[k] * y_changed_k * Util.dotProduct(V[k].X, V[i].X);
				}
					
				E_i = E_i + bias - y_changed_i;
					
				if((y_changed_i * E_i < -tol && alpha[i] < C) || (y_changed_i * E_i > tol && alpha[i] > 0)){
					// Selecting a random j != i
					do{
						j = rand.nextInt(N);
					}while(j == i);
						
					y_changed_j = V[j].cl.Y == 0? -1 : 1;
						
					E_j = 0;
						
					for(k = 0; k < N; k++){
						y_changed_k = V[k].cl.Y == 0? -1 : 1;
						E_j += alpha[k] * y_changed_k * Util.dotProduct(V[k].X, V[j].X);
					}
					E_j = E_j + bias - y_changed_j;
					old_alpha_i = alpha[i];
					old_alpha_j = alpha[j];
					if(y_changed_i != y_changed_j){
						L = Math.max(0, alpha[j] - alpha[i]);
						H = Math.min(C, C + alpha[j] - alpha[i]);
					}else if(y_changed_i == y_changed_j){
						L = Math.max(0, alpha[i] + alpha[j] - C);
						H = Math.min(C, alpha[i] + alpha[j]);
					}
					if(L == H)
						continue;
					eta = 2 * Util.dotProduct(V[i].X, V[j].X) 
						- Util.dotProduct(V[i].X, V[i].X) 
						- Util.dotProduct(V[j].X, V[j].X);
					if (eta >= 0)
						continue;
					alpha[j] = alpha[j] + y_changed_j * (E_j - E_i) / eta;
					if(alpha[j] > H)
						alpha[j] = H;
					else if(alpha[j] < L)
						alpha[j] = L;
					if(Math.abs(alpha[j] - old_alpha_j) < 1e-20)
						continue;
					alpha[i] = alpha[i] - y_changed_i * y_changed_j * (alpha[j] - old_alpha_j);
					b1 	= bias - E_i 
						- y_changed_i * (alpha[i] - old_alpha_i) * Util.dotProduct(V[i].X,V[i].X) 
						- y_changed_j * (alpha[j] - old_alpha_j) * Util.dotProduct(V[j].X,V[i].X);
					b2 = bias - E_j 
						- y_changed_i * (alpha[i] - old_alpha_i) * Util.dotProduct(V[i].X,V[j].X) 
						- y_changed_j * (alpha[j] - old_alpha_j) * Util.dotProduct(V[j].X,V[j].X);
					if(0 < alpha[i] && alpha[i] < C)
						bias = b1;
					else if(0 < alpha[j] && alpha[j] < C)
						bias = b2;
					else
						bias = (b1 + b2)/2;
					num_changed_alphas = num_changed_alphas + 1;
				}
			}
			if(num_changed_alphas == 0)
				passes++;
			else passes = 0;
		}	
		return alpha;
	}


	

	/* 
	Interior Point algorithm
	
	*/
	
	public static OptimSolution InteriorPoint(float[][] Q, float[][] v, float[][] A, float[][] c, float[][] x, float[][] y, float[][] lambda){
		int m = A.length; 		// no. of rows
		int n = A[0].length; 	// no. of columns
		int maxIteration = 100;
		float[][] tau = new float[maxIteration][1];
		for(int i = 0; i < maxIteration; i++)
			tau[i][0] = 0f;
		float zstar;
		float[][] zstar1 = new float[maxIteration + 1][1];
		zstar1[0][0] = Util.findValue(Q, x);
		
		for(int k = 0; k < maxIteration; k++) {
			float[][] rd = new float[Q.length][1];
			float[][] rp = new float[A.length][1];
			rd = Util.computeRD(Q, x, A, lambda, v);
			rp = Util.computeRP(A, x, y, c);
			int dimension = Q.length + 2 * m;
			float[][] F = new float[dimension][dimension];
			
			// Populating the F matrix
			for(int i = 0; i < Q.length; i++)
				for(int j = 0; j < Q.length; j++)
					F[i][j] = Q[i][j];
			for(int i = 0; i < Q.length; i++)
				for(int j = Q.length; j < Q.length + m; j++)
					F[i][j] = 0;
			float[][] trMinusA = Util.transposeMatrix(A);
			int z = 0;
			for(int i = 0; i < Q.length; i++) {
				int yy = 0;
				for(int j = Q.length + m; j < dimension; j++) {
					F[i][j] = (-1) * trMinusA[z][yy];
					yy += 1;
				}
				z += 1;	
			}
			int zz = 0;
			for(int i = Q.length; i < Q.length + m; i++) {
				int yyy = 0;
				for(int j = 0; j < n; j++) {
					F[i][j] = A[zz][yyy];
					yyy += 1;
				}
				zz += 1;	
			}
			for(int i = Q.length; i < Q.length + m; i++) {
				for(int j = Q.length; j < Q.length + m; j++) 
					if(i == j)
						F[i][j] = -1;
					else F[i][j] = 0;
			}
			for(int i = Q.length; i < Q.length + m; i++)
				for(int j = Q.length + m; j < dimension; j++)
					F[i][j] = 0;
			for(int i = Q.length + m; i < dimension; i++)
				for(int j = 0; j < n; j++)
					F[i][j] = 0;
			int zzz = 0;
			float[][] diaglambda = Util.diag(lambda);
			for(int i = Q.length + m; i < dimension; i++) {
				for(int j = Q.length; j < Q.length + m; j++) {
					if(i - j == m)
						F[i][j] = diaglambda[zzz][zzz];
				}
				zzz += 1;
			}
			int zzzz = 0;
			float[][] diagy = Util.diag(y);
			for(int i = Q.length + m; i < dimension; i++) {
				for(int j = Q.length + m; j < dimension; j++) {
					if(i == j)
						F[i][j] = diagy[zzzz][zzzz];
				}
				zzzz += 1;
			}
			// Populating the equation_b matrix
			float[][] equation_b = new float[dimension][1];
			for(int i = 0; i < rd.length; i++)
				equation_b[i][0] = (-1) * rd[i][0];
			int mm = 0;
			for(int i = rd.length; i < rd.length + rp.length; i++) {
				equation_b[i][0] = (-1) * rp[mm][0];
				mm++;
			}
			
			float[][] ones = new float[m][1];
			for(int i = 0; i < m; i++)
				ones[i][0] = 1;
			float[][] multiplyDiagLbdYOnes = Util.multiplyMatrix(Util.multiplyMatrixMinus(diaglambda, diagy), ones);
			int mo = 0;
			for(int i = rd.length + rp.length; i < dimension; i++) {
				equation_b[i][0] = multiplyDiagLbdYOnes[mo][0];
				mo++;
			}
			
			float[] b1 = new float[equation_b.length];
			for(int i = 0; i < b1.length; i++)
				b1[i] = equation_b[i][0];
			float[] solution = GaussJordanElimination.test(F, b1);
			float[][] delta_x_aff = new float[n][1];
			float[][] delta_y_aff = new float[m][1];
			float[][] delta_lambda_aff = new float[m][1];
			// Populating delta_x_aff matrix
			for(int i = 0; i < n; i++)
				delta_x_aff[i][0] = 0;
			// Populating delta_y_aff matrix
			for(int i = 0; i < m; i++)
				delta_y_aff[i][0] = 0;
			// Populating delta_lambda_aff matrix
			for(int i = 0; i < m; i++)
				delta_lambda_aff[i][0] = 0;
			
			for(int kk = 0; kk < n; kk++)
				delta_x_aff[kk][0] = solution[kk];
			for(int kk = n; kk < n + m; kk++)
				delta_y_aff[kk - n][0] = solution[kk];
			for(int kk = n + m; kk < n + 2 * m; kk++)
				delta_lambda_aff[kk - n - m][0] = solution[kk];
		
			// Computing mu
			float[][] muM = Util.multiplyMatrix(Util.transposeMatrix(y), lambda);
			float mu = muM[0][0] / m;
			
			float alpha_aff = 1;
			float max1 = Util.maxSumLessThan0(y, alpha_aff, delta_y_aff);
			float max2 = Util.maxSumLessThan0(lambda, alpha_aff, delta_lambda_aff);
			boolean compare = Util.compareTwoMax(max1, max2);
			while(Util.compareTwoMax(max1, max2)) {
				alpha_aff -= 0.01f;
				if(alpha_aff <= 0)
					break;
				max1 = Util.maxSumLessThan0(y, alpha_aff, delta_y_aff);
				max2 = Util.maxSumLessThan0(lambda, alpha_aff, delta_lambda_aff);
			}
			
			float[][] sumyalphadelta = Util.computeSum(y, alpha_aff, delta_y_aff);
			float[][] sumlambdaalphadelta = Util.computeSum(lambda, alpha_aff, delta_lambda_aff);
			float[][] muA = Util.multiplyMatrix(Util.transposeMatrix(sumyalphadelta),sumlambdaalphadelta);
			float mu_aff = muA[0][0] / m;
			
			float sigma = (mu_aff / mu) * (mu_aff / mu) * (mu_aff / mu);
			for(int i = 0; i < rd.length; i++)
				equation_b[i][0] = (-1) * rd[i][0];
			mm = 0;
			for(int i = rd.length; i < rd.length + rp.length; i++) {
				equation_b[i][0] = (-1) * rp[mm][0];
				mm++;
			}
			float[][] diagLambdaAff = Util.diag(delta_lambda_aff);
			float[][] diagYAff = Util.diag(delta_y_aff);
			float[][] scalars = Util.ScalarsMatrix(sigma, mu, ones);
			float[][] multiplyDiagLbdAffYOnes = Util.multiplyMatrix(Util.multiplyMatrixMinus(diagLambdaAff, diagYAff), ones);
			float[][] addedMatrix = Util.addMatrix(Util.addMatrix(multiplyDiagLbdYOnes, multiplyDiagLbdAffYOnes),scalars);
			
			mo = 0;
			for(int i = rd.length + rp.length; i < dimension; i++) {
				equation_b[i][0] = addedMatrix[mo][0];
				mo++;
			}
			
			float[][] delta_x = new float[n][1];
			for(int i = 0; i < n; i++)
				delta_x[i][0] = 0;
			float[][] delta_y = new float[m][1];
			for(int i = 0; i < m; i++)
				delta_y[i][0] = 0;
			float[][] delta_lambda = new float[m][1];
			for(int i = 0; i < m; i++)
				delta_lambda[i][0] = 0;
			
			for(int i = 0; i < b1.length; i++)
				b1[i] = equation_b[i][0];
			solution = GaussJordanElimination.test(F, b1);
			
			for(int kk = 0; kk < n; kk++)
				delta_x[kk][0] = solution[kk];
			for(int kk = n; kk < n + m; kk++)
				delta_y[kk - n][0] = solution[kk];
			for(int kk = n + m; kk < n + 2 * m; kk++)
				delta_lambda[kk - n - m][0] = solution[kk];
			
			tau[k][0] = 0.6f;
			float alpha_tau_pri = 1;
			float[][] yAlphaTauDelta = Util.computeSum(y, alpha_tau_pri, delta_y);
			float[][] oneMinusY = Util.oneMinusTau(tau[k][0], y);
			float[][] temp_cond = Util.compareIneq(yAlphaTauDelta, oneMinusY);
			while(Util.sumCompareIneq(temp_cond)) {
				alpha_tau_pri -= 0.01f;
				if(alpha_tau_pri <= 0)
					break;
				yAlphaTauDelta = Util.computeSum(y, alpha_tau_pri, delta_y);
				temp_cond = Util.compareIneq(yAlphaTauDelta, oneMinusY);
			}
			float alpha_tau_dual = 1;
			float[][] lambdaAlphaDTauDelta = Util.computeSum(lambda, alpha_tau_dual, delta_lambda);
			float[][] oneMinusLambda = Util.oneMinusTau(tau[k][0], lambda);
			temp_cond = Util.compareIneq(lambdaAlphaDTauDelta, oneMinusLambda);
			while(Util.sumCompareIneq(temp_cond)) {
				alpha_tau_dual -= 0.01f;
				if(alpha_tau_dual <= 0)
					break;
				lambdaAlphaDTauDelta = Util.computeSum(lambda, alpha_tau_dual, delta_lambda);
				temp_cond = Util.compareIneq(lambdaAlphaDTauDelta, oneMinusLambda);
			}
			
			float alpha = Math.min(alpha_tau_pri, alpha_tau_dual);
			x = Util.computeX(x, alpha, delta_x);
			y = Util.computeX(y, alpha, delta_y);
			lambda = Util.computeX(lambda, alpha, delta_lambda);
			
			float[][] zstea = Util.multiplyMatrix(Util.multiplyMatrix(Util.transposeMatrix(x),Q), x);
			zstea[0][0] *= 0.5;
			zstar1[k+1][0] = zstea[0][0];
			
			if(Math.abs(zstar1[k+1][0] - zstar1[k][0]) < 1e-8f)
				break;
		}

		// Solution of the problem
		OptimSolution OS = new OptimSolution();
		OS.sol = new float[x.length];
		for(int i = 0; i < OS.sol.length; i++)
			OS.sol[i] = x[i][0];
		
		// Minimum value of the problem
		float[][] zstar2 = Util.multiplyMatrix(Util.multiplyMatrix(Util.transposeMatrix(x),Q), x);
		zstar = zstar2[0][0] * 0.5f;		
		OS.val = zstar;
				
		return OS;
	}



}




//===========================================================




class OptimSolution{
	
	public double val_d;
	public double[][] Sol_d;
	public double[] sol_d;
	
	public float val_f;
	public float[][] Sol_f;
	public float[] sol_f;	
	
	public float val;
	public float[] sol;	
	
}




// https://introcs.cs.princeton.edu/java/95linear/GaussJordanElimination.java.html

class GaussJordanElimination {
    private static final float EPSILON = 1e-8f;

    private final int n;      // n-by-n system
    private float[][] a;     // n-by-(n+1) augmented matrix

    // Gauss-Jordan elimination with partial pivoting
    /**
     * Solves the linear system of equations <em>Ax</em> = <em>b</em>,
     * where <em>A</em> is an <em>n</em>-by-<em>n</em> matrix and <em>b</em>
     * is a length <em>n</em> vector.
     *
     * @param  a2 the <em>n</em>-by-<em>n</em> constraint matrix
     * @param  b the length <em>n</em> right-hand-side vector
     */
    public GaussJordanElimination(float[][] a2, float[] b) {
        n = b.length;

        // build augmented matrix
        a = new float[n][n+n+1];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                a[i][j] = a2[i][j];

        // only needed if you want to find certificate of infeasibility (or compute inverse)
        for (int i = 0; i < n; i++)
            a[i][n+i] = 1.0f;

        for (int i = 0; i < n; i++)
            a[i][n+n] = b[i];

        solve();

        assert certifySolution(a2, b);
    }

    public void solve() {

        // Gauss-Jordan elimination
        for (int p = 0; p < n; p++) {
            // show();

            // find pivot row using partial pivoting
            int max = p;
            for (int i = p+1; i < n; i++) {
                if (Math.abs(a[i][p]) > Math.abs(a[max][p])) {
                    max = i;
                }
            }

            // exchange row p with row max
            swap(p, max);

            // singular or nearly singular
            if (Math.abs(a[p][p]) <= EPSILON) {
                continue;
                // throw new ArithmeticException("Matrix is singular or nearly singular");
            }

            // pivot
            pivot(p, p);
        }
        // show();
    }

    // swap row1 and row2
    public void swap(int row1, int row2) {
        float[] temp = a[row1];
        a[row1] = a[row2];
        a[row2] = temp;
    }

    // pivot on entry (p, q) using Gauss-Jordan elimination
    public void pivot(int p, int q) {

        // everything but row p and column q
        for (int i = 0; i < n; i++) {
            float alpha = a[i][q] / a[p][q];
            for (int j = 0; j <= n+n; j++) {
                if (i != p && j != q) a[i][j] -= alpha * a[p][j];
            }
        }

        // zero out column q
        for (int i = 0; i < n; i++)
            if (i != p) a[i][q] = 0.0f;

        // scale row p (ok to go from q+1 to n, but do this for consistency with simplex pivot)
        for (int j = 0; j <= n+n; j++)
            if (j != q) a[p][j] /= a[p][q];
        a[p][q] = 1.0f;
    }

    /**
     * Returns a solution to the linear system of equations <em>Ax</em> = <em>b</em>.
     *      
     * @return a solution <em>x</em> to the linear system of equations
     *         <em>Ax</em> = <em>b</em>; {@code null} if no such solution
     */
    public float[] primal() {
        float[] x = new float[n];
        for (int i = 0; i < n; i++) {
            if (Math.abs(a[i][i]) > EPSILON)
                x[i] = a[i][n+n] / a[i][i];
            else if (Math.abs(a[i][n+n]) > EPSILON)
                return null;
        }
        return x;
    }

    /**
     * Returns a solution to the linear system of equations <em>yA</em> = 0,
     * <em>yb</em> &ne; 0.
     *      
     * @return a solution <em>y</em> to the linear system of equations
     *         <em>yA</em> = 0, <em>yb</em> &ne; 0; {@code null} if no such solution
     */
    public float[] dual() {
        float[] y = new float[n];
        for (int i = 0; i < n; i++) {
            if ((Math.abs(a[i][i]) <= EPSILON) && (Math.abs(a[i][n+n]) > EPSILON)) {
                for (int j = 0; j < n; j++)
                    y[j] = a[i][n+j];
                return y;
            }
        }
        return null;
    }

    /**
     * Returns true if there exists a solution to the linear system of
     * equations <em>Ax</em> = <em>b</em>.
     *      
     * @return {@code true} if there exists a solution to the linear system
     *         of equations <em>Ax</em> = <em>b</em>; {@code false} otherwise
     */
    public boolean isFeasible() {
        return primal() != null;
    }

    // print the tableaux
    public void show() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                System.out.println("%8.3f "+ a[i][j]);
            }
            System.out.println("| ");
            for (int j = n; j < n+n; j++) {
            	System.out.println("%8.3f "+ a[i][j]);
            }
            System.out.println("| %8.3f\n" + a[i][n+n]);
        }
        System.out.println();
    }

    // check that Ax = b or yA = 0, yb != 0
    public boolean certifySolution(float[][] A, float[] b) {

        // check that Ax = b
        if (isFeasible()) {
            float[] x = primal();
            for (int i = 0; i < n; i++) {
                float sum = 0.0f;
                for (int j = 0; j < n; j++) {
                    sum += A[i][j] * x[j];
                }
                if (Math.abs(sum - b[i]) > EPSILON) {
                	System.out.println("not feasible");
                	System.out.println("b[%d] = %8.3f, sum = %8.3f\n" + i + b[i] + sum);
                    return false;
                }
            }
            return true;
        }

        // or that yA = 0, yb != 0
        else {
            float[] y = dual();
            for (int j = 0; j < n; j++) {
                float sum = 0.0f;
                for (int i = 0; i < n; i++) {
                    sum += A[i][j] * y[i];
                }
                if (Math.abs(sum) > EPSILON) {
                	System.out.println("invalid certificate of infeasibility");
                	System.out.println("sum = %8.3f\n" + sum);
                    return false;
                }
            }
            float sum = 0.0f;
            for (int i = 0; i < n; i++) {
                sum += y[i] * b[i];
            }
            if (Math.abs(sum) < EPSILON) {
            	System.out.println("invalid certificate of infeasibility");
            	System.out.println("yb  = %8.3f\n" + sum);
                return false;
            }
            return true;
        }
    }

    public static float[] test(float[][] A, float[] b) {
    	GaussJordanElimination gaussian = new GaussJordanElimination(A, b);
        if (gaussian.isFeasible()) {
        	//System.out.println("Solution to Ax = b");
            float[] x = gaussian.primal();
            return x;
        }
        else {
        	//System.out.println("Certificate of infeasibility");
            float[] y = gaussian.dual();
            return y;
        }
    }

    public static void test(String name, float[][] A, float[] b) {
    	System.out.println("----------------------------------------------------");
    	System.out.println(name);
    	System.out.println("----------------------------------------------------");
        GaussJordanElimination gaussian = new GaussJordanElimination(A, b);
        if (gaussian.isFeasible()) {
        	System.out.println("Solution to Ax = b");
            float[] x = gaussian.primal();
            for (int i = 0; i < x.length; i++) {
            	System.out.print(x[i] + " ");
            }
        }
        else {
        	System.out.println("Certificate of infeasibility");
            float[] y = gaussian.dual();
            for (int j = 0; j < y.length; j++) {
            	System.out.println(y[j] + " ");
            }
        }
        System.out.println();
        System.out.println();
    }
}	