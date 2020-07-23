package alg;
import svm.SVM;
import alg.*;
import io.*;
import java.util.*;

public class FlexMargin extends Algorithm{
	SVM svm;
	
	public FlexMargin(SVM svm){
		super(svm);
		this.svm = svm;
		if(svm.ind.V != null){
			name = "Flexible Margin";
			svm.outd.algorithm = name;
			svm.outd.max_stages_count = 1;
			svm.outd.showInputData();
		}
	}
	
	public float sqrNorm(float[] a){
		float s = 0;
		for(int i = 0; i<a.length; ++i)
			s+=a[i]*a[i];
		return s;
	}
	
	public float dotProduct(float[] a, float[] b){
		float s = 0;
		for(int i = 0; i < a.length; i++)
			s += a[i] * b[i];
		return s;
	}
	
	
	public void run(){
		t = System.currentTimeMillis();
		int k = 0;
		int i, j, y_changed_i;
		float C = svm.settings.C;
		float tol = 0.01f;
		float bias = 0;
		int num_changed_alphas;
		int y_changed_j, y_changed_k;
		float[] alpha = new float[N];
		float old_alpha_i = 0;
		float old_alpha_j = 0;
		float E_i = 0;
		float E_j = 0;
		float L = 0;
		float H = 0;
		float eta;
		float b1, b2;
		float[] w_2 = new float[dim];
		float[] w_3 = new float[dim+1];
		float[] w_m2 = new float[dim+1];
		float[] w_p2 = new float[dim+1];
		Random rand = new Random();

		for(i = 0; i < N; i++)
			alpha[i] = 0;
		
		for(int passes = 0; passes < 100; ++passes){
			passes++;
			num_changed_alphas = 0;
			
			for(i = 0; i < N; i++){
                E_i = 0;
				y_changed_i = svm.ind.V[i].cl.Y == 0? -1 : 1;
				for(k = 0; k < N; k++){
					y_changed_k = svm.ind.V[k].cl.Y == 0? -1 : 1;
					E_i += alpha[k] * y_changed_k * dotProduct(svm.ind.V[k].X, svm.ind.V[i].X);
				}
				
				E_i = E_i + bias - y_changed_i;
				
				if((y_changed_i * E_i < -tol && alpha[i] < C) || (y_changed_i * E_i > tol && alpha[i] > 0)){
					do{
						j = rand.nextInt(N);
					}while(j == i);
					
                    y_changed_j = svm.ind.V[j].cl.Y == 0? -1 : 1;
					
                    E_j = 0;
					
					for(k = 0; k < N; k++){
						y_changed_k = svm.ind.V[k].cl.Y == 0? -1 : 1;
						E_j += alpha[k] * y_changed_k * dotProduct(svm.ind.V[k].X, svm.ind.V[j].X);
					}
					E_j = E_j + bias - y_changed_j;
					old_alpha_i = alpha[i];
					old_alpha_j = alpha[j];
					if(y_changed_i != y_changed_j){
						L = Math.max(0, alpha[j] - alpha[i]);
						H = Math.min(C, C + alpha[j] - alpha[i]);
					}
					else if(y_changed_i == y_changed_j){
						L = Math.max(0, alpha[i] + alpha[j] - C);
						H = Math.min(C, alpha[i] + alpha[j]);
					}
					if(L == H)
						continue;
					eta = 2 * dotProduct(svm.ind.V[i].X, svm.ind.V[j].X) - sqrNorm(svm.ind.V[i].X) - sqrNorm(svm.ind.V[j].X);
					if (eta >= 0)
						continue;
					alpha[j] = alpha[j] + y_changed_j * (E_j - E_i) / eta;
					if(alpha[j] > H)
						alpha[j] = H;
					else if(alpha[j] < L)
						alpha[j] = L;
					if(Math.abs(alpha[j] - old_alpha_j) < 1e-3)
						continue;
					alpha[i] = alpha[i] - y_changed_i * y_changed_j * (alpha[j] - old_alpha_j);
					b1 = bias - E_i 
					- y_changed_i * (alpha[i] - old_alpha_i) * dotProduct(svm.ind.V[i].X,svm.ind.V[i].X) 
					- y_changed_j * (alpha[j] - old_alpha_j) * dotProduct(svm.ind.V[j].X,svm.ind.V[i].X);
					b2 = bias - E_j 
					- y_changed_i * (alpha[i] - old_alpha_i) * dotProduct(svm.ind.V[i].X,svm.ind.V[j].X) 
					- y_changed_j * (alpha[j] - old_alpha_j) * dotProduct(svm.ind.V[j].X,svm.ind.V[j].X);
				
					bias = (b1 + b2)/2;
					num_changed_alphas = num_changed_alphas + 1;
				}
			}
			if(num_changed_alphas == 0)
                passes++;
            else passes = 0;
		}
	
		
		ArrayList<Integer> vectoriSuport = new ArrayList<Integer>();
		for(i = 0; i < N; i++){
			if(alpha[i] == Float.NaN)
				continue;
			if(alpha[i] != 0)
				vectoriSuport.add(i);
		}
					
		for(i = 0; i < vectoriSuport.size(); i++){
			int index = vectoriSuport.get(i);
			int y_i = svm.ind.V[index].cl.Y == 0 ? -1 : 1;
			for(k = 0; k < dim; k++)
				w_2[k] += alpha[index] * svm.ind.V[index].X[k] * y_i;
		}
		
		for(i = 0; i < dim; i++)
			System.out.print(w_2[i] + " ");
		System.out.println();
		
		System.out.println("b = " + bias);
		
		for(i = 0; i < dim; i++)
			w_3[i] = w_2[i];
		
		 w_3[dim] = bias;
		
		for(i = 0; i < dim; i++){
				w_p2[i] = w_3[i];
				w_m2[i] = w_3[i];
			}
		w_p2[dim] = bias + 1;
		w_m2[dim] = bias - 1;
		
		svm.outd.computing_time = System.currentTimeMillis() - t;
		
		float[] b = Tools.translate(w_3,svm.ind.V);
		w_3[dim] = -b[2];
		float[] w0 = new float[dim+1];
		for(j=0; j<dim; j++) w0[j] = w_3[j];
		w0[dim] = -b[0];
		float[] w1 = new float[dim+1];
		for(j=0; j<dim; j++) w1[j] = w_3[j];
		w1[dim] = -b[1];			
		if(dim==2) svm.design.setPointsOfMaxLine(w_3,w0,w1);	
		
		svm.outd.w = w_3;
		svm.outd.accuracy = getAccuracy(w_3);
		svm.outd.margin = b[3];
		svm.outd.stages_count = 1;
		svm.outd.max_stages_count = 1;
		svm.outd.showInputData();
		svm.outd.showOutputData();
		svm.design.calculates = false;
		svm.design.repaint();
		svm.control.start.enable(false);		
	}

}
