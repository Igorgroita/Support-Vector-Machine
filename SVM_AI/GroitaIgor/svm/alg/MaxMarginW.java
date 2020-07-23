package alg;

import svm.SVM;
import io.*;

public class MaxMarginW extends Algorithm{

	public MaxMarginW(SVM svm){
		super(svm);
		if(svm.ind.V != null){
			name = "MaxMargin (Wolfe)";
			svm.outd.algorithm = name;
			svm.outd.max_stages_count = 1;
			svm.outd.showInputData();
		}
	}


	public float[] getOptimHyperplane(io.Vector[] V){
		double[][] Q = new double[N+1][N+1];

		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++){
				double s = 0.;
				for (int k=0; k<dim; k++)
					s += V[i].X[k] * V[j].X[k];
				int y = V[i].cl.Y == 0 ? -1 : 1;			
				s *= y;
				y = V[j].cl.Y == 0 ? -1 : 1;			
				s *= y;				
				Q[i+1][j+1] = s;				
			}	

		double[] v = new double[N+1];
		for(int i=0; i<N; i++) v[i] = -1.;
			
		double[] alpha1 = new double[N+1];
		double[] alpha2 = new double[N+1];
		
		double[][] A1 = new double[2][N+1];
		for(int j=0; j<N; j++){
			int y = V[j].cl.Y == 0 ? -1 : 1;				
			A1[1][j+1] = -(double)y;
		}
		double[][] A2 = new double[2][N+1];
		for(int j=1; j<=N; j++) A2[1][j] = -A1[1][j];
		double[] c = {0.,0.};
	
		int status1 = Optimization.Wolfe(N, 1, Q, v, A1, c, 10000, alpha1);
		int status2 = Optimization.Wolfe(N, 1, Q, v, A2, c, 10000, alpha2);
		
		float[] w1 = new float[dim+1];
		float[] w2 = new float[dim+1];
		
		if (status1 == 0){
			for(int k=0; k<dim; k++){
				double s = 0.;
				for(int i=0; i<N; i++){
					int y = V[i].cl.Y == 0 ? -1 : 1;	
					s += alpha1[i+1]*y*V[i].X[k];
				}
				w1[k] = (float)s;
			}			
			if (status2 == 0){				
				for(int k=0; k<dim; k++){
					double s = 0.;
					for(int i=0; i<N; i++){
						int y = V[i].cl.Y == 0 ? -1 : 1;	
						s += alpha2[i+1]*y*V[i].X[k];
					}
					w2[k] = (float)s;
				}
				float[] b1 = Tools.translate(w1, V);
				float[] b2 = Tools.translate(w2, V);				
				if(b1[3] >= b2[3]) return w1;
				else return w2;
			}
			return w1;
		}else{
			if (status2 == 0){				
				for(int k=0; k<dim; k++){
					double s = 0.;
					for(int i=0; i<N; i++){
						int y = V[i].cl.Y == 0 ? -1 : 1;	
						s += alpha2[i+1]*y*V[i].X[k];
					}
					w2[k] = (float)s;
				}
				return w2;
			} else return null;
		}
	}	


	public void run(){
		t = System.currentTimeMillis();	
		
		float[] w = getOptimHyperplane(svm.ind.V);
		
		svm.outd.computing_time = System.currentTimeMillis() - t;			
		
		if(w != null){
			float[] b = Tools.translate(w, svm.ind.V);
			w[dim] = -b[2];
			float[] w0 = new float[dim+1];
			for(int j=0; j<dim; j++) w0[j] = w[j];
			w0[dim] = -b[0];
			float[] w1 = new float[dim+1];
			for(int j=0; j<dim; j++) w1[j] = w[j];
			w1[dim] = -b[1];
			
			if(dim==2) svm.design.setPointsOfMaxLine(w,w0,w1);
			svm.outd.accuracy = getAccuracy(w);
			svm.outd.w = w;
			svm.outd.stages_count = 1;
			svm.outd.max_stages_count = 1;
			svm.outd.margin = b[3];
			svm.outd.showInputData();
			svm.outd.showOutputData();
		}		
		svm.design.calculates = false;
		svm.design.repaint();
		svm.control.start.enable(false);

	}
	
}