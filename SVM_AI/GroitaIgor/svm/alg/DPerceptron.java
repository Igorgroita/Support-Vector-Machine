package alg;

import svm.SVM;
import io.*;

public class DPerceptron extends Algorithm{

	public DPerceptron(SVM svm){
		super(svm);
		if(svm.ind.V != null){
			name = "Perceptron Dual";
			svm.outd.algorithm = name;
			svm.outd.showInputData();
		}
	}
	
	public void run(){
		t = System.currentTimeMillis();
		boolean flag = false;
		int[] alpha = new int[N];
		for(int j = 0; j < N; j++) alpha[j] = 0;
		float b = 0;
		float[] w = new float[dim+1]; 
		
		for(long p = 1; p <= P; p++){
			boolean erori = false;
			for(int i = 0; i < N; i++) {
				w = alpha2ponders(alpha, b);
				float s = 0;
				for(int k = 0; k < dim; k++) s += w[k]*svm.ind.V[i].X[k];
				s += -b;
				int y = s < 0 ? 0 : 1; 
				int e = svm.ind.V[i].cl.Y - y;
				if(e != 0){
					erori = true;
					alpha[i]++;
					b = b - eta*e;
				}					
			}
			if(!erori) {
				svm.outd.stages_count = p;
				svm.outd.computing_time = System.currentTimeMillis() - t;
				svm.outd.alpha = alpha;
				svm.outd.b = b;
				svm.outd.w = alpha2ponders(alpha, b);
				svm.outd.accuracy = getAccuracy(svm.outd.w);
				svm.outd.margin = Tools.translate(w, svm.ind.V)[3];
				svm.outd.showInputData();
				svm.outd.showOutputData();
				svm.design.calculates = false;
				svm.design.repaint();
				flag = true;
				break;
			}
		}
		if(!flag) 
			System.out.println(P + " stages have passed. Increase the number of stages and reloaded.");		
		else{		
			for(int j = 0; j < N; j++) alpha[j] = 0;
			b = 0;
			for(int p = 1; p <= P; p++){
				boolean erori = false;
				for(int i = 0; i < N; i++) {
					w = alpha2ponders(alpha, b);
					float s = 0;
					for(int k = 0; k < dim; k++) s += w[k]*svm.ind.V[i].X[k];
					s += w[dim];
					int y = s < 0 ? 0 : 1; 
					int e = svm.ind.V[i].cl.Y - y;
					if(e != 0){
						erori = true;
						alpha[i]++;
						b = b - eta*e;					
					}					
				}
				svm.control.ta.append("Stage " + p + "\n");
				String s = "";
				for(int j = 0; j < w.length; j++) s += "w["+j+"] = " + w[j] + "; ";
				svm.control.ta.append(s + "\n");

				if(dim==2) svm.design.setPointsOfLine(w);
				try{Thread.sleep(250);}
				catch(InterruptedException ex){}			
				
				if(!erori) {
					svm.outd.alpha = alpha;
					svm.outd.b = b;					
					svm.outd.w = alpha2ponders(alpha, b);
					svm.outd.accuracy = getAccuracy(svm.outd.w);
					svm.outd.margin = Tools.translate(w, svm.ind.V)[3];
					svm.outd.showInputData();
					svm.outd.showOutputData();
					svm.control.start.setLabel("Start Simulation");
					svm.design.repaint();
					break;
				}
			}
		}
		svm.control.start.enable(false);		
	}
	
	public float[] alpha2ponders(int[] alpha, float b){
		float[] w = new float[dim+1]; 
		for(int k = 0; k < dim; k++) w[k] = 0;
		w[dim] = -b;
		for(int k = 0; k < dim; k++){
			for(int j = 0; j < N; j++){
				int y = svm.ind.V[j].cl.Y;
				if(y == 0) y = -1;
				w[k] += alpha[j]*y*svm.ind.V[j].X[k];
			}
		}
		
		return w;
	}
	
}