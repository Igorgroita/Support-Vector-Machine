package alg;

import svm.SVM;
import io.Vector;

public class MMPerceptron extends Algorithm{

	public MMPerceptron(SVM svm){
		super(svm);
		if(svm.ind.V != null){
			name = "2xMedian > Perceptron";
			svm.outd.algorithm = name;
			svm.outd.showInputData();
		}
	}
		
	public float[] median(Vector[] V){
		int dim = V[0].getDimension();
		int N = V.length;
		float[] w = new float[dim+1];	
		float[] M0 = new float[dim];
		float[] M1 = new float[dim];
		int k0 = 0, k1 = 0;		
		for(int i = 0; i < N; i++)
			if(V[i].cl.Y == 0) k0++; else k1++;
		for(int j = 0; j < dim; j++)
			for(int i = 0; i < N; i++)
				if(V[i].cl.Y == 0) 
					M0[j] += V[i].X[j];
				else 
					M1[j] += V[i].X[j];
		for(int j = 0; j < dim; j++){
			M0[j] /= k0;
			M1[j] /= k1;
		}
		float[] X0 = new float[dim];
		for(int j = 0; j < dim; j++){
			X0[j] = (M0[j] + M1[j])/2;
			w[j] =  M1[j] - M0[j];
			w[dim] -= w[j] * X0[j];
		}
		return w;		
	}
	
	public float[] getCloserHyperplan(Vector[] V){
		float min0 = Float.MAX_VALUE;
		float min1 = Float.MAX_VALUE;
		int j = -1, k = -1;
		
		float[] w = median(V);
		
		for(int i = 0; i < V.length; i++){
			float d =  distFromHiperplanToVector(w, V[i]);
			if(V[i].cl.Y == 0)
				if(d < min0){min0 = d; j = i;}
			else
				if(d < min1){min1 = d; k = i;}
		}
		
		float b0 = 0, b1 = 0;
		for(int i = 0; i < w.length-1; i++) b0 += w[i]*V[j].X[i];
		for(int i = 0; i < w.length-1; i++) b1 += w[i]*V[k].X[i];
		
		if(min0 < min1) w[w.length-1] = -b1;
		else w[w.length-1] = -b0;
		return w;
	}	
	
	public int getCloserClass(float[] w, Vector[] V){
		int cl = 0;
		float min0 = Float.MAX_VALUE;
		float min1 = Float.MAX_VALUE;
		for(int i = 0; i < V.length; i++){
			float d =  distFromHiperplanToVector(w, V[i]);
			if(V[i].cl.Y == 0) 
				if(min0 < d) min0 = d;
			else 
				if(min1 < d) min1 = d;
		}
		if(min1 < min0) cl = 1;
		return cl;
	}
	
	public float distFromHiperplanToVector(float[] w, Vector V){
		float dist = 0;
		float norm = 0;
		for(int j = 0; j < w.length-1; j++){
			dist += w[j]*V.X[j];
			norm += w[j]*w[j];
		}
		dist += w[w.length-1];
		dist = Math.abs(dist);
		norm = (float)Math.sqrt(norm);
		dist /= norm;
		return dist;
	}
	
	public int getOrientation(float[] w, Vector[] V){
		int cl = getCloserClass(w, V) == 0 ? 1 : 0;
		int nrm = 0, nrp = 0;
		for(int i = 0; i < V.length; i++){
			if(V[i].cl.Y == cl){
				float h = 0;
				for(int j = 0; j < w.length-1; j++)
					h += w[j]*V[i].X[j];
				h += w[w.length-1];
				if(h > 0) nrp++;
				else nrm++;
			}
		}
		return nrm < nrp ? 1 : -1;
	}

	public Vector[] removeSideVectors(float[] w, Vector[] V){
		int orientation = getOrientation(w, V);
		int cl = getCloserClass(w, V);
		Vector[] W = new Vector[V.length];
		int k = 0;
		for(int i = 0; i < V.length; i++){
			float h = 0;
			for(int j = 0; j < w.length-1; j++)
				h += w[j]*V[i].X[j];
			h += w[w.length-1];	
			if((orientation == 1 && h > 0) || (orientation == -1 && h <= 0)){
				W[k] = V[i];
				k++;
			}
		}
		Vector[] WW = new Vector[k];
		System.arraycopy(W, 0, WW, 0, k);
		return WW;
	}
	
	public float[] Median(Vector[] V){
		float[] w = median(V);
		V = removeSideVectors(w, V);
		w = median(V);
		V = removeSideVectors(w, V);
		w = median(V);
		V = removeSideVectors(w, V);
		w = median(V);		
		return w;
	}	
	
	public void run(){
		t = System.currentTimeMillis();
		boolean flag = false;
		
		//float[] w = median(svm.ind.V);
		//float[] w = Median(svm.ind.V);
		float[] w = getCloserHyperplan(svm.ind.V);

		float[] w1 = new float[dim+1];
		for(int j = 0; j < w.length; j++) w1[j] = w[j];
		
		for(long p = 1; p <= P; p++){
			boolean erori = false;
			for(int i = 0; i < N; i++) {
				float s = 0;
				for(int j = 0; j < dim; j++)
					s += w[j]*svm.ind.V[i].X[j];
				s += w[dim];
				int y = s < 0 ? 0 : 1;
				int e = svm.ind.V[i].cl.Y - y;
				if(e != 0){
					erori = true;
					for(int j = 0; j < dim; j++)
						w[j] += eta*svm.ind.V[i].X[j]*e;
					w[dim] += eta*e;					
				}
			}
			if(!erori) {
				svm.outd.stages_count = p;
				svm.outd.computing_time = System.currentTimeMillis() - t;
				svm.outd.w = w;
				svm.outd.accuracy = getAccuracy(w);
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
			for(int j = 0; j < w.length; j++) w[j] = w1[j];
			for(int p = 1; p <= P; p++){
				boolean erori = false;
				for(int i = 0; i < N; i++) {
					float s = 0;
					for(int j = 0; j < dim; j++)
						s += w[j]*svm.ind.V[i].X[j];
					s += w[dim];
					int y = s < 0 ? 0 : 1;
					int e = svm.ind.V[i].cl.Y - y;
					if(e != 0){
						erori = true;
						for(int j = 0; j < dim; j++)
							w[j] += eta*svm.ind.V[i].X[j]*e;
						w[dim] += eta*e;					
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
					svm.outd.w = w;
					svm.outd.accuracy = getAccuracy(w);
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
	
}