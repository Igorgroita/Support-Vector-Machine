package alg.stud;

import java.util.Arrays;
import java.util.Random;
import java.util.Vector;
import svm.SVM;
import alg.*;
import io.*;

public class MaxMarginSMOSorin extends Algorithm {
	SVM svm;
	float[] alphas, X2, errors, w;
	float C, tol, E2, b, a2, eps;
	Random r;
	int y2;

	public MaxMarginSMOSorin(SVM svm) {
		super(svm);
		this.svm = svm;
		if (svm.ind.V != null) {
			name = "MaxMargin SMO Sorin";
			svm.outd.algorithm = name;
			svm.outd.max_stages_count = 1;
			svm.outd.showInputData();
		}
		alphas = new float[N];
		errors = new float[N];
		w = new float[dim];
		r = new Random(System.currentTimeMillis());
		C = 1e+20f;
		b = 0;
		tol = 1e-2f;
		eps = 1e-10f;
	}

	void init_alg() {
		for (int i = 0; i < N; i++) {
			alphas[i] = 0;
			errors[i] = 0;
		}
		for (int i = 0; i < dim; i++)
			w[i] = 0;
		b = 0;
		E2 = 0;
		a2 = 2;
	}

	float dotProduct(float[] a, float b[]) {
		float s = 0;
		for (int i = 0; i < a.length; i++)
			s += a[i] * b[i];
		return s;
	}

	Vector<Integer> get_non_bound_indices() {
		Vector<Integer> non_bound_indeces = new Vector<>();
		for (int i = 0; i < N; i++)
			if (0 < alphas[i] && alphas[i] < C)
				non_bound_indeces.add(i);
		return non_bound_indeces;
	}

	int first_heuristic() {
		int num_changed = 0;
		Vector<Integer> non_bound_indeces = get_non_bound_indices();
		for (int i = 0; i < non_bound_indeces.size(); i++)
			num_changed += examine_exemple(non_bound_indeces.get(i));
		return num_changed;
	}

	int second_heuristic(Vector<Integer> non_bound_indices) {
		int i1 = -1;
		if (non_bound_indices.size() > 1) {
			float max = 0f;
			for (int j = 0; j < non_bound_indices.size(); j++) {
				int yj = svm.ind.V[j].cl.Y == 0 ? -1 : 1;
				float E1 = errors[j] - yj;
				float step = Math.abs(E1 - E2);
				if (step > max) {
					max = step;
					i1 = j;
				}
			}
		}
		return i1;
	}

	float compute_b(float E1, float a1, float a1_new, float a2_new, float k11, float k12, float k22, int y1) {
		float b1 = E1 + y1 * (a1_new - a1) * k11 + y2 * (a2_new - a2) * k12 + b;
		float b2 = E2 + y1 * (a1_new - a1) * k12 + y2 * (a2_new - a2) * k22 + b;
		float new_b;
		if (0 < a1_new && a1_new < C)
			new_b = b1;
		else if (0 < a2_new && a2_new < C)
			new_b = b2;
		else
			new_b = (b1 + b2) / 2f;
		return new_b;
	}

	float output(int i) {
		return dotProduct(w, svm.ind.V[i].X) - b;
	}

	float error(int i2) {
		int y2 = svm.ind.V[i2].cl.Y == 0 ? -1 : 1;
		return output(i2) - y2;
	}

	float get_error(int i1) {
		if (0 < alphas[i1] && alphas[i1] < C)
			return errors[i1];
		else
			return error(i1);
	}

	float accuracy() {
		int error = 0;
		for (int i = 0; i < N; i++) {
			int yi = svm.ind.V[i].cl.Y == 0 ? -1 : 1;
			if (yi * (dotProduct(w, svm.ind.V[i].X) - b) < 0)
				error++;
		}
		return (float) (error) / N;
	}

	boolean take_step(int i1, int i2) {
		if (i1 == i2)
			return false;
		float a1 = alphas[i1];
		int y1 = svm.ind.V[i1].cl.Y == 0 ? -1 : 1;
		float[] X1 = svm.ind.V[i1].X;
		float E1 = get_error(i1);
		float k11, k12, k22;
		int s = y1 * y2;

		float L, H;
		if (y1 != y2) {
			L = Math.max(0, a2 - a1);
			H = Math.min(C, C + a2 - a1);
		} else {
			L = Math.max(0, a2 + a1 - C);
			H = Math.min(C, a2 + a1);
		}
		if (L == H)
			return false;
		k11 = dotProduct(X1, X1);
		k12 = dotProduct(X1, svm.ind.V[i2].X);
		k22 = dotProduct(svm.ind.V[i2].X, svm.ind.V[i2].X);

		float eta = k11 + k22 - 2 * k12;

		float a2_new;
		if (eta > 0) {
			a2_new = a2 + y2 * (E1 - E2) / eta;

			if (a2_new < L)
				a2_new = L;
			else if (a2_new > H)
				a2_new = H;
		} else {
			float f1 = y1 * (E1 + b) - a1 * k11 - s * a2 * k12;
			float f2 = y2 * (E2 + b) - s * a1 * k12 - a2 * k22;
			float L1 = a1 + s * (a2 - L);
			float H1 = a1 + s * (a2 - H);
			float Lobj = L1 * f1 + L * f2 + 0.5f * (L1 * L1) * k11 + 0.5f * (L * L) * k22 + s * L * L1 * k12;
			float Hobj = H1 * f1 + H * f2 + 0.5f * (H1 * H1) * k11 + 0.5f * (H * H) * k22 + s * H * H1 * k12;

			if (Lobj < Hobj - eps)
				a2_new = L;
			else if (Lobj > Hobj + eps)
				a2_new = H;
			else
				a2_new = a2;
		}

		if (Math.abs(a2_new - a2) < eps * (a2_new + a2 + eps))
			return false;
		float a1_new = a1 + s * (a2 - a2_new);
		float new_b = compute_b(E1, a1, a1_new, a2_new, k11, k12, k22, y1);

		float delta_b = new_b - b;
		b = new_b;
		for (int k = 0; k < dim; k++)
			w[k] += y1 * (a1_new - a1) * X1[k] + y2 * (a2_new - a2) * X2[k];

		float delta1 = y1 * (a1_new - a1);
		float delta2 = y2 * (a2_new - a2);

		for (int i = 0; i < N; i++)
			if (0 < alphas[i] && alphas[i] > C)
				errors[i] += delta1 * dotProduct(X1, svm.ind.V[i].X) + delta2 * dotProduct(X2, svm.ind.V[i].X)
						- delta_b;
		errors[i1] = 0f;
		errors[i2] = 0f;
		alphas[i1] = a1_new;
		alphas[i2] = a2_new;
		return true;
	}

	int examine_exemple(int i2) {
		y2 = svm.ind.V[i2].cl.Y == 0 ? -1 : 1;
		a2 = alphas[i2];
		X2 = svm.ind.V[i2].X;
		E2 = get_error(i2);

		float r2 = E2 * y2;

		if (!((r2 < -tol && a2 < C) || (r2 > tol && a2 > 0)))
			return 0;
		Vector<Integer> non_bound_indeces = get_non_bound_indices();
		int i1 = second_heuristic(non_bound_indeces);

		if (i1 >= 0 && take_step(i1, i2))
			return 1;

		if (non_bound_indeces.size() > 0) {
			int rand_i = r.nextInt(non_bound_indeces.size());
			for (int i = 0; i < non_bound_indeces.size(); i++)
				if (take_step((rand_i + i) % non_bound_indeces.size(), i2))
					return 1;
		}

		int rand_i = r.nextInt(N);
		for (int i = 0; i < N; i++)
			if (take_step((rand_i + i) % N, i2))
				return 1;
		return 0;
	}

	void SMO_alg() {

		int num_changed = 0;
		boolean examine_all = true;

		while (num_changed > 0 || examine_all) {
			num_changed = 0;
			if (examine_all)
				for (int i = 0; i < N; i++)
					num_changed += examine_exemple(i);
			else
				num_changed += first_heuristic();
			if (examine_all)
				examine_all = false;
			else if (num_changed == 0)
				examine_all = true;

		}
	}

	public void run() {
		float max_accuracy = 1000;
		float max_margin = 1000;
		float[] best_w = new float[dim + 1];
		t = System.currentTimeMillis();
		int max_steps = 10000;
		for (int step = 0; step < max_steps; step++) {
			init_alg();
			SMO_alg();
			float[] final_w = new float[dim + 1];
			float margin = dotProduct(w, w);
			for (int i = 0; i < dim; i++)
				final_w[i] = w[i];
			final_w[dim] = -b;
			float accuracy = accuracy();
			if (accuracy < max_accuracy) {
				best_w = Arrays.copyOf(final_w, dim + 1);
				max_margin = margin;
				max_accuracy = accuracy;
				svm.outd.w = final_w;
				svm.design.setPointsOfLine(final_w);
				svm.design.repaint();
			} else if (Math.abs(accuracy - max_accuracy) <= 1e-10f) {
				if (margin < max_margin) {
					max_margin = margin;
					best_w = Arrays.copyOf(final_w, dim + 1);
					svm.outd.w = final_w;
					svm.design.setPointsOfLine(final_w);
					svm.design.repaint();
				}
			}
		}
		
		svm.outd.computing_time = System.currentTimeMillis() - t;
		
		float[] b = Tools.translate(best_w, svm.ind.V);
		best_w[dim] = -b[2];
		float[] w0 = new float[dim+1];
		for(int j=0; j<dim; j++) w0[j] = best_w[j];
		w0[dim] = -b[0];
		float[] w1 = new float[dim+1];
		for(int j=0; j<dim; j++) w1[j] = best_w[j];
		w1[dim] = -b[1];			
		if(dim==2) svm.design.setPointsOfMaxLine(best_w,w0,w1);		
		
		svm.outd.w = best_w;
		//svm.design.setPointsOfLine(best_w);
		svm.outd.accuracy = getAccuracy(svm.outd.w);
		svm.outd.margin = b[3];
		svm.outd.showInputData();
		svm.outd.showOutputData();
		svm.design.calculates = false;
		svm.control.start.enable(false);
		svm.design.repaint();
	}
	
}