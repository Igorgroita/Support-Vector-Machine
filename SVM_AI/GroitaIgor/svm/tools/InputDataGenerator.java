package tools;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import svm.SVM;

public class InputDataGenerator extends Dialog{
	SVM svm;
	TextField attributes_count, vectors_count, min, max, classes_count, margin, abs_error; 
	Checkbox liniar, errors;
	Label attributes_count_label, vectors_count_label, min_label, max_label, classes_count_label, liniar_label, errors_label, margin_label, abs_error_label;
	Button generate, save;
	public TextArea ta;
	String dir = ".\\svm\\data", path;
	public InputDataGenerator(SVM svm){
		super(svm, "Input Data Generator", true);
		this.svm = svm;
		setBackground(svm.settings.background_color_default);
		setResizable(false);
		resize(640,480);
		move((svm.res.width-640)/2,(svm.res.height-480)/2);			
		setLayout(null);
		
		int x1 = 10, x2 = 380;
		int y = 35;
		attributes_count_label = new Label("Attributes Count:");
		attributes_count_label.setBounds(x1,y,150,20);
		attributes_count_label.setForeground(Color.white);
		add(attributes_count_label);
		attributes_count = new TextField("2");
		attributes_count.setBounds(x1+150,y,100,20);
		add(attributes_count);
		
		vectors_count_label = new Label("Vectors Count:");
		vectors_count_label.setBounds(x2,y,150,20);
		vectors_count_label.setForeground(Color.white);
		add(vectors_count_label);
		vectors_count = new TextField("1000");
		vectors_count.setBounds(x2+150,y,100,20);
		add(vectors_count);
		
		y+=30;
		min_label = new Label("Minimum Coordinates:");
		min_label.setBounds(x1,y,150,20);
		min_label.setForeground(Color.white);
		add(min_label);
		min = new TextField("-1000");
		min.setBounds(x1+150,y,100,20);
		add(min);	

		max_label = new Label("Maximum Coordinates:");
		max_label.setBounds(x2,y,150,20);
		max_label.setForeground(Color.white);
		add(max_label);
		max = new TextField("1000");
		max.setBounds(x2+150,y,100,20);
		add(max);
		
		y+=30;
		classes_count_label = new Label("Classes Count:");
		classes_count_label.setBounds(x1,y,150,20);
		classes_count_label.setForeground(Color.white);
		add(classes_count_label);
		classes_count = new TextField("2");
		classes_count.setBounds(x1+150,y,100,20);
		add(classes_count);	
		classes_count.enable(false);
		
		y+=30;
		liniar_label = new Label("Liniar separated:");
		liniar_label.setBounds(x1,y,150,20);
		liniar_label.setForeground(Color.white);
		add(liniar_label);
		liniar = new Checkbox("");
		liniar.setBounds(x1+150,y,20,20);
		liniar.setState(true);
		add(liniar);	
		margin_label = new Label("Margin:");
		margin_label.setBounds(x2,y,150,20);
		margin_label.setForeground(Color.white);
		add(margin_label);
		margin = new TextField("3");
		margin.setBounds(x2+150,y,100,20);
		add(margin);		
		
		y+=30;
		
		errors_label = new Label("Include errors:");
		errors_label.setBounds(x1,y,150,20);
		errors_label.setForeground(Color.white);
		add(errors_label);
		errors = new Checkbox("");
		errors.setBounds(x1+150,y,20,20);
		errors.setState(false);
		errors.enable(liniar.getState());
		add(errors);
		
		
		abs_error_label = new Label("Max absolute error:");
		abs_error_label.setBounds(x2,y,150, 20);
		abs_error_label.setForeground(Color.gray);
		add(abs_error_label);
		abs_error = new TextField("200");
		abs_error.setBounds(x2+150,y,100,20);
		abs_error.enable(errors.getState());
		add(abs_error);
		y+=30;
		
		generate = new Button("Generate");
		generate.setBounds(x1,y,250,20);
		generate.setBackground(svm.settings.button_color_default);
		generate.setForeground(svm.settings.button_label_default);	
		add(generate);	
		
		save = new Button("Save");
		save.setBounds(x2,y,250,20);
		save.setBackground(svm.settings.button_color_default);
		save.setForeground(svm.settings.button_label_default);	
		add(save);
		
		y+=30;
		ta = new TextArea("");
		ta.setBounds(x1,y,size().width-2*x1,size().height-y-x1);
		ta.setBackground(svm.settings.button_color_default);
		ta.setForeground(svm.settings.string_color_default);
		add(ta);		

		show();
	}
	
	public boolean handleEvent(Event e){
		if(e.id==Event.WINDOW_DESTROY){
			dispose();
		}else if(e.id==Event.ACTION_EVENT && e.target == generate){
			generateData();
           	return true;
		}else if(e.id==Event.ACTION_EVENT && e.target == save){
			saveGeneratedData();
           	return true;
		}else if(e.id==Event.ACTION_EVENT && e.target == liniar){
			if(liniar.getState()){ 
				classes_count.setText("2");
				errors.enable(true);
				errors_label.setForeground(Color.white);
				margin_label.setForeground(Color.white);
			}
			else{				
				errors.setState(false);
				errors.enable(false);
				errors_label.setForeground(Color.gray);
				abs_error_label.setForeground(Color.gray);
				abs_error.enable(false);
				margin_label.setForeground(Color.gray);


			}
			classes_count.enable(!liniar.getState());
			margin.enable(liniar.getState());
           	return true;			
		}else if(e.id==Event.ACTION_EVENT && e.target==errors){
			if(errors.getState()){
				abs_error_label.setForeground(Color.white);
				abs_error.enable(true);
			}else{
				abs_error_label.setForeground(Color.gray);
				abs_error.enable(false);

			}
			
		}else return false;	
		return super.handleEvent(e);
	}

	void generateData(){
		ta.setText("");
		if(liniar.getState()){
			float marg = 3f; 
				try{
					float mg = Float.parseFloat(margin.getText());
					if (mg>0) marg = mg;
				}
				catch(NumberFormatException e){}
			if(errors.getState()){
				float error_abs = 50f;
				try{
					float ea = Float.parseFloat(abs_error.getText());
					if (ea>0) error_abs = ea;
				}
				catch(NumberFormatException e){}
				generateLinearDataWithErrors(error_abs,marg);
			}
			else{
				generateLiniarData(marg);
			}
			//generateMajoritarData();
		}else
			generateRandomData();		
	}
	
	void generateRandomData(){
		int N = Integer.parseInt(vectors_count.getText());
		int n = Integer.parseInt(attributes_count.getText());
		int MIN = Integer.parseInt(min.getText());
		int MAX = Integer.parseInt(max.getText());
		int C = Integer.parseInt(classes_count.getText()); 
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%\n");
		ta.append("% attributes count = " + n + "\n");
		ta.append("% vectors count = " + N + "\n");
		ta.append("%\n");
		ta.append("% Edit the attribute names (attrib_1, attrib_2, ...) in this text box and then save.\n");
		ta.append("%\n");
		ta.append("@relation xxxxxxxxxxx\n");
		for(int i=1; i<=n; i++)
			ta.append("@attribute attrib_" + i + " numeric\n");
		String ss = "";
		for(int i=0; i<C-1; i++)
			ss += i + ", ";
		ss += (C-1);
		ta.append("@attribute class {" + ss + "}\n");
		ta.append("@data\n");		
		if(N > 1 && n > 1 && MIN < MAX && C > 1){
			for(int k=0; k<N; k++){
				String s = "";
				for(int i=0; i<n; i++){
					s += (MIN + (float)Math.random()*(MAX-MIN)) + ",";
				}
				s += (int)(Math.random()*C) + "\n";
				ta.append(s);
			}
		}else mesaj();
	}

	void generateLinearDataWithErrors(float error_abs, float marg){
		int N = Integer.parseInt(vectors_count.getText());
		int n = Integer.parseInt(attributes_count.getText());
		int MIN = Integer.parseInt(min.getText());
		int MAX = Integer.parseInt(max.getText());
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%\n");
		ta.append("% attributes count = " + n + "\n");
		ta.append("% vectors count = " + N + "\n");
		ta.append("%\n");
		ta.append("% Edit the attribute names (attrib_1, attrib_2, ...) in this text box and then save.\n");
		ta.append("%\n");
		ta.append("@relation xxxxxxxxxxx\n");
		for(int i=1; i<=n; i++)
			ta.append("@attribute attrib_" + i + " numeric\n");
		ta.append("@attribute class {0, 1}\n");
		ta.append("@data\n");
		if(N > 1 && n > 1 && MIN < MAX){
			float[] w = new float[n+1];
			float max = Float.MIN_VALUE;
			for(int i=0; i<n; i++) {
				w[i] = MIN + (float)Math.random()*(MAX-MIN);
				if(w[i] > max) max = w[i];
			}
			w[n] = max*(MIN + (float)Math.random()*(MAX-MIN));
			for(int k=0; k<N; k++){
				while(true){
					String s = "";
					float[] x = new float[n];
					for(int i=0; i<n; i++)
						x[i] = MIN + (float)Math.random()*(MAX-MIN);
					float z = 0;
					for(int i=0; i<n; i++) 
						z += w[i]*x[i];
					z += w[n];
					int y = -1;
					if(z <= w[n]) y = 0;
					else if(z >= w[n]) y = 1;
					
					float comp_error =  (float)error_abs/(float)Math.sqrt((float)n);
					
					for(int i=0; i<n; i++)
						x[i] +=Math.random()*2*comp_error-comp_error;
					for(int i=0; i<n; i++) s += x[i] + ",";				
					s += y + "\n";	
					ta.append(s);
					break;
				}
			}
		}else mesaj();
	}

	void generateLiniarData(float marg){ //the margin: marg>=0
		int N = Integer.parseInt(vectors_count.getText());
		int n = Integer.parseInt(attributes_count.getText());
		int MIN = Integer.parseInt(min.getText());
		int MAX = Integer.parseInt(max.getText());
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%\n");
		ta.append("% attributes count = " + n + "\n");
		ta.append("% vectors count = " + N + "\n");
		ta.append("%\n");
		ta.append("% Edit the attribute names (attrib_1, attrib_2, ...) in this text box and then save.\n");
		ta.append("%\n");
		ta.append("@relation xxxxxxxxxxx\n");
		for(int i=1; i<=n; i++)
			ta.append("@attribute attrib_" + i + " numeric\n");
		ta.append("@attribute class {0, 1}\n");
		ta.append("@data\n");
		if(N > 1 && n > 1 && MIN < MAX){
			float[] w = new float[n+1];
			float max = Float.MIN_VALUE;
			for(int i=0; i<n; i++) {
				w[i] = MIN + (float)Math.random()*(MAX-MIN);
				if(w[i] > max) max = w[i];
			}
			w[n] = max*(MIN + (float)Math.random()*(MAX-MIN));
			for(int k=0; k<N; k++){
				while(true){
					String s = "";
					float[] x = new float[n];
					for(int i=0; i<n; i++)
						x[i] = MIN + (float)Math.random()*(MAX-MIN);
					float z = 0;
					for(int i=0; i<n; i++) 
						z += w[i]*x[i];
					z += w[n];
					int y = -1;
					if(z <= -max*marg) y = 0;
					else if(z >= max*marg) y = 1;
					else continue;
					for(int i=0; i<n; i++) s += x[i] + ",";				
					s += y + "\n";	
					ta.append(s);
					break;
				}
			}
		}else mesaj();
	}
	
	
	void generateMajoritarData(){
		int N = Integer.parseInt(vectors_count.getText());
		int n = Integer.parseInt(attributes_count.getText());
		int MIN = Integer.parseInt(min.getText());
		int MAX = Integer.parseInt(max.getText());
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%............comments............\n");
		ta.append("%\n");
		ta.append("% attributes count = " + n + "\n");
		ta.append("% vectors count = " + N + "\n");
		ta.append("%\n");
		ta.append("% Edit the attribute names (attrib_1, attrib_2, ...) in this text box and then save.\n");
		ta.append("%\n");
		ta.append("@relation xxxxxxxxxxx\n");
		for(int i=1; i<=n; i++)
			ta.append("@attribute attrib_" + i + " numeric\n");
		ta.append("@attribute class {0, 1}\n");
		ta.append("@data\n");
		if(N > 1 && n > 1 && MIN < MAX){			
			for(int k=0; k<N; k++){
				String s = "";
				float[] x = new float[n];
				int p = 0;
				for(int i=0; i<n; i++){
					x[i] = MIN + (float)Math.random()*(MAX-MIN);
					s += x[i] + ",";
					if(x[i]>0) p++;
				}
				int y = 1;
				if(p<n/2) y = 0;
				s += y + "\n";	
				ta.append(s);
			}
		}else mesaj();
	}	
	
	
	void saveGeneratedData(){
		if(!ta.getText().equals("")){
			try{ 	
				FileDialog fd=new FileDialog(this, "Save Generated Input Data", 1);
				if(dir!=null) fd.setDirectory(dir);
				fd.setFile("*.csv");
				fd.setVisible(true);
				if(fd.getFile() != null) {
					dir = fd.getDirectory();
					String fisier = fd.getFile();
					path = dir + fisier;
					File file = new File(path); 
					BufferedWriter bw = null;
					if(file.exists()) file.delete();
					try{
						bw = new BufferedWriter(new FileWriter(file));
						bw.write(ta.getText());
						bw.close();
					}
					catch(IOException e){e.printStackTrace();}				
				}
			}
			catch(Exception e) {e.printStackTrace();}	
		}		
	}
		
	void mesaj(){
		System.out.println("Vectors Dimension must be > 1.");
		System.out.println("Vectors Count must be > 1.");
		System.out.println("It is necessary that Minimum Coordonates < Maximum Coordonates.");		
	}
	
}