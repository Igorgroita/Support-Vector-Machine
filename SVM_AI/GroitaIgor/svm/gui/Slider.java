package gui;

import java.awt.*;
import svm.SVM;

public class Slider extends Panel implements Runnable {
    About about;
    Thread t;
    Image im;
    Graphics img;
    int w1, h1;
    int delay = 50;
    String[] students;
    Font fnt = new Font("Arial", 0, 12);
    int y, w = 10;
    boolean init = true, control;
    int ww, hh;

    public Slider(About about, int ww, int hh) {
    	this.about = about;
		this.ww = ww;
		this.hh = hh;    		
    	setStudents();
		setFont(fnt); 
		start();
    }

	public void setStudents(){
		String[] student = new String[16];
		student[0] = "Alexa Drago\u015F-Constantin";
		student[1] = "Anton Ioana-Andreea";
		student[2] = "Butnaru Gheorghi\u021B\u0103";
		student[3] = "Cimpoe\u0219u Tudor-Andrei";
		student[4] = "David Andreea";
		student[5] = "Dalnicenco Elizaveta";
		student[6] = "Dr\u0103gu\u015F Teodora-Andreea";
		student[7] = "G\u00EEnga Raluca-Andreea";
		student[8] = "Groi\u021B\u0103 Igor";
		student[9] = "Iacob Eduard";
		student[10] = "Mitrea Gabriel-Claudiu";
		student[11] = "Negru \u0218tefan";
		student[12] = "Pardo Santi Romero";
		student[13] = "Purice Crina-Otilia";
		student[14] = "Simco Sorin";
		student[15] = "";
		
		students = new String[student.length*100];
		for(int i = 0; i < students.length; i++)
			students[i] = student[i % student.length];
	}
	
    public void start() {
    	if(t == null){
    		t = new Thread(this); 
    		t.start();
	        try{Thread.sleep(1000);}
			catch(InterruptedException e) { }    			
    	}
    }
    	
    public void stop() {if(t != null){ t.stop(); t = null;}}	
	
    public void run() {
	    do {
	        repaint();
			try {Thread.sleep(delay);}
			catch(InterruptedException e) {return;}
	    } while(true);
    }  	
	
	public void reset(){
		y = hh + 10;	
		repaint();
		stop();
	}

    public final void paint(Graphics g) {
    	if(init){
			im = createImage(ww, hh);
			img = im.getGraphics();	
			for(int i = 0; i < students.length; i++) {
				FontMetrics fm = img.getFontMetrics(fnt);
				h1 += fm.getHeight();
				if(fm.stringWidth(students[i]) > w1) w1 = fm.stringWidth(students[i]);
			}
			y = hh + 10; 	
			init = false;
		}
		Color color = null;
		for(int l = 0; l < students.length; l++) {
			float f = (float)hh / 4.0F;
			float f1 = 1.0F, f2 = 1.0F, f3 = 1.0F;
			int i1 = y + (int)(1.5 * getFont().getSize() * l);
			if(i1 >= 0 && i1 <= hh) {
				float ff = 0;
				if((float)i1 <= f)
					ff = (float)i1 / f;
				else if((float)i1 >= (float)hh - f)
					ff = ((float)hh - (float)i1) / f;
				else
					ff = 1.0F;
				color = new Color((int)((float)255 * ff), (int)((float)255 * ff), (int)((float)255 * ff));
			}else color = new Color(0, 0, 0);
			img.setColor(color);
			img.setFont(fnt);
			img.drawString(students[l], w, i1);
	    }
	    g.drawImage(im, 0, 0, this);
    }

    public final void update(Graphics g) {   
    	if(img!=null){	
			img.setColor(Color.black);
			img.fillRect(0,0,ww,hh);		
		}
		if(y < -(int)((float)h1*0.75f))
			y = hh + 10;
		else
			y --;
		paint(g);
	}		

}
