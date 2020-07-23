package gui;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import svm.SVM;

public class Design extends Panel implements MouseListener, MouseMotionListener{
	SVM svm;
	Image im;
	Graphics img;
	int ww, hh;
	int Ox, Oy, cx, cy, ccx, ccy;
	boolean init = true;
	String coords = "";
	public boolean show_coords = false;
	public boolean show_line = false;
	public boolean show_lines = false;
	public boolean calculates = false;
	public int x1, y1, x2, y2;
	public int x01, y01, x02, y02;
	public int x11, y11, x12, y12;
	
	
	public Design(SVM svm){
		this.svm = svm;
		addMouseListener(this);
		addMouseMotionListener(this);		
	}
	
	public void initO(){
		Ox = ww/2; Oy = hh/2;
		cx = 0; cy = 0;
		repaint();
	}
	
	public void paint(Graphics g){update(g);}
	
	public void update(Graphics g){
		if(init){
			ww = size().width;
			hh = size().height;
			im = createImage(ww, hh);
			img = im.getGraphics();	
			initO();
			init = false;			
		}
		
		img.setColor(svm.settings.background_color);
		img.fillRect(0,0,ww,hh);	
		if(svm.ind.V != null && svm.ind.V[0].getDimension() == 2){
			if(show_lines){
				img.setColor(svm.settings.background_color.darker());
				int[] x = {Ox+x01,Ox+x11,Ox+x12,Ox+x02};
				int[] y = {Oy-y01,Oy-y11,Oy-y12,Oy-y02}; 
				img.fillPolygon(x, y, 4);			
			}				
			
			drawAxis(img);
			
			if(show_line){
				img.setColor(svm.settings.line_color);
				img.drawLine(Ox+x1,Oy-y1,Ox+x2,Oy-y2);
			}
			
			if(show_lines){
				img.setColor(svm.settings.line_color);
				img.drawLine(Ox+x1,Oy-y1,Ox+x2,Oy-y2);
				
				img.setColor(svm.settings.line_color.darker());
				img.drawLine(Ox+x01,Oy-y01,Ox+x02,Oy-y02);	

				img.setColor(svm.settings.line_color.darker());
				img.drawLine(Ox+x11,Oy-y11,Ox+x12,Oy-y12);				
			}			
			
			for(int i = 0; i < svm.ind.V.length; i++){
				Point p = new Point(Ox + (int)(svm.ind.V[i].X[0]+0.5), Oy - (int)(svm.ind.V[i].X[1]+0.5));
				img.setColor(svm.ind.V[i].cl.color);
				int r = svm.settings.point_radius;
				img.fillOval(p.x-r,p.y-r,2*r,2*r);
				img.setColor(Color.black);
				img.drawOval(p.x-r,p.y-r,2*r,2*r);			
			}
			
			if(show_coords){
				img.setColor(svm.settings.string_color);
				img.drawString(coords, ccx+15, ccy+30);
			}

		}
		
		if(calculates) img.drawImage(svm.calculates, (ww-svm.calculates.getWidth(this))/2, (hh-svm.calculates.getHeight(this))/2, this);
		
		g.drawImage(im,0,0,this);
	}

	public void drawAxis(Graphics g){
		//deseneaza gridul
		if(svm.settings.grid){
			g.setColor(svm.settings.grid_color);	
			for(int i=Ox+svm.settings.axis_min; i<=Ox+svm.settings.axis_max; i+=svm.settings.grid_size) g.drawLine(i, Oy+svm.settings.axis_min, i, Oy+svm.settings.axis_max);
			for(int j=Oy+svm.settings.axis_min; j<=Oy+svm.settings.axis_max; j+=svm.settings.grid_size) g.drawLine(Ox+svm.settings.axis_min, j, Ox+svm.settings.axis_max, j);
		}
		//deseneaza axele
		if(svm.settings.axis){
			g.setColor(svm.settings.axis_color);		
			g.drawLine(Ox+svm.settings.axis_min, Oy, Ox+svm.settings.axis_max, Oy);
			g.drawLine(Ox, Oy+svm.settings.axis_min, Ox, Oy+svm.settings.axis_max);	
			//deseneaza gradatiile
			if(svm.settings.gradations){
				for(int i=Ox; i<=Ox+svm.settings.axis_max; i+=svm.settings.axis_gradations) g.drawLine(i, Oy-2, i, Oy+2);
				for(int i=Ox; i>=Ox+svm.settings.axis_min; i-=svm.settings.axis_gradations) g.drawLine(i, Oy-2, i, Oy+2);
				for(int j=Oy; j<=Oy+svm.settings.axis_max; j+=svm.settings.axis_gradations) g.drawLine(Ox-2, j, Ox+2, j);
				for(int j=Oy; j>=Oy+svm.settings.axis_min; j-=svm.settings.axis_gradations) g.drawLine(Ox-2, j, Ox+2, j);
			}
		}
	}

	public void setPointsOfLine(float[] w){
		show_line = true;
		if(Math.abs(w[0]) < Math.abs(w[1])){
			x1 = svm.settings.axis_min;
			y1 = (int)((-w[2]-w[0]*x1)/w[1]+0.5);
			x2 = svm.settings.axis_max;
			y2 = (int)((-w[2]-w[0]*x2)/w[1]+0.5);
		}else{
			y1 = svm.settings.axis_min;
			x1 = (int)((-w[2]-w[1]*y1)/w[0]+0.5);
			y2 = svm.settings.axis_max;
			x2 = (int)((-w[2]-w[1]*y2)/w[0]+0.5);			
		}
		repaint();
	}	
	
	public void setPointsOfMaxLine(float[] w, float[] w0, float[] w1){
		show_lines = true;
		if(Math.abs(w[0]) < Math.abs(w[1])){
			x1 = svm.settings.axis_min;
			y1 = (int)((-w[2]-w[0]*x1)/w[1]+0.5);
			x2 = svm.settings.axis_max;
			y2 = (int)((-w[2]-w[0]*x2)/w[1]+0.5);
		}else{
			y1 = svm.settings.axis_min;
			x1 = (int)((-w[2]-w[1]*y1)/w[0]+0.5);
			y2 = svm.settings.axis_max;
			x2 = (int)((-w[2]-w[1]*y2)/w[0]+0.5);			
		}
		
		if(Math.abs(w0[0]) < Math.abs(w0[1])){
			x01 = svm.settings.axis_min;
			y01 = (int)((-w0[2]-w0[0]*x01)/w0[1]+0.5);
			x02 = svm.settings.axis_max;
			y02 = (int)((-w0[2]-w0[0]*x02)/w0[1]+0.5);
		}else{
			y01 = svm.settings.axis_min;
			x01 = (int)((-w0[2]-w0[1]*y01)/w0[0]+0.5);
			y02 = svm.settings.axis_max;
			x02 = (int)((-w0[2]-w0[1]*y02)/w0[0]+0.5);			
		}		
		
		if(Math.abs(w1[0]) < Math.abs(w1[1])){
			x11 = svm.settings.axis_min;
			y11 = (int)((-w1[2]-w1[0]*x11)/w1[1]+0.5);
			x12 = svm.settings.axis_max;
			y12 = (int)((-w1[2]-w1[0]*x12)/w1[1]+0.5);
		}else{
			y11 = svm.settings.axis_min;
			x11 = (int)((-w1[2]-w1[1]*y11)/w1[0]+0.5);
			y12 = svm.settings.axis_max;
			x12 = (int)((-w1[2]-w1[1]*y12)/w1[0]+0.5);			
		}		
		
		repaint();
	}		
	
	public void mouseClicked(MouseEvent me) {initO();}

	public void mouseEntered(MouseEvent me) {}

	public void mouseExited(MouseEvent me) {}

	public void mouseMoved(MouseEvent me) {
		ccx = me.getX(); ccy = me.getY();
		coords = "(" + (ccx-Ox) + "," + (Oy-ccy) + ")";
		if(ccx<=2 || ccx >= ww-5 || ccy<=5 || ccy>=hh-5) coords = "";
		repaint();
	} 
		
	public void mousePressed(MouseEvent me) { 
		cx = me.getX(); cy = me.getY();  
		coords = "";		
	}
			
	public void mouseDragged(MouseEvent me) {
		int x = me.getX(), y = me.getY();
		if(svm.ind.V != null && svm.ind.V[0].getDimension() == 2){
			cx = x - cx; cy = y - cy;  
			Ox += cx; Oy += cy;
			cx = x; cy = y;
			coords = "";
			repaint();	
		}
	}

	public void mouseReleased(MouseEvent me) {
		int x = me.getX(), y = me.getY(); 
	}	
	
}