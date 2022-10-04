package maximal_likelihood;

import java.util.Scanner;

public class Gradient_Descent_positioning {
	public static void main(String[] args) {
		Scanner scanner = new Scanner(System.in);
		int ap = scanner.nextInt();
		
		double[] x = new double[ap+1];
		double[] y = new double[ap+1];
		double[] z = new double[ap+1];
		double[] r = new double[ap+1];

		double [] c = new double[ap+1];
		for(int i = 1; i <= ap; i++) {
			x[i] = scanner.nextDouble();
			y[i] = scanner.nextDouble();
			z[i] = scanner.nextDouble();
			c[i] = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
		}
		double[][] xx2 = new double[ap+1][ap+1];
		double[][] yy2 = new double[ap+1][ap+1];
		double[][] zz2 = new double[ap+1][ap+1];
		double[][] cc = new double[ap+1][ap+1];
		for(int i = 1; i <= ap; i++) {
			for(int j = 1; j <= ap; j++) {
				xx2[i][j] = (x[i] - x[j])*2;
				yy2[i][j] = (y[i] - y[j])*2;
				zz2[i][j] = (z[i] - z[j])*2;
				cc[i][j] = c[i] - c[j];
			}
		}
		
		for(int i = 1; i <= ap; i++) {
			r[i] = scanner.nextDouble();
		}
		double small = 0.000001;
		int limit = 1000000;
		double rate = 0.0;
		double tx = 0, ty = 0, tz = 0;
		double gx = 0, gy = 0, gz = 0, gl = 1;
		double e = 0;
		final int window=1000,interval=10000;
		double[] eh1 = new double[window];
		double[] eh2 = new double[window];
		int t = 0;
		do {
			rate = 1/gl*Math.log(e/10+1);
			tx -= gx*rate;
			ty -= gy*rate;
			tz -= gz*rate;
			gx = 0;
			gy = 0;
			gz = 0;
			e = 0;
			for(int i = 1; i <= ap; i++) {
				double raw = tx*tx - 2*x[i]*tx + ty*ty - 2*y[i]*ty + tz*tz - 2*z[i]*tz + (c[i] - r[i]*r[i]);
				double share = 2*raw/(raw*raw + 1);
				gx += share*2*(tx - x[i]);
				gy += share*2*(ty - y[i]);
				gz += share*2*(tz - z[i]);
				e += Math.log(raw*raw + 1);
				for(int j = 1; j <= ap; j++) {
					if(i < j) {
						raw = xx2[i][j]*tx + yy2[i][j]*ty + zz2[i][j]*tz + (cc[j][i] + r[i]*r[i] - r[j]*r[j]);
						share = 2*raw/(raw*raw + 1);
						gx += share*xx2[i][j];
						gy += share*yy2[i][j];
						gz += share*zz2[i][j];
						e += Math.log(raw*raw + 1);
					}
				}
			}
			e /= (ap + ap*(ap-1)/2);
			gl = Math.sqrt(gx*gx + gy*gy + gz*gz);
			t++;
			eh1[t%window] = e;
			if(t%interval == 0) {
				double m1=0,m2=0;
				for(int i=0;i<window;i++) {
					m1+=eh1[i];
					m2+=eh2[i];
				}
				m1/=window;
				m2/=window;
				if(Math.abs(m2-m1) < small) {
					break;
				}
				for(int i=0;i<window;i++) {
					eh2[i] = eh1[i];
				}
			}
		}while(t < limit);
		System.out.println(tx + " " + ty + " " + tz);
	}
}
