import java.awt.Point;
import java.math.BigDecimal;
import java.util.ArrayList;


public class Trajectory {
	private double x,y,dydx,d2ydx2,A,C,a,b,step;
	//private double a,b;
	//private BigDecimal x,y,dydx,d2ydx2,A,C;
	private Complex a1,a2,a3,a4,b1,b2,b3,b4;
	public Trajectory(double X, double Y, double Velocity, double Acceleration){
		x = X;
		y = Y;
		dydx = Velocity;
		d2ydx2 = Acceleration;
		//A = new BigDecimal(-1).multiply(dydx.pow(4)).subtract(new BigDecimal(2).multiply(y.multiply(dydx.pow(2).multiply(d2ydx2)))).subtract(d2ydx2.pow(2).multiply(y.pow(2))).subtract(dydx.pow(6)).subtract(new BigDecimal(2).multiply(y).multiply(dydx.pow(4)).multiply(d2ydx2)).subtract(new BigDecimal(3).multiply(y.pow(2)).multiply(dydx.pow(2).multiply(d2ydx2.pow(2)))).subtract(d2ydx2.multiply(y).multiply(dydx.pow(4))).subtract(d2ydx2.pow(3).multiply(y.pow(3)));
		//C = (y.pow(2).multiply(dydx.pow(2))).subtract(new BigDecimal(2).multiply(x).multiply(y).multiply((dydx.pow(2)).add(d2ydx2.multiply(y))).multiply(dydx)).add(x.pow(2).multiply((dydx.pow(4)).add(new BigDecimal(2).multiply(y).multiply(dydx.pow(2)).multiply(d2ydx2)).add(d2ydx2.pow(2).multiply(y.pow(2)))));
		A = -(dydx*dydx*dydx*dydx)-(2*y*dydx*dydx*d2ydx2)-(d2ydx2*d2ydx2*y*y)-(dydx*dydx*dydx*dydx*dydx*dydx)-(2*y*dydx*dydx*dydx*dydx*d2ydx2)-(3*y*y*dydx*dydx*d2ydx2*d2ydx2)-(d2ydx2*y*dydx*dydx*dydx*dydx)-(d2ydx2*d2ydx2*d2ydx2*y*y*y);
		C = (y*y*dydx*dydx)-(2*x*y*((dydx*dydx)+(d2ydx2*y))*dydx)+(x*x*((dydx*dydx*dydx*dydx)+(2*y*dydx*dydx*d2ydx2)+(d2ydx2*d2ydx2*y*y)));
		//System.out.println(A);
		//System.out.println(C);
		a1 = (new Complex(-1, 0).power(0.25).times(new Complex (C,0).integerRoot(4))).divides(new Complex (A,0).integerRoot(4)).times(new Complex(-1,0));
		a2 = (new Complex(-1, 0).power(0.25).times(new Complex (C,0).integerRoot(4))).divides(new Complex (A,0).integerRoot(4));
		a3 = (new Complex(-1, 0).power(0.75).times(new Complex (C,0).integerRoot(4))).divides(new Complex (A,0).integerRoot(4)).times(new Complex(-1,0));
		a4 = (new Complex(-1, 0).power(0.75).times(new Complex (C,0).integerRoot(4))).divides(new Complex (A,0).integerRoot(4));

		//System.out.println(a1);
		//System.out.println(a2);
		//System.out.println(a3);
		//System.out.println(a4);
		//if (a1.im() == 0 && a1.re() > 0) a = a1.re()*a1.re();
		a = a2.re()*a2.re();
		//if (a2.im() == 0 && a2.re() > 0) a = a2.re()*a2.re();
		//if (a3.im() == 0 && a3.re() > 0) a = a3.re()*a3.re();
		//if (a4.im() == 0 && a4.re() > 0) a = a4.re()*a4.re();
		//b = (new BigDecimal(-1).multiply(new BigDecimal(a).pow(2)).multiply(dydx.pow(2))).subtract(d2ydx2.multiply(y).multiply(new BigDecimal(a).pow(2))).doubleValue();
		b = Math.sqrt((-1*(a*a*dydx*dydx))-(d2ydx2*y*a*a));
		//System.out.println("A = "+a+"\r\nB = "+b);
	}
	public void setPrecision(double Precision){
		step = Precision;
	}
	public ArrayList <Point> getPoints(){
		ArrayList <Point> j = new ArrayList<>();
		double x1 = 0;
		while (x1 >0)
			x1 += step;
		return null;
		
	}
	public String toString() {
		return "A = "+a+"\r\nB = "+b;
		
	}
	private class Complex {
	    private final double re;   // the real part
	    private final double im;   // the imaginary part
	    public boolean infinite;

	    // create a new object with the given real and imaginary parts
	    public Complex(BigDecimal c, double imag) {
	        re = c.doubleValue();
	        im = imag;
	    }
	 // create a new object with the given real and imaginary parts
	    public Complex(double c, double imag) {
	        re = c;
	        im = imag;
	    }
	    public Complex log() {
	        double modulus = Math.sqrt(re*re + im*im);
	        double arg = Math.atan2(im,re);
	        return new Complex(Math.log(modulus), arg);
	      }
	    // return a string representation of the invoking Complex object
	    public String toString() {
	        if (im == 0) return re + "";
	        if (re == 0) return im + "i";
	        if (im <  0) return re + " - " + (-im) + "i";
	        return re + " + " + im + "i";
	    }
	    /**
	     * Returns the absolute value, "r" in polar coordinates, of this.
	     * @return the square root of (real part squared plus imaginary part squared)
	     */
	    public double r() {
	      return Math.sqrt(re*re + im*im);
	    }
	    
	    /**
	     * Returns arg(this), the angular polar coordinate of this complex number, in the range -pi to pi.
	     * The return value is simply Math.atan2(imaginary part, real part).
	     */
	    public double theta() {
	      return Math.atan2(im,re);
	    }
	    /**
	     * Returns a complex k-th root of this complex number.  The root that is returned is 
	     * the one with the smallest positive arg.
	     * (If k is 0, the return value is 1.  If k is negative, the value is 1/integerRoot(-k).)
	     */
	    public Complex integerRoot(int k) {
	      double a,b;
	      boolean neg = false;
	      if (k < 0) {
	        k = -k;
	        neg = true;
	      }
	      if (k == 0) { 
	        a = 1;
	        b = 0;
	      }
	      else if (k == 1) {
	        a = re;
	        b = im;
	      }
	      else {
	        double length = r();
	        double angle = theta();
	        if (angle < 0)
	          angle += Math.PI*2;
	        length = Math.pow(length,1.0/k);
	        angle = angle / k;
	        a = length*Math.cos(angle);
	        b = length*Math.sin(angle);
	      }
	      if (neg) {
	        double denom = a*a + b*b;
	        a = a/denom;
	        b = -b/denom;
	      }
	      return new Complex(a,b);
	    }
	    // return abs/modulus/magnitude and angle/phase/argument
	    public double abs()   { return Math.hypot(re, im); }  // Math.sqrt(re*re + im*im)
	    public double phase() { return Math.atan2(im, re); }  // between -pi and pi
	    public Complex power(double x) {
	        double modulus = Math.sqrt(re*re + im*im);
	        double arg = Math.atan2(im,re);
	        double log_re = Math.log(modulus);
	        double log_im = arg;
	        double x_log_re = x * log_re;
	        double x_log_im = x * log_im;
	        double modulus_ans = Math.exp(x_log_re);
	        return new Complex(modulus_ans*Math.cos(x_log_im), modulus_ans*Math.sin(x_log_im));
	      }
	    // return a new Complex object whose value is (this + b)
	    public Complex plus(Complex b) {
	        Complex a = this;             // invoking object
	        double real = a.re + b.re;
	        double imag = a.im + b.im;
	        return new Complex(real, imag);
	    }

	    // return a new Complex object whose value is (this - b)
	    public Complex minus(Complex b) {
	        Complex a = this;
	        double real = a.re - b.re;
	        double imag = a.im - b.im;
	        return new Complex(real, imag);
	    }

	    // return a new Complex object whose value is (this * b)
	    public Complex times(Complex b) {
	        Complex a = this;
	        double real = a.re * b.re - a.im * b.im;
	        double imag = a.re * b.im + a.im * b.re;
	        return new Complex(real, imag);
	    }

	    // scalar multiplication
	    // return a new object whose value is (this * alpha)
	    public Complex times(double alpha) {
	        return new Complex(alpha * re, alpha * im);
	    }

	    // return a new Complex object whose value is the conjugate of this
	    public Complex conjugate() {  return new Complex(re, -im); }

	    // return a new Complex object whose value is the reciprocal of this
	    public Complex reciprocal() {
	        double scale = re*re + im*im;
	        return new Complex(re / scale, -im / scale);
	    }

	    // return the real or imaginary part
	    public double re() { return re; }
	    public double im() { return im; }

	    // return a / b
	    public Complex divides(Complex b) {
	        Complex a = this;
	        return a.times(b.reciprocal());
	    }

	    // return a new Complex object whose value is the complex exponential of this
	    public Complex exp() {
	        return new Complex(Math.exp(re) * Math.cos(im), Math.exp(re) * Math.sin(im));
	    }

	    // return a new Complex object whose value is the complex sine of this
	    public Complex sin() {
	        return new Complex(Math.sin(re) * Math.cosh(im), Math.cos(re) * Math.sinh(im));
	    }

	    // return a new Complex object whose value is the complex cosine of this
	    public Complex cos() {
	        return new Complex(Math.cos(re) * Math.cosh(im), -Math.sin(re) * Math.sinh(im));
	    }

	    // return a new Complex object whose value is the complex tangent of this
	    public Complex tan() {
	        return sin().divides(cos());
	    }
	    


	    // a static version of plus
	    public Complex plus(Complex a, Complex b) {
	        double real = a.re + b.re;
	        double imag = a.im + b.im;
	        Complex sum = new Complex(real, imag);
	        return sum;
	    }
	}
}
