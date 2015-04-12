
public class Main {
	public static void main(String[] args) {
		double a = 502;
		double b = 500;
		double x = Math.sqrt((a*a)-(b*b));
		double y = b;
		for (int i = 0; i < 6; i++) {
			a++;
			x = Math.sqrt((a*a)-(b*b));
			y = b;			
			System.out.println("@ ("+x+","+y+"):"); 
			System.out.println("A really = "+a);
			System.out.println("B really = "+b);
			
			Trajectory l = new Trajectory(x, y, getVelocity(a,b,x,y), getAcceleration(a,b,x,y));
			System.out.println(l);
			
		}
	}
	public static double getVelocity(double a, double b,double x,double y) {
		return (b*b*(Math.sqrt((a*a)-(b*b)) - x))/(a*a*y);
	}
	public static double getAcceleration(double a, double b, double x, double y) {
		System.out.println("Velocity: "+getVelocity(a,b,x,y)+"\nAcceleration: "+(((-b*b)/(a*a))-(getVelocity(a,b,x,y)*getVelocity(a,b,x,y))/y));
		return (((-1*b*b)/(a*a))-(getVelocity(a,b,x,y)*getVelocity(a,b,x,y)))/y;
	}
}
 