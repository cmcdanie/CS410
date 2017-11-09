
public class Light {
	
	public double[] location;
	public double W;
	public double[] emission;
	
	public Light() {
		location = new double[3];
		emission = new double[3];
	}
	
	public void setLocation(double x, double y, double z) {
		location[0] = x;
		location[1] = y;
		location[2] = z;
	}
	
	public void setW(double w) {
		W = w;
	}
	
	public void setEmission(double Kr, double Kg, double Kb) {
		emission[0] = Kr;
		emission[1] = Kg;
		emission[2] = Kb;
	}
	
	public double[] getLocation() {
		return location;
	}
	
	public double getW() {
		return W;
	}
	
	public double[] getEmission() {
		return emission;
	}
	
}

