
public class ModelSphere {
	
	public double[] c;
	public double r;
	
	public ModelSphere() {
		c = new double[3];
		//setCenter(cx, cy, cz);
		//r = radius;
	}
	
	public void setCenter(double cx, double cy, double cz) {
		c[0] = cx;
		c[1] = cy;
		c[2] = cz;
	}
	
	public void setRadius(double radius) {
		r = radius;
	}
	
	public double[] getCenter() {
		return c;
	}
	
	public double getRadius() {
		return r;
	}
}
