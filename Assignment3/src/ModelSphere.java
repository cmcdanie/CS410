
public class ModelSphere {
	
	public double[] center;
	public double r;
	
	public double[] ambient;
	public double[] diffuse;
	public double[] specular;
	public double[] attenuation;
	
	public ModelSphere() {
		center = new double[3];
		ambient = new double[3];
		diffuse = new double[3];
		specular = new double[3];
		attenuation = new double[3];
		//setCenter(cx, cy, cz);
		//r = radius;
	}
	
	public void setCenter(double cx, double cy, double cz) {
		center[0] = cx;
		center[1] = cy;
		center[2] = cz;
	}
	
	public void setRadius(double radius) {
		r = radius;
	}
	
	public void setAmbient(double Kar, double Kag, double Kab) {
		ambient[0] = Kar;
		ambient[1] = Kag;
		ambient[2] = Kab;
	}
	
	public void setDiffuse(double Kdr, double Kdg, double Kdb) {
		diffuse[0] = Kdr;
		diffuse[1] = Kdg;
		diffuse[2] = Kdb;
	}
	
	public void setSpecular(double Ksr, double Ksg, double Ksb) {
		specular[0] = Ksr;
		specular[1] = Ksg;
		specular[2] = Ksb;
	}
	
	public void setAttenuation(double Krr, double Krg, double Krb) {
		attenuation[0] = Krr;
		attenuation[1] = Krg;
		attenuation[2] = Krb;
	}
	
	public double[] getCenter() {
		return center;
	}
	
	public double getRadius() {
		return r;
	}
	
	public double[] getAmbient() {
		return ambient;
	}
	
	public double[] getDiffuse() {
		return diffuse;
	}
	
	public double[] getSpecular() {
		return specular;
	}
	
	public double[] getAttenuation() {
		return attenuation;
	}
}
