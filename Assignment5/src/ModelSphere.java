import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class ModelSphere {
	
	private RealVector center;
	private double r;
	
	private RealVector ambient;
	private RealVector diffuse;
	private RealVector specular;
	private RealVector attenuation;
	private RealVector opacity;
	private double eta;
	
	public ModelSphere() {
		center = new ArrayRealVector(3);
		ambient = new ArrayRealVector(3);
		diffuse = new ArrayRealVector(3);
		specular = new ArrayRealVector(3);
		attenuation = new ArrayRealVector(3);
		opacity = new ArrayRealVector(3);
	}
	
	public void setCenter(RealVector cntr) {
		center = cntr;
	}
	
	public void setRadius(double radius) {
		r = radius;
	}
	
	public void setAmbient(RealVector Ka) {
		ambient = Ka;
	}
	
	public void setDiffuse(RealVector Kd) {
		diffuse = Kd;
	}
	
	public void setSpecular(RealVector Ks) {
		specular = Ks;
	}
	
	public void setAttenuation(RealVector Kr) {
		attenuation = Kr;
	}
	
	public void setOpacity(RealVector Ko){
		opacity = Ko;
	}
	
	public void setEta(double indexOfReflect){
		eta = indexOfReflect;
	}
	
	public RealVector getCenter() {
		return center;
	}
	
	public double getRadius() {
		return r;
	}
	
	public RealVector getAmbient() {
		return ambient;
	}
	
	public RealVector getDiffuse() {
		return diffuse;
	}
	
	public RealVector getSpecular() {
		return specular;
	}
	
	public RealVector getAttenuation() {
		return attenuation;
	}
	
	public RealVector getOpacity(){
		return opacity;
	}
	
	public double getEta(){
		return eta;
	}
}
