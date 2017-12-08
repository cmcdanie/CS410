import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class MaterialObject {
	
	private String materialName;
	private RealVector ambient;
	private RealVector diffuse;
	private RealVector specular;
	private RealVector attenuation;
	private double Ns;
	private double Ko;
	private double eta;
	
	MaterialObject(String n){
		ambient = new ArrayRealVector(3);
		diffuse = new ArrayRealVector(3);
		specular = new ArrayRealVector(3);
		attenuation = new ArrayRealVector(new double[]{1, 1, 1,}) ;
		materialName = n;
		Ns = 0;
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
	
	public void setNs(double ns) {
		Ns = ns;
	}
	
	public void setOpacity(double d) {
		Ko = d;
	}
	
	public void setEta(double indexOfRefraction){
		eta = indexOfRefraction;
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
	
	public double getNs() {
		return Ns;
	}
	
	public double getOpacity() {
		return Ko;
	}
	
	public double getEta(){
		return eta;
	}
	
	public String getName() {
		return materialName;
	}
	
}
