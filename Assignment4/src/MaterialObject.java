import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class MaterialObject {
	
	public String materialName;
	public RealVector ambient;
	public RealVector diffuse;
	public RealVector specular;
	public RealVector reflection;
	public double Ns;
	
	MaterialObject(String n){
		ambient = new ArrayRealVector(3);
		diffuse = new ArrayRealVector(3);
		specular = new ArrayRealVector(3);
		reflection = new ArrayRealVector(new double[]{1, 1, 1,}) ;
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
	
	public void setReflect(RealVector Kr) {
		reflection = Kr;
	}
	
	public void setNs(double ns) {
		Ns = ns;
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
	
	public RealVector getReflect() {
		return reflection;
	}
	
	public double getNs() {
		return Ns;
	}
	
	public String getName() {
		return materialName;
	}
	
}
