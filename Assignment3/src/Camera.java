import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class Camera {

	public static RealVector eye;
	public static RealVector look;
	public static RealVector up;
	public static double d;
	public static double[] bounds;
	public static int[] res;
	public static RealMatrix cameraAxis;
	
	
	public Camera() {
		//System.out.println("Initializing Camera");
		eye = new ArrayRealVector(new double[3]);
		look = new ArrayRealVector(new double[3]);
		up = new ArrayRealVector(new double[3]);
		d = 0.0;
		bounds = new double[4];
		res = new int[2];
		cameraAxis = MatrixUtils.createRealMatrix(4, 4);
	}
	
	public void setEye(double x, double y, double z) {
		//System.out.println("eye : {" + x + ", " + y + ", " + z + "}");
		eye.setEntry(0, x);
		eye.setEntry(1, y);
		eye.setEntry(2, z);
	}
	
	public void setLook(double x, double y, double z) {
		look.setEntry(0, x);
		look.setEntry(1, y);
		look.setEntry(2, z);
	}
	
	public void setUp(double x, double y, double z) {
		up.setEntry(0, x);
		up.setEntry(1, y);
		up.setEntry(2, z);
	}
	
	public void setDistance(double distance) {
		d = distance;
	}
	
	public void setBounds(double left, double bottom, double right, double top) {
		bounds[0] = left;
		bounds[1] = bottom;
		bounds[2] = right;
		bounds[3] = top;
	}
	
	public void setRes(int horiz, int vert) {
		res[0] = horiz;
		res[1] = vert;
	}
	
	public void setCameraAxis(RealMatrix camAxis) {
		cameraAxis = camAxis;
	}
	
	public RealVector getEye() {
		return eye;
	}
	
	public RealVector getLook() {
		return look;
	}
	
	public RealVector getUp() {
		return up;
	}
	
	public double getDistance() {
		return d;
	}
	
	public double[] getBounds() {
		return bounds;
	}
	
	public int[] getRes() {
		return res;
	}
	
	public RealVector getGaze() {
		//double[] gaze = {eye[0] - look[0], eye[1] - look[1], eye[2] - look[2]};
		return eye.subtract(look);
	}
	
	public RealMatrix getCameraAxis() {
		return cameraAxis;
	}
	
	public int getWidth() {
		return res[0];
	}
	
	public int getHeight() {
		return res[1];
	}
	
	public double getLeft() {
		return bounds[0];
	}
	
	public double getBottom() {
		return bounds[1];
	}
	
	public double getRight() {
		return bounds[2];
	}
	
	public double getTop() {
		return bounds[3];
	}
	
	public RealVector getU() {
		return cameraAxis.getRowVector(0);
	}
	
	public RealVector getV() {
		return cameraAxis.getRowVector(1);
	}

	public RealVector getW() {
		return cameraAxis.getRowVector(2);
	}
	
}
