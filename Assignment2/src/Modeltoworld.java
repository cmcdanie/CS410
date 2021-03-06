//Colin McDaniel
//CS 410 Assignment #1
//CSU ID: 830293766



import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Scanner;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
public class Modeltoworld {
	
	public LinkedList<ModelObject> modelObjects;
	public LinkedList<ModelSphere> modelSpheres;
	public Camera cam;
	public int[][] pixelR;
	public int[][] pixelG;
	public int[][] pixelB;
	public double[][] tValues;
	public double tMin;
	public double tMax;
	public double betaT;
	public double gammaT;
	
	public Modeltoworld(){
		modelObjects = new LinkedList<ModelObject>();
		modelSpheres = new LinkedList<ModelSphere>();
	}
	
	//Read the Driver file line by line
	public void readDriver(String filename) throws FileNotFoundException{
			Scanner scan = new Scanner(new File(filename));
			String line;
			//modelObjects = new LinkedList<ModelObject>();
			//modelSpheres = new LinkedList<ModelSphere>();
			
			//Driver variables
			String type;
			double wx, wy, wz, theta, sc, tx, ty, tz;
			String model;
			
			cam = new Camera();
			while(scan.hasNextLine()){
								
				//Read transformations from Driver File
				line = scan.nextLine();
				String[] words = line.split("\\s");
				
				/*
				for(int i = 0; i < words.length; i++){
					System.out.println(words[i]);
				}
				*/
				
				type = words[0];
				
				switch(type) {
					case "model":
						wx = Double.parseDouble(words[1]);
						wy = Double.parseDouble(words[2]);
						wz = Double.parseDouble(words[3]);
						theta = Double.parseDouble(words[4]);
						sc = Double.parseDouble(words[5]);
						tx = Double.parseDouble(words[6]);
						ty = Double.parseDouble(words[7]);
						tz = Double.parseDouble(words[8]);
						model = words[9];
						
						//Apply Transformations to Model Object
						ModelObject obj = new ModelObject(model);
						rotateObj(wx,wy,wz,theta,obj);
						scaleObj(sc, obj);
						translateObj(tx, ty, tz, obj);
						
						modelObjects.add(obj);
						
						
						//To-Do remove this section
						
						//Save obj
						//exportObj(obj, model, filename);
						break;
						
					case "eye":
						cam.setEye(Double.parseDouble(words[1]), Double.parseDouble(words[2]), Double.parseDouble(words[3]));
						break;
					case "look":
						cam.setLook(Double.parseDouble(words[1]), Double.parseDouble(words[2]), Double.parseDouble(words[3]));
						break;
					case "up":
						cam.setUp(Double.parseDouble(words[1]), Double.parseDouble(words[2]), Double.parseDouble(words[3]));
						break;
					case "d":
						cam.setDistance(Double.parseDouble(words[1]));
						break;
					case "bounds":
						cam.setBounds(Double.parseDouble(words[1]), Double.parseDouble(words[2]), Double.parseDouble(words[3]), Double.parseDouble(words[4]));
						break;
					case "res":
						cam.setRes(Integer.parseInt(words[1]), Integer.parseInt(words[2]));
						break;
					case "sphere":
						ModelSphere sphere = new ModelSphere();
						sphere.setCenter(Double.parseDouble(words[1]), Double.parseDouble(words[2]), Double.parseDouble(words[3]));
						sphere.setRadius(Double.parseDouble(words[4]));
						modelSpheres.add(sphere);
						break;
					case "#":
						break;
					
				}
				
				
			}
			scan.close();
			setupCamera();
		}
	
	//Generate Camera Axis and other parameters
	public void setupCamera() {
		
		RealVector W = cam.getGaze();
		W.unitize();
		
		RealVector up = cam.getUp();
		RealVector U = cross3(up, W);
		U.unitize();
		
		RealVector V = cross3(W, U);
		
		W = W.append(0);
		U = U.append(0);
		V = V.append(0);
		
		//System.out.println("W: " + W);
		//System.out.println("U: " + U);
		//System.out.println("V: " + V);
		
		double[][] camAxisArray = {U.toArray(), V.toArray(), W.toArray(), {0,0,0,1}};
		RealMatrix camAxis = MatrixUtils.createRealMatrix(camAxisArray);
		
		cam.setCameraAxis(camAxis);
		
		//System.out.println("Camera Axis:");
		//printMatrix(cam.getCameraAxis());
		
		
		
	}
	
	public void rayTrace() {
		
		pixelR = new int[cam.getHeight()][cam.getWidth()];
		pixelG = new int[cam.getHeight()][cam.getWidth()];
		pixelB = new int[cam.getHeight()][cam.getWidth()];
		tValues = new double[cam.getHeight()][cam.getWidth()];
		tMin = Double.POSITIVE_INFINITY;
		tMax = 0.0;
		
		
		double pastPercent = 0.0;
		//Iterate through each pixel
		for(int i = 0; i < cam.getHeight(); i++) {
			for(int j = 0; j < cam.getWidth(); j++) {
				
				
				//System.out.println("Pixel(" + i + ", " + j + ")");
				double percent = (double)i / (double)(cam.getWidth() - 1);
				percent = percent * 100;
				DecimalFormat df = new DecimalFormat("#.##");
				if(!df.format(percent).equals(df.format(pastPercent))){
					System.out.println("Ray Trace: " + df.format(percent) + "% Completed"); 
				}
				
				pastPercent = percent;
				//System.out.println("Pixel(" + i + ", " + j + ")");
				//System.out.println("Right: " + cam.getRight());
				//System.out.println("Left: " + cam.getLeft());
				
				double pixelX = (double)j / (double)(cam.getWidth() - 1) * (cam.getRight() - cam.getLeft()) + cam.getLeft();
				double pixelY = (double)i / (double)(cam.getHeight() - 1) * (cam.getBottom() - cam.getTop()) + cam.getTop();
				
				
				//System.out.println("pixelX: " + pixelX);
				//System.out.println("pixelY: " + pixelY);
				
				//Code used for Class Slides
				
				//Generate pixel
				RealVector pixelPt = cam.getEye();										//Eye
				pixelPt = pixelPt.append(0);
				pixelPt = pixelPt.add(cam.getW().mapMultiply(-1 * cam.getDistance()));	//Eye + (W * -d)
				pixelPt = pixelPt.add(cam.getU().mapMultiply(pixelX));					//Eye + (W * -d) + (U * pixelX)
				pixelPt = pixelPt.add(cam.getV().mapMultiply(pixelY));					//Eye + (W * -d) + (U * pixelX) + (V * pixelY)
				pixelPt = pixelPt.getSubVector(0, 3);
				//System.out.println("Pixel Pt: " + pixelPt);
				
				//Generate Ray
				RealVector ray = pixelPt.subtract(cam.getEye());
				ray.unitize();
				//System.out.println("Ray: " + ray);
				
				double t = Double.POSITIVE_INFINITY;
				
				//Iterate through all objects
				for(int k = 0; k < modelObjects.size(); k++) {
					//System.out.println("Current Object: " + modelObjects.get(k).getName());
					//Iterate through faces
					for(int l = 0; l < modelObjects.get(k).getNumOfFaces(); l++) {
						
						if(i == 102 && j == 102){
							System.out.println("Current Face: " + (l + 1));
						}
						
						
						int[] face = modelObjects.get(k).getFace(l);
						//System.out.println(face[0]);
						double[] ptA = modelObjects.get(k).getVertex(face[0] - 1);
						double[] ptB = modelObjects.get(k).getVertex(face[1] - 1);
						double[] ptC = modelObjects.get(k).getVertex(face[2] - 1);
						
						//System.out.println("ptA: {" + ptA[0] + "," + ptA[1] + "," + ptA[2] + "}");
						//System.out.println("ptB: {" + ptB[0] + "," + ptB[1] + "," + ptB[2] + "}");
						//System.out.println("ptC: {" + ptC[0] + "," + ptC[1] + "," + ptC[2] + "}");
						
						//System.out.println("ray:dx: " + ray.getEntry(0));
						double[][] M = {{ptA[0] - ptB[0], ptA[0] - ptC[0], ray.getEntry(0)},
										{ptA[1] - ptB[1], ptA[1] - ptC[1], ray.getEntry(1)},
										{ptA[2] - ptB[2], ptA[2] - ptC[2], ray.getEntry(2)}};
						RealMatrix coefficients = new Array2DRowRealMatrix(M, false);
						DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
						
						//System.out.println("eyex: " + cam.getEye().getEntry(0));
						double[] Y = 	{ptA[0] - pixelPt.getEntry(0), 
										 ptA[1] - pixelPt.getEntry(1), 
										 ptA[2] - pixelPt.getEntry(2)};
						
						RealVector constants = new ArrayRealVector(Y, false);
						RealVector solution = solver.solve(constants);
						double beta = solution.getEntry(0);
						double gamma = solution.getEntry(1);
						double t1 = solution.getEntry(2);
						if((beta < 0) && (beta > -0.00000000001)){
							beta = 0;
						}
						if((gamma < 0) && (gamma > -0.00000000001)){
							gamma = 0;
						}
						if(i == 102 && j == 102){
							System.out.println("Beta = " + beta);
							System.out.println("Gamma = " + gamma);
							System.out.println("t = " + t);
							System.out.println("t1 = " + t1);
						}
						
						
						
						//Test if intersection is on face
						if((beta >= 0) && (gamma >= 0)){
							
							if(i == 102 && j == 102){
								System.out.println("yes1");
								System.out.println("beta + gamma: " + (beta + gamma));
							}
							
							if((beta + gamma) <= 1){
								
								if(i == 102 && j == 102){
									System.out.println("yes2");
								}
								
								if(t1 > 0){
									if(i == 102 && j == 102){
										System.out.println("yes3");
									}
									//System.out.println("*****************************************");
									//System.out.println("Intersection!!!");
									//System.out.println("*****************************************");
									
									//System.out.println("point of intersection:");
									//RealVector intersection = ray.mapMultiply(t1);
									//System.out.println("{" + intersection.getEntry(0) + "," + intersection.getEntry(1) + "," + intersection.getEntry(2) + "}");
									
									if(t1 < t){
										betaT = beta;
										gammaT = gamma;
										//System.out.println("Current Face: " + (l + 1));
										t = t1;
									}
								}
							}
						}
						
					}
				}
				//System.out.println("\n");
				
				
				//Iterate through spheres
				for(int k = 0; k < modelSpheres.size(); k++){
					double[] cntr = modelSpheres.get(k).getCenter();
					double radius = modelSpheres.get(k).getRadius();
					RealVector center = new ArrayRealVector(cntr);
					RealVector C = center.subtract(pixelPt);
					
					double V = C.dotProduct(ray);
					double csq = C.dotProduct(C);
					double disc = (radius * radius) - (csq - (V * V));
					if(disc >= 0){
						double d = Math.sqrt(disc);
						RealVector Q = ray.mapMultiply(V - d);
						Q = Q.add(pixelPt);
						double t1 = V - d;
						if(t1 < t){
							t = t1;
						}
					}
				}
				//System.out.println("t: " + t);
				//System.out.println("tMaxPre: " + tMax);
				if((t > tMax) && !(t == Double.POSITIVE_INFINITY)){
					tMax = t;
				}
				
				if(t < tMin){
					tMin = t;
				}
				//System.out.println("Pixel(" + i + ", " + j + ")");
				//System.out.println("tMax: " + tMax);
				//System.out.println("beta: " + betaT);
				//System.out.println("gamma: " + gammaT);
	
				
				//System.out.println("tMin: " + tMin);
				tValues[i][j] = t;
				
				if(t == Double.POSITIVE_INFINITY){
					pixelR[i][j] = 239;
					pixelG[i][j] = 239;
					pixelB[i][j] = 239;
				}
			}
		}
	}
	
	public void rotateObj(double wx, double wy, double wz, double theta, ModelObject obj){
		
		//Normalize rotation (W) axis
		double[] w = {wx,wy,wz};
		RealVector W = new ArrayRealVector(w);
		W.unitize();
				//System.out.println("W: " + W);
		
		//Create axis M
		int minIndex = W.getMinIndex();
		RealVector M = W.copy();
		M.setEntry(minIndex, 1);
		M.unitize();
				//System.out.println("M: " + M);
				//System.out.println("minIndex = " + minIndex);
		
		//Create U axis (U = W x M)
		RealVector U = cross3(W,M);
		if(U.getL1Norm() == 0.0) {
			switch (minIndex) {
				case 0:	M.setEntry(1, 1);
						break;
				case 1:	M.setEntry(2, 1);
						break;
				case 2:	M.setEntry(0, 1);
						break;
			}
			M.unitize();
			U = cross3(W,M);
		}
		U.unitize();
				//System.out.println("U: " + U);
		
		//Create -V axis (-V = W x U) or (V = U x W)
		RealVector V = cross3(W,U);
				//System.out.println("V: " + V);
		
		//Create Rotation Matrix
		U = U.append(0);
		V = V.append(0);
		W = W.append(0);
				//System.out.println("U: " + U);
				//System.out.println("V: " + V);
				//System.out.println("W: " + W);
		double[][] rot = {U.toArray(), V.toArray(), W.toArray(), {0,0,0,1}};
		RealMatrix rotMatrix = MatrixUtils.createRealMatrix(rot);
		
		
		//Create Z axis rotation matrix
		double[][] z = {
						{Math.cos(Math.toRadians(theta)),	-1 * Math.sin(Math.toRadians(theta)),	0,	0},
						{Math.sin(Math.toRadians(theta)),	Math.cos(Math.toRadians(theta)),		0,	0},
						{0,									0,										1,	0},
						{0,									0,										0,	1}
										};
		RealMatrix zRot = MatrixUtils.createRealMatrix(z);
				
				//System.out.println("rotMatrix: ");
				//printMatrix(rotMatrix);
				//System.out.println("zRot: ");
				//printMatrix(zRot);
		
		//Compute rotation
		for(int i = 0; i < obj.numOfVert(); i++){
			RealMatrix vertex = MatrixUtils.createRealMatrix(obj.getVertex(i).length + 1, 1);
			double[] adjust = Arrays.copyOf(obj.getVertex(i), obj.getVertex(i).length + 1);
			adjust[adjust.length - 1] = 1;
			vertex.setColumn(0, adjust);
			
					//System.out.println("vertex: ");
					//printMatrix(vertex);
			
			//Rw * P
			RealMatrix temp = rotMatrix.multiply(vertex);
					//System.out.println("Rw: ");
					//printMatrix(temp);
			
			//Rz * P
			temp = zRot.multiply(temp);
					//System.out.println("Rz: ");
					//printMatrix(temp);
			
			//RwT * P
			
			RealMatrix newVertex = rotMatrix.transpose().multiply(temp);
					//System.out.println("newVertex: ");
					//printMatrix(newVertex);
			
			double[] adjust2 = newVertex.getColumn(0);
			adjust2 = Arrays.copyOf(adjust2, adjust2.length - 1);
			
			
			obj.setVertex(i, adjust2);
			
			
		}
		
	}
	
	public void scaleObj(double sc, ModelObject obj){
		//Scale Matrix
		double[][] s = {{sc,0,0,0}, {0,sc,0,0}, {0,0,sc,0}, {0,0,0,1}};
		RealMatrix scaleMatrix = MatrixUtils.createRealMatrix(s);
		
		for(int i = 0; i < obj.numOfVert(); i++){
			RealMatrix vertex = MatrixUtils.createRealMatrix(obj.getVertex(i).length + 1, 1);
			double[] adjust = Arrays.copyOf(obj.getVertex(i), obj.getVertex(i).length + 1);
			adjust[adjust.length - 1] = 1;
			vertex.setColumn(0, adjust);
			
					//stem.out.println("vertex: ");
					//printMatrix(vertex);
					//System.out.println("Scale: " + sc);
			RealMatrix newVertex = scaleMatrix.multiply(vertex);
			
			double[] adjust2 = newVertex.getColumn(0);
			adjust2 = Arrays.copyOf(adjust2, adjust2.length - 1);
			
			
			obj.setVertex(i, adjust2);
					//System.out.println("newVertex: ");
					//printMatrix(newVertex);
		}
			
	}
	
	public void translateObj(double tx, double ty, double tz, ModelObject obj){
		//Translation Matrix
		double[][] t = {{1,0,0,tx}, {0,1,0,ty}, {0,0,1,tz}, {0,0,0,1}};
		RealMatrix tranMatrix = MatrixUtils.createRealMatrix(t);
		
		for(int i = 0; i < obj.numOfVert(); i++){
			RealMatrix vertex = MatrixUtils.createRealMatrix(obj.getVertex(i).length + 1, 1);
			double[] adjust = Arrays.copyOf(obj.getVertex(i), obj.getVertex(i).length + 1);
			adjust[adjust.length - 1] = 1;
			vertex.setColumn(0, adjust);
			
			
			RealMatrix newVertex = tranMatrix.multiply(vertex);
			
			double[] adjust2 = newVertex.getColumn(0);
			adjust2 = Arrays.copyOf(adjust2, adjust2.length - 1);
			
			
			obj.setVertex(i, adjust2);
			//System.out.println("newVertex: ");
			//printMatrix(newVertex);
		}
	}

	public RealVector cross3(RealVector u, RealVector v){
		double[] c = { (u.getEntry(1) * v.getEntry(2)) - (u.getEntry(2) * v.getEntry(1)) 
						, (u.getEntry(2) * v.getEntry(0)) - (u.getEntry(0) * v.getEntry(2)) 
							, (u.getEntry(0) * v.getEntry(1)) - (u.getEntry(1) * v.getEntry(0))};
		
		RealVector cross = new ArrayRealVector(c);
		return cross;
	}

	public void exportObj(ModelObject obj, String model, String driverName){
		//Create driver folder
		//System.out.println(driverName);
		String[] a = driverName.split("\\.");
		new File(a[0]).mkdir();
		//System.out.println(a[0]);
		
		//file name
		int suffix = 0;
		String[] modelFile = model.split("\\.");
		//System.out.println(a[0] + "\\" + modelFile[0] + "_mw" + String.format("%02d", suffix) + ".obj");
		File testFile = new File(a[0] + "/" + modelFile[0] + "_mw" + String.format("%02d", suffix) + ".obj");
		boolean exists = testFile.exists();
		while(exists){
			suffix++;
			testFile = new File(a[0] + "/" + modelFile[0] + "_mw" + String.format("%02d", suffix) + ".obj");
			exists = testFile.exists();
		}
		
		try {
			FileWriter fw = new FileWriter(a[0] + "/" + modelFile[0] + "_mw" + String.format("%02d", suffix) + ".obj");
			int j = 0;
			//loop for each line of text
			for(int i = 0; i < obj.getObjTxt().size(); i++){
				String line = obj.getObjTxt(i);
				String[] words = line.split("\\s");
				//if line is vertex, use newer data
				if(words[0].equals("v")){
					double[] temp = obj.getVertex(j);
					j++;
					String points = "";
					for(int k = 0; k < temp.length; k++){
						points += String.format("%.7f", temp[k]) + " ";
					}
					fw.write("v " + points + "\n");
				}
				else if(words[0].equals("vn")){
					
				}
				else{
					fw.write(obj.getObjTxt(i) + "\n");
				}
			}
			
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void render(){
		
		//System.out.println("tMax: " + tMax);
		//System.out.println("tMin: " + tMin);
		
		for(int i = 0; i < cam.getHeight(); i++) {
			for(int j = 0; j < cam.getWidth(); j++) {
				if((pixelR[i][j] == 239) && (pixelG[i][j] == 239) && (pixelB[i][j] == 239)){
					
				}
				else{
					double t = tValues[i][j];
					double ratio = 2 * (t - tMin) / (tMax - tMin);
					pixelR[i][j] = (int) Math.max(0, 255 * (1 - ratio));
					pixelB[i][j] = (int) Math.max(0, 255 * (ratio - 1));
					pixelG[i][j] = 255 - pixelB[i][j] - pixelR[i][j];
				}
			}
		}
	}
	
	public void exportPPM(String filename){
		try {
			FileWriter fw = new FileWriter(new File(filename));
			fw.write("P3\n");
			fw.write(cam.getWidth() + " " + cam.getHeight() + " 255\n");
			for(int i = 0; i < cam.getHeight(); i++) {
				for(int j = 0; j < cam.getWidth(); j++) {
					fw.write(pixelR[i][j] + " " + pixelG[i][j] + " " + pixelB[i][j] + "\n");
				}
			}
			
			
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public void printMatrix(RealMatrix a){
		DecimalFormat df = new DecimalFormat("#.####");
		for(int i = 0; i < a.getRowDimension(); i++){
			for(int j = 0; j < a.getColumnDimension(); j++){
				System.out.print(df.format(a.getEntry(i, j)) + "\t");
			}
			System.out.println(" ");
		}
				
	}
	
	//Main
	public static void main(String[] args) throws FileNotFoundException {
		Modeltoworld a = new Modeltoworld();
		a.readDriver(args[0]);
		System.out.println("Completed Reading Driver file");
		a.rayTrace();
		System.out.println("Completed Ray Tracing");
		a.render();
		System.out.println("Completed Storing RGB values");
		a.exportPPM(args[1]);
		System.out.println("Exported Image");
		System.out.println("Done!");
		
	}

}
