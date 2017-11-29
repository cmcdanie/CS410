//Colin McDaniel
//CS 410 Assignment #2
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
public class Raytracer {
	
	//Scene objects
	public LinkedList<ModelObject> modelObjects;
	public LinkedList<ModelSphere> modelSpheres;
	public Camera cam;
	
	//Output
	public int[][] pixelR;
	public int[][] pixelG;
	public int[][] pixelB;
	//public double[][] tValues;
	public double tMin;
	public double tMax;
	public double betaT;
	public double gammaT;
	
	//Lighting
	public RealVector ambientLight;
	public LinkedList<Light> lights;
	public RealVector[][] iValues;
	public double phongConstant = 19;
	public double[] background;
	public int level;
	
	public Raytracer(){
		modelObjects = new LinkedList<ModelObject>();
		modelSpheres = new LinkedList<ModelSphere>();
		lights = new LinkedList<Light>();	
		background = new double[] {0.0,0.0,0.0};
	}
	
	//Read the Driver file line by line
	public void readDriver(String filename) throws FileNotFoundException{
			Scanner scan = new Scanner(new File(filename));
			String line;
			
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
					//Model related
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
						break;
					
					
					case "sphere":
						ModelSphere sphere = new ModelSphere();
						sphere.setCenter(new ArrayRealVector(new double[] {Double.parseDouble(words[1]), Double.parseDouble(words[2]), Double.parseDouble(words[3])}));
						sphere.setRadius(Double.parseDouble(words[4]));
						sphere.setAmbient(new ArrayRealVector(new double[] {Double.parseDouble(words[5]), Double.parseDouble(words[6]), Double.parseDouble(words[7])}));
						sphere.setDiffuse(new ArrayRealVector(new double[] {Double.parseDouble(words[8]), Double.parseDouble(words[9]), Double.parseDouble(words[10])}));
						sphere.setSpecular(new ArrayRealVector(new double[] {Double.parseDouble(words[11]), Double.parseDouble(words[12]), Double.parseDouble(words[13])}));
						sphere.setAttenuation(new ArrayRealVector(new double[] {Double.parseDouble(words[14]), Double.parseDouble(words[15]), Double.parseDouble(words[16])}));
						modelSpheres.add(sphere);
						break;
					
					//Camera related
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
					
					//Lighting Related
					case "ambient":
						double[] Ba = new double[]{Double.parseDouble(words[1]), Double.parseDouble(words[2]), Double.parseDouble(words[3])};
						ambientLight = new ArrayRealVector(Ba);
						break;
					case "light":
						Light light = new Light();
						light.setLocation(Double.parseDouble(words[1]), Double.parseDouble(words[2]), Double.parseDouble(words[3]));
						light.setW(Double.parseDouble(words[4]));
						light.setEmission(Double.parseDouble(words[5]), Double.parseDouble(words[6]), Double.parseDouble(words[7]));
						lights.add(light);
						break;
						
					case "recursionLevel":
						level = Integer.parseInt(words[1]);
					
					//Comments
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
	
	//Rotates and Object using axis-angle
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
	
	//Scales and Object
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
	
	//Translates an Object
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

	//Returns Cross Product RealVector of two RealVectors
	public RealVector cross3(RealVector u, RealVector v){
		double[] c = { (u.getEntry(1) * v.getEntry(2)) - (u.getEntry(2) * v.getEntry(1)) 
						, (u.getEntry(2) * v.getEntry(0)) - (u.getEntry(0) * v.getEntry(2)) 
							, (u.getEntry(0) * v.getEntry(1)) - (u.getEntry(1) * v.getEntry(0))};
		
		RealVector cross = new ArrayRealVector(c);
		return cross;
	}

	//Saves Object to new file after transformations applied
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
	
	//
	public void render() {
		
		pixelR = new int[cam.getHeight()][cam.getWidth()];
		pixelG = new int[cam.getHeight()][cam.getWidth()];
		pixelB = new int[cam.getHeight()][cam.getWidth()];
		iValues = new RealVector[cam.getHeight()][cam.getWidth()];
		tMin = Double.POSITIVE_INFINITY;
		tMax = 0.0;
		
		double pastPercent = 0.0;
		
		//Iterate through each pixel
		for(int i = 0; i < cam.getHeight(); i++) {
			for(int j = 0; j < cam.getWidth(); j++) {
				
				//System.out.println("Pixel(" + i + ", " + j + ")");
				
				//Display percentage of completed pixels
				double percent = (double)i / (double)(cam.getHeight() - 1);
				percent = percent * 100;
				DecimalFormat df = new DecimalFormat("#.##");
				if(!df.format(percent).equals(df.format(pastPercent))){
					if((percent - pastPercent) >= 1 | percent == 100) {
						System.out.println("Ray Trace: " + df.format(percent) + "% Completed"); 
						pastPercent = percent;
					}
				}
				
				
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
				
				//Illumination Array
				RealVector Illum = new ArrayRealVector(new double[]{0.0, 0.0, 0.0});
				RealVector reffatt = new ArrayRealVector(new double[] {0.0, 0.0, 0.0});
				if(i == 278 & j == 233) {
					System.out.println("\tHere!");
				}
				iValues[i][j] = rayTrace(i, j, ray, pixelPt, Illum, reffatt, level);
				
			}
		}
	}

	//
	public RealVector rayTrace(int i, int j, RealVector ray, RealVector originPt, RealVector accum, RealVector reffatt, int level) {
		double t = Double.POSITIVE_INFINITY;
		
		String closestType = "";
		ModelObject closestObject = null;
		RealVector normalVect = null;
		int faceIndex = -1;
		RealVector intersectPt = null;
		ModelSphere closestSphere = null;
		int sphereIndex = -1;
		
		//Iterate through all objects
		for(int k = 0; k < modelObjects.size(); k++) {
			//System.out.println("Current Object: " + modelObjects.get(k).getName());
			
			//Iterate through faces
			for(int l = 0; l < modelObjects.get(k).getNumOfFaces(); l++) {
				
				int[] face = modelObjects.get(k).getFace(l);
				
				//System.out.println("Current Face: " + face[0] + ", " + face[1] + ", " + face[2]);
				
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
				double[] Y = 	{ptA[0] - originPt.getEntry(0), 
								 ptA[1] - originPt.getEntry(1), 
								 ptA[2] - originPt.getEntry(2)};
				
				RealVector constants = new ArrayRealVector(Y, false);
				RealVector solution = solver.solve(constants);
				double beta = solution.getEntry(0);
				double gamma = solution.getEntry(1);
				double t1 = solution.getEntry(2);
				
				if((beta < 0) && (beta > -0.00001)){
					beta = 0;
				}
				if((gamma < 0) && (gamma > -0.00001)){
					gamma = 0;
				}
				
				
				
				//Test if intersection is on face
				
				if((beta >= 0) && (gamma >= 0)){
					if((beta + gamma) <= 1){
						if((t1 > 0.00001) & (t1 <= t)){
							betaT = beta;
							gammaT = gamma;
							t = t1;
							closestType = "object";
							closestObject = modelObjects.get(k);
							
							//Calculate faceNormal
							//A - B
							double[] ab = {ptB[0] - ptA[0], ptB[1] - ptA[1], ptB[2] - ptA[2]};
							RealVector AB = new ArrayRealVector(ab);
							//A - C
							double[] ac = {ptC[0] - ptA[0], ptC[1] - ptA[1], ptC[2] - ptA[2]};
							RealVector AC = new ArrayRealVector(ac);
							
							normalVect = cross3(AB, AC);
							normalVect.unitize();
							
							//Calculate IntersectionPt
							RealVector ptAv = new ArrayRealVector(ptA);
							intersectPt = ptAv.add(AB.mapMultiply(beta).add(AC.mapMultiply(gamma)));
							faceIndex = l;
						}
					}
				}
			}
		}
		//System.out.println("\n");
		
		
		//Iterate through spheres
		for(int k = 0; k < modelSpheres.size(); k++){
			double radius = modelSpheres.get(k).getRadius();
			RealVector center = modelSpheres.get(k).getCenter();
			RealVector C = center.subtract(originPt);
			
			double V = C.dotProduct(ray);
			double csq = C.dotProduct(C);
			double disc = (radius * radius) - (csq - (V * V));
			if(disc >= 0){
				double d = Math.sqrt(disc);
				RealVector Q = ray.mapMultiply(V - d);
				Q = Q.add(originPt);
				double t1 = V - d;
				if((t1 < t) & (t1 > 0.00001)){
					t = t1;
					closestType = "sphere";
					closestSphere = modelSpheres.get(k);
					intersectPt = Q;
					sphereIndex = k;
					
					//Calculate sphereNormal
					normalVect = Q.subtract(center);
					normalVect.unitize();
					
				}
			}
		}
		if(i == 278 & j == 233) {
			System.out.println("current level: " + level);
			System.out.println("t: " + t);
			System.out.println("intersectPt: " + intersectPt);
		}
		
		//Calculate Illumination
		RealVector illumination = new ArrayRealVector(3);
		RealVector materialKr = new ArrayRealVector(3);
		
		//If no objects intersected, set illumination to background value
		if(t == Double.POSITIVE_INFINITY){
			illumination = new ArrayRealVector(background);
			accum = accum.add(illumination);
			//System.out.println("\tival: " + illumination);
		}
		else {
			//Lighting for Objects
			if(closestType.equals("object")) {
				illumination = objectIllumination(closestObject, faceIndex, normalVect, originPt, intersectPt);
				materialKr = closestObject.getFaceMaterial(faceIndex).getReflect();
			}
			
			//Lighting for Spheres
			else if(closestType.equals("sphere")) {
				illumination = sphereIllumination(closestSphere, sphereIndex, normalVect, originPt, intersectPt);
				materialKr = closestSphere.getAttenuation();
			}
			
			accum = accum.add(illumination);
			
			//Check recursion level
			if(level > 0) {
				//Calculate reflection ray from incoming ray
				RealVector rayInv = ray.mapMultiply(-1);
				double rayInvDotN = rayInv.dotProduct(normalVect);
				rayInvDotN = rayInvDotN * 2;
				RealVector Rv = normalVect.mapMultiply(rayInvDotN);
				Rv = Rv.subtract(rayInv);
				Rv.unitize();
				
				//Recursive call
				if(i == 278 & j == 233) {
					System.out.println("Rv: " + Rv);
					System.out.println("originPt: " + originPt);
					System.out.println("illumination: " + illumination);
				}
				//System.out.println("\t\taccum at level: " + level + ": " + accum);
				
				accum = rayTrace(i, j, Rv, intersectPt, accum, materialKr.ebeMultiply(reffatt), level - 1);
				
				if(i == 278 & j == 233) {
					System.out.println("accum: " + accum);
				}
			}
			
		}
		
		
		return accum;
	}
	
	//Returns true if no object or sphere in path from intersectPt to light
	public boolean pathToLightTest(RealVector ray, RealVector originPt, ModelObject currObj, int faceIndex) {
		
		double t = Double.POSITIVE_INFINITY;
		
		//Iterate through all objects
		for(int k = 0; k < modelObjects.size(); k++) {
			//System.out.println("Current Object: " + modelObjects.get(k).getName());
			
			//Iterate through faces
			for(int l = 0; l < modelObjects.get(k).getNumOfFaces(); l++) {
				
				int[] face = modelObjects.get(k).getFace(l);
				
				//System.out.println("Current Face: " + face[0] + ", " + face[1] + ", " + face[2]);
				
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
				double[] Y = 	{ptA[0] - originPt.getEntry(0), 
								 ptA[1] - originPt.getEntry(1), 
								 ptA[2] - originPt.getEntry(2)};
				
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
				
				
				
				//Test if intersection is on different face
				//if(modelObjects.get(k).getFace(l) != currObj.getFace(faceIndex)) {
					if((beta >= 0) && (gamma >= 0)){
						if((beta + gamma) <= 1){
							if((t1 > 0.00001) & (t1 <= t)){
								betaT = beta;
								gammaT = gamma;
								t = t1;
								
							}
						}
					}
				//}
				
				
				if(t < Double.POSITIVE_INFINITY) {
					//System.out.println("Path to light returned false: Object");
					return false;
				}
			}
		}
		//System.out.println("\n");
		
		
		//Iterate through spheres
		for(int k = 0; k < modelSpheres.size(); k++){
			double radius = modelSpheres.get(k).getRadius();
			RealVector center = modelSpheres.get(k).getCenter();
			RealVector C = center.subtract(originPt);
			
			double V = C.dotProduct(ray);
			double csq = C.dotProduct(C);
			double disc = (radius * radius) - (csq - (V * V));
			if(disc >= 0){
				double d = Math.sqrt(disc);
				RealVector Q = ray.mapMultiply(V - d);
				Q = Q.add(originPt);
				double t1 = V - d;
				if((t1 < t) & (t1 > 0.00001)){
					t = t1;
					
					
				}
			}
			
			if(t < Double.POSITIVE_INFINITY) {
				//System.out.println("Path to light returned false: Sphere");
				return false;
			}
			
		}
		return true;
	}
	
	public boolean pathToLightTest(RealVector ray, RealVector originPt, ModelSphere currSphere) {
		
		double t = Double.POSITIVE_INFINITY;
		
		//Iterate through all objects
		for(int k = 0; k < modelObjects.size(); k++) {
			//System.out.println("Current Object: " + modelObjects.get(k).getName());
			
			//Iterate through faces
			for(int l = 0; l < modelObjects.get(k).getNumOfFaces(); l++) {
				
				int[] face = modelObjects.get(k).getFace(l);
				
				//System.out.println("Current Face: " + face[0] + ", " + face[1] + ", " + face[2]);
				
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
				double[] Y = 	{ptA[0] - originPt.getEntry(0), 
								 ptA[1] - originPt.getEntry(1), 
								 ptA[2] - originPt.getEntry(2)};
				
				RealVector constants = new ArrayRealVector(Y, false);
				RealVector solution = solver.solve(constants);
				double beta = solution.getEntry(0);
				double gamma = solution.getEntry(1);
				double t1 = solution.getEntry(2);
				
				if((beta < 0) && (beta > -0.00001)){
					beta = 0;
				}
				if((gamma < 0) && (gamma > -0.00001)){
					gamma = 0;
				}
				
				
				
				//Test if intersection is on face
				if((beta >= 0) && (gamma >= 0)){
					if((beta + gamma) <= 1){
						if((t1 > 0.00001) & (t1 <= t)){
							betaT = beta;
							gammaT = gamma;
							t = t1;
							
						}
					}
				}
				
				
				if(t < Double.POSITIVE_INFINITY) {
					//System.out.println("Path to light returned false: Object");
					return false;
				}
			}
		}
		//System.out.println("\n");
		
		
		//Iterate through spheres
		for(int k = 0; k < modelSpheres.size(); k++){
			if(modelSpheres.get(k) != currSphere) {
				double radius = modelSpheres.get(k).getRadius();
				RealVector center = modelSpheres.get(k).getCenter();
				RealVector C = center.subtract(originPt);
				
				double V = C.dotProduct(ray);
				double csq = C.dotProduct(C);
				double disc = (radius * radius) - (csq - (V * V));
				if(disc >= 0){
					double d = Math.sqrt(disc);
					RealVector Q = ray.mapMultiply(V - d);
					Q = Q.add(originPt);
					double t1 = V - d;
					if((t1 < t) & (t1 > 0.00001)){
						t = t1;
						
						
					}
				}
			}
			
			if(t < Double.POSITIVE_INFINITY) {
				//System.out.println("Path to light returned false: Sphere");
				return false;
			}
			
		}
		return true;
	}
	
	//Calculates Illumination value for an object
	public RealVector objectIllumination(ModelObject closestObject, int faceIndex, RealVector faceNormal, RealVector originPt, RealVector intersectPt) {
		//Ambient
		//Ba = ambient light from scene; Ka = ambient light from material
		RealVector Ba = ambientLight;
		RealVector Ka = closestObject.getFaceMaterial(faceIndex).getAmbient();
		RealVector Ia = Ba.ebeMultiply(Ka);
		
		//Diffuse
		RealVector Id = new ArrayRealVector(3);
		for(int m = 0; m < lights.size(); m++) {
			
			//Create LightVector 
			//Light location - surface point location
			RealVector L = new ArrayRealVector(lights.get(m).getLocation());
			L = L.subtract(intersectPt);
			L.unitize();
			
			if(pathToLightTest(L, intersectPt, closestObject, faceIndex)) {
				//B = brightness from light source; Kd = diffuse value from material
				RealVector B = new ArrayRealVector(lights.get(m).getEmission());
				RealVector Kd = closestObject.getFaceMaterial(faceIndex).getDiffuse();
				
				//I += Bd * Kd * (N dot L)
				//System.out.println("FaceNormal: " + faceNormal);
				double NdotL = faceNormal.dotProduct(L);
				
				if(NdotL < 0) {
					NdotL = -1 * NdotL;
				}
				
				//System.out.println("\tNdotLv: " + NdotLv);
				Id = Id.add(B.ebeMultiply(Kd).mapMultiply(NdotL));
			}
			
		}
		
		//Specular
		RealVector Is = new ArrayRealVector(3);
		for(int m = 0; m < lights.size(); m++) {
			
			//Create LightVector 
			//Light location - surface point location
			RealVector L = new ArrayRealVector(lights.get(m).getLocation());
			L = L.subtract(intersectPt);
			L.unitize();
			
			if(pathToLightTest(L, intersectPt, closestObject, faceIndex)) {
				//B = brightness from light source; Kd = diffuse value from material
				RealVector B = new ArrayRealVector(lights.get(m).getEmission());
				RealVector Ks = closestObject.getFaceMaterial(faceIndex).getSpecular();
				
				//Create Reflection Vector
				double NdotL = faceNormal.dotProduct(L);
				if(NdotL < 0) {
					faceNormal.mapMultiplyToSelf(-1);
					NdotL = faceNormal.dotProduct(L);
					//NdotL = -1 * NdotL;
				}
				
				RealVector R = faceNormal.mapMultiply(2 * NdotL).subtract(L);
				
				
				//Create View Vector
				RealVector V = originPt.subtract(intersectPt);
				V.unitize();
				double VdotR = V.dotProduct(R);
				
				//System.out.println("\tNdotL: " + NdotL);
				if(VdotR < 0) {
					VdotR = 0;
				}
				Is = Is.add(B.ebeMultiply(Ks).mapMultiply(Math.pow(VdotR, closestObject.getFaceMaterial(faceIndex).getNs())));
			}
			
				
		}
		return Ia.add(Id.add(Is));
		
	}
	
	//Calculates Illumination value for spheres
	public RealVector sphereIllumination(ModelSphere closestSphere, int sphereIndex, RealVector sphereNormal, RealVector originPt, RealVector intersectPt) {
		//Ambient
		//Ba = ambient light from scene; Ka = ambient light from material
		RealVector Ba = ambientLight;
		RealVector Ka = new ArrayRealVector(closestSphere.getAmbient());
		RealVector Ia = Ka.ebeMultiply(Ba);
		
		
		//Diffuse
		RealVector Id = new ArrayRealVector(3);
		for(int m = 0; m < lights.size(); m++) {
			
			//Create LightVector 
			//Light location - surface point location
			RealVector L = new ArrayRealVector(lights.get(m).getLocation());
			L = L.subtract(intersectPt);
			L.unitize();
			
			if(pathToLightTest(L, intersectPt, closestSphere)) {
				//B = brightness from light source; Kd = diffuse value from material
				RealVector B = new ArrayRealVector(lights.get(m).getEmission());
				//System.out.println("B: " + B);
				RealVector Kd = new ArrayRealVector(closestSphere.getDiffuse());
				//System.out.println("Kd: " + Kd);
				
				//I += Bd * Kd * (N dot L)
				//System.out.println("\tSN: " + sphereNormal);
				//System.out.println("\tL: " + L);
				double NdotL = sphereNormal.dotProduct(L);
				
				if(NdotL < 0) {
					NdotL = 0;
				}
				
				//System.out.println("\tNdotL: " + NdotL);
				Id = Id.add(B.ebeMultiply(Kd).mapMultiply(NdotL));
			}
			
		}
		
		
		//Specular
		RealVector Is = new ArrayRealVector(3);
		for(int m = 0; m < lights.size(); m++) {
			
			//Create LightVector 
			//Light location - surface point location
			RealVector L = new ArrayRealVector(lights.get(m).getLocation());
			L = L.subtract(intersectPt);
			L.unitize();
			
			if(pathToLightTest(L, intersectPt, closestSphere)) {
				//B = brightness from light source; Kd = diffuse value from material
				RealVector B = new ArrayRealVector(lights.get(m).getEmission());
				RealVector Ks = new ArrayRealVector(closestSphere.getSpecular());
				
				//I += Bd * Kd * (N dot L)
				//System.out.println("\tSN: " + sphereNormal);
				//System.out.println("\tL: " + L);
				
				//Create Reflection Vector
				double NdotL = sphereNormal.dotProduct(L);
				if(NdotL < 0) {
					NdotL = 0;
				}
				RealVector R = sphereNormal.mapMultiply(2 * NdotL).subtract(L);
				
				//Create View Vector
				RealVector V = originPt.subtract(intersectPt);
				V.unitize();
				double VdotR = V.dotProduct(R);
				
				//System.out.println("\tNdotL: " + NdotL);
				if(VdotR < 0) {
					VdotR = 0;
				}
				Is = Is.add(B.ebeMultiply(Ks).mapMultiply(Math.pow(VdotR, phongConstant)));
			}
			
		}
		return Ia.add(Id.add(Is));
	}
	
	//Converts illumination values to RGB (0-255)
	public void illumToRGB() {
		//System.out.println("tMax: " + tMax);
		//System.out.println("tMin: " + tMin);
		
		for(int i = 0; i < cam.getHeight(); i++) {
			for(int j = 0; j < cam.getWidth(); j++) {
				int Red = (int) (iValues[i][j].getEntry(0) * 255);
				if(Red > 255) {
					Red = 255;
				}
				int Green = (int) (iValues[i][j].getEntry(1) * 255);
				if(Green > 255) {
					Green = 255;
				}
				int Blue = (int) (iValues[i][j].getEntry(2) * 255);
				if(Blue > 255) {
					Blue = 255;
				}
				pixelR[i][j] = Red;
				pixelG[i][j] = Green;
				pixelB[i][j] = Blue;
			}
		}
		
		
	}
	
	//Saves RGB values as a ppm image file
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
	
	//Prints a Matrix to the console
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
		Raytracer a = new Raytracer();
		
		/*
		double[] a1 = {0.1, 0.2, 0.3};
		RealVector Ba = new ArrayRealVector(a1);
		double[] a2 = {0.4, 0.5, 0.6};
		RealVector Ka = new ArrayRealVector(a2);
		RealVector ill = Ba.ebeMultiply(Ka);
		System.out.println(ill);
		*/
		
		//a.readDriver(args[0]);
		
		//a.readDriver("driver00.txt");
		//a.readDriver("driver01.txt");
		//a.readDriver("driver02.txt");
		
		//a.readDriver("drivers_models\driver_models\driver00.txt");
		a.readDriver("simface.txt");
		
		System.out.println("Completed Reading Driver file");
		a.render();
		System.out.println("Completed Ray Tracing");
		a.illumToRGB();
		System.out.println("Completed Storing RGB values");
		a.exportPPM(args[1]);
		System.out.println("Exported Image");
		System.out.println("Done!");
		
	}

}
