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
import java.util.Scanner;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.ArrayRealVector;
public class Modeltoworld {
	
	
	
	public static void rotateObj(double wx, double wy, double wz, double theta, ModelObject obj){
		
		//Normalize rotation (W) axis
		double[] w = {wx,wy,wz};
		RealVector W = new ArrayRealVector(w);
		W.unitize();
		System.out.println("W: " + W);
		
		//Create axis M
		int minIndex = W.getMinIndex();
		System.out.println(minIndex);
		RealVector M = W.copy();
		M.setEntry(minIndex, 1);
		M.unitize();
		System.out.println("M: " + M);
		
		//Create U axis (U = W x M)
		RealVector U = cross3(W,M);
		U.unitize();
		System.out.println("U: " + U);
		
		//Create -V axis (-V = W x U) or (V = U x W)
		RealVector V = cross3(W,U);
		System.out.println("V: " + V);
		
		//Create Rotation Matrix
		U = U.append(0);
		V = V.append(0);
		W = W.append(0);
		System.out.println("U: " + U);
		System.out.println("V: " + V);
		System.out.println("W: " + W);
		double[][] rot = {U.toArray(), V.toArray(), W.toArray(), {0,0,0,1}};
		RealMatrix rotMatrix = MatrixUtils.createRealMatrix(rot);
		
		
		//Create Z axis rotation matrix
		double[][] z = {{Math.cos(Math.toRadians(theta)), -1 * Math.sin(Math.toRadians(theta)), 0, 0},
							{Math.sin(Math.toRadians(theta)), Math.cos(Math.toRadians(theta)), 0, 0},
								{0, 0, 1, 0},
									{0, 0, 0, 1}};
		RealMatrix zRot = MatrixUtils.createRealMatrix(z);
		
		//Compute rotation
		for(int i = 0; i < obj.numOfVert(); i++){
			RealMatrix vertex = MatrixUtils.createRealMatrix(obj.getVertex(i).length + 1, 1);
			double[] adjust = Arrays.copyOf(obj.getVertex(i), obj.getVertex(i).length + 1);
			adjust[adjust.length - 1] = 1;
			vertex.setColumn(0, adjust);
			
			
			
			
			//System.out.println("rotMatrix: ");
			//printMatrix(rotMatrix);
			//System.out.println("vertex: ");
			//printMatrix(vertex);
			
			//Rw * P
			RealMatrix temp = rotMatrix.multiply(vertex);
			//System.out.println("temp: ");
			//printMatrix(temp);
			
			//Rz * P
			temp = zRot.multiply(temp);
			
			//RwT * P
			RealMatrix newVertex = rotMatrix.transpose().multiply(temp);
			//System.out.println("newVertex: ");
			//printMatrix(newVertex);
			
			double[] adjust2 = newVertex.getColumn(0);
			adjust2 = Arrays.copyOf(adjust2, adjust2.length - 1);
			
			
			obj.setVertex(i, adjust2);
			
			
		}
		
	}
	
	public static void scaleObj(double sc, ModelObject obj){
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
	
	public static void translateObj(double tx, double ty, double tz, ModelObject obj){
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
	
	
	
	public static RealVector cross3(RealVector u, RealVector v){
		double[] c = { (u.getEntry(1) * v.getEntry(2)) - (u.getEntry(2) * v.getEntry(1)) 
						, (u.getEntry(2) * v.getEntry(0)) - (u.getEntry(0) * v.getEntry(2)) 
							, (u.getEntry(0) * v.getEntry(1)) - (u.getEntry(1) * v.getEntry(0))};
		
		RealVector cross = new ArrayRealVector(c);
		return cross;
	}
	
	//Read the Driver file line by line
	public static void readDriver(String filename) throws FileNotFoundException{
		Scanner scan = new Scanner(new File(filename));
		String line;
		
		//Driver variables
		String type;
		double wx, wy, wz, theta, sc, tx, ty, tz;
		String model;
		
		
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
			
			//Save obj
			exportObj(obj, model, filename);
			
		}
		scan.close();
	}
	
	public static void exportObj(ModelObject obj, String model, String driverName){
		//Create driver folder
		//System.out.println(driverName);
		String[] a = driverName.split("\\.");
		new File(a[0]).mkdir();
		//System.out.println(a[0]);
		
		//file name
		int suffix = 0;
		String[] modelFile = model.split("\\.");
		//System.out.println(a[0] + "\\" + modelFile[0] + "_mw" + String.format("%02d", suffix) + ".obj");
		File testFile = new File(a[0] + "\\" + modelFile[0] + "_mw" + String.format("%02d", suffix) + ".obj");
		boolean exists = testFile.exists();
		while(exists){
			suffix++;
			testFile = new File(a[0] + "\\" + modelFile[0] + "_mw" + String.format("%02d", suffix) + ".obj");
			exists = testFile.exists();
		}
		
		try {
			FileWriter fw = new FileWriter(a[0] + "\\" + modelFile[0] + "_mw" + String.format("%02d", suffix) + ".obj");
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
	
	public static void printMatrix(RealMatrix a){
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
		readDriver(args[0]);
	}

}
