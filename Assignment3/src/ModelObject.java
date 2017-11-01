//Colin McDaniel
//CS 410 Assignment #1
//CSU ID: 830293766

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Scanner;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class ModelObject {

	public String modelName;				//Name of model
	public ArrayList<String> objTxt;		//List of lines of modelObject file
	public ArrayList<double[]> vertices;	//List of vertices
	public ArrayList<int[]> faces;			//List of faces
	public String materialFile;
	public LinkedList<MaterialObject> materials;
	
	
	
	public ModelObject(String n) throws FileNotFoundException{
		modelName = n;
		readFile(n);
		
	}
	
	public void readFile(String fileName) throws FileNotFoundException{
		Scanner scan = new Scanner(new File(fileName));
		String line;
		objTxt = new ArrayList<String>();
		vertices = new ArrayList<double[]>();
		faces = new ArrayList<int[]>();
		
		//Scan file
		while(scan.hasNextLine()){
			line = scan.nextLine();
			objTxt.add(line);
			String[] words = line.split("\\s");
			//If line stars with "v", add vertex to list
			if(words[0].equals("v")){
				double[] token = new double[words.length - 1];
				for(int i = 0; i < words.length - 1; i++){
					token[i] = Double.parseDouble(words[i + 1]);
				}
				
				vertices.add(token);
			}
			else if(words[0].equals("f")) {
				int[] token = new int[words.length - 1];
				for(int i = 0; i < words.length - 1; i++) {
					//System.out.println(words[i + 1]);
					String[] letters = words[i + 1].split("//");
					token[i] = Integer.parseInt(letters[0]);
				}
				faces.add(token);
			}
			else if(words[0].equals("mtllib")) {
				readMaterialFile(words[1]);
			}
			
			else if(words[0].equals("Ka")) {
				double[] token = new double[words.length - 1];
				for(int i = 0; i < words.length - 1; i++){
					token[i] = Double.parseDouble(words[i + 1]);
				}
				materials.get(materials.size()).setAmbient(new ArrayRealVector(token));
			}
			
			else if(words[0].equals("Kd")) {
				double[] token = new double[words.length - 1];
				for(int i = 0; i < words.length - 1; i++){
					token[i] = Double.parseDouble(words[i + 1]);
				}
				materials.get(materials.size()).setDiffuse(new ArrayRealVector(token));
			}
			
			else if(words[0].equals("Ks")) {
				double[] token = new double[words.length - 1];
				for(int i = 0; i < words.length - 1; i++){
					token[i] = Double.parseDouble(words[i + 1]);
				}
				materials.get(materials.size()).setSpecular(new ArrayRealVector(token));
			}
			
		}
		/*
		for(int i = 0; i < faces.size(); i++){
			for(int j = 0; j < 3; j++){
				System.out.print(faces.get(i)[j] + " ");
			}
			System.out.println(" ");
		}
		*/
		scan.close();
	}
	
	public void readMaterialFile(String filename) throws FileNotFoundException {
		Scanner scan = new Scanner(new File(filename));
		String line;
		materials = new LinkedList<MaterialObject>();
		
		while(scan.hasNextLine()){
			line = scan.nextLine();
			String[] words = line.split("\\s");
			
			
			if(words[0].equals("newmtl")) {
				MaterialObject mtl = new MaterialObject(words[1]);
				materials.add(mtl);
			}
			
			else if(words[0].equals("Ka")) {
				double[] token = new double[words.length - 1];
				for(int i = 0; i < words.length - 1; i++){
					token[i] = Double.parseDouble(words[i + 1]);
				}
				RealVector Ka = new ArrayRealVector(token);
				materials.getLast().setAmbient(Ka);
				
			}
			
			else if(words[0].equals("Kd")) {
				double[] token = new double[words.length - 1];
				for(int i = 0; i < words.length - 1; i++){
					token[i] = Double.parseDouble(words[i + 1]);
				}
				RealVector Kd = new ArrayRealVector(token);
				materials.getLast().setDiffuse(Kd);
			}
			
			else if(words[0].equals("Ks")) {
				double[] token = new double[words.length - 1];
				for(int i = 0; i < words.length - 1; i++){
					token[i] = Double.parseDouble(words[i + 1]);
				}
				RealVector Ks = new ArrayRealVector(token);
				materials.getLast().setSpecular(Ks);
			}
		}
		
		scan.close();
	}
	
	public double[] getVertex(int index){
		double[] result = vertices.get(index);
		return result;
	}
	
	public void setVertex(int index, double[] vertex){
		vertices.set(index, vertex);
	}
	
	public ArrayList<String> getObjTxt(){
		return objTxt;
	}
	
	public String getObjTxt(int index){
		return objTxt.get(index);
	}
	
	public int numOfVert(){
		return vertices.size();
	}
	
	public int[] getFace(int index) {
		int[] result = faces.get(index);
		return result;
	}
	
	public int getNumOfFaces() {
		return faces.size();
	}
	
	public String getName() {
		return modelName;
	}
	
	public MaterialObject getMaterial(int index) {
		return materials.get(index);
	}

	public static void main(String[] args) throws FileNotFoundException {
		ModelObject a = new ModelObject("cube_centered.obj");
		//a.readFile("cube_centered.obj");

	}

}
