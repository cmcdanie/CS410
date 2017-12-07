//Colin McDaniel
//CS 410 Assignment #1
//CSU ID: 830293766

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class ModelObject {

	private String modelName;				//Name of model
	private ArrayList<String> objTxt;		//List of lines of modelObject file
	private ArrayList<double[]> vertices;	//List of vertices
	private ArrayList<int[]> faces;			//List of faces
	private String materialFile;					//File name of material file
	private Map<String, MaterialObject> materials;			//List of materials found in object file
	private Map<Integer, MaterialObject> faceMaterial;		//Map of face index to material name
	private MaterialObject defaultMaterial;
	
	
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
		faceMaterial = new HashMap<Integer, MaterialObject>();
		String currentMtl = "";
		
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
				if(currentMtl.equals("")) {
					currentMtl = defaultMaterial.getName();
				}
				faceMaterial.put(faces.size() - 1, materials.get(currentMtl));
			}
			else if(words[0].equals("mtllib")) {
				readMaterialFile(words[1]);
			}
			
			else if(words[0].equals("usemtl")) {
				currentMtl = words[1];
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
		materials = new HashMap<String, MaterialObject>();
		String currentMtl = "";
		
		while(scan.hasNextLine()){
			line = scan.nextLine();
			String[] words = line.split("\\s");
			
			
			if(words[0].equals("newmtl")) {
				currentMtl = words[1];
				MaterialObject mtl = new MaterialObject(currentMtl);
				materials.put(currentMtl, mtl);	
			}
			
			else if(words[0].equals("Ns")) {
				double token = Double.parseDouble(words[1]);
				materials.get(currentMtl).setNs(token);
			}
			
			else if(words[0].equals("Ka")) {
				double[] token = new double[words.length - 1];
				for(int i = 0; i < words.length - 1; i++){
					token[i] = Double.parseDouble(words[i + 1]);
				}
				RealVector Ka = new ArrayRealVector(token);
				materials.get(currentMtl).setAmbient(Ka);
				
			}
			
			else if(words[0].equals("Kd")) {
				double[] token = new double[words.length - 1];
				for(int i = 0; i < words.length - 1; i++){
					token[i] = Double.parseDouble(words[i + 1]);
				}
				RealVector Kd = new ArrayRealVector(token);
				materials.get(currentMtl).setDiffuse(Kd);
			}
			
			else if(words[0].equals("Ks")) {
				double[] token = new double[words.length - 1];
				for(int i = 0; i < words.length - 1; i++){
					token[i] = Double.parseDouble(words[i + 1]);
				}
				RealVector Ks = new ArrayRealVector(token);
				materials.get(currentMtl).setSpecular(Ks);
			}
			
			else if(words[0].equals("Kr")) {
				double[] token = new double[words.length - 1];
				for(int i = 0; i < words.length - 1; i++){
					token[i] = Double.parseDouble(words[i + 1]);
				}
				RealVector Kr = new ArrayRealVector(token);
				materials.get(currentMtl).setAttenuation(Kr);;
			}
		}
		
		defaultMaterial = materials.get(currentMtl);
		
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
	
	public MaterialObject getFaceMaterial(int faceIndex) {
		return faceMaterial.get(faceIndex);
	}

	public static void main(String[] args) throws FileNotFoundException {
		ModelObject a = new ModelObject("simface.obj");
		a.readFile("cube_centered.obj");

	}

}
