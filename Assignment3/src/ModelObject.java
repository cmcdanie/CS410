//Colin McDaniel
//CS 410 Assignment #1
//CSU ID: 830293766

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

public class ModelObject {

	public String modelName;				//Name of model
	public ArrayList<String> objTxt;		//List of lines of modelObject file
	public ArrayList<double[]> vertices;	//List of vertices
	public ArrayList<int[]> faces;			//List of faces
	
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
			if(words[0].equals("f")) {
				int[] token = new int[words.length - 1];
				for(int i = 0; i < words.length - 1; i++) {
					//System.out.println(words[i + 1]);
					String[] letters = words[i + 1].split("//");
					token[i] = Integer.parseInt(letters[0]);
				}
				faces.add(token);
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

	public static void main(String[] args) throws FileNotFoundException {
		//readFile("cube.obj");

	}

}
