//Colin McDaniel
//CS 410 Assignment #1
//CSU ID: 830293766

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class ModelObject {

	public static String modelName;
	public static ArrayList<String> objTxt;
	public static ArrayList<double[]> vertices;
	
	public ModelObject(String n) throws FileNotFoundException{
		modelName = n;
		readFile(n);
		
	}
	
	public static void readFile(String fileName) throws FileNotFoundException{
		Scanner scan = new Scanner(new File(fileName));
		String line;
		objTxt = new ArrayList<String>();
		vertices = new ArrayList<double[]>();
		while(scan.hasNextLine()){
			line = scan.nextLine();
			objTxt.add(line);
			String[] words = line.split("\\s");
			if(words[0].equals("v")){
				double[] token = new double[words.length - 1];
				for(int i = 0; i < words.length - 1; i++){
					token[i] = Double.parseDouble(words[i + 1]);
				}
				
				vertices.add(token);
			}
		}
		/*
		for(int i = 0; i < vertices.size(); i++){
			for(int j = 0; j < 3; j++){
				System.out.print(vertices.get(i)[j] + " ");
			}
			System.out.println(" ");
		}
		*/	
		
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
	
	

	public static void main(String[] args) throws FileNotFoundException {
		readFile("cube.obj");

	}

}
