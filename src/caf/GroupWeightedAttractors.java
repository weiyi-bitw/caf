package caf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import obj.DataFile;
import obj.ValIdx;

public class GroupWeightedAttractors {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		String path = "/home/weiyi/workspace/javaworks/caf/output/629/";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		DataFile ma = DataFile.parse("/home/weiyi/workspace/data/brca/gse2034/ge.13271x286.var.txt");
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		ArrayList<String> genes = ma.getProbes();
		
		BufferedReader br = new BufferedReader(new FileReader(path + "attractors.gwt"));
		BufferedReader br2 = new BufferedReader(new FileReader(path + "attractees.gwt"));
		String line = br.readLine();
		String line2 = br2.readLine();
		PrintWriter pw = new PrintWriter(new FileWriter(path + "attractors.top10.gwt"));
		
		while(line != null){
			ArrayList<ValIdx> vec = new ArrayList<ValIdx>();
			String[] tokens = line.split("\t");
			int numBasins = line2.split("\t").length - 1;
			String name = tokens[0];
			pw.print(name + "\t" + numBasins);
			float sum = 0;
			for(int i = 0; i < m; i++){
				float f= Float.parseFloat(tokens[i+1]);
				float w = (float) Math.exp(5 * Math.log(f));
				vec.add(new ValIdx(i, w));
				
				sum += w;
			}
			Collections.sort(vec);
			for(int i = 0; i < 10; i++){
				pw.print("\t" + genes.get(vec.get(i).idx) + "(" + vec.get(i).val/sum + ")");
			}pw.println();
			
			line2 = br2.readLine();
			line = br.readLine();
		}
		br.close();
		pw.close();
		System.out.println("Done.");
	}

}
