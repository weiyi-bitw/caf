package caf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import obj.Annotations;
import obj.DataFile;
import obj.Genome;
import obj.ValIdx;

public class GroupWeightedAttractor {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		String path = "/home/weiyi/workspace/javaworks/caf/output/weighted/ov.gse9891.jetset.mean";
		if(path.endsWith("/")){
			path = path.substring(0, path.length()-1);
		}
		System.out.println("Loading files...");
		int IDX = 2;
		int outgenes = 30;
		int quantile = 30;
		
		String[] dataFiles = {
				"/home/weiyi/workspace/data/brca/gse2034/ge.12160x286.jetset.mean.txt",
				"/home/weiyi/workspace/data/coad/gse14333/ge.19189x290.jetset.mean.txt",
				"/home/weiyi/workspace/data/ov/gse9891/ge.19189x285.jetset.mean.txt",
				"/home/weiyi/workspace/data/brca/tcga/ge/ge.17814x536.knn.txt",
				"/home/weiyi/workspace/data/coad/tcga/ge/ge.17814x154.knn.txt",
				"/home/weiyi/workspace/data/ov/tcga/ge/ge.12042x582.txt"
		};
		
		
		/*String[] annots = {
				"/home/weiyi/workspace/data/annot/affy/u133a/annot.jetset.csv",
				"/home/weiyi/workspace/data/annot/tcga/4502a073/annot.csv",
				"/home/weiyi/workspace/data/annot/affy/u133p2/annot.jetset.csv",
				"/home/weiyi/workspace/data/annot/tcga/4502a073/annot.csv",
				"/home/weiyi/workspace/data/annot/affy/u133p2/annot.jetset.csv",
				"/home/weiyi/workspace/data/annot/affy/u133a/annot.jetset.csv"
		};*/
		
		
		
		//Annotations annot = Annotations.parseAnnotations(annots[IDX]);
		DataFile ma = DataFile.parse(dataFiles[IDX]);
		
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		ArrayList<String> genes = ma.getProbes();
		
		BufferedReader br = new BufferedReader(new FileReader(path + "/attractors.gwt"));
		BufferedReader br2 = new BufferedReader(new FileReader(path + "/attractees.gwt"));
		String line = br.readLine();
		String line2 = br2.readLine();
		/*String[] tokens = line.split("\t");
		int nt = tokens.length;
		String[] allgenes = new String[nt-2];
		int mg = allgenes.length;
		System.arraycopy(tokens, 2, allgenes, 0, nt-2);*/
		
		new File(path + "/../mergeroom").mkdir();
		String outFileName = path.substring(path.lastIndexOf("/"));
		PrintWriter pw = new PrintWriter(new FileWriter(path + "/../mergeroom/" + outFileName));
		PrintWriter pw2 = new PrintWriter(new FileWriter(path + "/attractees.decoded.gwt"));
		
		//line = br.readLine();
		
		while(line != null){
			ArrayList<ValIdx> vec = new ArrayList<ValIdx>();
			String[] tokens = line.split("\t");
			int numBasins = line2.split("\t").length - 2;
			String name = tokens[0];
			String[] t2 = line2.split("\t");
			pw.print(name + "\t" + numBasins);
			pw2.print(name + "\t" + numBasins);
			float sum = 0;
			for(int i = 0; i < m; i++){
				//float f= Float.parseFloat(tokens[i+2]);
				//float w = (float) Math.exp(5 * Math.log(f));
				float w = Float.parseFloat(tokens[i+2]);
				vec.add(new ValIdx(i, w));
					
				sum += w;
			}
			for(int i = 0; i < numBasins; i++){
				pw2.print("\t" + genes.get(Integer.parseInt(t2[i+2])) );
			}
			Collections.sort(vec);
			for(int i = 0; i < outgenes; i++){
				String g = genes.get(vec.get(i).idx);
				pw.print("\t" + g);
			}
			pw.print("\t" + vec.get(quantile-1).val);
			pw.println();
			pw2.println();
			
			line2 = br2.readLine();
			line = br.readLine();
		}
		br.close();
		pw.close();
		pw2.close();
		System.out.println("Done.");
	}

}
