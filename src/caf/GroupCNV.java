package caf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import obj.Genome;
import obj.IntPair;

public class GroupCNV {

	static class CNV implements Comparable<CNV>{
		static Genome gn;
		static int quantile = 5;
		static int winSize = 100;
		
		String basin;
		String chrArm;
		String[] genes;
		float score;
		IntPair range;
		
		CNV(String basin, String chrArm, String[] genes, float score, IntPair range){
			this.basin = basin;
			this.chrArm = chrArm;
			this.genes = genes;
			this.score = score;
			this.range = range;
		}
		
		boolean ovlapWith(CNV other){
			if(!this.chrArm.equals(other.chrArm)){
				return false;
			}
			return this.range.overlapWith(other.range);
		}
		
		static void linkToGenome(Genome gn){
			CNV.gn = gn;
		}
		
		public int compareTo(CNV other) {
			return -Double.compare(this.score, other.score);
		}
		
		public String toString(){
			String s = gn.getChrBand(genes[0]);
			for(String g : genes){
				s += "\t" + g;
			}
			s += "\t" + score + "\t" + range.x + "-" + range.y;
			return s;
		}
		
		static void setWindowSize(int winSize){
			CNV.winSize = winSize;
		}
		
		static CNV parseCNV(String line, int n){
			String[] tokens = line.split("\t");
			int nt = tokens.length;
			if(nt-2 < n){
				return null;
			}
			String basin = tokens[0];
			String chrArm = tokens[1];
			String[] genes = new String[n];
			float[] mis = new float[n];
			float score = 0;
			int x= Integer.MAX_VALUE;
			int y = -1;
			
			int x5 = Integer.MAX_VALUE;
			int y5 = -1;
			for(int i = 0; i < n; i++){
				String[] t2 = tokens[i+2].split(":");
				int idx = gn.getIdx(t2[0]);
				
				if(idx < x) x = idx;
				if(idx > y) y = idx;
				if(i < quantile){
					if(idx < x5) x5 = idx;
					if(idx > y5) y5 = idx;
				}
				genes[i] = t2[0];
				mis[i] = Float.parseFloat(t2[1]);
				
				if(i == quantile-1){
					score = Float.parseFloat(t2[1]);
				}
			}
			
			if(mis[0] - mis[1] > 0.4 || y - x > winSize){
				return null;
			}
			
			
			return new CNV(basin, chrArm, genes, score, new IntPair(x, y));
		}
		
		
		
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String[] dataFiles = {
				//"/home/weiyi/workspace/data/brca/gse2034/ge.12160x286.jetset.mean.txt",
				//"/home/weiyi/workspace/data/coad/gse14333/ge.19189x290.jetset.mean.txt",
				"/home/weiyi/workspace/data/ov/gse9891/ge.19189x285.jetset.mean.txt",
				//"/home/weiyi/workspace/data/brca/tcga/ge/ge.17814x536.knn.txt",
				//"/home/weiyi/workspace/data/coad/tcga/ge/ge.17814x154.knn.txt"
				"/home/weiyi/workspace/data/ov/tcga/ge/ge.12042x582.txt"
		};
		
		final String[] outputDirs={
				//"brca.gse2034.jetset.mean",
				//"coad.gse14333.jetset.mean",
				"ov.gse9891.jetset.mean",
				//"brca.tcga",
				//"coad.tcga"
				"ov.tcga"
		};
		
		String outPath = "/home/weiyi/workspace/javaworks/caf/output/cnv/";
		if(!outPath.endsWith("/")){
			outPath = outPath + "/";
		}
		int excludeSize = 25;
		int quantile = 5;
		int loadIn = 10;
		int windowSize = 300;
		
		System.out.println("Loading files...");
		final String geneLocFile = "/home/weiyi/workspace/data/annot/affy/u133p2/gene.location4";
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		CNV.linkToGenome(gn);
		CNV.setWindowSize(windowSize);
		
		for(int qq = 0; qq < dataFiles.length; qq++){
			ArrayList<CNV> allCNVs = new ArrayList<CNV>();
			
			BufferedReader br = new BufferedReader(new FileReader(outPath + outputDirs[qq] + "/caf.txt"));
			String line = br.readLine();
			while(line != null){
				CNV c = CNV.parseCNV(line, loadIn);
				if(c != null){
					allCNVs.add(c);
				}
				line = br.readLine();
			}
			br.close();
			
			System.out.println(allCNVs.size() + " CNVs loaded.");
			/*Collections.sort(allCNVs);
			for(int i = allCNVs.size() -1; i > 0; i--){
				CNV c = allCNVs.get(i);
				for(int j = 0 ; j < i; j++){
					if(c.ovlapWith(allCNVs.get(j))){
						allCNVs.remove(i);
						break;
					}
				}
			}
			System.out.println(allCNVs.size() + " CNVs remained after filtering.");*/
			
			new File(outPath + "mergeroom").mkdir();
			PrintWriter pw = new PrintWriter(new FileWriter(outPath + "mergeroom/" + outputDirs[qq]));
			for(CNV c : allCNVs){
				pw.println(c);
			}
			pw.close();
			
		} // END qq iteration
		System.out.println("Done.");
	}

}
