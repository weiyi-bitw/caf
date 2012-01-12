package caf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import obj.DataFile;

public class GroupAttractors2 {
	static class Attractor{
		HashSet<Integer> geneIdx;
		int sz;
		int strength; // number of attractees
		HashSet<Attractor> child = null;
		Attractor parent = null;
		
		Attractor(HashSet<Integer> geneIdx, int sz, int strength){
			this.geneIdx = geneIdx;
			this.sz = sz;
			this.strength = strength;
		}
		Attractor(Attractor a){
			this.geneIdx = new HashSet<Integer>();
			this.geneIdx.addAll(a.geneIdx);
			this.sz = a.sz;
			this.strength = a.strength;
			this.addChild(a);
		}
		void addChild(Attractor a){
			if(this.child==null){
				this.child = new HashSet<Attractor>();
			}
			this.child.add(a);
		}
		void setParent(Attractor a){
			this.parent = a;
		}
		
		int overlap(Attractor a){
			int ovlp = 0;	
			for(int i : a.geneIdx){
				if(this.geneIdx.contains(i)){
					ovlp++;
				}
			}
			return ovlp;
		}
		
	}
	private static ArrayList<Attractor> parseAttractorInOneFile(
			String file, HashMap<String, Integer> rowmap) throws IOException {
		
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = br.readLine();
		ArrayList<Attractor> allAttractors = new ArrayList<Attractor>();
		while(line != null){ // each line is an attractor
			String[] tokens = line.split("\t");
			String name = tokens[0];
			int nt = tokens.length;
			int numChild = Integer.parseInt(tokens[1].split(":")[1]);
			HashSet<Integer> gidx = new HashSet<Integer>();
			for(int j = 2; j < nt; j++){
				String[] t2 = tokens[j].split(":");
				int i = rowmap.get(t2[0]);
				gidx.add(i);
			}
			allAttractors.add(new Attractor(gidx, nt-2, numChild));
			line = br.readLine();
		}
		br.close();
		return allAttractors;
	}
	
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		String path = "/home/weiyi/workspace/javaworks/caf/output/freeRun/coad.gse17536.rownorm.z7";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		int minSize = 10;
		
		System.out.println("Loading files...");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/brca/gse2034/ge.13271x286.var.txt");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/brca/tcga/ge/ge.17814x536.knn.txt");
		DataFile ma = DataFile.parse("/home/weiyi/workspace/data/coad/gse17536/ge.20765x177.var.txt");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/coad/tcga/ge/ge.17814x154.knn.txt");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/ov/gse9891/ge.20765x285.var.txt");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/ov/tcga/ge/ge.12042x582.txt");
		//ma.normalizeRows();
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		HashMap<String, Integer> rowmap = ma.getRows();
		
		ArrayList<String> probeNames = ma.getProbes();
		System.out.print("Loading gene sets...");
		
		ArrayList<Attractor> allGeneSet = parseAttractorInOneFile(path + "attractors.gwt", rowmap);
		int N = allGeneSet.size();
		System.out.println(N + " gene sets are loaded.");
		
		
		System.out.println("Merging gene sets...");
		int progress = 0;
		int numGroups = 0;
		for(int i = 0; i < N; i++){
			Attractor aa = allGeneSet.get(i);
			for(int j = 0; j < N; j++){
				if(j == i) continue;
				Attractor ab = allGeneSet.get(j);
				int ovlp = aa.overlap(ab);
				if(ovlp > 0){
					
				}
			}
		}
		
	}

	
}
