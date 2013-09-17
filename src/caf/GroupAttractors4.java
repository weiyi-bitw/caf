package caf;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import obj.Attractor;
import obj.Chromosome;
import obj.DataFile;
import obj.GeneSet;
import obj.ValIdx;
import worker.Converger;
public class GroupAttractors4 {
	static class DistPair implements Comparable<DistPair>{
		static int n = 1000000;
		String x;
		String y;
		float sim; // similarities
		
		DistPair(String a, String b, float sim){
			if(a.compareTo(b) < 0){ // x must always be smaller than y
				this.x = a;
				this.y = b;
			}else{
				this.x = b;
				this.y = a;
			}
			this.sim = sim;
		}
		boolean contains(String k){
			return this.x.equals(k) || this.y.equals(k);
		}
		public boolean equals(Object other){
			boolean result = false;
	        if (other instanceof DistPair) {
	        	DistPair that = (DistPair) other;
	        	result = this.x.equals(that.x) && this.y.equals(that.y);
	        }
	        return result;
		}
		public int hashCode(){
			return(x.hashCode() * n + y.hashCode());
		}
		public int compareTo(DistPair other) {
			return -Double.compare(this.sim, other.sim);
		}
		static void setTotalIdx(int n){
			DistPair.n = n;
		}
		public String toString(){
			String s = this.x + "\t" + this.y + ":" + this.sim;
			return s;
		}
	}
	
	private static ArrayList<Attractor> parseAttractorInOneFile(
			String file, HashMap<String, Integer> rowmap, int minSize) throws IOException {
		
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = br.readLine();
		ArrayList<Attractor> allAttractors = new ArrayList<Attractor>();
		while(line != null){ // each line is an attractor
			String[] tokens = line.split("\t");
			String name = tokens[0];
			int nt = tokens.length;
			if(nt - 2 < minSize) {
				line = br.readLine();
				continue;
			}
			float numChild = Float.parseFloat(tokens[1].split(":")[1]);
			float z = 0;
			HashMap<String, Float> zMap = new HashMap<String, Float>();
			ArrayList<String> genes = new ArrayList<String>();
			for(int j = 2; j < nt; j++){
				String[] t2 = tokens[j].split(":");
				if(j==2){
					z = Float.parseFloat(t2[2]);
				}
				genes.add(t2[0]);
				zMap.put(t2[0], z = Float.parseFloat(t2[2]));
			}
			allAttractors.add(new Attractor(name, genes, zMap, rowmap, numChild, z));
			line = br.readLine();
		}
		br.close();
		return allAttractors;
	}
	
	private static float[] getMetaGene(float[][] data, Attractor a, int n){
		int m = a.size();
		float[] out = new float[n];
		for(int j = 0; j < n; j++){
			for(Integer i : a.geneIndices()){
				out[j] += data[i][j];
			}
			out[j] /= m;
		}
		return out;
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] arg) throws Exception {
		// TODO Auto-generated method stub
		String path = "/home/weiyi/workspace/javaworks/caf/output/ov.tcga.agi.minorm.z8/";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		int minSize = 10;
		float zScore = 8f;
		float reconvergeTh = 0.8f; // reconverge threshold
		float ovlpTh = 0.5f; // overlap threshold
		boolean rowNormalization = false;
		boolean MINormalization = true;
		
		System.out.println("Loading files...");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/brca/gse2034/ge.13271x286.var.txt");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/brca/tcga/ge/ge.17814x536.knn.txt");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/coad/gse14333/ge.20765x290.var.txt");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/coad/tcga/ge/ge.17814x154.knn.txt");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/ov/gse9891/ge.20765x285.var.txt");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/ov/tcga/ge/ge.12042x582.txt");
		DataFile ma = DataFile.parse("/home/weiyi/workspace/data/ov/tcga/ge/ge.17814x584.knn.txt");
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		if(rowNormalization) ma.normalizeRows();
		float[][] data = ma.getData();
		
		HashMap<String, Integer> rowmap = ma.getRows();
		
		ArrayList<String> probeNames = ma.getProbes();
		System.out.print("Loading gene sets...");
		
		ArrayList<Attractor> allGeneSet = parseAttractorInOneFile(path + "attractors.gwt", rowmap, minSize);
		int N = allGeneSet.size();
		
		System.out.println(N + " gene sets are loaded.");
		Collections.sort(allGeneSet);
		/*PrintWriter pw2 = new PrintWriter(new FileWriter(path + "SortedAttractors.gwt"));
		for(Attractor a: allGeneSet){
			pw2.println(a);
		}
		pw2.close();
		*/
		// Merging using UPGMA
		ArrayList<DistPair> allDist = new ArrayList<DistPair>();
		for(int i = 0; i < N-1; i++){
			Attractor ai = allGeneSet.get(i);
			for(int j = i+1; j < N; j++){
				Attractor aj = allGeneSet.get(j);
				float ovlpRatio = ai.similarities(aj);
				if(ovlpRatio > ovlpTh){
					allDist.add(new DistPair(ai.name(), aj.name(), ovlpRatio));
				}
			}
		}
		Collections.sort(allDist);
		
		Converger cvg = new Converger(0, 1, System.currentTimeMillis(), "ZSCORE", 100, false);
		cvg.miNormalization(MINormalization);
		cvg.setZThreshold(zScore);
		cvg.setAttractorSize(minSize);
		cvg.setMIParameter(6, 3);
		GeneSet.setOverlapRatio(ovlpTh);
		
		while(allDist.size() > 0){
			DistPair target = allDist.get(0); // merge the most similar pair
			System.out.print(allDist.size() + "\t" + target.toString() + "...");
			
			int xIdx = allGeneSet.indexOf(new Attractor(target.x));
			int yIdx = allGeneSet.indexOf(new Attractor(target.y));
			if(xIdx > yIdx){
				int tmp = yIdx;
				yIdx = xIdx;
				xIdx = tmp;
			}
			
			//System.out.println(target.x + "(" + xIdx + ")" + "\t" + target.y + "(" + yIdx + ")");
			
		// 1. merge gene set, put in new key
			if(target.sim == 1){
				allGeneSet.get(xIdx).merge(allGeneSet.get(yIdx));
			}else if(target.sim >= reconvergeTh){
				allGeneSet.get(xIdx).fight(allGeneSet.get(yIdx));
			}else{
				allGeneSet.get(xIdx).merge(allGeneSet.get(yIdx));
				float[] vec = getMetaGene(data,allGeneSet.get(xIdx),n );
				ArrayList<ValIdx> metaIdx = cvg.findAttractor(data, vec);
				allGeneSet.get(xIdx).reset(metaIdx, probeNames);
			}
			String mergedName = allGeneSet.get(xIdx).name();
			System.out.println(" into " + mergedName + " (" + xIdx + ")." + allGeneSet.get(xIdx).strength());
		// 2. remove merged gene set, remove keys
			allGeneSet.remove(yIdx);
			
		// 3. remove distances regarding these two gene set
			for(int i = allDist.size()-1; i >= 0; i--){
				if(allDist.get(i).contains(target.x) || allDist.get(i).contains(target.y)){
					allDist.remove(i);
				}
			}
			
			if(allGeneSet.get(xIdx).size() < minSize){
				allGeneSet.remove(xIdx);
				continue;
			}
		// 4. add new distances
			for(int i = 0; i < allGeneSet.size(); i++){
				if(i != xIdx){
					Attractor ai = allGeneSet.get(i);
					float ovlpRatio = ai.similarities(allGeneSet.get(xIdx));
					if(ovlpRatio > ovlpTh){
						allDist.add(new DistPair(ai.name(), mergedName, ovlpRatio));
					}
				}
			}
			Collections.sort(allDist);
		}
		N = allGeneSet.size();
		
		Collections.sort(allGeneSet);
		for(int i = N-1; i >=0; i--){
			Attractor a = allGeneSet.get(i);
			if(a.size() < minSize){
				allGeneSet.remove(i);
				continue;
			}
			for(int j = 0; j < i; j++){
				if(a.similarities(allGeneSet.get(j)) > ovlpTh){
					allGeneSet.remove(i);
					break;
				}
			}
		}
		
		N = allGeneSet.size();
		System.out.println(N + " meta-attractors left.");
		PrintWriter pw2 = new PrintWriter(new FileWriter(path + "MetaAttractors.gwt"));
		for(int i = 0; i < N; i++ ){
			
			Attractor a = allGeneSet.get(i);
			//System.out.println(a.name());
			pw2.println(a.toStringGenesOnly());
		}
		pw2.close();
		
		System.out.println(N + " Done.");
		
		/*System.out.println("Reconverging...");
		
		ArrayList<GeneSet> convergedGeneSets = new ArrayList<GeneSet>();
		for(int i = 0; i < N; i++){
			Attractor a = allGeneSet.get(i);
			int k = 0;
			for(Integer ii : a.geneIdx){
				k = ii; break;
			}
			Chromosome chr = chromap.get(k);
			System.out.println(i + "\t" + a.getName() + "...");
			float[] vec = getMetaGene(data, a.geneIdx, n);
			ArrayList<ValIdx> metaIdx = CNV? cvg.findCNV(data, vec, chr):cvg.findAttractor(data, vec);
			GeneSet gs;
			if(CNV){
				gs = new GeneSet(a.getName() + "_" + chr.name(), metaIdx.toArray(new ValIdx[0]), a.strength);
			}else{
				gs = new GeneSet(a.getName(), metaIdx.toArray(new ValIdx[0]), a.strength);
			}
			if(gs.size() < minSize){
				System.out.println("Too small, filtered.");
				continue;
			}
			if(gs.overlapWith(convergedGeneSets)){
				System.out.println("Overlap, filtered.");
				continue;
			}
			convergedGeneSets.add(gs);
			if(convergedGeneSets.size() >= 20){
				break;
			}
		}
		
		PrintWriter pw = new PrintWriter(new FileWriter(path + "ConvergedMetaAttractors.gwt"));
		GeneSet.setProbeNames(probeNames);
		for(GeneSet gs : convergedGeneSets){
			pw.print(gs.getName() + "\t" + gs.size() + ":" + gs.getNumChild() + "\t");
			pw.println(gs.toProbes());
		}
		pw.close();*/
	}
}
