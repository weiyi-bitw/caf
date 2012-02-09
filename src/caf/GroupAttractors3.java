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

import obj.Chromosome;
import obj.DataFile;
import obj.GeneSet;
import worker.Converger;
import worker.Converger.ValIdx;

public class GroupAttractors3 {
	static class DistPair implements Comparable<DistPair>{
		static int n = 10000;
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
	
	static class Attractor implements Comparable<Attractor>{
		static int counter = 0;
		static ArrayList<String> geneNames;
		String name;
		HashSet<Integer> geneIdx;
		int sz;
		int strength; // number of attractees
		HashSet<Attractor> child = null;
		HashSet<String> stomach = null;
		Attractor parent = null;
		
		
		Attractor(String name, HashSet<Integer> geneIdx, int sz, int strength){
			this.name = name;
			this.geneIdx = geneIdx;
			this.sz = sz;
			this.strength = strength;
			this.stomach = new HashSet<String>();
			this.stomach.add(this.name);
		}
		Attractor(String name){
			this.name = name;
		}
		Attractor(Attractor a){
			this.name = "MetaAttractor" + counter;
			counter++;
			this.sz = a.sz;
			this.child = new HashSet<Attractor>();
			this.addChild(a);
		}
		void combine(Attractor a){ //combine the children of two attractors 
			for(Attractor aa : a.child){
				this.addChild(aa);
			}
		}
		void merge(Attractor a){
			this.name = "MetaAttractor" + counter;
			counter++;
			this.geneIdx.addAll(a.geneIdx);
			this.sz = this.geneIdx.size();
			this.strength += a.strength;
			this.stomach.add(a.name);
		}
		
		void addChild(Attractor a){
			if(this.child==null){
				this.child = new HashSet<Attractor>();
			}
			a.parent = this;
			this.child.add(a);
		}
		void setParent(Attractor a){
			this.parent = a;
		}
		String getName(){
			return name;
		}
		public String toString(){
			String s = this.name + "\t" + this.sz + ":" + this.strength;
			for(int i : this.geneIdx){
				s += ("\t" + geneNames.get(i));
			}
			return s;
		}
		public boolean equals(Object other){
			boolean result = false;
	        if (other instanceof Attractor) {
	        	Attractor that = (Attractor) other;
	        	result = this.name.equals(that.name);
	        }
	        return result;
		}
		public int hashCode(){
			return(this.name.hashCode());
		}
		void cronus(){
			this.geneIdx = new HashSet<Integer>();
			this.strength = 0;
			for(Attractor a : this.child){
				this.geneIdx.addAll(a.geneIdx);
				this.strength += a.strength;
			}
			this.sz = this.geneIdx.size();
		}
		int getOverlap(Attractor a){
			int ovlp = 0;	
			for(int i : a.geneIdx){
				if(this.geneIdx.contains(i)){
					ovlp++;
				}
			}
			return ovlp;
		}
		boolean overlap(Attractor a){
			int th = Math.min(this.sz, a.sz)/2;
			int ovlp = 0;
			for(int i : a.geneIdx){
				if(this.geneIdx.contains(i)){
					ovlp++;
					if(ovlp > th)return true;
				}
			}
			return false;
		}
		float similarities(Attractor a){
			float ovlp = 0;
			float minSz = (float) Math.min(this.sz, a.sz);
			for(int i : a.geneIdx){
				if(this.geneIdx.contains(i)){
					ovlp++;
				}
			}
			return ovlp / minSz;
		}
		static void setGeneNames(ArrayList<String> geneNames){
			Attractor.geneNames = geneNames;
		}
		public int compareTo(Attractor other) {
			return -Double.compare(this.strength, other.strength);
		}
	}
	static ArrayList<Chromosome> parseChromGenes(String file, ArrayList<String> genes, HashMap<String, Integer> geneMap, HashMap<Integer, Chromosome> chromap) throws IOException{
		ArrayList<Chromosome> chrs = new ArrayList<Chromosome>();
		Chromosome.setGeneMap(geneMap);
		Chromosome.setGeneNames(genes);
		BufferedReader br = new BufferedReader(new FileReader(file));
		br.readLine(); // first line header
		
		String line = br.readLine();
		while(line != null){
			//System.out.println(line);
			String[] tokens = line.split("\t");
			int nt = tokens.length;
			Chromosome chr = new Chromosome(tokens[nt-1]);
			float coord = Float.parseFloat(tokens[5]);
			boolean strand = tokens[2].equals("(+)");
			if(chrs.contains(chr)){
				chrs.get(chrs.indexOf(chr)).addGene(tokens[0], coord, strand);
			}else{
				chr.addGene(tokens[0], coord, strand);
				chrs.add(chr);
				
			}
			if(geneMap.get(tokens[0]) != null){
				chromap.put(geneMap.get(tokens[0]), chrs.get(chrs.indexOf(chr)));
			}
			line = br.readLine();
		}
		return chrs;
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
			/*if(nt - 2 < minSize) {
				line = br.readLine();
				continue;
			}*/
			int numChild = Integer.parseInt(tokens[1].split(":")[1]);
			HashSet<Integer> gidx = new HashSet<Integer>();
			for(int j = 2; j < nt; j++){
				String[] t2 = tokens[j].split(":");
				int i = rowmap.get(t2[0]);
				gidx.add(i);
			}
			allAttractors.add(new Attractor(name, gidx, nt-2, numChild));
			line = br.readLine();
		}
		br.close();
		return allAttractors;
	}
	
	private static float[] getMetaGene(float[][] data, HashSet<Integer> idx, int n){
		int m = idx.size();
		float[] out = new float[n];
		for(int j = 0; j < n; j++){
			for(Integer i : idx){
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
		String path = "/home/weiyi/workspace/javaworks/caf/output/caf/ov.gse9891.rownorm.z8/";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		int minSize = 10;
		float zScore = 8f;
		float ovlpTh = 0.5f; // overlap threshold
		boolean rowNormalization = true;
		boolean MINormalization = false;
		boolean CNV = false;
		
		System.out.println("Loading files...");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/brca/gse2034/ge.13271x286.var.txt");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/brca/tcga/ge/ge.17814x536.knn.txt");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/coad/gse14333/ge.20765x290.var.txt");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/coad/tcga/ge/ge.17814x154.knn.txt");
		DataFile ma = DataFile.parse("/home/weiyi/workspace/data/ov/gse9891/ge.20765x285.var.txt");
		//DataFile ma = DataFile.parse("/home/weiyi/workspace/data/ov/tcga/ge/ge.12042x582.txt");
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		if(rowNormalization) ma.normalizeRows();
		float[][] data = ma.getData();
		
		HashMap<String, Integer> rowmap = ma.getRows();
		HashMap<Integer, Chromosome> chromap = new HashMap<Integer, Chromosome>();
		ArrayList<Chromosome> chrs = new ArrayList<Chromosome>();
		
		if(CNV){
				chrs = parseChromGenes("/home/weiyi/workspace/data/annot/affy/u133p2/gene.location3", 
				ma.getProbes(), ma.getRows(), chromap);
		}
		
		ArrayList<String> probeNames = ma.getProbes();
		Attractor.setGeneNames(probeNames);
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
				if(ovlpRatio >= ovlpTh){
					allDist.add(new DistPair(ai.name, aj.name, ovlpRatio));
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
			System.out.print(target.toString());
			
			int xIdx = allGeneSet.indexOf(new Attractor(target.x));
			int yIdx = allGeneSet.indexOf(new Attractor(target.y));
			if(xIdx > yIdx){
				int tmp = yIdx;
				yIdx = xIdx;
				xIdx = tmp;
			}
			
		// 1. merge gene set, put in new key
			allGeneSet.get(xIdx).merge(allGeneSet.get(yIdx));
			String mergedName = allGeneSet.get(xIdx).name;
			System.out.println(" into " + mergedName + " (" + xIdx + ").");
			
			
			
		// 2. remove merged gene set, remove keys
			allGeneSet.remove(yIdx);
			
		// 3. remove distances regarding these two gene set
			for(int i = allDist.size()-1; i >= 0; i--){
				if(allDist.get(i).contains(target.x) || allDist.get(i).contains(target.y)){
					allDist.remove(i);
				}
			}
		// 4. add new distances
			for(int i = 0; i < allGeneSet.size(); i++){
				if(i != xIdx){
					Attractor ai = allGeneSet.get(i);
					float ovlpRatio = ai.similarities(allGeneSet.get(xIdx));
					if(ovlpRatio >= ovlpTh){
						allDist.add(new DistPair(ai.name, mergedName, ovlpRatio));
					}
				}
			}
			Collections.sort(allDist);
		}
		N = allGeneSet.size();
		
		Collections.sort(allGeneSet);
		PrintWriter pw2 = new PrintWriter(new FileWriter(path + "MetaAttractors.gwt"));
		for(int i = N-1; i >= 0; i-- ){
			Attractor a = allGeneSet.get(i);
			if(a.sz < minSize){
				allGeneSet.remove(i);
			}else{
				pw2.println(a);
			}
		}
		pw2.close();
		
		N = allGeneSet.size();
		System.out.println(N + " meta-attractors left.");
		
		System.out.println("Reconverging...");
		
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
			/*if(convergedGeneSets.size() >= 10){
				break;
			}*/
		}
		
		PrintWriter pw = new PrintWriter(new FileWriter(path + "ConvergedMetaAttractors.gwt"));
		GeneSet.setProbeNames(probeNames);
		for(GeneSet gs : convergedGeneSets){
			pw.print(gs.getName() + "\t" + gs.size() + ":" + gs.getNumChild() + "\t");
			pw.println(gs.toProbes());
		}
		pw.close();
	}

}
