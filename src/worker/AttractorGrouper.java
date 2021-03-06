package worker;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;


import obj.Chromosome;
import obj.DataFile;
import obj.GeneSet;
import obj.ValIdx;

public class AttractorGrouper extends DistributedWorker {
	static float ovlpTh = 0.5f;
	
	
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
		
		void intersect(Attractor a){
			this.name = "MetaAttractor" + counter;
			counter++;
			HashSet<Integer> newGeneIdx = new HashSet<Integer>();
			for(Integer i : a.geneIdx){
				if(this.geneIdx.contains(i)){
					newGeneIdx.add(i);
				}
			}
			
			this.geneIdx = newGeneIdx;
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
	
	public AttractorGrouper(int segment, int numSegments, long jobID) {
		id = segment;
		totalComputers = numSegments;
		AttractorGrouper.jobID = jobID;
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
			if(nt-2 < minSize){
				line = br.readLine();
				continue;
			}
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
	private static ArrayList<GeneSet> parseGeneSetInOneFile(
			String file, HashMap<String, Integer> rowmap) throws IOException {
		
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = br.readLine();
		ArrayList<GeneSet> allAttractors = new ArrayList<GeneSet>();
		while(line != null){ // each line is an attractor
			String[] tokens = line.split("\t");
			String name = tokens[0];
			int nt = tokens.length;
			int numChild = Integer.parseInt(tokens[1].split(":")[1]);
			ArrayList<ValIdx> gidx = new ArrayList<ValIdx>();
			for(int j = 2; j < nt; j++){
				String[] t2 = tokens[j].split(":");
				int i = rowmap.get(t2[0]);
				float f = Float.parseFloat(t2[1]);
				gidx.add(new ValIdx(i, f));
			}
			allAttractors.add(new GeneSet(name, gidx.toArray(new ValIdx[0]), numChild));
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
	
	public void mergeAndReconverge(String path, DataFile ma, Converger cvg, int minSize) throws Exception{
		ArrayList<String> probeNames = ma.getProbes();
		Attractor.setGeneNames(probeNames);
		System.out.print("Loading gene sets...");
		
		ArrayList<Attractor> allGeneSet = parseAttractorInOneFile(path + "attractors.gwt", ma.getRows(), minSize);
		int N = allGeneSet.size();
		
		System.out.println(N + " gene sets are loaded.");
		Collections.sort(allGeneSet);
		
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
		
		while(allDist.size() > 0){
			DistPair target = allDist.get(0); // merge the most similar pair
			//System.out.print(target.toString());
			
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
			//System.out.println(" into " + mergedName + " (" + xIdx + ").");
			
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
		
		if(this.id == 0){
			new File("output").mkdir();
			new File("output/" + jobID).mkdir();
			PrintWriter pw2 = new PrintWriter(new FileWriter("output/" + jobID + "/MetaAttractors.gwt"));
			for(int i = N-1; i >= 0; i-- ){
				Attractor a = allGeneSet.get(i);
				pw2.println(a);
			}
			pw2.close();
		}
		
		N = allGeneSet.size();
		System.out.println(N + " meta-attractors left.");
		
		System.out.println("Reconverging...");
		
		
		ArrayList<GeneSet> convergedGeneSets = new ArrayList<GeneSet>();
		float[][] data = ma.getData();
		int start = N * id / totalComputers;
		int end = N * (id + 1) / totalComputers;
		
		for(int i = 0; i < N; i++){
			if(i >= start && i < end){
				Attractor a = allGeneSet.get(i);
				System.out.println(i + "\t" + a.getName() + "...");
				float[] vec = getMetaGene(data, a.geneIdx, ma.getNumCols());
				ArrayList<ValIdx> metaIdx = cvg.findAttractor(data, vec);
				GeneSet gs;
				gs = new GeneSet(a.getName(), metaIdx.toArray(new ValIdx[0]), a.strength);
				convergedGeneSets.add(gs);
			}
		}
		prepare("geneset");
		PrintWriter pw = new PrintWriter(new FileWriter("tmp/" + jobID + "/geneset/caf." + String.format("%05d", id)+".txt"));
		GeneSet.setProbeNames(probeNames);
		for(GeneSet gs : convergedGeneSets){
			pw.print(gs.getName() + "\t" + gs.size() + ":" + gs.getNumChild() + "\t");
			pw.println(gs.toProbes());
		}
		pw.close();
	}
	
	public void outputConvergedAttractors(String path, DataFile ma, int minSize) throws Exception{
		ArrayList<GeneSet> allAttractors = new ArrayList<GeneSet>();
		String[] files = new File(path).list();
		Arrays.sort(files);
		
		for(String s : files){
			allAttractors.addAll(parseGeneSetInOneFile(path + s, ma.getRows()));
		}
		int N = allAttractors.size();
		
		System.out.println(N + " attractors loaded.");
		System.out.println("Examine for overlapping...");
		Collections.sort(allAttractors);
		for(int i = N-1; i >= 0; i--){
			GeneSet gs = allAttractors.get(i);
			if(gs.size() < minSize){
				allAttractors.remove(i);
				continue;
			}
			for(int j = 0; j < i; j++){
				GeneSet gs2 = allAttractors.get(j);
				
				int th = (int) (Math.min(gs2.size(), gs.size()) * ovlpTh);
				int k = gs2.overlapWith(gs);
				if(k > th){
					allAttractors.remove(i);
					break;
				}
			}
		}
		
		N = allAttractors.size();
		System.out.println(N + " attractors left.");
		
		new File("output").mkdir();
		new File("output/" + jobID).mkdir();
		PrintWriter pw = new PrintWriter(new FileWriter("output/" + jobID + "/ConvergedAttractors.gwt"));
		GeneSet.setProbeNames(ma.getProbes());
		for(GeneSet gs : allAttractors){
			pw.print(gs.getName() + "\t" + gs.size() + ":" + gs.getNumChild() + "\t");
			pw.println(gs.toProbes());
		}
		pw.close();
	}
	public static void setOvlpTh(float ovlpTh){
		AttractorGrouper.ovlpTh = ovlpTh;
	}
	
}
