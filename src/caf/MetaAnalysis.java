package caf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import util.StatOps;

public class MetaAnalysis {
	static class Gene implements Comparable<Gene>{
		String name;
		String chr;
		float idx;
		
		Gene(String name){
			this.name = name;
		}
		Gene(String name, String chr, float idx){
			this.name = name;
			this.chr = chr;
			this.idx = idx;
		}
		public int hashCode(){
			return name.hashCode();
		}
		public boolean equals(Object other){
			boolean result = false;
	        if (other instanceof Gene) {
	        	Gene that = (Gene) other;
	            result = (this.name.equals(that.name));
	        }
	        return result;
		}
		public int compareTo(Gene other) {
			return Double.compare(this.idx, other.idx);
		}
		public void resetIdx(float i){
			this.idx = i;
		}
	}
	
	static class Genome{
		ArrayList<Gene> genes;
		HashMap<String, String> chrMap;
		HashMap<String, Integer> idxMap;
		
		static Genome parseGeneLocation(String file)throws Exception{
			BufferedReader br = new BufferedReader(new FileReader(file));
			ArrayList<Gene> genes = new ArrayList<Gene>();
			br.readLine(); // first row header
			String line = br.readLine();
			int lncnt = 0;
			while(line != null){
				String[] tokens = line.split("\t");
				String name = tokens[0];
				String chr = tokens[6];
				genes.add(new Gene(name, chr, lncnt));
				lncnt++;
				line = br.readLine();
			}
			br.close();
			
			System.out.println("Genome contains " + lncnt + " genes.");
			
			return new Genome(genes);
		}
		
		Genome(ArrayList<Gene> genes){
			this.genes = genes;
			this.chrMap = new HashMap<String, String>();
			this.idxMap = new HashMap<String, Integer>();
			for(Gene g : genes){
				chrMap.put(g.name, g.chr);
				idxMap.put(g.name, (int) g.idx);
			}
		}
		
		ArrayList<Gene> getAllGenesInChr(String chr){
			ArrayList<Gene> out = new ArrayList<Gene>();
			for(Gene g: genes){
				if (g.chr.equalsIgnoreCase(chr)){
					out.add(g);
				}
			}
			return out;
		}
		
		int getIdx(String gene){
			return idxMap.get(gene);
		}
		
		String getChr(String gene){
			return chrMap.get(gene);
		}
		
		boolean contains(String gene){
			return genes.contains(new Gene(gene));
		}
	}
	
	
	static class Attractor implements Comparable<Attractor>{
		static Genome gn;
		static int minSize;
		static boolean CNV;
		static float CNVTh;
		String name;
		ArrayList<String> genes;
		HashMap<String, Float> zMap;
		float idxStd;
		float strength;
		int size;
		String chr;
		
		
		Attractor(String name, ArrayList<String> genes, float strength){
			this.name = name;
			this.genes = genes;
			this.size = genes.size();
			this.strength = strength;
		}
		
		Attractor(String name, ArrayList<String> genes, float idxStd, float strength){
			this.name = name;
			this.genes = genes;
			this.idxStd = idxStd;
			this.size = genes.size();
			this.strength = strength;
		}
		Attractor(String name, ArrayList<String> genes, HashMap<String, Float> zMap, float idxStd, float strength, String chr){
			this.name = name;
			this.genes = genes;
			this.idxStd = idxStd;
			this.size = genes.size();
			this.strength = strength;
			this.chr = chr;
		}
		static Attractor parseAttractor(String line){
			HashMap<String, Float> zMap = new HashMap<String, Float>();
			String[] tokens = line.split("\t");
			int nt = tokens.length;
			float strength = Float.parseFloat(tokens[1].split(":")[1]);
			ArrayList<String> genes = new ArrayList<String>();
			String chr = "NA";
			float[] indices = new float[nt-2];
			for(int i = 2; i < nt; i++){
				String[] t2 = tokens[i].split(":");
				genes.add(t2[0]);
				zMap.put(t2[0], Float.parseFloat(t2[1]));
				if(i==2){
					if(gn.contains(t2[0])){
						chr = gn.getChr(t2[0]);
					}
				}
				if(gn.contains(t2[0])){
					indices[i-2] = (float)gn.getIdx(t2[0]);
				}else{
					indices[i-2] = Float.NaN;
				}
			}
			float idxStd = StatOps.std(indices);
			
			return new Attractor(tokens[0], genes, zMap,idxStd, strength, chr);
		}
		ArrayList<String> getOvlp(Attractor other){
			ArrayList<String> ovlp = new ArrayList<String>();
			ArrayList<String> genesOther = other.genes;
			for(String g : this.genes){
				if(genesOther.contains(g)){
					ovlp.add(g);
				}
			}
			return ovlp;
		}
		float getZ(String s){
			return zMap.get(s);
		}
		public int compareTo(Attractor other) {
			return -Double.compare(this.strength, other.strength);
		}
		static public void setGenome(Genome gn){
			Attractor.gn = gn;
		}
		static public void setMinSize(int minSize){
			Attractor.minSize = minSize;
		}
		static public void CNV(boolean CNV){
			Attractor.CNV = CNV;
		}
		static public void setCNVTh(float CNVTh){
			Attractor.CNVTh = CNVTh;
		}
	}
	static ArrayList<Attractor> parseAttractors(String filename) throws Exception{
		ArrayList<Attractor> out = new ArrayList<Attractor>();
		BufferedReader br = new BufferedReader(new FileReader(filename));
		br.readLine();
		String line = br.readLine();
		while(line != null){
			out.add(Attractor.parseAttractor(line));
			line = br.readLine();
		}
		br.close();
		Collections.sort(out);
		return out;
	}
	static String[] parseTarget(String filename) throws Exception{
		ArrayList<String> targets = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String line = br.readLine();
		while(line != null){
			targets.add(line);
			line = br.readLine();
		}
		br.close();
		return (targets.toArray(new String[0]));
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{
		boolean CNV = true;
		float CNVTh = 150;
		int minSize = 10;
		
		String path = "/home/weiyi/workspace/javaworks/caf/output/cnv/metaana/";
		Genome gn = Genome.parseGeneLocation("/home/weiyi/workspace/data/annot/affy/u133p2/gene.location3");
		
		String[] targets = parseTarget(path + "targets.txt");
		int nf = targets.length;
		
		System.out.println(nf + " target files are specified.");
		
		Attractor.setMinSize(minSize);
		Attractor.setCNVTh(CNVTh);
		Attractor.setGenome(gn);
		Attractor.CNV(CNV);
		
		ArrayList<ArrayList<Attractor>> results = new ArrayList<ArrayList<Attractor>>(nf);
		for(String t : targets){
			results.add(parseAttractors(t));
		}
		ArrayList<Attractor> origin = results.get(0);
		ArrayList<Attractor> encounter = results.get(1);
		PrintWriter pw = new PrintWriter(new FileWriter(path + "matchedAttractors.txt"));
		for(int i = 0; i < origin.size(); i++){
			Attractor a = origin.get(i);
			ArrayList<String> best = new ArrayList<String>();
			Attractor bestb = null;
			for(Attractor b: encounter){
				ArrayList<String> c = a.getOvlp(b);
				if(c.size() > best.size()){
					best = c;
					bestb = b;
				}
			}
			if(best.size() > 0){
				encounter.remove(bestb);
				pw.print(a.chr);
				for(String g : best){
					pw.print("\t" + g);
				}
				pw.println();
			}
		}
		pw.close();
		
		System.out.println("Done.");
	}

}
