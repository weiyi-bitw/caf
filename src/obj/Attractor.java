package obj;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class Attractor implements Comparable<Attractor>{
	static class GeneIdxVal extends ValIdx{
		String g;
		public GeneIdxVal(String g, int i, float v) {
			super(i, v);
			this.g = g;
		}
		
	}
		static int counter = 0;
		static Genome gn;
		static int minSize;
		static boolean CNV;
		static float CNVTh;
		String name;
		ArrayList<String> genes;
		HashMap<String, Float> zMap;
		HashMap<String, Integer> idxMap;
		float strength;
		float z; // representative z-score
		float range = Float.NaN;
		int size;
		String chr = "NA";
		
		public Attractor(String name){
			this.name = name;
		}
		
		public Attractor(String name, ArrayList<String> genes, HashMap<String, Float> zMap, HashMap<String, Integer> idxMap, float strength){
			this.name = name;
			this.genes = genes;
			this.zMap = zMap;
			this.idxMap = idxMap;
			this.size = genes.size();
			this.strength = strength;
		}
		public Attractor(String name, ArrayList<String> genes, HashMap<String, Float> zMap, HashMap<String, Integer> idxMap, float strength, float z){
			this.name = name;
			this.genes = genes;
			this.z = z;
			this.zMap = zMap;
			this.idxMap = idxMap;
			this.size = genes.size();
			this.strength = strength;
		}
		Attractor(String name, ArrayList<String> genes, HashMap<String, Float> zMap, HashMap<String, Integer> idxMap, float range, float strength, String chr){
			this.name = name;
			this.genes = genes;
			this.zMap = zMap;
			this.idxMap = idxMap;
			this.size = genes.size();
			this.range = range;
			this.strength = strength;
			this.chr = chr;
		}
		public static Attractor parseAttractor(String line){
			HashMap<String, Float> zMap = new HashMap<String, Float>();
			HashMap<String, Integer> idxMap = new HashMap<String, Integer>();
			String[] tokens = line.split("\t");
			int nt = tokens.length;
			
			ArrayList<String> genes = new ArrayList<String>();
			String chr = "NA";
			float idxMin = Float.MAX_VALUE;
			float idxMax = -1;
			for(int i = 2; i < nt; i++){
				String[] t2 = tokens[i].split(":");
				genes.add(t2[0]);
				zMap.put(t2[0], Float.parseFloat(t2[1]));
				if(i==2){
					if(CNV && gn.contains(t2[0])){
						chr = gn.getChr(t2[0]);
					}
				}
				if(CNV && gn.contains(t2[0])){
					float idx = gn.getIdx(t2[0]);
					idxMap.put(t2[0], (int) idx);
					if(idx < idxMin){
						idxMin = idx;
					}
					if(idx > idxMax){
						idxMax = idx;
					}
				}
			}
			float range = CNV? idxMax - idxMin : Float.NaN;
			float strength = CNV? 1/range : Float.parseFloat(tokens[1].split(":")[1]);
			return new Attractor(tokens[0], genes, zMap,idxMap, range, strength, chr);
		}
		
		public int compareTo(Attractor other) {
			return -Double.compare(this.strength, other.strength);
		}
		public int hashCode(){
			return name.hashCode();
		}
		public boolean equals(Object other){
			boolean result = false;
	        if (other instanceof Attractor) {
	        	Attractor that = (Attractor) other;
	            result = (this.name.equals(that.name));
	        }
	        return result;
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
		public void aggregatingMode(){
			for(String g : genes){
				this.zMap.put(g, strength/range);
			}
		}
		public void aggregate(Attractor a){ // merge a and add the strength to the genes
			for(String s : a.genes){
				if(!genes.contains(s)){
					genes.add(s);
					zMap.put(s, a.zMap.get(s));
					idxMap.put(s, a.idxMap.get(s));
				}else{
					zMap.put(s, this.zMap.get(s) + a.zMap.get(s));
				}
			}
			this.size = genes.size();
			this.strength += a.strength;
		}
		public void merge(Attractor a){
			this.name = "MetaAttractor" + counter;
			counter++;
			/*if(size < a.size){
				zMap = a.zMap;
			}*/
			for(String s : a.genes){
				if(!genes.contains(s)){
					genes.add(s);
					zMap.put(s, a.zMap.get(s));
				}
			}
			this.size = genes.size();
			this.strength += a.strength;
		}
		public void intersect(Attractor a){
			this.name = "MetaAttractor" + counter;
			counter++;
			HashSet<String> genesHs = new HashSet<String>();
			ArrayList<String> newGenes = new ArrayList<String>();
			genesHs.addAll(this.genes);
			for(String s : a.genes){
				if(genesHs.contains(s)){
					newGenes.add(s);
				}
			}
			
			this.genes = newGenes;
			this.size = newGenes.size();
			if(size > a.size) {this.zMap = a.zMap;}
			this.strength += a.strength;
		}
		public void fight(Attractor a){
			if(a.strength > this.strength){
				this.name = "MetaAttractor" + counter;
				counter ++;
				this.genes = a.genes;
				this.size = a.size;
				this.zMap = a.zMap;
				this.z = a.z;
			}
			this.strength += a.strength;
			
		}
		
		public void reset(ArrayList<ValIdx> metaIdx, ArrayList<String> probeNames){
			ArrayList<String> newGenes = new ArrayList<String>();
			HashMap<String, Float> newZMap = new HashMap<String, Float>();
			for(ValIdx vi : metaIdx){
				String g = probeNames.get(vi.idx());
				newGenes.add(g);
				newZMap.put(g, vi.val());
			}
			this.genes = newGenes;
			this.size = newGenes.size();
			this.zMap =newZMap;
		}
		public String name(){
			return this.name;
		}
		public ArrayList<String> getOvlp(Attractor other){
			ArrayList<String> ovlp = new ArrayList<String>();
			ArrayList<String> genesOther = other.genes;
			for(String g : this.genes){
				if(genesOther.contains(g)){
					ovlp.add(g);
				}
			}
			return ovlp;
		}
		public float getZ(String s){
			return zMap.get(s);
		}
		public float similarities(Attractor a){
			HashSet<String> myGenes = new HashSet<String>();
			myGenes.addAll(this.genes);
			float ovlp = 0;
			for(String s : a.genes){
				if(myGenes.contains(s)){
					ovlp++;
				}
			}
			return ovlp / (float)Math.min(this.size, a.size);
		}
		public ArrayList<String> genes(){
			return this.genes;
		}
		public ArrayList<Integer> geneIndices(){
			ArrayList<Integer> geneIdx = new ArrayList<Integer>();
			for(String s: genes){
				geneIdx.add(idxMap.get(s));
			}
			return geneIdx;
		}
		public int size(){
			return size;
		}
		public String chr(){
			return this.chr;
		}
		public float range(){
			return this.range;
		}
		public float strength(){
			return this.strength;
		}
		public String toStringInDetail(){
			ArrayList<GeneIdxVal> genesVal = new ArrayList<GeneIdxVal>();
			for(String g : genes){
				genesVal.add(new GeneIdxVal(g, idxMap.get(g), zMap.get(g)));
			}
			Collections.sort(genesVal);
			String s = name + "\t" + chr + "\t" + strength + "\t" + range;
			for(GeneIdxVal giv : genesVal){
				float z = zMap.get(giv.g);
				z = Math.round(z*1000)/(float)1000;
				s = s + "\t" + giv.g + "(" + idxMap.get(giv.g) + ")" + ":" + z;
			}
			return s;
		}
		public String toString(){
			ArrayList<GeneIdxVal> genesVal = new ArrayList<GeneIdxVal>();
			for(String g : genes){
				genesVal.add(new GeneIdxVal(g, idxMap.get(g), zMap.get(g)));
			}
			Collections.sort(genesVal);
			String s = name + "\t" + chr + "\t" + strength + "\t" + range;
			for(GeneIdxVal giv : genesVal){
				s = s + "\t" + giv.g;
			}
			return s;
		}
		public String toStringGenesOnly(){
			ArrayList<GeneIdxVal> genesVal = new ArrayList<GeneIdxVal>();
			for(String g : genes){
				/*System.out.println(g);
				System.out.println(idxMap.get(g));
				System.out.println(zMap.get(g));*/
				
				genesVal.add(new GeneIdxVal(g, idxMap.get(g), zMap.get(g)));
			}
			Collections.sort(genesVal);
			String s = name + "\t" + size + ":" + (int)strength;
			for(GeneIdxVal giv : genesVal){
				s = s + "\t" + giv.g + ":" + giv.val ;
			}
			return s;
		}
}
