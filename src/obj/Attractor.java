package obj;

import java.util.ArrayList;
import java.util.HashMap;

public class Attractor implements Comparable<Attractor>{
		static Genome gn;
		static int minSize;
		static boolean CNV;
		static float CNVTh;
		String name;
		ArrayList<String> genes;
		HashMap<String, Float> zMap;
		HashMap<String, Integer> idxMap;
		float strength;
		float range;
		int size;
		String chr;
		
		
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
			float strength = CNV? 1/(range) : Float.parseFloat(tokens[1].split(":")[1]);
			return new Attractor(tokens[0], genes, zMap,idxMap, range, strength, chr);
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
			String s = name + "\t" + chr + "\t" + strength + "\t" + range;
			for(String ss : genes){
				s = s + "\t" + ss + "(" + idxMap.get(ss) + ")" + ":" + zMap.get(ss);
			}
			return s;
		}
		public String toString(){
			String s = name + "\t" + chr + "\t" + strength + "\t" + range;
			for(String ss : genes){
				s = s + "\t" + ss ;
			}
			return s;
		}
}
