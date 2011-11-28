package obj;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;



public class GeneSet implements Comparable<GeneSet>{
	static class ValIdx implements Comparable<ValIdx>{
		float val;
		int idx;
		ValIdx(int i, float v){
			this.idx = i;
			this.val = v;
		}
		
		public int compareTo(ValIdx other) {
			return -Double.compare(this.val, other.val);
		}
	}
	static ArrayList<String> probeNames;
	static Annotations annot;
	static int lowestMergeFold = 2; 
	
	HashSet<Integer> attractees;
	HashMap<Integer, Float> weightMap; 
	ArrayList<String> geneNames;
	HashMap<String, Float> geneWeightMap;
	int[] geneIdx;
	String seed;
	int sz;
	int minIdx; // used as an temporary id
	int numChild;
	
	public GeneSet(int[] idx){
		Arrays.sort(idx);
		this.geneIdx = idx;
		this.sz = geneIdx.length;
		weightMap = new HashMap<Integer, Float>();
		for(Integer i : idx){
			weightMap.put(i, 1f);
		}
		numChild = 1;
	}
	
	public GeneSet(HashSet<Integer> attractees, int[] idx){
		Arrays.sort(idx);
		this.geneIdx = idx;
		this.sz = geneIdx.length;
		this.attractees = attractees;
		weightMap = new HashMap<Integer, Float>();
		for(Integer i : idx){
			weightMap.put(i, 1f);
		}
		numChild = 1;
	}
	
	public GeneSet(HashSet<Integer> attractees, int[] idx, float[] wts, int numChild){
		Arrays.sort(idx);
		this.geneIdx = idx;
		this.sz = geneIdx.length;
		this.attractees = attractees;
		weightMap = new HashMap<Integer, Float>();
		for(int i = 0; i < sz; i++){
			weightMap.put(geneIdx[i], wts[i]);
		}
		this.numChild = numChild;
	}
	
	
	public static void setProbeNames(ArrayList<String> probeNames){
		GeneSet.probeNames = probeNames;
	}
	public static void setAnnotations(Annotations annot){
		GeneSet.annot = annot;
	}
	public static void setMergeThreshold(int th){
		GeneSet.lowestMergeFold = th;
	}
	public int compareTo(GeneSet other) {
		return Double.compare(this.minIdx, other.minIdx);
	}
	public boolean equals(Object other){
		boolean result = false;
        if (other instanceof GeneSet) {
        	GeneSet that = (GeneSet) other;
        	if(this.sz != that.sz){
        		return false;
        	}else{
        		result = Arrays.equals(this.geneIdx, that.geneIdx);
        	}
        }
        return result;
	}
	
	/*// Merging the gene sets using their intersection	
	public boolean merge(GeneSet other){
		ArrayList<Integer> newGeneIdx = new ArrayList<Integer>();
		int[] otherGeneIdx = other.geneIdx;
		int cnt = 0;
		for(int i : this.geneIdx){
			if(Arrays.binarySearch(otherGeneIdx, i) >= 0){
				cnt ++;
				newGeneIdx.add(i);
			}
		}
		int mergeThreshold = Math.min(this.sz, other.sz)/lowestMergeFold;
		if(cnt < mergeThreshold || cnt <= 1){
			return false;
		}else{
			int[] ngIdx = new int[cnt];
			for(int i = 0; i < cnt; i++){
				ngIdx[i] = newGeneIdx.get(i);
			}
			for(int j : other.attractees){
				this.attractees.add(j);
			}
			this.geneIdx = ngIdx;
			this.sz = cnt;
			return true;
		}
	}*/
	
	// Merging the gene sets using their union	
		public boolean merge(GeneSet other){
			HashSet<Integer> newGeneIdx = new HashSet<Integer>();
			HashMap<Integer, Float> newWeightMap = new HashMap<Integer, Float>();
			int[] otherGeneIdx = other.geneIdx;
			int cnt = 0;
			for(int i : this.geneIdx){
				newGeneIdx.add(i);
				float w = weightMap.get(i);
				if(other.weightMap.get(i) != null){
					cnt++;
					w += other.weightMap.get(i);
				}
				newWeightMap.put(i, w);
			}
			for(int i : otherGeneIdx){
				newGeneIdx.add(i);
				float w = other.weightMap.get(i);
				if(weightMap.get(i)==null){
					newWeightMap.put(i, w);
				}
			}
			int mergeThreshold = Math.min(this.sz, other.sz)/lowestMergeFold;
			if(cnt <= mergeThreshold || cnt <= 1){
				// did not pass merge threshold, skip
				return false;
			}else{
				int[] ngIdx = new int[newGeneIdx.size()];
				cnt = 0;
				for(Integer i : newGeneIdx){
					ngIdx[cnt] = i;
					cnt++;
				}
				for(Integer j : other.attractees){
					this.attractees.add(j);
				}
				this.geneIdx = ngIdx;
				this.sz = cnt;
				this.numChild += other.numChild;
				this.weightMap = newWeightMap;
				return true;
			}
		}
	
	public String toString(){ // to a string of attractees, numChild, attractors, all in index
		String s = "";
		
		boolean first = true;
		for(Integer i : attractees){
			if(first){
				s = s + i;
				first = false;
			}else{
				s = s + "," + i;
			}
		}
		
		s = s + "\t" + numChild;
		for(int i : geneIdx){
				s = s + "\t" + i + "," + weightMap.get(i);
		}
		return s;
	}
	
	public int size(){
		return sz;
	}
	
	public String toGenes(){
		if(annot == null){
			String s = "";
			boolean first = true;
			for(int i : geneIdx){
				if(first){
					s = probeNames.get(i);
					first = false;
				}else{
					s = s + "\t" + probeNames.get(i);
				}
			}
			return s;
		
		}else{
			ArrayList<String> geneNames = new ArrayList<String>();
			String s;
			for(int i : geneIdx){
				s = probeNames.get(i);
				if(!geneNames.contains(s)){
					geneNames.add(s);
				}
			}
			Collections.sort(geneNames);
			s = "";
			boolean first = true;
			for(String ss : geneNames){
				if(first){
					s = ss;
					first = false;
				}else{
					s = s + "\t" + ss;
				}
			}
			return s;
		}
	}
	
	public String getAttractees(){
		if(annot == null){
			String s = "";
			boolean first = true;
			for(Integer i : attractees){
				if(first){
					s = probeNames.get(i);
					first = false;
				}else{
					s = s + "\t" + probeNames.get(i);
				}
			}
			return s;
		
		}else{
			ArrayList<String> geneNames = new ArrayList<String>();
			String s;
			for(Integer i : attractees){
				s = annot.getGene(probeNames.get(i));
				if(!geneNames.contains(s)){
					geneNames.add(s);
				}
			}
			Collections.sort(geneNames);
			s = "";
			boolean first = true;
			for(String ss : geneNames){
				if(first){
					s = ss;
					first = false;
				}else{
					s = s + "\t" + ss;
				}
			}
			return s;
		}
	}
	
	public void sort(){ // sort the indices according to their weights decreasingly
		ValIdx[] vec = new ValIdx[sz];
		for(int i = 0; i < sz; i++){
			vec[i] = new ValIdx(geneIdx[i], weightMap.get(geneIdx[i]));
		}
		Arrays.sort(vec);
		for(int i = 0; i < sz; i++){
			geneIdx[i] = vec[i].idx;
		}
	}
	
	public String getWeight(){
		geneNames = new ArrayList<String>();
		geneWeightMap = new HashMap<String, Float>();
		if(annot == null){
			String s = "";
			boolean first = true;
			for(Integer i : geneIdx){
				geneNames.add(probeNames.get(i));
				float w = weightMap.get(i)/numChild;
				geneWeightMap.put(probeNames.get(i), w);
				if(first){
					s = s+w;
					first = false;
				}else{
					s = s + "\t" + w;
				}
			}
			return s;
		
		}else{
			String s;
			for(Integer i : geneIdx){
				s = annot.getGene(probeNames.get(i));
				if(!geneNames.contains(s)){
					geneNames.add(s);
					geneWeightMap.put(s, weightMap.get(i));
				}else{
					float w = weightMap.get(s);
					if(w < weightMap.get(i)){
						geneWeightMap.put(s, weightMap.get(i));
					}
				}
			}
			
			s = "";
			boolean first = true;
			for(String ss : geneNames){
				if(first){
					s = s+(geneWeightMap.get(ss)/numChild);
					first = false;
				}else{
					s = s + "\t" + (geneWeightMap.get(ss)/numChild);
				}
			}
			return s;
		}
	}
	public ArrayList<String> getGeneNames(){
		return geneNames;
	}
	public HashMap<String, Float> getGeneWeightMap(){
		return geneWeightMap;
	}
	public int getAttracteeSize(){
		return attractees.size();
	}
}