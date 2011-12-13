package obj;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import worker.Converger.ValIdx;


public class GeneSet implements Comparable<GeneSet>{
	static ArrayList<String> probeNames;
	static Annotations annot;
	static int lowestMergeFold = 2; 
	
	HashSet<Integer> attractees;
	ArrayList<String> geneNames;
	ValIdx[] geneIdx;
	String seed;
	int sz;
	int minIdx; // used as an temporary id
	int numChild;
	
	public GeneSet(ValIdx[] idx){
		Arrays.sort(idx);
		this.geneIdx = idx;
		this.sz = geneIdx.length;
		numChild = 1;
	}
	
	public GeneSet(HashSet<Integer> attractees, ValIdx[] idx){
		Arrays.sort(idx);
		this.geneIdx = idx;
		this.sz = geneIdx.length;
		this.attractees = attractees;
		numChild = 1;
	}
	
	public GeneSet(HashSet<Integer> attractees, ValIdx[] idx, int numChild){
		this.geneIdx = idx;
		this.sz = geneIdx.length;
		this.attractees = attractees;
		this.numChild = numChild;
		Arrays.sort(this.geneIdx);
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
	
	// merge two gene set when they are identical
	public void merge(GeneSet other){
		for(Integer i : other.attractees){
			this.attractees.add(i);
		}
		numChild += other.numChild;
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
	
	/*// Merging the gene sets using their union	
		public boolean merge(GeneSet other){
			HashSet<Integer> newGeneIdx = new HashSet<Integer>();
			HashMap<Integer, Float> newWeightMap = new HashMap<Integer, Float>();
			ValIdx[] otherGeneIdx = other.geneIdx;
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
		}*/
	
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
		for(ValIdx vi : geneIdx){
				s = s + "\t" + vi.idx() + "," + vi.val();
		}
		return s;
	}
	
	public int size(){
		return sz;
	}
	
	public String toProbes(){
		String s = "";
		boolean first = true;
		for(ValIdx vi : geneIdx){
			if(first){
				s = probeNames.get(vi.idx()) + ":" + vi.val();
				first = false;
			}else{
				s = s + "\t" + probeNames.get(vi.idx()) + ":" + vi.val();
			}
		}
		return s;
	}
	
	public String toGenes(){
		ArrayList<String> geneNames = new ArrayList<String>();
		ArrayList<String> output = new ArrayList<String>();
		
		String s;
		for(ValIdx vi : geneIdx){
			s = annot.getGene(probeNames.get(vi.idx()));
			if(!geneNames.contains(s)){
				geneNames.add(s);
				output.add(s);
			}
		}
		s = "";
		boolean first = true;
		for(String ss : output){
			if(first){
				s = ss;
				first = false;
			}else{
				s = s + "\t" + ss;
			}
		}
		return s;
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
		Arrays.sort(geneIdx);
	}
	
	/*public void calcWeight(){
		geneNames = new ArrayList<String>();
		geneWeightMap = new HashMap<String, Float>();
		if(annot == null){
			for(Integer i : geneIdx){
				geneNames.add(probeNames.get(i));
				float w = weightMap.get(i)/numChild;
				geneWeightMap.put(probeNames.get(i), w);
			}
		}else{
			for(Integer i : geneIdx){
				String s = annot.getGene(probeNames.get(i));
				if(!geneNames.contains(s)){
					geneNames.add(s);
					geneWeightMap.put(s, weightMap.get(i));
				}else{
					float w = weightMap.get(i);
					if(w < weightMap.get(i)){
						geneWeightMap.put(s, weightMap.get(i));
					}
				}
			}
		}
	}*/
	
	/*public String getWeight(){
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
					float w = weightMap.get(i);
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
	}*/
	public ArrayList<String> getGeneNames(){
		return geneNames;
	}
	/*public HashMap<String, Float> getGeneWeightMap(){
		return geneWeightMap;
	}*/
	public int getAttracteeSize(){
		return attractees.size();
	}
	public static boolean hasAnnot(){
		return annot != null;
	}
	public ValIdx[] getGeneIdx(){
		return geneIdx;
	}
	/*public float getOneWeight(int i ){
		return weightMap.get(i)/numChild;
	}*/
	public String getOnePair(int i){
		String s = probeNames.get(geneIdx[i].idx());
		if(hasAnnot()){
			s +=  "\t" + annot.getGene(probeNames.get(geneIdx[i].idx()));
		}
		s += "\t" + geneIdx[i].val();
		
		return s;
	}
}