package obj;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

public class GeneSet implements Comparable<GeneSet>{
	static ArrayList<String> probeNames;
	static Annotations annot;
	static int lowestMergeFold; 
	
	int[] geneIdx;
	String seed;
	int sz;
	int minIdx; // used as an temporary id
	
	public GeneSet(int[] idx){
		Arrays.sort(idx);
		this.geneIdx = idx;
		this.sz = geneIdx.length;
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
        		result = this.geneIdx.equals(that.geneIdx);
        	}
        }
        return result;
	}
	
	// Merging the gene sets using their intersection	
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
		if(cnt < mergeThreshold || cnt == 0){
			return false;
		}else{
			int[] ngIdx = new int[cnt];
			for(int i = 0; i < cnt; i++){
				ngIdx[i] = newGeneIdx.get(i);
			}
			this.geneIdx = ngIdx;
			this.sz = cnt;
			return true;
		}
	}
	
	public String toString(){
		String s = "";
		boolean first = true;
		for(Integer i : geneIdx){
			if(first){
				s = i.toString();
				first = false;
			}else{
				s = s + "\t" + i.toString();
			}
		}
		return s;
	}
	
	public String toGenes(){
		if(annot == null){
			String s = "";
			boolean first = true;
			for(Integer i : geneIdx){
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
			for(Integer i : geneIdx){
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
	
	
}