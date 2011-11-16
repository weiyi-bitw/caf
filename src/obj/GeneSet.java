package obj;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

public class GeneSet implements Comparable<GeneSet>{
	static ArrayList<String> probeNames;
	static Annotations annot;
	static int mergeThreshold; 
	
	HashSet<Integer> geneIdx;
	String seed;
	int sz;
	int minIdx; // used as an temporary id
	
	public GeneSet(HashSet<Integer> idx){
		this.geneIdx = idx;
		this.sz = geneIdx.size();
		mergeThreshold = sz/2;
	}
	
	public static void setProbeNames(ArrayList<String> probeNames){
		GeneSet.probeNames = probeNames;
	}
	public static void setAnnotations(Annotations annot){
		GeneSet.annot = annot;
	}
	public static void setMergeThreshold(int th){
		GeneSet.mergeThreshold = th;
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
		HashSet<Integer> newGeneIdx = new HashSet<Integer>();
		HashSet<Integer> otherGeneIdx = other.geneIdx;
		int cnt = 0;
		for(Integer i : this.geneIdx){
			if(otherGeneIdx.contains(i)){
				cnt ++;
				newGeneIdx.add(i);
			}
		}
		if(cnt < mergeThreshold){
			return false;
		}else{
			this.geneIdx = newGeneIdx;
			this.sz = newGeneIdx.size();
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
			ArrayList<Integer> idx = new ArrayList<Integer>();
			for(Integer i : geneIdx){
				idx.add(i);
			}
			Collections.sort(idx);
			String s = "";
			boolean first = true;
			for(Integer i : idx){
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