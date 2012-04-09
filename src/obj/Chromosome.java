package obj;

import java.util.ArrayList;
import java.util.HashMap;

import obj.ValIdx;

public class Chromosome{
	String name;
	ArrayList<ValIdx> geneIdx;
	static ArrayList<String> genes;
	static HashMap<String, Integer> geneMap;
	ArrayList<Boolean> plusStrand;
	
	public Chromosome(String name){
		this.name = name;
		geneIdx = new ArrayList<ValIdx>();
		plusStrand = new ArrayList<Boolean>();
	}
	
	public int hashCode(){
		return name.hashCode();
	}
	public boolean equals(Object other){
		boolean result = false;
        if (other instanceof Chromosome) {
        	Chromosome that = (Chromosome) other;
            result = (this.name.equals(that.name));
        }
        return result;
	}
	
	static public void setGeneNames(ArrayList<String> genes){
		Chromosome.genes = genes;
	}
	
	static public void setGeneMap(HashMap<String, Integer> geneMap){
		Chromosome.geneMap = geneMap;
	}
	
	public void addGene(String g, float coord, boolean plus){
		if(geneMap.get(g) != null){
			int idx = geneMap.get(g);
			geneIdx.add(new ValIdx(idx, coord));
			plusStrand.add(plus);
		}
	}
	public int size(){
		return geneIdx.size();
	}
	public ArrayList<ValIdx> geneIdx(){
		return geneIdx;
	}
	public String name(){
		return name;
	}
}
