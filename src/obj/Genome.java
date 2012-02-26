package obj;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

public class Genome {
	
	static class Gene implements Comparable<Gene>{
		static boolean useChrCoord = false;
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
		static void useChrCoord(boolean ucc){
			Gene.useChrCoord = ucc;
		}
	}
	
	ArrayList<Gene> genes;
	HashMap<String, String> chrMap;
	HashMap<String, Integer> idxMap;
	
	public static Genome parseGeneLocation(String file, boolean useChrCoord)throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(file));
		ArrayList<Gene> genes = new ArrayList<Gene>();
		br.readLine(); // first row header
		String line = br.readLine();
		int lncnt = 0;
		while(line != null){
			String[] tokens = line.split("\t");
			String name = tokens[0];
			String chr = tokens[6];
			float chrCoord = useChrCoord? Float.parseFloat(tokens[5]) : lncnt;
			genes.add(new Gene(name, chr, chrCoord));
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
