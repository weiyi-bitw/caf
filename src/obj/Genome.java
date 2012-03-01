package obj;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class Genome {
	static class IntPair{
		int x;
		int y;
		
		IntPair(int x){
			this.x = x;
		}
		
		IntPair(int x, int y){
			this.x = x;
			this.y = y;
		}
		
		void setX(int x){
			this.x = x;
		}
		
		void setY(int y){
			this.y = y;
		}
		
		public String toString(){
			return x + "\t" + y;
			
		}
	}
	static class Gene implements Comparable<Gene>{
		String name;
		String chr;
		float coord;
		
		Gene(String name){
			this.name = name;
		}
		Gene(String name, String chr, float coord){
			this.name = name;
			this.chr = chr;
			this.coord = coord;
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
			int k = this.chr.compareTo(other.chr);
			if( k == 0){
				return Double.compare(this.coord, other.coord);
			}else{
				return k;
			}
			
		}
	}
	
	ArrayList<Gene> genes;
	HashMap<String, String> chrMap;
	HashMap<String, Integer> idxMap;
	HashMap<String, Float> coordMap;
	HashMap<String, IntPair> chrIdxRangeMap;
	HashMap<String, IntPair> chrCoordRangeMap;
	
	public static Genome parseGeneLocation(String file)throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(file));
		ArrayList<Gene> genes = new ArrayList<Gene>();
		HashMap<String, IntPair> chrCoordRangeMap = new HashMap<String, IntPair>();
		
		br.readLine(); // first row header
		String line = br.readLine();
		int lncnt = 0;
		String preChr = null;
		int preMaxCoord = -1;
		
		while(line != null){
			String[] tokens = line.split("\t");
			String name = tokens[0];
			String chr = tokens[6];
			int startCoord = Integer.parseInt(tokens[3]);
			int endCoord = Integer.parseInt(tokens[4]);
			
			if(chrCoordRangeMap.get(chr) == null){ // new chr coming up!
				if(preChr != null){ // previous chr, push the end value
					chrCoordRangeMap.get(preChr).setY(preMaxCoord);
				}
				chrCoordRangeMap.put(chr, new IntPair(startCoord));
				preChr = chr;
			}
			
			float chrCoord = Float.parseFloat(tokens[5]);
			genes.add(new Gene(name, chr, chrCoord));
			preMaxCoord = endCoord;
			lncnt++;
			line = br.readLine();
		}
		br.close();
		chrCoordRangeMap.get(preChr).setY(preMaxCoord);
		
		System.out.println("Genome contains " + lncnt + " genes.");
		
		return new Genome(genes, chrCoordRangeMap);
	}
	
	Genome(ArrayList<Gene> genes, HashMap<String, IntPair> chrCoordRangeMap){
		this.genes = genes;
		Collections.sort(genes);
		this.chrMap = new HashMap<String, String>();
		this.idxMap = new HashMap<String, Integer>();
		this.coordMap = new HashMap<String, Float>();
		this.chrIdxRangeMap = new HashMap<String, IntPair>();
		this.chrCoordRangeMap = chrCoordRangeMap;
		int cnt = 0;
		String preChr = null;
		for(Gene g : genes){
			String chr = g.chr;
			if(chrIdxRangeMap.get(chr) == null){ // new chr coming up!
				if(preChr != null){ // previous chr, push the end value
					chrIdxRangeMap.get(preChr).setY(cnt);
				}
				chrIdxRangeMap.put(chr, new IntPair(cnt));
				preChr = chr;
				
			}
			
			chrMap.put(g.name, chr);
			idxMap.put(g.name, cnt);
			coordMap.put(g.name, g.coord);
			cnt++;
		}
		chrIdxRangeMap.get(preChr).setY(cnt);
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
	
	public float getCoord(String gene){
		return coordMap.get(gene);
	}
	
	public float getChrCoordRange(String chr){
		IntPair ip = chrCoordRangeMap.get(chr);
		return (ip.y - ip.x);
	}
	
// filtering genes to those contained in the DataFile
	public void linkToDataFile(DataFile ma){ 
		ArrayList<String> maGenes = ma.getProbes();
		int k = genes.size();
		for(int i = k-1; i >= 0; i--){
			String g = this.genes.get(i).name;
			if(!maGenes.contains(g)){
				chrMap.remove(g);
				genes.remove(i);
			}
		}
		k = genes.size();
		this.chrIdxRangeMap = new HashMap<String, IntPair>();
		idxMap = new HashMap<String, Integer>();
		String preChr = null;
		for(int i = 0; i < k; i++){
			Gene g = genes.get(i);
			String chr = g.chr;
			if(chrIdxRangeMap.get(chr) == null){ // new chr coming up!
				if(preChr != null){ // previous chr, push the end value
					chrIdxRangeMap.get(preChr).setY(i);
				}
				chrIdxRangeMap.put(chr, new IntPair(i));
				preChr = chr;
				
			}
			
			chrMap.put(g.name, chr);
			idxMap.put(g.name, i);
			coordMap.put(g.name, g.coord);
		}
		chrIdxRangeMap.get(preChr).setY(k);
		
		/*for(String c : chrCoordRangeMap.keySet()){
			System.out.println(c + "\t" + chrCoordRangeMap.get(c));
		}*/
		System.out.println(k + " genes linked.");
		
	}
	
// return a gene set within the windowsize centered at the gene, -1 for the whole chromosome
	public String[] getNeighbors(String gene, int windowsize){
		int idx = genes.indexOf(new Gene(gene));
		String chr = getChr(gene);
		ArrayList<String> outGenes = new ArrayList<String>();
		IntPair ip = chrIdxRangeMap.get(chr);
		if(windowsize == -1){
			for(int i = ip.x; i < ip.y; i++){
				outGenes.add(genes.get(i).name);
			}
		}else{
			if(windowsize > (ip.y - ip.x)){
				return null;
			}
			int start = idx - windowsize/2;
			int end = idx + windowsize/2;
			
			//if(start < ip.x) {end = end + ip.x - start; start = ip.x; }
			//if(end >= ip.y) {start = start - end + ip.y - 1; end = ip.y-1; }
			System.out.println(idx + "\t" + start + "\t" + end);
			for(int i = start; i <= end; i++){
				outGenes.add(genes.get(i).name);
			}
		}
		return outGenes.toArray(new String[0]);
	}
	
	public String getChr(String gene){
		return chrMap.get(gene);
	}
	
	public boolean contains(String gene){
		return genes.contains(new Gene(gene));
	}
	
}
