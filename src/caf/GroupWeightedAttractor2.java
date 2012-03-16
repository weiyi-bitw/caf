package caf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class GroupWeightedAttractor2 {
	static class StringIntPair implements Comparable<StringIntPair>{
		String s;
		int i;
		
		StringIntPair(String s){
			this.s = s;
			this.i = 1;
		}
		StringIntPair(String s, int i){
			this.s = s;
			this.i = i;
		}
		
		public boolean equals(Object other){
			boolean result = false;
	        if (other instanceof StringIntPair) {
	        	StringIntPair that = (StringIntPair) other;
	        	result = this.s.equals(that.s);
	        }
	        return result;
		}
		
		public int hashCode(){
			return this.s.hashCode();
		}
		
		public void incre(){
			this.i += 1;
		}
		
		public String toString(){
			String out = s + ":" + i;
			return out;
		}
		
		public int compareTo(StringIntPair other) {
			return -Double.compare(this.i, other.i);
		}
		
	}
	
	static class WtdAttractor implements Comparable<WtdAttractor>{
		String name;
		int source;
		ArrayList<String> geneNames;
		int basins;
		float val;
		
		WtdAttractor(String name, int basins, ArrayList<String> geneNames, float val,int source){
			this.name = name;
			this.basins = basins;
			this.geneNames = geneNames;
			this.val = val;
			this.source = source;
		}
		static WtdAttractor parseWtdAttractor(String line, int source){
			String[] tokens = line.split("\t");
			int nt = tokens.length;
			String name = tokens[0];
			int basins = Integer.parseInt(tokens[1]);
			float val = Float.parseFloat(tokens[nt-1]);
			ArrayList<String> genes = new ArrayList<String>();
			for(int i = 2; i < nt-2; i++){
				String g = tokens[i];
				genes.add(g);
			}
			return new WtdAttractor(name, basins, genes, val, source);
		}
		public boolean overlapWith(WtdAttractor other){
			for(String g : other.geneNames){
				if(this.geneNames.contains(g)){
					return true;
				}
			}
			return false;
			
		}
		
		public int compareTo(WtdAttractor other) {
			return -Double.compare(this.basins, other.basins);
		}
		
		public String toString(){
			String s = name + "\t" + basins;
			for(String g : geneNames){
				s += "\t" + g;
			}
			s += "\t" + val;
			return s;
		}
	}
	
	static class WtdAttractorSet implements Comparable<WtdAttractorSet>{
		static String[] names;
		static int k;
		WtdAttractor[] content;
		int matchNumber;
		ArrayList<StringIntPair> allGenes;
		float avgMI;
		int numCommonGenes;
		
		static void setNames(String[] names){
			WtdAttractorSet.names = names;
			WtdAttractorSet.k = names.length;
		}
		
		WtdAttractorSet(WtdAttractor[] content){
			this.content = content;
			this.allGenes = new ArrayList<StringIntPair>();
			this.avgMI = 0;
			this.matchNumber = 0;
			
			for(WtdAttractor wa : content){
				if(wa != null){
					matchNumber++;
					avgMI += wa.val;
					for(String g : wa.geneNames){
						StringIntPair sip = new StringIntPair(g);
						if(allGenes.contains(sip)){
							allGenes.get(allGenes.indexOf(sip)).incre();
						}else{
							allGenes.add(sip);
						}
					}
				}
			}
			avgMI /= matchNumber;
			this.numCommonGenes = 0;
			Collections.sort(allGenes);
			for(StringIntPair sip: allGenes){
				if(sip.i == k){
					numCommonGenes++;
				}else{
					break;
				}
			}
		}
		
		public int compareTo(WtdAttractorSet other) {
			if(this.matchNumber == other.matchNumber){
				return Double.compare(this.allGenes.size(), other.allGenes.size());
			}else{
				return -Double.compare(this.matchNumber, other.matchNumber);
			}
		}
		public String toString(){
			String s = matchNumber + "\t" + allGenes.size() + "\t" + avgMI + "\t" + numCommonGenes + "\n";
			for(int i = 0; i < k; i++){
				s += names[i];
				s += "\t" + content[i] + "\n";
			}
			s += "Top genes:";
			for(int i = 0; i < 20; i++){
				s += "\t" + allGenes.get(i);
			}
			s += "\n";
			return s;
		}
		public int getNumGenes(){
			return allGenes.size();
		}
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String inPath = "/home/weiyi/workspace/javaworks/caf/output/weighted.6/";
		if(!inPath.endsWith("/")){
			inPath += "/";
		}
		
		String[] files = new File(inPath + "mergeroom").list();
		Arrays.sort(files);
		int nf = files.length;
		System.out.println(nf + " files in the directory.");
		
		ArrayList<WtdAttractor> allWtdAttractors = new ArrayList<WtdAttractor>();
		
		int cnt = 0;
		for(String f : files){
			System.out.println("Loading file " + f + "...");
			BufferedReader br = new BufferedReader(new FileReader(inPath + "mergeroom/" + f));
			String line = br.readLine();
			while(line != null){
				allWtdAttractors.add(WtdAttractor.parseWtdAttractor(line, cnt));
				line = br.readLine();
			}
			cnt++;
			br.close();
		}
		
		
		System.out.println(allWtdAttractors.size() + " CNV windows were loaded.");
		
		System.out.println("Making matches...");
		
		Collections.sort(allWtdAttractors);
		ArrayList<WtdAttractorSet> out = new ArrayList<WtdAttractorSet>();
		WtdAttractorSet.setNames(files);
		
		while(allWtdAttractors.size() > 0){
			WtdAttractor wa = allWtdAttractors.get(0);
			allWtdAttractors.remove(0);
			WtdAttractor[] row = new WtdAttractor[nf];
			row[wa.source] = wa;
			int matchNumber = 1;
			for(int j = allWtdAttractors.size()-1; j >= 0; j --){
				WtdAttractor wa2 = allWtdAttractors.get(j);
				if(wa.source != wa2.source ){
					if(wa.overlapWith(wa2)){
						if(row[wa2.source] == null) matchNumber++;
						row[wa2.source] = wa2;
						allWtdAttractors.remove(j);
					}
				}
			}
			if(matchNumber > 1){
				out.add(new WtdAttractorSet(row));
			}
			
			Collections.sort(out);
			
			
		}
		System.out.println("Output to file...");
		PrintWriter pw = new PrintWriter(new FileWriter(inPath + "/matchTable." + ".txt"));
		int ii = 1;
		
		for(WtdAttractorSet was : out){
			pw.print(ii + "\t");
			pw.println(was);
			ii++;
		}
		
		pw.close();

		System.out.println("Done.");
		
	}

}
