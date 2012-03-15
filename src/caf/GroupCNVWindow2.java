package caf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import obj.Genome;

public class GroupCNVWindow2 {
	static class IntPair{
		int x;
		int y;
		
		IntPair(int x, int y){
			if(x <= y){
				this.x = x;
				this.y = y;
			}else{
				this.x = y;
				this.y = x;
			}
		}
		
		public boolean overlapWith(IntPair ip2){
			if(this.x > ip2.y && this.y > ip2.y){
				return false;
			}else if(this.x < ip2.x && this.y < ip2.x){
				return false;
			}else{
				return true;
			}
		}
		
		void setX(int x){
			this.x = x;
		}
		
		void setY(int y){
			this.y = y;
		}
		
	}
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
	
	static class CNVWindowSet implements Comparable<CNVWindowSet>{
		static String[] names;
		static int k;
		String chr;
		CNVWindow[] content;
		int matchNumber;
		ArrayList<StringIntPair> allGenes;
		float avgMI;
		int numCommonGenes;
		
		static void setNames(String[] names){
			CNVWindowSet.names = names;
			CNVWindowSet.k = names.length;
		}
		
		CNVWindowSet(String chr, CNVWindow[] content){
			this.chr = chr;
			this.content = content;
			this.allGenes = new ArrayList<StringIntPair>();
			this.avgMI = 0;
			this.matchNumber = 0;
			for(CNVWindow cnvw : content){
				if(cnvw != null){
					matchNumber++;
					avgMI += cnvw.val;
					for(String g : cnvw.geneNames){
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
		
		public int compareTo(CNVWindowSet other) {
			if(this.matchNumber == other.matchNumber){
				return Double.compare(this.getNumGenes(), other.getNumGenes());
			}
			return -Double.compare(this.matchNumber, other.matchNumber);
		}
		
		public String toString(){
			String s = chr + "\t" + matchNumber + "\t" + allGenes.size() + "\t" + avgMI + "\t" + numCommonGenes + "\n";
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
	
	static class CNVWindow implements Comparable<CNVWindow>{
		static Genome gn;
		String name;
		String chr;
		String chrArm;
		int source;
		ArrayList<String> geneNames;
		IntPair range;
		
		float val;
		
		static void linkGenome(Genome gn){
			CNVWindow.gn = gn;
		}
		
		CNVWindow(String name, String chr, String chrArm, ArrayList<String> geneNames, IntPair range, float val, int source){
			this.name = name;
			this.chr = chr;
			this.chrArm = chrArm;
			this.geneNames = geneNames;
			this.range = range;
			this.val = val;
			this.source = source;
		}
		
		static CNVWindow parseCNVWindow(String line, int source){
			String[] tokens = line.split("\t");
			String name = tokens[0];
			String chr = tokens[1];
			int nt = tokens.length;
			String chrArm = tokens[nt-2];
			float val = Float.parseFloat(tokens[nt-1]);
			ArrayList<String> genes = new ArrayList<String>();
			int x = Integer.MAX_VALUE;
			int y = -1;
			for(int i = 2; i < nt-2; i++){
				String g = tokens[i];
				genes.add(g);
				int j = gn.getIdx(g);
				if(j < x){
					x = j;
				}
				if(j > y){
					y = j;
				}
			}
			
			return new CNVWindow(name, chr, chrArm, genes, new IntPair(x, y), val, source);
			
		}
		
		public boolean overlapWith(CNVWindow other){
			if(!this.chr.equals(other.chr)){
				return false;
			}
			if(this.range.overlapWith(other.range)){
				return true;
			}else{
				return false;
			}
		}
		
		public int compareTo(CNVWindow other) {
			return -Double.compare(this.val, other.val);
		}
		
		public String toString(){
			String s = name + "\t" + chr + "\t" + chrArm;
			for(String g : geneNames){
				s += "\t" + g;
			}
			s += "\t" + val;
			return s;
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		String inPath = "/home/weiyi/workspace/javaworks/caf/output/window51/";
		int loadIn = 500;
		
		String[] files = new File(inPath + "mergeroom").list();
		Arrays.sort(files);
		int nf = files.length;
		System.out.println(nf + " files in the directory.");
		
		final String geneLocFile = "/home/weiyi/workspace/data/annot/affy/u133p2/gene.location3";
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		CNVWindow.linkGenome(gn);
		
		ArrayList<CNVWindow> allCNVWindows = new ArrayList<CNVWindow>();
		int cnt = 0;
		for(String f : files){
			System.out.println("Loading file " + f + "...");
			BufferedReader br = new BufferedReader(new FileReader(inPath + "mergeroom/" + f));
			String line = br.readLine();
			int cnt2 = 0;
			while(line != null && cnt2 < loadIn){
				allCNVWindows.add(CNVWindow.parseCNVWindow(line, cnt));
				cnt2++;
				line = br.readLine();
			}
			cnt++;
			br.close();
			
		}
		
		System.out.println(allCNVWindows.size() + " CNV windows were loaded.");
		
		System.out.println("Making matches...");
		
		Collections.sort(allCNVWindows);
		ArrayList<CNVWindowSet> out = new ArrayList<CNVWindowSet>();
		CNVWindowSet.setNames(files);
		
		while(allCNVWindows.size() > 0){
			CNVWindow cnvw = allCNVWindows.get(0);
			allCNVWindows.remove(0);
			CNVWindow[] row = new CNVWindow[nf];
			String chr = cnvw.chr;
			row[cnvw.source] = cnvw;
			int matchNumber = 1;
			for(int j = allCNVWindows.size() - 1; j >= 0; j--){
				CNVWindow cnvw2 = allCNVWindows.get(j);
				if(cnvw.source != cnvw2.source && cnvw.chr.equals(cnvw2.chr)){
					if(cnvw.overlapWith(cnvw2)){
						if(row[cnvw2.source] == null) matchNumber++;
						row[cnvw2.source] = cnvw2;
						allCNVWindows.remove(j);
						
					}
				}
			}
			if(matchNumber >1){
				out.add(new CNVWindowSet(chr, row));
			}
		}
		
		Collections.sort(out);
		
		System.out.println("Output to file...");
		PrintWriter pw = new PrintWriter(new FileWriter(inPath + "/matchTable." + loadIn + ".txt"));
		int ii = 1;
		for(CNVWindowSet cnvws : out){
			pw.print(ii + "\t");
			pw.println(cnvws);
			ii++;
		}
		
		pw.close();

		System.out.println("Done.");
		
	}

}
