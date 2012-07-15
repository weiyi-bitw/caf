package caf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import obj.Genome;
import obj.IntPair;
import obj.ValString;

public class GroupCNVWindow2 {
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
		double avgMI;
		double minMI;
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
			this.minMI = Float.MAX_VALUE;
			this.matchNumber = 0;
			for(CNVWindow cnvw : content){
				if(cnvw != null){
					matchNumber++;
					avgMI += cnvw.val;
					if(cnvw.val < minMI){
						minMI = cnvw.val;
					}
					for(ValString g : cnvw.geneNames){
						StringIntPair sip = new StringIntPair(g.s);
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
				return -Double.compare(this.minMI, other.minMI);
			}
			return -Double.compare(this.matchNumber, other.matchNumber);
		}
		
		public String toString(){
			String s = chr + "\t" + matchNumber + "\t" + minMI + "\n";
			for(int i = 0; i < k; i++){
				s += names[i];
				s += "\t" + content[i] + "\n";
			}
			s += "Top genes:";
			for(int i = 0; i < 10; i++){
				s += "\t" + allGenes.get(i);
			}
			s += "\n";
			return s;
		}
		public int getNumGenes(){
			return allGenes.size();
		}
		
		public ArrayList<ValString> calAverageMI(){
			ArrayList<ValString> out = new ArrayList<ValString>();
			for(int i = 0; i < k; i++){
				if(content[i] != null){
					for(ValString g : content[i].geneNames){
						if(out.contains(g)){
							out.get(out.indexOf(g)).val += g.val;
						}else{
							out.add(g);
						}
					}
				}
			}
			
			for(int i = 0; i < out.size(); i++){
				out.get(i).val /= this.matchNumber;
			}
			Collections.sort(out);
			
			return out;
		}
	}
	
	static class CNVWindow implements Comparable<CNVWindow>{
		static Genome gn;
		String name;
		String chr;
		String chrArm;
		int source;
		ArrayList<ValString> geneNames;
		IntPair range;
		
		double val;
		
		static void linkGenome(Genome gn){
			CNVWindow.gn = gn;
		}
		
		CNVWindow(String name, String chr, String chrArm, ArrayList<ValString> geneNames, IntPair range, double val, int source){
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
			double val = Double.parseDouble(tokens[nt-1]);
			ArrayList<ValString> genes = new ArrayList<ValString>();
			int x = Integer.MAX_VALUE;
			int y = -1;
			for(int i = 2; i < nt-2; i++){
				String[] t2 = tokens[i].split(":");
				String g = t2[0];
				Double v = Double.parseDouble(t2[1]);
				genes.add(new ValString(g, v));
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
			String a = gn.getChrBand(gn.getGene(range.x));
			String b = gn.getChrBand(gn.getGene(range.y));
			
			String s = gn.getChr(geneNames.get(0).s);
			if(a.equals(b)){
				s += a;
			}else{
				int idx = 0;
				if(b.contains("p")) idx = b.lastIndexOf("p");
				if(b.contains("q")) idx = b.lastIndexOf("q");
				
				int idxa = a.length();
				if(a.contains("-")) idxa = a.indexOf("-");
				if(a.contains("|")) idxa = a.indexOf("|");
				s += a.substring(0, idxa)+ "-" +b.substring(idx) ;
			}
			
			for(ValString g : geneNames){
				s += "\t" + g.s;
			}
			s += "\t" + val;
			return s;
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		String inPath = "/home/weiyi/workspace/javaworks/caf/output/window/";
		int loadIn = 500;
		
		String[] files = new File(inPath + "mergeroom").list();
		Arrays.sort(files);
		int nf = files.length;
		System.out.println(nf + " files in the directory.");
		
		final String geneLocFile = "/home/weiyi/workspace/data/annot/ncbi/gene.location.ncbi";
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
		
		PrintWriter pw2 = new PrintWriter(new FileWriter(inPath + "/consensus." + loadIn + ".txt"));
		DecimalFormat df = new DecimalFormat("0.0000"); 
		ii = 1;
		for(CNVWindowSet cnvws : out){
			ArrayList<ValString> consensus = cnvws.calAverageMI();
			pw2.println();
			pw2.println(ii + ". " + gn.getChr(consensus.get(0).s) +  gn.getChrBand(consensus.get(0).s) + "\t\t\tMinimum Strength: " 
			+ df.format(cnvws.minMI) + "\t\t\tCommon Dataset: " + cnvws.matchNumber);
			pw2.print("Gene");
			for(int i = 0; i < 10; i++){
				ValString vs = consensus.get(i);
				pw2.print("\t" + vs.s);
			}pw2.println();
			pw2.print("Avg MI");
			for(int i = 0; i < 10; i++){
				ValString vs = consensus.get(i);
				pw2.print("\t" + df.format(vs.val));
			}pw2.println();
			ii++;
		}
		pw2.close();
		
		System.out.println("Done.");
		
	}

}
