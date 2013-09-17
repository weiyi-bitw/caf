package caf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;

import obj.ValString;

public class GroupWeightedAttractor2 {
	public static class DistPair implements Comparable<DistPair>{
		static int n = 1000000;
		WtdAttractorSet x;
		WtdAttractorSet y;
		double sim; // similarities
		
		DistPair(WtdAttractorSet a, WtdAttractorSet b, double sim){
			if(a.id < b.id){ // x must always be smaller than y
				this.x = a;
				this.y = b;
			}else{
				this.x = b;
				this.y = a;
			}
			this.sim = sim;
		}
		boolean contains(WtdAttractorSet k){
			return this.x.equals(k) || this.y.equals(k);
		}
		public boolean equals(Object other){
			boolean result = false;
	        if (other instanceof DistPair) {
	        	DistPair that = (DistPair) other;
	        	result = this.x.equals(that.x) && this.y.equals(that.y);
	        }
	        return result;
		}
		public int hashCode(){
			return(x.hashCode() * n + y.hashCode());
		}
		public int compareTo(DistPair other) {
			return -Double.compare(this.sim, other.sim);
		}
		static void setTotalIdx(int n){
			DistPair.n = n;
		}
		public String toString(){
			String s = this.x.id + "\t" + this.y.id + "\t" + this.sim;
			return s;
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
	
	public static class WtdAttractor implements Comparable<WtdAttractor>{
		static int quantile = 10;
		static long randseed = System.currentTimeMillis();
		static Random randMachine;
		String name;
		int source;
		ArrayList<ValString> genes;
		int basins;
		double val;
		
		WtdAttractor(String name, int basins, ArrayList<ValString> genes,int source){
			this.name = name;
			this.basins = basins;
			this.genes = genes;
			this.val = genes.get(quantile-1).val;
			this.source = source;
		}
		static WtdAttractor parseWtdAttractor(String line, int source){
			String[] tokens = line.split("\t");
			int nt = tokens.length;
			String name = tokens[0];
			int basins = Integer.parseInt(tokens[1]);
			ArrayList<ValString> genes = new ArrayList<ValString>();
			for(int i = 2; i < 50+2; i++){
				//System.out.println(tokens[i]);
				String[] t2 = tokens[i].split(":");
				genes.add(new ValString(t2[0], Double.parseDouble(t2[1])));
			}
			return new WtdAttractor(name, basins, genes, source);
		}
		static void setSeed(long s){
			WtdAttractor.randseed = s;
			WtdAttractor.randMachine = new Random(s);
		}
		static WtdAttractor generateRandomAttractor(int id, int source, int geneSize){
			HashSet<Integer> hs = new HashSet<Integer>();
			while(hs.size() < quantile){
				int k = randMachine.nextInt(geneSize);
				hs.add(k);
			}
			String name = "Attractor" + String.format("%05d", id);
			int basins = 0;
			ArrayList<ValString> genes = new ArrayList<ValString>();
			for(int k : hs){
				genes.add(new ValString("G" + k, 0.0));
			}
			return new WtdAttractor(name, basins, genes, source);
			
		}
		
		
		
		public ArrayList<ValString> getOvlp(WtdAttractor other){
			if(other == null) return null;
			ArrayList<ValString> al = new ArrayList<ValString>();
			for(ValString g : other.genes){
				if(this.genes.contains(g)){
					int idx = this.genes.indexOf(g);
					ValString gb = this.genes.get(idx);
					al.add(new ValString(g.s, 2));
				}
			}
			return al;
		}
		public ArrayList<ValString> getOvlp(ArrayList<ValString> in){
			ArrayList<ValString> al = new ArrayList<ValString>();
			for(ValString g : in){
				if(this.genes.contains(g)){
					al.add(new ValString(g.s, g.val+1));
				}
			}
			return al;
		}
		
		public double overlapWith(WtdAttractor other){
			if(other == null) return 0;
			double cnt = 0;
			for(ValString g : other.genes){
				if(this.genes.contains(g)){
					int i = this.genes.indexOf(g);
					//cnt += (g.val + this.genes.get(i).val);
					cnt ++;
				}
			}
			return cnt;
			
		}
		
		public int compareTo(WtdAttractor other) {
			if(this.basins == other.basins){
				return -Double.compare(this.val, other.val);
			}
			return -Double.compare(this.basins, other.basins);
		}
		
		public boolean equals(Object other){
			boolean result = false;
	        if (other instanceof WtdAttractor) {
	        	WtdAttractor that = (WtdAttractor) other;
	        	result = (this.source == that.source) && (this.name.equals(that.name));
	        }
	        return result;
		}
		
		public String toString(){
			String s ="";
			for(ValString g : genes){
				if(g.s.startsWith("+AF8AXwBf")) continue;
				s += g.s + "\t";
			}
			//s += val;
			return s;
		}
		public static void setquantile(int k){
			WtdAttractor.quantile = k;
		}
	}
	
	public static class WtdAttractorSet implements Comparable<WtdAttractorSet>{
		static String[] names;
		static int k;
		static int count = 0;
		int id;
		WtdAttractor[] content;
		int matchNumber;
		ArrayList<StringIntPair> allGenes;
		double avgMI;
		int numCommonGenes;
		double minMI;
		
		static void setNames(String[] names){
			WtdAttractorSet.names = names;
			WtdAttractorSet.k = names.length;
		}
		
		WtdAttractorSet(WtdAttractor wa){
			this.id = count;
			count++;
			this.content = new WtdAttractor[k];
			this.content[wa.source] = wa;
			this.allGenes = new ArrayList<StringIntPair>();
			this.avgMI = 0;
			this.minMI = wa.val;
			this.matchNumber = 1;
			this.numCommonGenes = 0;
			
			for(ValString g : wa.genes){
				allGenes.add(new StringIntPair(g.s,1));
			}
			
			avgMI = wa.val;
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
					for(ValString g : wa.genes){
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
		
		public int compareTo(WtdAttractorSet other) {
			if(this.matchNumber == other.matchNumber){
				return -Double.compare(this.minMI, other.minMI);
			}else{
				return -Double.compare(this.matchNumber, other.matchNumber);
			}
		}
		public String toString(){
			String s = "ID" + id + "\t" + matchNumber + "\n";//"\t" + minMI + "\n";
			for(int i = 0; i < k; i++){
				s += names[i] + "\t" + content[i] + "\n";
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
		public int ovlapWith(WtdAttractor wa){
			int nn = 0;
			int sz = wa.genes.size();
			for(int i = 0; i < sz; i++){
				StringIntPair sip = new StringIntPair(wa.genes.get(i).s, 1);
				if(allGenes.contains(sip)){
					nn += i;
				}
			}
			return nn;
		}
		public double ovlapWith(WtdAttractorSet was){
			double nn = 0;
			for(int i = 0; i < k; i++){
				if(content[i] == null) continue;
				WtdAttractor wa = content[i];
				for(int j = 0; j < k; j++){
					if(i != j){
						nn += wa.overlapWith(was.content[j]);
					}
				}
			}
			return nn;
		}
		public boolean merge(WtdAttractorSet was){
			int ovlpcnt = 0;
			ArrayList<Integer> conflictIdx = new ArrayList<Integer>();
			for(int i = 0; i < k; i++){
				if (this.content[i] != null && was.content[i] != null){
					conflictIdx.add(i);
					ovlpcnt++;
					if(ovlpcnt > 0){
						return false;
					}
				}
			}
			
			for(int i = 0; i < k; i++){
				if(this.content[i] == null && was.content[i] != null){
					this.content[i] = was.content[i];
					matchNumber++;
				}
			}
			for(int i = 0; i < k; i++){
				if(this.content[i] != null && was.content[i] != null){
					WtdAttractor wa = was.content[i];
					int wasCnt = 0;
					int thisCnt = 0;
					for(int j = 0; j < k; j++){
						if(!conflictIdx.contains(j)){
							wasCnt += wa.overlapWith(this.content[j]);
							thisCnt += this.content[i].overlapWith(this.content[j]);
						}
					}
					
					if(wasCnt > thisCnt){
					//if(was.content[i].basins > this.content[i].basins){
						//WtdAttractor wa = this.content[i];
						//System.out.println(wa.source + "\t" + wa.name + " replaced.");
						this.content[i] = was.content[i];
					}
				}
			}
			allGenes = new ArrayList<StringIntPair>();
			numCommonGenes = 0;
			avgMI = 0;
			minMI = Double.MAX_VALUE;
			for(WtdAttractor wa2 : content){
				if(wa2 == null) continue;
				avgMI += wa2.val;
				if(wa2.val < minMI) minMI = wa2.val;
				for(ValString s : wa2.genes){
					StringIntPair sip = new StringIntPair(s.s, 1);
					if(!allGenes.contains(sip)){
						allGenes.add(sip);
					}else{
						sip = allGenes.get(allGenes.indexOf(sip));
						sip.incre();
						if(sip.i == k){
							numCommonGenes++;
						}
					}
				}
			}
			avgMI /= matchNumber;
			Collections.sort(allGenes);
			return true;
		}
		public void addWtd(WtdAttractor wa){
			if(content[wa.source] == null){
				content[wa.source] = wa;
				if(wa.val < minMI) minMI = wa.val;
				avgMI = avgMI * matchNumber + wa.val;
				matchNumber ++ ;
				avgMI /= matchNumber;
				for(ValString vs : wa.genes){
					StringIntPair sip = new StringIntPair(vs.s, 1);
					if(!allGenes.contains(sip)){
						allGenes.add(sip);
					}else{
						sip = allGenes.get(allGenes.indexOf(sip));
						sip.incre();
						if(sip.i == k){
							numCommonGenes++;
						}
					}
				}
			}
			else if(wa.basins > content[wa.source].basins){
				WtdAttractor preWa = content[wa.source];
				content[wa.source] = wa;
				if(wa.val < minMI) minMI = wa.val;
				avgMI = avgMI + (wa.val - preWa.val)/matchNumber;
				allGenes = new ArrayList<StringIntPair>();
				for(WtdAttractor wa2 : content){
					if(wa2 == null) continue;
					for(ValString s : wa2.genes){
						StringIntPair sip = new StringIntPair(s.s, 1);
						if(!allGenes.contains(sip)){
							allGenes.add(sip);
						}else{
							sip = allGenes.get(allGenes.indexOf(sip));
							sip.incre();
							if(sip.i == k){
								numCommonGenes++;
							}
						}
					}
				}
			}
			Collections.sort(allGenes);
			
		}
		public boolean contains(WtdAttractor a){
			if(content[a.source] == null){
				return false;
			}else{
				return content[a.source].equals(a);
			}
		}
		public boolean equals(WtdAttractorSet was){
			for(int i = 0; i < k; i++){
				if(this.content[i] != null && was.content[i] != null){
					if(!this.content[i].equals(was.content[i])){
						return false;
					}
				}else if((this.content[i] != null && was.content[i] == null) || (this.content[i] == null && was.content[i] != null)){
					return false;
				}
			}
			return true;
		}
		public int hashCode(){
			return id;
		}
		
		public ArrayList<ValString> calAverageMI(){
			ArrayList<ValString> out = new ArrayList<ValString>();
			for(int i = 0; i < k; i++){
				if(content[i] != null){
					for(ValString g : content[i].genes){
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
	
	/*static class WaCoreSet implements Comparable<WaCoreSet>{
		static String[] names;
		static int k;
		static int count = 0;
		int id;
		WtdAttractor[] content;
		int matchNumber;
		ArrayList<ValString> coreGenes;
		double sumMI;
		int numCommonGenes;
		double minMI;
		
		static void setNames(String[] names){
			WaCoreSet.names = names;
			WaCoreSet.k = names.length;
		}
		
		WaCoreSet(WtdAttractor wa, WtdAttractor wb){
			this.id = count;
			count++;
			this.content = new WtdAttractor[k];
			this.content[wa.source] = wa;
			this.content[wb.source] = wb;
			this.coreGenes = wa.getOvlp(wb);
			this.sumMI = wa.val + wb.val;
			this.minMI = Math.min(wa.val, wb.val);
			this.matchNumber = 2;
			this.numCommonGenes = coreGenes.size();
		}
		WtdAttractor fillInBest(ArrayList<WtdAttractor> in ){
			double score = 0;
			WtdAttractor winner = null;
			for(WtdAttractor wa : in ){
				if(content[wa.source] != null) continue;
				double ovlp = this.ovlapWith(wa);
				if(ovlp > score){
					score = ovlp;
					winner = wa;
				}
			}
			if(winner != null){
				content[winner.source] = winner;
				coreGenes = winner.getOvlp(coreGenes);
				matchNumber++;
			}
			return winner;
		}
		public double ovlapWith(WtdAttractor wa){
			double nn = 0;
			for(ValString vs : wa.genes){
				if(coreGenes.contains(vs)){
					//nn += (coreGenes.get(coreGenes.indexOf(vs)).val + vs.val);
					nn++;
				}
			}
			return nn;
		}
		
		public boolean contains(WtdAttractor a){
			if(content[a.source] == null){
				return false;
			}else{
				return content[a.source].equals(a);
			}
		}
		public int compareTo(WaCoreSet o) {
			if(this.matchNumber == o.matchNumber){
				return -Double.compare(this.minMI, o.minMI);
			}else{
				return -Double.compare(this.matchNumber, o.matchNumber);
			}
		}
		public boolean equals(WaCoreSet was){
			for(int i = 0; i < k; i++){
				if(this.content[i] != null && was.content[i] != null){
					if(!this.content[i].equals(was.content[i])){
						return false;
					}
				}else if((this.content[i] != null && was.content[i] == null) || (this.content[i] == null && was.content[i] != null)){
					return false;
				}
			}
			return true;
		}
		public int hashCode(){
			return id;
		}
		
		public String toString(){
			String s = "ID" + id + "\t" + matchNumber + "\t" + coreGenes.size() + "\t" + minMI + "\t" + numCommonGenes + "\n";
			for(int i = 0; i < k; i++){
				s += names[i];
				s += "\t" + content[i] + "\n";
			}
			s += "Top genes:";
			for(ValString vs : coreGenes){
				s += "\t" + vs.toString();
			}
			s += "\n";
			return s;
		}
		
		public boolean full(){
			return matchNumber == k || coreGenes.size() == 0;
		}
	}
	*/
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		//String inPath = "/home/weiyi/workspace/javaworks/caf/output/brca/weighted/";
		String inPath = "/home/weiyi/workspace/data/caf_clusterValid/pca/";
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
			System.out.print("Loading file " + f + "...");
			ArrayList<WtdAttractor> waInThisFile = new ArrayList<WtdAttractor>();
			BufferedReader br = new BufferedReader(new FileReader(inPath + "mergeroom/" + f));
			String line = br.readLine();
			while(line != null){
				WtdAttractor wa = WtdAttractor.parseWtdAttractor(line, cnt);
				if(wa.basins >= 0) waInThisFile.add(wa);
				line = br.readLine();
			}
			allWtdAttractors.addAll(waInThisFile);
			System.out.println(" (" + waInThisFile.size() + ") ");
			cnt++;
			br.close();
		}
		
		int n = allWtdAttractors.size();
		System.out.println(n + " attractors were loaded.");
		DistPair.setTotalIdx(n);
		Collections.sort(allWtdAttractors);
		
		ArrayList<DistPair> allDistPairs = new ArrayList<DistPair>();
		
		WtdAttractorSet.setNames(files);
		ArrayList<WtdAttractorSet> out = new ArrayList<WtdAttractorSet>();
		System.out.println("Calculating distance...");
		for(int i = 0; i < n; i++){
			WtdAttractorSet a = new WtdAttractorSet(allWtdAttractors.get(i));
			out.add(a);
		}
		
		for(int i = 0; i < n; i++){
			WtdAttractorSet a = out.get(i);
			
			for(int j = i+1; j < n; j++){
				WtdAttractorSet b = out.get(j);
				double ovlp = a.ovlapWith(b);
				if(ovlp > 0){
					allDistPairs.add(new DistPair(a, b , ovlp));
				}
			}
		}
		Collections.sort(allDistPairs);
		System.out.println(allDistPairs.size() + " pairs of distances have been added.");
		
		while(allDistPairs.size() > 0){
			DistPair dp = allDistPairs.get(0);
			allDistPairs.remove(dp);
			
		
			System.out.println(dp);
			WtdAttractorSet x = dp.x;
			boolean mergeable = x.merge(dp.y);
			if(mergeable){
				out.remove(dp.y);
				for(int i = allDistPairs.size() - 1; i >=0; i--){
					DistPair dpp = allDistPairs.get(i);
					if(dpp.contains(dp.x) || dpp.contains(dp.y)){
						allDistPairs.remove(i);
					}
				}
				int xi = out.indexOf(dp.x);
				for(int i = 0; i < out.size(); i++){
					if(i != xi){
						WtdAttractorSet b = out.get(i);
						double ovlp = x.ovlapWith(b);
						if(ovlp > 0){
							allDistPairs.add(new DistPair(x, b , ovlp));
						}
					}
				}
				Collections.sort(allDistPairs);
			}
		}
		
		System.out.println("Making matches...");
		
		
		Collections.sort(out);
		System.out.println("Output to file...");
		PrintWriter pw = new PrintWriter(new FileWriter(inPath + "/matchTable.txt"));
		int ii = 1;
		
		for(WtdAttractorSet was : out){
			pw.print(ii + "\t");
			pw.println(was);
			ii++;
		}
		
		pw.close();
		DecimalFormat df = new DecimalFormat("0.0000"); 
		PrintWriter pw1 = new PrintWriter(new FileWriter(inPath + "/consensus.txt"));
		int ii1 = 0;
		for(WtdAttractorSet was : out){
			if(ii1 >= 50) break;
			pw1.println();
			pw1.println((ii1+1) + ".\t\t\tMinimum Strength: " + was.minMI + "\t\t\tCommon Datasets:" + was.matchNumber);
			ArrayList<ValString> consensus = was.calAverageMI();
			pw1.print("Gene");
			for(int i = 0; i < 50; i++){
				ValString vs = consensus.get(i);
				if(vs.s.startsWith("___")) continue;
				pw1.print("\t" + vs.s);
			}pw1.println();
			pw1.print("Avg MI");
			for(int i = 0; i < 50; i++){
				ValString vs = consensus.get(i);
				pw1.print("\t" + df.format(vs.val));
			}pw1.println();
			ii1++;
		}
		pw1.close();
		System.out.println("Done.");
		
	}

}
