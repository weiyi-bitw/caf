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
	static class DistPair implements Comparable<DistPair>{
		static int n = 1000000;
		WtdAttractorSet x;
		WtdAttractorSet y;
		float sim; // similarities
		
		DistPair(WtdAttractorSet a, WtdAttractorSet b, float sim){
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
			for(int i = 2; i < nt-1; i++){
				String g = tokens[i];
				genes.add(g);
			}
			return new WtdAttractor(name, basins, genes, val, source);
		}
		public int overlapWith(WtdAttractor other){
			if(other == null) return 0;
			int cnt = 0;
			for(String g : other.geneNames){
				if(this.geneNames.contains(g)){
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
		static int count = 0;
		int id;
		WtdAttractor[] content;
		int matchNumber;
		ArrayList<StringIntPair> allGenes;
		float avgMI;
		int numCommonGenes;
		float minMI;
		
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
			
			for(String g : wa.geneNames){
				allGenes.add(new StringIntPair(g,1));
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
				return -Double.compare(this.minMI, other.minMI);
			}else{
				return -Double.compare(this.matchNumber, other.matchNumber);
			}
		}
		public String toString(){
			String s = "ID" + id + "\t" + matchNumber + "\t" + allGenes.size() + "\t" + minMI + "\t" + numCommonGenes + "\n";
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
		public int ovlapWith(WtdAttractor wa){
			int nn = 0;
			int sz = wa.geneNames.size();
			for(int i = 0; i < sz; i++){
				StringIntPair sip = new StringIntPair(wa.geneNames.get(i), 1);
				if(allGenes.contains(sip)){
					nn += i;
				}
			}
			return nn;
		}
		public int ovlapWith(WtdAttractorSet was){
			int nn = 0;
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
					if(ovlpcnt > 2){
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
			minMI = Float.MAX_VALUE;
			for(WtdAttractor wa2 : content){
				if(wa2 == null) continue;
				avgMI += wa2.val;
				if(wa2.val < minMI) minMI = wa2.val;
				for(String s : wa2.geneNames){
					StringIntPair sip = new StringIntPair(s, 1);
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
				for(String s : wa.geneNames){
					StringIntPair sip = new StringIntPair(s, 1);
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
					for(String s : wa2.geneNames){
						StringIntPair sip = new StringIntPair(s, 1);
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
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String inPath = "/home/weiyi/workspace/javaworks/caf/output/weighted/";
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
			ArrayList<WtdAttractor> waInThisFile = new ArrayList<WtdAttractor>();
			BufferedReader br = new BufferedReader(new FileReader(inPath + "mergeroom/" + f));
			String line = br.readLine();
			while(line != null){
				waInThisFile.add(WtdAttractor.parseWtdAttractor(line, cnt));
				line = br.readLine();
			}
			/*Collections.sort(waInThisFile);
			for(int i = 0; i < waInThisFile.size(); i++){
				WtdAttractor w = waInThisFile.get(i);
				for(int j = waInThisFile.size() - 1; j > i; j--){
					if(w.overlapWith(waInThisFile.get(j))>0){
						waInThisFile.remove(j);
					}
				}
			}*/
			
			allWtdAttractors.addAll(waInThisFile);
			cnt++;
			br.close();
		}
		
		int n = allWtdAttractors.size();
		System.out.println(n + " attractors were loaded.");
		DistPair.setTotalIdx(n);
		Collections.sort(allWtdAttractors);
		
		WtdAttractorSet.setNames(files);
		ArrayList<DistPair> allDistPairs = new ArrayList<DistPair>();
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
				int ovlp = a.ovlapWith(b);
				if(ovlp > 0){
					allDistPairs.add(new DistPair(a, b , ovlp));
				}
			}
		}
		Collections.sort(allDistPairs);
		System.out.println(allDistPairs.size() + " pairs of distances have been added.");
		
		
		
		int cnt2 = 0;
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
						int ovlp = x.ovlapWith(b);
						if(ovlp > 0){
							allDistPairs.add(new DistPair(x, b , ovlp));
						}
					}
				}
				Collections.sort(allDistPairs);
			}
			cnt2++;
		}
		
		/*System.out.println("Making matches...");
		
		while(allWtdAttractors.size() > 0){
			WtdAttractor wa = allWtdAttractors.get(0);
			WtdAttractorSet was = new WtdAttractorSet(wa);
			//System.out.println(wa.name + "\t" + files[wa.source] + "\n");
			allWtdAttractors.remove(0);
			int[] ovlpCnt = new int[nf];
			ovlpCnt[wa.source] = Integer.MAX_VALUE;
			for(int j = allWtdAttractors.size()-1; j > 0; j --){
				WtdAttractor wa2 = allWtdAttractors.get(j);
				if(wa.source != wa2.source ){
					int ovlp = was.ovlapWith(wa2);
					if(ovlp > ovlpCnt[wa2.source]){
						//System.out.println(wa2.name + "\t" + files[wa2.source]);
						was.addWtd(wa2);
						ovlpCnt[wa2.source] = ovlp;
						allWtdAttractors.remove(j);
					}
				}
			}
			//System.out.println();
			for(WtdAttractor wa2 : row){
				if(wa2 == null) System.out.print("null\t");
				else System.out.print(wa2.name + "\t");
				allWtdAttractors.remove(wa2);
			}System.out.println();
			for(int gg : ovlpCnt){
				System.out.print(gg + "\t");
			}System.out.println();
			if(was.matchNumber > 1){
				out.add(was);
			}
		}*/
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

		System.out.println("Done.");
		
	}

}
