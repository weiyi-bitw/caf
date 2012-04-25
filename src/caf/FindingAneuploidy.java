package caf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import obj.Genome;
import obj.IntPair;
import obj.ValString;

public class FindingAneuploidy {

	static class AneuploidyAttractor implements Comparable<AneuploidyAttractor>{
		static Genome gn;
		static int quantile = 50;
		String source;
		int basin;
		String chr;
		ValString[] genes;
		IntPair range;
		
		AneuploidyAttractor(String source, int basin, String chr, ValString[] genes, IntPair range){
			this.source = source;
			this.basin = basin;
			this.chr = chr;
			this.genes = genes;
			this.range = range;
		}
		
		public static AneuploidyAttractor parseAA(String line, String source){
			String[] tokens = line.split("\t");
			int basin = Integer.parseInt(tokens[1]);
			if(basin < 3){
				return null;
			}
			int nt = Math.min(tokens.length-2, quantile);
			ValString genes[] = new ValString[nt];
			String chr = null;
			int x = Integer.MAX_VALUE;
			int y = -1;
			for(int i = 2; i < nt+2; i++){
				String[] t2 = tokens[i].split(":");
				String g = t2[0];
				if(gn.contains(g)){
					int idx = gn.getIdx(g);
					if(idx < x) x = idx;
					if(idx > y) y = idx;
					
					if(chr != null){
						//if(!gn.getChr(g).equals(chr)) return null;
					}else{
						chr = gn.getChr(g);
					}
				}
				float mi = Float.parseFloat(t2[1]);
				genes[i-2] = new ValString(g, mi);
			}
			Arrays.sort(genes);
			IntPair range = new IntPair(x, y);
			if(range.range() < 150){
				//return null;
			}
			return new AneuploidyAttractor(source, basin, chr, genes, range);
		}
		
		public int compareTo(AneuploidyAttractor o) {	
			if(this.chr.equals(o.chr)){
				return -Double.compare(this.range.range(), o.range.range());
			}
			return this.chr.compareTo(o.chr);
		}
		
		static void linkGenome(Genome gn){
			AneuploidyAttractor.gn = gn;
		}
		
		static void setQuantile(int quantile){
			AneuploidyAttractor.quantile = quantile;
		}
		public String toString(){
			String s = source + "_" + chr + "\t" + "Strength: " + genes[quantile-1].val + "\tRange: " 
					 + range + "\n";
			int n = genes.length;
			s += "Genes";
			for(int i = 0; i < n; i++){
				s += "\t" + genes[i].s;
			}
			s += "\nChromosomal Band";
			for(int i = 0; i < n; i++){
				String g = genes[i].s;
				if(gn.contains(g)){
					s += "\t" + gn.getChrBand(g);
				}else{
					s += "\t" + "NA";
				}
			}
			s += "\n";
			return s;
			
		}
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		final String[] outs={
				//"brca.gse2034",
				"coad.gse14333",
				//"ov.gse9891",
				//"brca.tcga",
				//"coad.tcga",
				//"ov.tcga"
		};
		
		final String geneLocFile = "/home/weiyi/workspace/data/annot/affy/u133p2/gene.location4";
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		
		String path = "/home/weiyi/workspace/javaworks/caf/output/weighted/";
		if(!path.endsWith("/")){
			path += "/";
		}
		
		AneuploidyAttractor.linkGenome(gn);
		AneuploidyAttractor.setQuantile(50);
		ArrayList<AneuploidyAttractor> aaList = new ArrayList<AneuploidyAttractor>();
		
		
		for(int qq = 0; qq < outs.length; qq++){
			System.out.println("Processing " + outs[qq] + "...");
			BufferedReader br = new BufferedReader(new FileReader(path + "/mergeroom.300/" + outs[qq]));
			String line = br.readLine();
			while(line != null){
				AneuploidyAttractor aa = AneuploidyAttractor.parseAA(line, outs[qq]);
				if(aa != null){
					aaList.add(aa);
				}
				line = br.readLine();
			}
			br.close();
		}
		Collections.sort(aaList);
		System.out.println(aaList.size() + " aneuploidy attractor loaded.");
		
		PrintWriter pw = new PrintWriter(new FileWriter(path + "aneuploidy.txt"));
		for(AneuploidyAttractor aa : aaList){
			pw.println(aa);
		}
		pw.close();
		
	}

}
