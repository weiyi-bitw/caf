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
import java.util.HashMap;

import caf.GroupCNVWindow2.CNVWindow;

import obj.DataFile;
import obj.DataFileD;
import obj.Genome;
import obj.IntPair;
import obj.ValIdx;
import obj.ValIdxD;
import worker.Converger;
import worker.ITComputer;

public class GroupCNVWindow {
	static class Window implements Comparable<Window>{
		static ArrayList<String> genes;
		static Genome gn;
		static int quantile = 5;
		int startIdx;
		int centerIdx;
		ValIdx[] mis;
		String name;
		IntPair range;
		
		Window(String name, int startIdx, int centerIdx, ValIdx[] mis, IntPair range){
			this.name = name;
			this.startIdx = startIdx;
			this.centerIdx = centerIdx;
			this.mis = mis;
			this.range = range;
		}
		
		static Window parseWindow(String line, int startIdx, int readIn){
			String[] tokens = line.split("\t");
			String name = "Attractor_" + tokens[0];
			int nt = tokens.length;
			if(nt <= 3){
				return null;
			}
			ValIdx[] mis = new ValIdx[nt-2];
			for(int i = 2; i < nt; i++){
				String[] t2 = tokens[i].split(":");
				int idx = Integer.parseInt(t2[0]);
				mis[i-2] = new ValIdx(idx, Float.parseFloat(t2[1]));
			}
			Arrays.sort(mis);
			int x = Integer.MAX_VALUE;
			int y = -1;
			int k = Math.min(readIn, nt-2);
			for(int i = 0; i < k; i++){
				int idx = mis[i].idx;
				if(idx < x) x = idx;
				if(idx > y) y = idx;
			}
			return new Window(name, startIdx, gn.getIdx(genes.get(mis[0].idx)), mis, new IntPair(x, y));
		}
		
		public int hashCode(){
			return startIdx;
		}
		
		public String toString(int n){
			String s = name + "\t" + gn.getChr(genes.get(mis[0].idx)) + gn.getChrBand(genes.get(mis[0].idx));
			int k = Math.min(mis.length, n);
			for(int i = 0; i < k; i++){
				s += "\t" + genes.get(mis[i].idx) + ":" + mis[i].val;
			}
			s += "\t" + mis[quantile-1].val;
			return s;
			
		}
		
		public String toString(){
			String s = name + "\t" + gn.getChr(genes.get(startIdx));
			int n = mis.length;
			for(int i = 0; i < n; i++){
				s += "\t" + genes.get(mis[i].idx) + ":" + mis[i].val;
			}
			s += "\t" + gn.getChrBand(genes.get(startIdx)) + "\t" + mis[quantile-1].val;
			return s;
			
		}
		
		boolean ovlpWith(Window w){
			return this.range.overlapWith(w.range);
		}
		
		static void linkGenes(ArrayList<String> genes){
			Window.genes = genes;
		}
		static void linkGenome(Genome gn){
			Window.gn = gn;
		}
		static void setQuantile(int i){
			Window.quantile = i;
		}

		public int compareTo(Window other) {
			return -Double.compare(this.mis[quantile-1].val, other.mis[quantile-1].val);
		}
	}
	
	
	private static double[] getMetaGene(double[][] data, ArrayList<Integer> idx, int n){
		int m = idx.size();
		double[] out = new double[n];
		for(int j = 0; j < n; j++){
			for(Integer i : idx){
				out[j] += data[i][j];
			}
			out[j] /= m;
		}
		return out;
	}
	private static ArrayList<Integer> getAttractorIdx(double[][] data, double[][] val, ArrayList<Integer> inputIdx, int m, int n, ITComputer itc, int attractorSize, int maxIter) throws Exception{
		ArrayList<Integer> metaIdx = inputIdx;
		ArrayList<Integer> preMetaIdx = new ArrayList<Integer>();
		ArrayList<Integer> cycMetaIdx = new ArrayList<Integer>();
		
		preMetaIdx.addAll(inputIdx);
		cycMetaIdx.addAll(inputIdx);
		
		double[] metagene = getMetaGene(data, metaIdx, n);
		ValIdxD[] vec = new ValIdxD[m];
		int cnt = 0;
		
		
		while(cnt < maxIter){
			double[] mi = itc.getAllDoubleMIWith(metagene, val);
			for(int i = 0; i < m; i++){
				vec[i] = new ValIdxD(i, mi[i]);
			}
			Arrays.sort(vec);
			
			System.out.println(cnt + "\t" + vec[0].idx() + "\t" + vec[attractorSize/2].idx() + "\t" 
					+ "\t" + vec[attractorSize-1].idx());
			
			metaIdx.clear();
			for(int i = 0; i < attractorSize; i++){
				metaIdx.add(vec[i].idx());
			}
			if(metaIdx.equals(preMetaIdx) || metaIdx.equals(cycMetaIdx)){
				break;
			}
			cycMetaIdx.clear();
			cycMetaIdx.addAll(preMetaIdx);
			preMetaIdx.clear();
			preMetaIdx.addAll(metaIdx);
			cnt++;
			metagene = getMetaGene(data, metaIdx, n);
			
		}
		
		return metaIdx;
	}
	/*private static ArrayList<String> slidingWindowSelector(String inFileName, int winSize, int topAttractors) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(inFileName));
		ArrayList<String> genesAl = new ArrayList<String>();
		ArrayList<ValIdx> scoresAl = new ArrayList<ValIdx>();
		String line = br.readLine();
		int cnt = 0;
		
		while(line != null){
			String[] tokens = line.split("\t");
			genesAl.add(tokens[0]);
			scoresAl.add(new ValIdx(cnt,Float.parseFloat(tokens[1])));
			line = br.readLine();
			cnt ++;
		}
		br.close();
		
		int k = genesAl.size();
		Collections.sort(scoresAl);
		ArrayList<String> outG = new ArrayList<String>();
		boolean[] eliminate = new boolean[k];
		
		for(ValIdx vi : scoresAl){
			int idx = vi.idx;
			if(!eliminate[idx]){
				for(int i = (idx - winSize + 1); i <= (idx + winSize -1); i++){
					if(i >= 0 && i < k){
						eliminate[i] = true;
					}
				}
				outG.add(genesAl.get(idx));
				if(topAttractors > 0){
					if(outG.size() >= topAttractors){
						break;
					}
				}
			}
		}
		
		
		return outG;
	}*/
	
	private static ArrayList<Window> slidingWindowSelector(String inFileName, int excludeSize, int readIn) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(inFileName));
		ArrayList<Window> out = new ArrayList<Window>();
		
		String line = br.readLine();
		int cnt = 0;
		
		while(line != null){
			Window w = Window.parseWindow(line, cnt, readIn);
			if(w != null) out.add(w);
			line = br.readLine();
			cnt ++;
		}
		br.close();
		
		int k = out.size();
		
		Collections.sort(out);
		
		PrintWriter pw = new PrintWriter(new FileWriter("debug.txt"));
		for(Window w : out){
			pw.println(w.name + "\t" + w.centerIdx + "\t" + w.mis[4].val);
		}
		
		pw.close();
		
		for(int i = k -1; i >=0; i--){
			Window w = out.get(i);
			for(int j = 0; j < i; j++){
				//if(Math.abs(w.centerIdx - out.get(j).centerIdx) < excludeSize || w.ovlpWith(out.get(j))){
				if(w.ovlpWith(out.get(j))){
					out.remove(i);
					break;
				}
			}
		}
		return out;
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		final int IDX = 0;
		
		String[] dataFiles = {
				"/home/weiyi/workspace/data/brca/gse3494/ge.12160x251.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/brca/gse32646/ge.19190x115.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/brca/gse36771/ge.19190x107.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/brca/gse31448/ge.19190x353.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/brca/gse2034/ge.12160x286.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/brca/tcga/ge/ge.17475x536.ncbi.txt",
				"/home/weiyi/workspace/data/dream7/preTraining/train/ge.24940x500.mean.txt",
				"/home/weiyi/workspace/data/coad/gse14333/ge.19190x290.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/ov/gse9891/ge.19190x285.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/coad/tcga/ge/ge.17475x154.ncbi.txt",
				"/home/weiyi/workspace/data/ov/tcga/ge/ge.11963x582.ncbi.txt",
				//"/home/weiyi/workspace/data/gbm/tcga/ge/ge.12042x545.txt",
				"/home/weiyi/workspace/data/dream7/preTraining/train/cnv.21533x500.txt"
		};
		
		final String[] outputDirs={
				"brca.gse3494",
				"brca.gse32646",
				"brca.gse36771",
				"brca.gse31448",
				"brca.gse2034",
				"brca.tcga",
				"dream7",
				"coad.gse14333",
				"ov.gse9891",
				"coad.tcga",
				"ov.tcga",
				//"gbm.tcga",
				"dream7.cnv"
				
		};
		
		
		String outPath = "/home/weiyi/workspace/javaworks/caf/output/window/";
		
		final String geneLocFile = "/home/weiyi/workspace/data/annot/ncbi/gene.location.ncbi";
		int excludeSize = 25;
		int quantile = 5;
		int readIn = 15;
		
		ArrayList<Integer> allSizes = new ArrayList<Integer>();
		
		for(int qq = 0; qq < dataFiles.length; qq++)
		{
		
		System.out.println("Loading file " + dataFiles[qq]);
		DataFileD ma = DataFileD.parse(dataFiles[qq]);
		
		//ma.normalizeRows();
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		double[][] data = ma.getData();
		
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		//gn.linkToDataFile(ma);
		ma = ma.getSubProbes(gn.getAllGenes());
		
		Window.linkGenes(ma.getProbes());
		Window.linkGenome(gn);
		Window.setQuantile(quantile);
		
		System.out.println("Parsing windows...");
		ArrayList<Window> out = slidingWindowSelector(outPath + outputDirs[qq] + "/basinScores.txt", excludeSize, readIn);
		for(Window w : out){
			allSizes.add(w.mis.length);
		}
		System.out.println(out.size() + " CNVWindow selected.");
		
		String outFileName = outputDirs[qq];
		new File(outPath + "mergeroom").mkdirs();
		PrintWriter pw = new PrintWriter(new FileWriter(outPath + "mergeroom/" + outFileName));
		for(Window w : out){
			pw.println(w.toString(readIn));
		}
		
		pw.close();
		
				
		}
		System.out.println("Done.");
	}

}
