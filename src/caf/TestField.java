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
import java.util.HashSet;

import org.apache.commons.math.distribution.NormalDistributionImpl;

import obj.Annotations;
import obj.DataFile;
import obj.Genome;
import obj.ValIdx;

import util.StatOps;
import worker.Converger;
import worker.ITComputer;

public class TestField {
	private static float[] getMetaGene(float[][] data, ArrayList<Integer> idx, int n){
		int m = idx.size();
		float[] out = new float[n];
		for(int j = 0; j < n; j++){
			for(Integer i : idx){
				out[j] += data[i][j];
			}
			out[j] /= m;
		}
		return out;
	}
	private static ArrayList<Integer> getAttractorIdx(float[][] data, float[][] val, ArrayList<Integer> inputIdx, int m, int n, ITComputer itc, int attractorSize, int maxIter) throws Exception{
		ArrayList<Integer> metaIdx = inputIdx;
		ArrayList<Integer> preMetaIdx = new ArrayList<Integer>();
		ArrayList<Integer> cycMetaIdx = new ArrayList<Integer>();
		
		preMetaIdx.addAll(inputIdx);
		cycMetaIdx.addAll(inputIdx);
		
		float[] metagene = getMetaGene(data, metaIdx, n);
		ValIdx[] vec = new ValIdx[m];
		int cnt = 0;
		
		
		while(cnt < maxIter){
			float[] mi = itc.getAllMIWith(metagene, val);
			for(int i = 0; i < m; i++){
				vec[i] = new ValIdx(i, mi[i]);
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
	private static ArrayList<String> slidingWindowSelector(String inFileName, int winSize) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(inFileName));
		ArrayList<String> genesAl = new ArrayList<String>();
		ArrayList<Float> scoresAl = new ArrayList<Float>();
		String line = br.readLine();
		while(line != null){
			String[] tokens = line.split("\t");
			genesAl.add(tokens[0]);
			scoresAl.add(Float.parseFloat(tokens[1]));
			line = br.readLine();
		}
		
		int k = genesAl.size();
		ArrayList<String> outG = new ArrayList<String>();
		for(int i = 0; i <= k-winSize; i+=20){
			int maxIdx = -1;
			float maxF = -1;
			for(int j = 0; j < winSize; j++){
				if(scoresAl.get(i+j) > maxF){
					maxF = scoresAl.get(i+j);
					maxIdx = (i+j);
				}
			}
			if(maxIdx == (i + winSize-1)){
				continue;
			}
			String g = genesAl.get(maxIdx);
			if(!outG.contains(g)){
				outG.add(g);
			}
		}
		br.close();
		return outG;
	}
	
	
	
/*	static class ValIdx implements Comparable<ValIdx>{
		float val;
		int idx;
		ValIdx(int i, float v){
			this.idx = i;
			this.val = v;
		}
		
		public int compareTo(ValIdx other) {
			return -Double.compare(this.val, other.val);
		}
	}*/
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String path = "/home/weiyi/workspace/data/ov/tcga/ge";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		System.out.println("Loading files...");
		
		//String outPath = "/home/weiyi/workspace/javaworks/caf/output/656/";
		String outPath = "/home/weiyi/workspace/javaworks/caf/tmp/";
		if(!outPath.endsWith("/")){
			outPath = outPath + "/";
		}
		DataFile ma = DataFile.parse(path + "ge.12042x582.txt");
		ArrayList<String> genes = ma.getProbes();
		
		final String geneLocFile = "/home/weiyi/workspace/data/annot/affy/u133p2/gene.location4";
		//final String geneLocFile = "/home/weiyi/workspace/javaworks/caf/output/639/gene.location3";
		
		//ma.normalizeRows();
		
		ArrayList<String> gs = new ArrayList<String>();
		gs.add("CYC1");
		gs.add("EXOSC4");
		gs.add("PSMD12");
		gs.add("AEBP1");
		
		long jobID = System.currentTimeMillis();
		
		//String annotPath = "/home/weiyi/workspace/data/annot/affy/u133p2/annot.csv";
		//Annotations annot = Annotations.parseAnnotations(annotPath);
		
		Converger cvg = new Converger(0, 1, jobID);
		ITComputer itc = new ITComputer(6, 3, 0, 1, true);
		cvg.linkITComputer(itc);
		
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		//if(command.equals("CNV")) gn.linkToDataFile(ma);
		for(String g : gs){
			if(genes.contains(g)){
				System.out.println("Processing " + g + " (" + gn.getChrArm(g) + ")" + "...");
				String[] neighbors = gn.getAllGenesInChrArm(gn.getChrArm(g));
				if(neighbors == null){
					System.out.println("No neighbors :(");
					break;
				}
				
				DataFile ma2 = ma.getSubProbes(neighbors);
				HashMap<String, Integer> rowmap = ma2.getRows();
				ArrayList<String> genes2 = ma2.getProbes();
				int idx = rowmap.get(g);
				ValIdx[] out = cvg.findWeightedAttractorOptimizePower(ma2, idx, 1f, 6, 0.5f, 5);
				
				if(out == null){
					System.out.println("No legit attractor.");
					continue;
				}
				
				new File("tmp").mkdir();
				PrintWriter pw = new PrintWriter(new FileWriter("tmp/CNV_" + g + ".txt"));
				for(int i = 0; i < out.length; i++){
					pw.println(genes2.get(out[i].idx) + "\t" + out[i].val);
				}
				pw.close();
			}else{
				System.out.println("Does not contain gene " + g);
			}
			
		}
		
		System.out.println("Done.");
	}

}
