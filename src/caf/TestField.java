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
		String path = "/home/weiyi/workspace/data/brca/gse2034/";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		System.out.println("Loading files...");
		
		//String outPath = "/home/weiyi/workspace/javaworks/caf/output/656/";
		String outPath = "/home/weiyi/workspace/javaworks/caf/tmp/";
		if(!outPath.endsWith("/")){
			outPath = outPath + "/";
		}
		/*DataFile ma = DataFile.parse(path + "ge.12160x286.jetset.mean.txt");
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();*/
		
		final String geneLocFile = "/home/weiyi/workspace/data/annot/affy/u133p2/gene.location4";
		//final String geneLocFile = "/home/weiyi/workspace/javaworks/caf/output/639/gene.location3";
		
		String command = "CNV";
		float power = 5f;
		boolean excludeTop = false;
		boolean miDecay = false;
		int winSize = 51;
		
		//ma.normalizeRows();
		
		ArrayList<String> gs = new ArrayList<String>();
		/*BufferedReader br = new BufferedReader(new FileReader("COL11A1_50"));
		br.readLine();
		String line = br.readLine();
		while(line != null){
			String[] tokens = line.split("\t");
			gs.add(tokens[0]);
			line = br.readLine();
		}
		br.close();*/
		
		long jobID = System.currentTimeMillis();
		
		//String annotPath = "/home/weiyi/workspace/data/annot/affy/u133p2/annot.csv";
		//Annotations annot = Annotations.parseAnnotations(annotPath);
		Annotations annot = null;
		
		Converger cvg = new Converger(0, 1, jobID);
		
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		String[] chr17q = gn.getAllGenesInChrArm("chr17q");
		
		for(String s : chr17q){
			System.out.println(s);
		}
		//if(command.equals("CNV")) gn.linkToDataFile(ma);
		
		
		
		System.out.println("Done.");
	}

}
