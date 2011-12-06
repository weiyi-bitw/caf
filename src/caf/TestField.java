package caf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import org.apache.commons.math.distribution.NormalDistributionImpl;

import obj.DataFile;

import util.StatOps;
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
			
			System.out.println(cnt + "\t" + vec[0].idx + "\t" + vec[attractorSize/2].idx + "\t" 
					+ "\t" + vec[attractorSize-1].idx);
			
			metaIdx.clear();
			for(int i = 0; i < attractorSize; i++){
				metaIdx.add(vec[i].idx);
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
	static class ValIdx implements Comparable<ValIdx>{
		float val;
		int idx;
		ValIdx(int i, float v){
			this.idx = i;
			this.val = v;
		}
		
		public int compareTo(ValIdx other) {
			return -Double.compare(this.val, other.val);
		}
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String path = "/home/weiyi/workspace/data/ov/tcga/mergeroom";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		System.out.println("Loading files...");
		DataFile ma = DataFile.parse(path + "ge.17814x514.common.txt");
		//ma.normalizeRows();
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		float[][] val = new float[m][n];
		ArrayList<Integer> metaIdx = new ArrayList<Integer>();
		for(int i = 0; i < m; i++){
			System.arraycopy(data[i], 0, val[i], 0, n);
			metaIdx.add(i);
		}
		ITComputer itc = new ITComputer(7, 3, 0, 1);
		itc.negateMI(true);
		float[] allMeta = getMetaGene(data, metaIdx,n);
		float[] mi = itc.getAllMIWith(allMeta, val);
		ValIdx[] vec = new ValIdx[m];
		for(int i = 0; i < m; i++){
			vec[i] = new ValIdx(i, mi[i]);
		}
		Arrays.sort(vec);
		int idx = 16831;
		metaIdx.clear();
		metaIdx.add(idx);
		int maxIter = 100;
		ArrayList<Integer> topMeta = new ArrayList<Integer>();
		topMeta.addAll(getAttractorIdx(data, val, metaIdx, m, n, itc, m/4, maxIter));
		
		/*for(int i = 0; i < m; i++){
			val[i] = StatOps.rank(data[i]);
		}*/
		
		
		ArrayList<String> probeNames = ma.getProbes();
		System.out.println("Top meta genes");
		for(int i = 0; i < 20; i++){
			System.out.println(topMeta.get(i) + "\t" + probeNames.get(topMeta.get(i)));
		}
		
	}

}
