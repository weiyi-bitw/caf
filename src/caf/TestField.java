package caf;

import java.util.ArrayList;
import java.util.HashSet;

import org.apache.commons.math.distribution.NormalDistributionImpl;

import obj.DataFile;

import util.StatOps;
import worker.Converger;
import worker.ITComputer;

public class TestField {
	private static float[] getMetaGene(float[][] data, HashSet<Integer> idx, int n){
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
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String path = "/home/weiyi/workspace/data/gbm/tcga/ge_mir_meth";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		System.out.println("Loading files...");
		DataFile ma = DataFile.parse(path + "ge.17814x278.knn.txt");
		ma.normalizeRows();
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		float[][] val = new float[m][n];
		for(int i = 0; i < m; i++){
			val[i] = StatOps.rank(data[i]);
		}
		
		int idx = ma.getRows().get("C3");
		ITComputer itc = new ITComputer(7, 3, 0, 1);
		itc.negateMI(true);
		NormalDistributionImpl norm = new NormalDistributionImpl();
		double pth = 0.05/m;
		float zThreshold = (float) -norm.inverseCumulativeProbability(pth);
		
		float[] mi = itc.getAllMIWith(val[idx], val);
		float[] z = StatOps.xToZ(mi, m);
		HashSet<Integer> metaIdx = new HashSet<Integer>();
		for(int i = 0; i < m; i++){
			if(z[i] > zThreshold){
				metaIdx.add(i);
			}
		}
		int cnt = 0;
		HashSet<Integer> preMetaIdx = new HashSet<Integer>();
		for(Integer i : metaIdx){
			preMetaIdx.add(i);
		}
		System.out.println("Initial gene set size " + metaIdx.size() );
		
		/*
		 * Step 2: Calculate metagene, find the genes that have correlation exceeding the 
		 *         threshold as the new metagene
		 */
		int maxIter = 100;
		while(cnt < maxIter){
			
			// cannot find significant associated genes, exit.
			
			if(metaIdx.size() == 0){
				//System.out.println("Empty set, exit.");
				break;
			}
			//System.out.print("Iteration " + cnt + "...");
			float[] metaGene = getMetaGene(data,metaIdx, n);
			metaGene = StatOps.rank(metaGene);
			mi = itc.getAllMIWith(metaGene, val);
			z = StatOps.xToZ(mi, m);
			metaIdx = new HashSet<Integer>();
			for(int i = 0; i < m; i++){
				/*if(r[i] > corrThreshold){
					metaIdx.add(i);
				}*/
				/*if(padj[i] < fdrThreshold){
					metaIdx.add(i);
				}*/
				if(z[i] > zThreshold){
					metaIdx.add(i);
				}
			}
			if(preMetaIdx.equals(metaIdx)){
				System.out.println("Converged."); 
				System.out.println("Gene Set Size: " + metaIdx.size());
				break;
			}else{
				preMetaIdx = metaIdx;
				System.out.println("Gene Set Size: " + metaIdx.size());
				cnt++;
			}
			
		}
		
		ArrayList<String> probeNames = ma.getProbes();
		for(Integer i : metaIdx){
			System.out.println(probeNames.get(i));
		}
		
		
	}

}
