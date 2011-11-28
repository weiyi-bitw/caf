package worker;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashSet;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import util.StatOps;

public class Converger extends DistributedWorker{
	private static double fdrThreshold = 0.05;
	private static float zThreshold = 3;
	private static int maxIter = 100;
	private static float corrThreshold = 0.7f;
	private static boolean rankBased = true;
	
	public Converger(int id, int totalComputers, long jobID){
		super(id, totalComputers, jobID);
	}
	public Converger(int id, int totalComputers, long jobID, double fdrTh, int maxIter, float corrTh, boolean rankBased){
		super(id, totalComputers, jobID);
		Converger.fdrThreshold = fdrTh;
		Converger.maxIter = maxIter;
		Converger.corrThreshold = corrTh;
		Converger.rankBased = rankBased;
	}
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
	public void findAttractor(float[][] val, float[][] data) throws Exception{
		int m = val.length;
		int n = val[0].length;
		
		int start = id * m / totalComputers;
		int end = (id+1) * m / totalComputers;
		
		System.out.println("Processing gene " + (start+1) + " to " + end);
		
		ITComputer itc = new ITComputer(7, 3, id, totalComputers);
		itc.negateMI(true);
		prepare("geneset");
		PrintWriter pw = new PrintWriter(new FileWriter("tmp/" + jobID + "/geneset/caf." + String.format("%05d", id)+".txt"));
		for(int idx = start; idx < end; idx++){
			
			/*
			 * Step 1: find the genes that are significantly associated with the seed gene 
			 *         as the initial metagene
			 */
			
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
			//System.out.println("Initial gene set size " + metaIdx.size() );
			
			/*
			 * Step 2: Calculate metagene, find the genes that have correlation exceeding the 
			 *         threshold as the new metagene
			 */
			
			while(cnt < maxIter){
				
				// cannot find significant associated genes, exit.
				
				if(metaIdx.size() == 0){
					//System.out.println("Empty set, exit.");
					break;
				}
				//System.out.print("Iteration " + cnt + "...");
				float[] metaGene = getMetaGene(data,metaIdx, n);
				if(rankBased){
					metaGene = StatOps.rank(metaGene);
				}
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
					/*System.out.println("Converged."); 
					System.out.println("Gene Set Size: " + metaIdx.size());
					*/break;
				}else{
					preMetaIdx = metaIdx;
					//System.out.println("Gene Set Size: " + metaIdx.size());
					cnt++;
				}
				
			}
			// first token: attractee index
			pw.print(idx);
			if(metaIdx.size() > 1){
				for(Integer i: metaIdx){
						pw.print("\t" + i);
				}
			}else{
				pw.print("\tNA");
			}
			pw.println();
		}
		pw.close();
		
	}
	
	public void setZThreshold(float z) throws MathException{
		Converger.zThreshold = z;
	}
	public void setZThreshold(int m) throws MathException{
		NormalDistributionImpl norm = new NormalDistributionImpl();
		double pth = 0.05/m;
		Converger.zThreshold = (float) -norm.inverseCumulativeProbability(pth);
	}
	public float getZThreshold(){
		return zThreshold;
	}
	
}
