package worker;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import util.StatOps;

public class Converger extends DistributedWorker{
	private static double fdrThreshold = 0.05;
	private static float zThreshold = 3;
	private static int maxIter = 100;
	private static float corrThreshold = 0.7f;
	private static boolean rankBased = false;
	private static int attractorSize = 20;
	private static String convergeMethod = "FIXEDSIZE";
	
	public static class ValIdx implements Comparable<ValIdx>{
		float val;
		int idx;
		public ValIdx(int i, float v){
			this.idx = i;
			this.val = v;
		}
		public int hashCode(){
			return idx;
		}
		
		public int compareTo(ValIdx other) {
			return -Double.compare(this.val, other.val);
		}
		
		public int idx(){
			return idx;
		}
		public float val(){
			return val;
		}
		public boolean equals(Object other){
			boolean result = false;
	        if (other instanceof ValIdx) {
	        	ValIdx that = (ValIdx) other;
	            result = (this.idx == that.idx);
	        }
	        return result;
		}
		
	}
	
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
	private static float[] getMetaGene(float[][] data, ArrayList<ValIdx> idx, int n){
		int m = idx.size();
		float[] out = new float[n];
		for(int j = 0; j < n; j++){
			for(ValIdx vi : idx){
				out[j] += data[vi.idx][j];
			}
			out[j] /= m;
		}
		return out;
	}
	public ArrayList<ValIdx> findAttractor(float[][] data, int idx, int size) throws Exception{
		int m = data.length;
		int n = data[0].length;
		
		ITComputer itc = new ITComputer(7, 3, id, totalComputers);
		float[] mi = itc.getAllMIWith(data[idx], data);
		
		ArrayList<ValIdx> metaIdx = new ArrayList<ValIdx>();
		ValIdx[] vec = new ValIdx[m];
		for(int i = 0; i < m; i++){
			vec[i] = new ValIdx(i, mi[i]);
		}
		Arrays.sort(vec);
		for(int i = 0; i < size; i++){
			metaIdx.add(vec[i]);
		}
		int cnt = 0;
		ArrayList<ValIdx> preMetaIdx = new ArrayList<ValIdx>();
		preMetaIdx.addAll(metaIdx);
		
		while(cnt < maxIter){
			
			// cannot find significant associated genes, exit.
			
			if(metaIdx.size() == 0){
				//System.out.println("Empty set, exit.");
				break;
			}
			//System.out.print("Iteration " + cnt + "...");
			float[] metaGene = getMetaGene(data,metaIdx, n);
			mi = itc.getAllMIWith(metaGene, data);
			metaIdx = new ArrayList<ValIdx>();	
			vec = new ValIdx[m];
			for(int i = 0; i < m; i++){
				vec[i] = new ValIdx(i, mi[i]);
			}
			Arrays.sort(vec);
			for(int i = 0; i < size; i++){
				metaIdx.add(vec[i]);
			}
			if(preMetaIdx.equals(metaIdx)){
				break;
			}else{
				preMetaIdx = metaIdx;
				//System.out.println("Gene Set Size: " + metaIdx.size());
				cnt++;
			}
			
		}
		if(cnt == maxIter){
			System.out.println("Not converged.");
		}
		return metaIdx;
		
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
			ArrayList<ValIdx> metaIdx = new ArrayList<ValIdx>();
			ValIdx[] vec = new ValIdx[m];
			for(int i = 0; i < m; i++){
				vec[i] = new ValIdx(i, mi[i]);
			}
			if(convergeMethod.equals("FIXEDSIZE")){
				Arrays.sort(vec);
				for(int i = 0; i < attractorSize; i++){
					metaIdx.add(vec[i]);
				}
			}else if(convergeMethod.equals("ZSCORE")){
				float[] z = StatOps.xToZ(mi, m);
				for(int i = 0; i < m; i++){
					if(z[i] > zThreshold){
						metaIdx.add(vec[i]);
					}
				}
			}
			
			int cnt = 0;
			ArrayList<ValIdx> preMetaIdx = new ArrayList<ValIdx>();
			preMetaIdx.addAll(metaIdx);
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
				metaIdx = new ArrayList<ValIdx>();
				vec = new ValIdx[m];
				for(int i = 0; i < m; i++){
					vec[i] = new ValIdx(i, mi[i]);
				}
				if(convergeMethod.equals("FIXEDSIZE")){
					Arrays.sort(vec);
					for(int i = 0; i < attractorSize; i++){
						metaIdx.add(vec[i]);
					}
				}else if(convergeMethod.equals("ZSCORE")){
					float[] z = StatOps.xToZ(mi, m);
					metaIdx = new ArrayList<ValIdx>();
					for(int i = 0; i < m; i++){
						/*if(r[i] > corrThreshold){
							metaIdx.add(i);
						}
						if(padj[i] < fdrThreshold){
							metaIdx.add(i);
						}*/
						if(z[i] > zThreshold){
							metaIdx.add(vec[i]);
						}
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
			pw.print("\t" + 1);
			if(metaIdx.size() > 1){
				for(ValIdx vi: metaIdx){
						pw.print("\t" + vi.idx + "," + vi.val);
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
	public void setAttractorSize(int sz){
		Converger.attractorSize = sz;
	}
	public float getZThreshold(){
		return zThreshold;
	}
	
}
