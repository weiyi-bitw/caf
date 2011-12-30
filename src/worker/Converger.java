package worker;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

import obj.Chromosome;
import obj.DataFile;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import util.StatOps;

public class Converger extends DistributedWorker{
	private static float zThreshold = 10;
	private static int maxIter = 100;
	private static boolean rankBased = false;
	private static int attractorSize = 10; // minimum size of an attractor
	private static String convergeMethod = "FIXEDSIZE";
	private static int bins = 7;
	private static int splineOrder = 3;
	
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
	// search ValIdx array by its INDICES!
	public static int biSearch(ValIdx[] x, ValIdx k){
		int out = -1;
		int n = x.length;
		int r = n;
		int l = 0;
		int m;
		while(r > l){
			System.out.println("(" + l + "\t" + r + ")");
			m = (r+l)/2;
			if(k.idx == x[m].idx){
				return m;
			}else if(k.idx > x[m].idx){
				l = m;
			}else if(k.idx < x[m].idx){
				r = m;
			}
		}
		return out;
	}
	
	// Sort ValIdx BY INDICES!!! Note: sorting by value can be directly applied Arrays.sort()
	public static ValIdx[] mergeSort(ValIdx[] x){
		int n = x.length;
		if(n <= 1) return x;
		
		ValIdx left[], right[],result[];
		int m = n/2;
		left = new ValIdx[m];
		right = new ValIdx[n - m];
		
		for(int i = 0; i < m; i++){
			left[i] = x[i];
		}
		for(int i = m; i < n; i++){
			right[i-m] = x[i];
		}
		left = mergeSort(left);
		right = mergeSort(right);
		
		result = merge(left, right);
		
		return result;
	}
	
	private static ValIdx[] merge(ValIdx[] left, ValIdx[] right){
		int nL = left.length;
		int nR = right.length;
		ValIdx[] result = new ValIdx[nL + nR];
		int iL = 0, iR = 0, iOut = 0;;
		while(iL < nL || iR < nR){
			if(iL < nL && iR < nR){
				if(left[iL].idx <= right[iR].idx){
					result[iOut] = left[iL];
					iL++; 
				}else{
					result[iOut] = right[iR];
					iR++;
				}
			}else if(iL < nL){
				result[iOut] = left[iL];
				iL++; 
			}else if(iR < nR){
				result[iOut] = right[iR];
				iR++;
			}
			iOut++;
		}
		return result;
	}
	
	public Converger(int id, int totalComputers, long jobID){
		super(id, totalComputers, jobID);
	}
	public Converger(int id, int totalComputers, long jobID, String method, int maxIter, boolean rankBased){
		super(id, totalComputers, jobID);
		Converger.maxIter = maxIter;
		Converger.rankBased = rankBased;
		Converger.convergeMethod = method;
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
	private static float[] getMetaGene(float[][] data, int start, int len, int n){
		float[] out = new float[n];
		for(int j = 0; j < n; j++){
			for(int l = start; l < (start + len); l++){
				out[j] += data[l][j];
			}
			out[j] /= len;
		}
		return out;
	}
	public ArrayList<ValIdx> findAttractor(float[][] data, int idx) throws Exception{
		int m = data.length;
		int n = data[0].length;
		
		ITComputer itc = new ITComputer(bins, splineOrder, id, totalComputers);
		float[] mi = itc.getAllMIWith(data[idx], data);
		ArrayList<ValIdx> metaIdx = new ArrayList<ValIdx>();
		ValIdx[] vec = new ValIdx[m];
		if(convergeMethod.equals("FIXEDSIZE")){
			
			for(int i = 0; i < m; i++){
				vec[i] = new ValIdx(i, mi[i]);
			}
			Arrays.sort(vec);
			for(int i = 0; i < attractorSize; i++){
				metaIdx.add(vec[i]);
			}
		}else if(convergeMethod.equals("ZSCORE")){
			float[] z = StatOps.xToZ(mi, m);
			for(int i = 0; i < m; i++){
				vec[i] = new ValIdx(i, mi[i]);
			}
			//Arrays.sort(vec);
			for(int i = 0; i < m; i++){
				if(z[i] > zThreshold){
					metaIdx.add(vec[i]);
				}
			}
			Collections.sort(metaIdx);
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
			System.out.print("Iteration " + cnt + "...");
			float[] metaGene = getMetaGene(data,metaIdx, n);
			mi = itc.getAllMIWith(metaGene, data);
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
				for(int i = 0; i < m; i++){
					vec[i] = new ValIdx(i, mi[i]);
				}
				//Arrays.sort(vec);
				for(int i = 0; i < m; i++){
					if(z[i] > zThreshold){
						metaIdx.add(vec[i]);
					}
				}
				Collections.sort(metaIdx);
			}
			if(preMetaIdx.equals(metaIdx)){
				System.out.println("Converged.");
				break;
			}else{
				preMetaIdx = metaIdx;
				System.out.println("Gene Set Size: " + metaIdx.size());
				cnt++;
			}
			
		}
		if(cnt == maxIter){
			System.out.println("Not converged.");
		}
		return metaIdx;
		
	}
	public void findCNV(DataFile ma, ArrayList<Chromosome> chrs, int winsize) throws Exception{
		int n = ma.getNumCols();
		
		int numTasks = 0;
		// Calculate total number of tasks
		for(Chromosome chr : chrs){
			int sz = chr.size() - winsize + 1;
			if (sz > 0){
				numTasks += sz;
			}
		}
		
		int start = id * numTasks / totalComputers;
		int end = (id+1) * numTasks / totalComputers;
		
		System.out.println("Processing task " + (start+1) + " to " + end);
		
		ITComputer itc = new ITComputer(bins, splineOrder, id, totalComputers);
		//itc.negateMI(true);
		prepare("geneset");
		PrintWriter pw = new PrintWriter(new FileWriter("tmp/" + jobID + "/geneset/caf." + String.format("%05d", id)+".txt"));
		int tt = 0;
		for(Chromosome chr : chrs){
			System.out.println("At chromosome " + chr.name() + "...");
			ArrayList<ValIdx> geneIdx = chr.geneIdx();
			// sort gene idx ascendantly
			Collections.sort(geneIdx);
			Collections.reverse(geneIdx);
			
			float[][] data = ma.getSubProbes(geneIdx).getData();
			int mm = data.length;
			
			System.out.println("Genes: " + mm);
			
			float[][] val = rankBased? new float[mm][n] : data;
			if(rankBased){
				for(int i = 0; i < mm; i++){
					val[i] = StatOps.rank(data[i]);
				}
			}
			
			for(int idx = 0; idx < mm-winsize+1; idx++){
				if(tt >= start && tt < end){
					System.out.println("Processing " + tt + "...");
					/*
					 * Step 1: find the genes that are significantly associated with the seed gene 
					 *         as the initial metagene
					 */
					float[] row = getMetaGene(data, idx, winsize, n);
					if(rankBased){
						row = StatOps.rank(row);
					}
					float[] mi = itc.getAllMIWith(row, val);
					ArrayList<ValIdx> metaIdx = new ArrayList<ValIdx>();
					
					float[] z = StatOps.xToZ(mi, mm);
					ValIdx[] vec = new ValIdx[mm];
					for(int i = 0; i < mm; i++){
						vec[i] = new ValIdx(i, z[i]);
					}
					Arrays.sort(vec);
					for(int i = 0; i < winsize; i++){
						metaIdx.add(vec[i]);
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
						vec = new ValIdx[mm];
						z = StatOps.xToZ(mi, mm);
						for(int i = 0; i < mm; i++){
							vec[i] = new ValIdx(i, z[i]);
						}
						Arrays.sort(vec);
						metaIdx = new ArrayList<ValIdx>();
						for(int i = 0; i < winsize; i++){
							metaIdx.add(vec[i]);
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
					pw.print(geneIdx.get(idx).idx());
					pw.print("\t" + 1);
					if(metaIdx.size() > 1){
						for(ValIdx vi: metaIdx){
								pw.print("\t" + geneIdx.get(vi.idx).idx() + "," + vi.val);
						}
					}else{
						pw.print("\tNA");
					}
					pw.println();
					
				}
				tt++;
			}
		}
		pw.close();
	}
		
	public void findAttractor(float[][] val, float[][] data) throws Exception{
		int m = val.length;
		int n = val[0].length;
		
		int start = id * m / totalComputers;
		int end = (id+1) * m / totalComputers;
		
		System.out.println("Processing gene " + (start+1) + " to " + end);
		
		ITComputer itc = new ITComputer(bins, splineOrder, id, totalComputers);
		//itc.negateMI(true);
		prepare("geneset");
		PrintWriter pw = new PrintWriter(new FileWriter("tmp/" + jobID + "/geneset/caf." + String.format("%05d", id)+".txt"));
		for(int idx = start; idx < end; idx++){
			System.out.println("Processing " + idx + "...");
			/*
			 * Step 1: find the genes that are significantly associated with the seed gene 
			 *         as the initial metagene
			 */
			
			float[] mi = itc.getAllMIWith(val[idx], val);
			ArrayList<ValIdx> metaIdx = new ArrayList<ValIdx>();
			ValIdx[] vec = new ValIdx[m];
			
			if(convergeMethod.equals("FIXEDSIZE")){
				for(int i = 0; i < m; i++){
					vec[i] = new ValIdx(i, mi[i]);
				}
				Arrays.sort(vec);
				for(int i = 0; i < attractorSize; i++){
					metaIdx.add(vec[i]);
				}
			}else if(convergeMethod.equals("ZSCORE")){
				float[] z = StatOps.xToZ(mi, m);
				for(int i = 0; i < m; i++){
					vec[i] = new ValIdx(i, z[i]);
				}
				Arrays.sort(vec);
				for(int i = 0; i < attractorSize; i++){
					metaIdx.add(vec[i]);
				}
				for(int i = attractorSize; i < m; i++){
					if(z[i] > zThreshold){
						metaIdx.add(vec[i]);
					}else{
						break;
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
				if(convergeMethod.equals("FIXEDSIZE")){
					for(int i = 0; i < m; i++){
						vec[i] = new ValIdx(i, mi[i]);
					}
					Arrays.sort(vec);
					for(int i = 0; i < attractorSize; i++){
						metaIdx.add(vec[i]);
					}
				}else if(convergeMethod.equals("ZSCORE")){
					float[] z = StatOps.xToZ(mi, m);
					metaIdx = new ArrayList<ValIdx>();
					for(int i = 0; i < m; i++){
						vec[i] = new ValIdx(i, z[i]);
					}
					Arrays.sort(vec);
					for(int i = 0; i < attractorSize; i++){
						metaIdx.add(vec[i]);
					}
					for(int i = attractorSize; i < m; i++){
						if(z[i] > zThreshold){
							metaIdx.add(vec[i]);
						}else{
							break;
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
	public void setConvergeMethos(String mthd){
		Converger.convergeMethod = mthd;
	}
	public void setMIParameter(int bins, int so){
		Converger.bins = bins;
		Converger.splineOrder = so;
	}
}
