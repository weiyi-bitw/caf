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
	public void findCNV(float[][] data, float[][] val, ArrayList<Chromosome> chrs, float zth) throws Exception{
		int n = data[0].length;
		int m = data.length;
		
		int mm = 0;
		for(Chromosome chr : chrs){
			mm += chr.size();
		}
		int start = id * mm / totalComputers;
		int end = (id+1) * mm / totalComputers;
		
		
		System.out.println("Processing task " + (start+1) + " to " + end);
		
		ITComputer itc = new ITComputer(bins, splineOrder, id, totalComputers);
		//itc.negateMI(true);
		prepare("geneset");
		PrintWriter pw = new PrintWriter(new FileWriter("tmp/" + jobID + "/geneset/caf." + String.format("%05d", id)+".txt"));
		int tt = 0;
		for(Chromosome chr : chrs){
			System.out.println("At chromosome " + chr.name() + "...");
			ArrayList<ValIdx> geneIdx = chr.geneIdx();
			/*// sort gene idx ascendantly
			Collections.sort(geneIdx);
			Collections.reverse(geneIdx);
			*/
			int k = geneIdx.size();
			
			for(int i = 0; i < k; i++){
				if(tt >= start && tt < end){
					int idx = geneIdx.get(i).idx;
					System.out.print("Processing " + tt + "...");
					/*
					 * Step 1: find the genes that are significantly associated with the seed gene 
					 *         as the initial metagene
					 */
					float[] mi = itc.getAllMIWith(val[idx], val);
					ArrayList<ValIdx> metaIdx = new ArrayList<ValIdx>();
					
					float[] z = StatOps.xToZ(mi, m);
					ValIdx[] vec = new ValIdx[m];
					for(int j = 0; j < m; j++){
						vec[j] = new ValIdx(j, z[j]);
					}
					Arrays.sort(vec);
					for(int j = 0; j < m; j++){
						if(vec[j].val > zth){
							if(geneIdx.contains(vec[j])){
								metaIdx.add(vec[j]);
							}
						}else{
							break;
						}
					}
					
					int cnt = 0;
					ArrayList<ValIdx> prepreMetaIdx = new ArrayList<ValIdx>();
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
						z = StatOps.xToZ(mi, m);
						for(int j = 0; j < m; j++){
							vec[j] = new ValIdx(j, z[j]);
						}
						Arrays.sort(vec);
						metaIdx = new ArrayList<ValIdx>();
						for(int j = 0; j < m; j++){
							if(vec[j].val > zth){
								if(geneIdx.contains(vec[j])){
									metaIdx.add(vec[j]);
								}
							}else{
								break;
							}
						}
						if(preMetaIdx.equals(metaIdx)){
							System.out.print("Converged. "); 
							System.out.println("Gene Set Size: " + metaIdx.size());
							break;
						}else if (prepreMetaIdx.equals(metaIdx)){
							System.out.println("Cycled.");
							if(metaIdx.size() >= preMetaIdx.size()){
								break;
							}else{
								metaIdx = preMetaIdx;
								break;
							}
						}
						else{
							prepreMetaIdx = preMetaIdx;
							preMetaIdx = metaIdx;
							//System.out.println("Gene Set Size: " + metaIdx.size());
							cnt++;
						}
						
					}
					if(cnt == maxIter){
						System.out.println("Not converged.");
					}
					// first token: attractee index
					pw.print(geneIdx.get(i).idx());
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
				tt++;
			}
		}
		pw.close();
	}
	public ArrayList<ValIdx> findAttractor(float[][] data, float[] vec)throws Exception{
		int m = data.length;
		int n = data[0].length;
		ITComputer itc = new ITComputer(bins, splineOrder, id, totalComputers);
			/*
			 * Step 1: find the genes that are significantly associated with the seed gene 
			 *         as the initial metagene
			 */
			
			float[] mi = itc.getAllMIWith(vec, data);
			ArrayList<ValIdx> metaIdx = new ArrayList<ValIdx>();
			ValIdx[] vecMI = new ValIdx[m];
			for(int i = 0; i < m; i++){
				vecMI[i] = new ValIdx(i, mi[i]);
			}
			ValIdx[] vecZ = new ValIdx[m];
			
			if(convergeMethod.equals("FIXEDSIZE")){
				Arrays.sort(vecMI);
				for(int i = 0; i < attractorSize; i++){
					metaIdx.add(vecMI[i]);
				}
			}else if(convergeMethod.equals("ZSCORE")){
				float[] z = StatOps.xToZ(mi, m);
				for(int i = 0; i < m; i++){
					vecZ[i] = new ValIdx(i, z[i]);
				}
				Arrays.sort(vecZ);
				/*for(int i = 0; i < attractorSize; i++){
					metaIdx.add(vecZ[i]);
				}*/
				for(int i = 0; i < m; i++){
					if(vecZ[i].val() > zThreshold){
						metaIdx.add(vecZ[i]);
					}else{
						break;
					}
				}
			}
			int cnt = 0;
			ArrayList<ValIdx> prepreMetaIdx = new ArrayList<ValIdx>();
			ArrayList<ValIdx> preMetaIdx = new ArrayList<ValIdx>();
			preMetaIdx.addAll(metaIdx);
			System.out.println("Initial gene set size " + metaIdx.size() );
			
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
				mi = itc.getAllMIWith(metaGene, data);
				vecMI = new ValIdx[m];
				for(int i = 0; i < m; i++){
					vecMI[i] = new ValIdx(i, mi[i]);
				}
				metaIdx = new ArrayList<ValIdx>();
				if(convergeMethod.equals("FIXEDSIZE")){
					Arrays.sort(vecMI);
					for(int i = 0; i < attractorSize; i++){
						metaIdx.add(vecMI[i]);
					}
				}else if(convergeMethod.equals("ZSCORE")){
					vecZ = new ValIdx[m];
					float[] z = StatOps.xToZ(mi, m);
					for(int i = 0; i < m; i++){
						vecZ[i] = new ValIdx(i, z[i]);
					}
					Arrays.sort(vecZ);
					/*for(int i = 0; i < attractorSize; i++){
						metaIdx.add(vecZ[i]);
					}*/
					for(int i = 0; i < m; i++){
						if(vecZ[i].val() > zThreshold){
							metaIdx.add(vecZ[i]);
						}else{
							break;
						}
					}
				}
				if(preMetaIdx.equals(metaIdx)){
					System.out.println("Converged."); 
					System.out.println("Gene Set Size: " + metaIdx.size());
					break;
				}else if (prepreMetaIdx.equals(metaIdx)){
					System.out.println("Cycled.");
					if(metaIdx.size() >= preMetaIdx.size()){
						break;
					}else{
						metaIdx = preMetaIdx;
						break;
					}
				}
				else{
					prepreMetaIdx = preMetaIdx;
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
			System.out.print("Processing " + idx + "...");
			/*
			 * Step 1: find the genes that are significantly associated with the seed gene 
			 *         as the initial metagene
			 */
			
			float[] mi = itc.getAllMIWith(val[idx], val);
			ArrayList<ValIdx> metaIdx = new ArrayList<ValIdx>();
			ValIdx[] vecMI = new ValIdx[m];
			for(int i = 0; i < m; i++){
				vecMI[i] = new ValIdx(i, mi[i]);
			}
			ValIdx[] vecZ = new ValIdx[m];
			
			if(convergeMethod.equals("FIXEDSIZE")){
				Arrays.sort(vecMI);
				for(int i = 0; i < attractorSize; i++){
					metaIdx.add(vecMI[i]);
				}
			}else if(convergeMethod.equals("ZSCORE")){
				float[] z = StatOps.xToZ(mi, m);
				for(int i = 0; i < m; i++){
					vecZ[i] = new ValIdx(i, z[i]);
				}
				Arrays.sort(vecZ);
				for(int i = 0; i < attractorSize; i++){
					metaIdx.add(vecZ[i]);
				}
				for(int i = attractorSize; i < m; i++){
					if(vecZ[i].val() > zThreshold){
						metaIdx.add(vecZ[i]);
					}else{
						break;
					}
				}
			}
			int cnt = 0;
			ArrayList<ValIdx> prepreMetaIdx = new ArrayList<ValIdx>();
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
				vecMI = new ValIdx[m];
				for(int i = 0; i < m; i++){
					vecMI[i] = new ValIdx(i, mi[i]);
				}
				metaIdx = new ArrayList<ValIdx>();
				if(convergeMethod.equals("FIXEDSIZE")){
					Arrays.sort(vecMI);
					for(int i = 0; i < attractorSize; i++){
						metaIdx.add(vecMI[i]);
					}
				}else if(convergeMethod.equals("ZSCORE")){
					vecZ = new ValIdx[m];
					float[] z = StatOps.xToZ(mi, m);
					for(int i = 0; i < m; i++){
						vecZ[i] = new ValIdx(i, z[i]);
					}
					Arrays.sort(vecZ);
					for(int i = 0; i < attractorSize; i++){
						metaIdx.add(vecZ[i]);
					}
					for(int i = attractorSize; i < m; i++){
						if(vecZ[i].val() > zThreshold){
							metaIdx.add(vecZ[i]);
						}else{
							break;
						}
					}
				}
				if(preMetaIdx.equals(metaIdx)){
					System.out.println("Converged."); 
					//System.out.println("Gene Set Size: " + metaIdx.size());
					break;
				}else if (prepreMetaIdx.equals(metaIdx)){
					System.out.println("Cycled.");
					if(metaIdx.size() >= preMetaIdx.size()){
						break;
					}else{
						metaIdx = preMetaIdx;
						break;
					}
				}
				else{
					prepreMetaIdx = preMetaIdx;
					preMetaIdx = metaIdx;
					//System.out.println("Gene Set Size: " + metaIdx.size());
					cnt++;
				}
				
			}
			if(cnt == maxIter){
				System.out.println("Not converged.");
			}
			// first token: attractee index
			pw.print(idx);
			pw.print("\t" + 1);
			if(metaIdx.size() > 1){
				if(convergeMethod.equals("ZSCORE")){
					for(ValIdx vi: metaIdx){
						pw.print("\t" + vi.idx + "," + vecMI[vi.idx].val + "," + vi.val);
					}
				}else{
					for(ValIdx vi: metaIdx){
						pw.print("\t" + vi.idx + "," + vi.val + ",NaN");
					}
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
