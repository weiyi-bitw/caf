package worker;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import obj.Annotations;
import obj.Chromosome;
import obj.DataFile;
import obj.DataFileD;
import obj.Genome;
import obj.InverseAnnotations;
import obj.ValIdx;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import util.StatOps;

public class Converger extends DistributedWorker{
	private static float zThreshold = 10;
	private static int maxIter = 100;
	private static boolean rankBased = false;
	private static int attractorSize = 10; // minimum size of an attractor
	private static String convergeMethod = "WEIGHTED";
	private static int bins = 6;
	private static int splineOrder = 3;
	private static boolean miNorm = false;
	private static double precision = (float) 1E-4;
	private static double epsilon = 5E-14;
	static ITComputer itc;
	
	
	/*public static class ValIdx implements Comparable<ValIdx>{
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
		
	}*/
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
	public Converger(int id, int totalComputers, long jobID, int maxIter, boolean rankBased){
		super(id, totalComputers, jobID);
		Converger.maxIter = maxIter;
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
	private static double sigmoid(double x, double a){
		return (1/(1 + Math.exp(-2 * Math.PI * a * x)));
	}
	private static float[] getWeightedMetaGene(float[][] data, float[] w, float power, int m, int n){
		float[] out = new float[n];
		double sum = 0;
		//float[] z = StatOps.xToZ(w, m);
		//float[] r = StatOps.rank(w);
		for(int i = 0; i < m; i++){
			//if( (z[i]) > 0){
			if(w[i] > 0){
				//double ww = Math.exp(power*Math.log(w[i]));
				double f = Math.exp(power*Math.log(w[i]));
				//double sig =  Math.exp(power * Math.log(sigmoid(2*w[i]-1, 1.0)));
				//double f = w[i] * sig;
				//double f = sigmoid(2*w[i]-1, power);
				//double f = sigmoid(2*r[i]/m - 1, power);
				//double f = Math.tan(w[i] * Math.PI/2);
				sum += f;
				for(int j = 0; j < n; j++){
					out[j] += data[i][j] * f;
				}
			}
		}
		for(int j = 0; j < n; j++){
			out[j] /= sum;
		}
		return out;
	}
	private static double[] getWeightedMetaGene(double[][] data, double[] w, double power, int m, int n){
		double[] out = new double[n];
		double sum = 0;
		for(int i = 0; i < m; i++){
			if(w[i] > 0){
				double f = Math.exp(power*Math.log(w[i]));
				sum += f;
				for(int j = 0; j < n; j++){
					out[j] += data[i][j] * f;
				}
			}
		}
		for(int j = 0; j < n; j++){
			out[j] /= sum;
		}
		return out;
	}
	
	public static float[] getWeightedMetaGene(float[][] data, double[] w, double power, int m, int n){
		double[] tmp = new double[n];
		double sum = 0;
		for(int i = 0; i < m; i++){
			if(w[i] > 0){
				double f = Math.exp(power*Math.log(w[i]));
				sum += f;
				for(int j = 0; j < n; j++){
					tmp[j] += (double)data[i][j] * f;
				}
			}
		}
		
		float[] out = new float[n];
		for(int j = 0; j < n; j++){
			out[j] = (float) (tmp[j]/sum);
		}
		return out;
	}
	private static float[] getChrWeightedMetaGene(float[][] data, float[] w, ArrayList<String> genes, 
			String chr, Genome gn, float power, int m, int n){
		float[] out = new float[n];
		double sum = 0;
		//float[] z = StatOps.xToZ(w, m);
		//float[] r = StatOps.rank(w);
		for(int i = 0; i < m; i++){
			String g = genes.get(i);
			if(gn.contains(g)){
				if(gn.getChr(g).equals(chr)){
				//if( (z[i]) > 8){
				//if(w[i] > 0){
					//double ww = Math.exp(power*Math.log(w[i]));
					double f = Math.exp(power*Math.log(w[i]));
					//double sig =  Math.exp(power * Math.log(sigmoid(2*w[i]-1, 1.0)));
					//double f = w[i] * sig;
					//double f = sigmoid(2*w[i]-1, power);
					//double f = sigmoid(2*r[i]/m - 1, power);
					//double f = Math.tan(w[i] * Math.PI/2);
					sum += f;
					for(int j = 0; j < n; j++){
						out[j] += data[i][j] * f;
					}
				}
			}
		}
		for(int j = 0; j < n; j++){
			out[j] /= sum;
		}
		return out;
	}
	private static float[] getSmoothedWeightedMetaGene(float[][] data, float[] w, float[] pw, float power, int m, int n){
		float[] out = new float[n];
		double sum = 0;
		float vw = StatOps.var(w, m);
		float a= vw / (vw + StatOps.var(pw, m));
		for(int i = 0; i < m; i++){
			if(w[i] > 0){
				w[i] = a * w[i] + (1-a)*pw[i];
				double f = (float) Math.exp(power*Math.log(w[i]));
				sum += f;
				for(int j = 0; j < n; j++){
					out[j] += data[i][j] * f;
				}
			}
		}
		for(int j = 0; j < n; j++){
			out[j] /= sum;
		}
		return out;
	}
	public static float calcMSE(float[] a, float[] b, int n){
		float err = 0;
		for(int i = 0; i < n; i++){
			err += (a[i] - b[i]) * (a[i] - b[i]);
		}
		return err / n;
	}
	public static double calcMSE(double[] a, double[] b, int n){
		double err = 0;
		for(int i = 0; i < n; i++){
			err += (a[i] - b[i]) * (a[i] - b[i]);
		}
		return err / n;
	}
	public static boolean equal(float[] a, float[] b, int n, float delta){
		for(int i = 0; i < n; i++){
			if(Math.abs(a[i] - b[i]) > delta){
				//System.out.println(Math.abs(a[i] - b[i]));
				return false;
			}
		}
		return true;
	}
	public static boolean identical(double[] w1, double[] w2, int n, double precision){
		for(int i = 0; i < n ;i++){
			if(Math.abs(w1[i] - w2[i]) > precision) return false;
		}
		return true;
	}
	
	public void findWeightedCNVCoef(DataFileD ma, Genome gn, int wstart, int wend, int delw, double pstart, double pend, double delp, int quantile) throws Exception{
		ma = ma.getSubProbes(gn.getAllGenes());
		
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		ArrayList<String> genes = ma.getProbes();
		HashMap<String, Integer> maMap = ma.getRows();
		
		int start = id * m / totalComputers;
		int end = (id+1) * m / totalComputers;
		
		System.out.println("Processing gene " + (start+1) + " to " + end);
		new File("output").mkdir();
		new File("output/" + jobID).mkdir();
		PrintWriter pw = new PrintWriter("output/" + jobID + "/basinScores." + String.format("%05d", id)+ ".txt");
		
		for(int idx = start; idx < end; idx++)
		{
			double bestScore = -1;
			int bestWinSize = -1;
			double bestExp = -1;
			ValIdx[] bestVec = null;
			String g = genes.get(idx);
			System.out.print("Processing " + g + "..."); 
			
			for(int winSize = wstart; winSize <= wend; winSize += delw)
			{
				
				String[] neighbors = gn.getNeighbors(g, winSize);
				if(neighbors == null){
					System.out.println("No neighbors :(");
					break;
				}
				
				DataFileD ma2 = ma.getSubProbes(neighbors);
				ArrayList<String> ma2Genes = ma2.getProbes();
				int m2 = ma2.getNumRows();
				if(m2 < quantile){
					continue;
				}
				double[][] data = ma2.getData();
				int idx2 = ma2.getRows().get(g);
				double[] vec = data[idx2];
				double convergeTh = precision * precision /m2;
				
				for(double power = pstart; power <= pend; power += delp)
				{
					
					double[] wVec = itc.getAllDoubleMIWith(vec, data);
					double[] preWVec = new double[m2];
					System.arraycopy(wVec, 0, preWVec, 0, m2);
					int c = 0;
					boolean converge = false;
					double score = -1;
					
					while(c < maxIter){
						double[] metaGene = getWeightedMetaGene(data, wVec, power,  m2, n);
						wVec = itc.getAllDoubleMIWith(metaGene, data);
						
						double err = calcMSE(wVec, preWVec, m2);
						System.arraycopy(wVec, 0, preWVec, 0, m2);
						
						if(err < convergeTh){
							Arrays.sort(preWVec);
							score = preWVec[m2 - quantile];
							converge = true;
							break;
						}
						
						c++;
					}
					
					if(converge && score > bestScore){
						bestScore = score;
						bestWinSize = winSize;
						bestVec = new ValIdx[m2];
						for(int i = 0; i < m2; i++){
							bestVec[i] = new ValIdx(maMap.get(ma2Genes.get(i)), (float)wVec[i]);
						}
					}
					
					
				} // END power iteration
			
			} // END winSize iteration
			
			System.out.println(bestScore);
			
			if(bestScore < 0){
				continue;
			}
			
			String chr = gn.getChr(g);
			
			pw.print(g + "\t" + chr);
			for(int i = 0; i < bestVec.length; i++){
				pw.print("\t" + bestVec[i].idx + ":" + bestVec[i].val);
			}pw.println();
			
		}// END idx iteration
		pw.close();
		
	}

	public double[] findWeightedAttractorDouble(DataFile ma, float[] vec, double power) throws Exception{
		float[][] data = ma.getData();
		int m = data.length;
		int n = data[0].length;
		//ArrayList<String> genes = ma.getProbes();
		
		double[] wVec = itc.getAllDoubleMIWith(vec, data);
		
		double[] preWVec = new double[m];
		System.arraycopy(wVec, 0, preWVec, 0, m);
		int c = 0;
		double convergeTh = epsilon;
		
		while(c < maxIter){
			float[] metaGene = getWeightedMetaGene(data, wVec, power,  m, n);
			wVec = itc.getAllDoubleMIWith(metaGene, data);
			
			double err = calcMSE(wVec, preWVec, m);
			System.out.println(err);
			if(err < convergeTh){
				System.out.println("Converged.");
				return wVec;
			}
			System.arraycopy(wVec, 0, preWVec, 0, m);
			c++;
		}
		System.out.println("Not converged.");
		wVec[0] = -1;
		return wVec;
	}
	
	public void findWeightedAttractor(DataFileD ma, double power) throws Exception{
		double[][] data = ma.getData();
		int m = data.length;
		int n = data[0].length;
		
		int start = id * m / totalComputers;
		int end = (id+1) * m / totalComputers;
		
		System.out.println("Processing gene " + (start+1) + " to " + end);
		
		ArrayList<double[]> wVecs = new ArrayList<double[]>();
		ArrayList<ArrayList<Integer>> basins = new ArrayList<ArrayList<Integer>>();
		
		for(int idx = start; idx < end; idx++){
			System.out.print("Processing " + ma.getProbes().get(idx) + "( " + idx + " )" + "...");
			double[] wVec = itc.getAllDoubleMIWith(data[idx], data);
			//float[] wVec = StatOps.pearsonCorr(vec, data, m, n);
			//float[] wVec = StatOps.cov(vec, data, m, n);
			double[] preWVec = new double[m];
			System.arraycopy(wVec, 0, preWVec, 0, m);
			int c = 0;
			boolean converge = false;
			while(c < maxIter){
				double[] metaGene = getWeightedMetaGene(data, wVec, power,  m, n);
				wVec = itc.getAllDoubleMIWith(metaGene, data);
				
				double err = calcMSE(wVec, preWVec, m);
				System.arraycopy(wVec, 0, preWVec, 0, m);
				if(err < epsilon){
					Arrays.sort(preWVec);
					if(preWVec[m-1] - preWVec[m-2] > 0.5){
						System.out.println("Top dominated.");
						converge=false;
					}else{
						System.out.println("Converged.");
						converge=true;
					}
					break;
				}
				
				c++;
			}
			if(converge){
				boolean newOne = true;
				for(int i = 0; i < wVecs.size(); i++){
					double[] fs = wVecs.get(i);
					if(identical(fs, wVec, m, precision)){ 
						newOne = false;
						basins.get(i).add(idx);
						break;
					}
				}
				if(newOne){
					wVecs.add(wVec);
					ArrayList<Integer> basin = new ArrayList<Integer>();
					basin.add(idx);
					basins.add(basin);
				}
			}
			
		}
		
		prepare("geneset");
		PrintWriter pw = new PrintWriter(new FileWriter("tmp/" + jobID + "/geneset/caf." + String.format("%05d", id)+".txt"));
		for(int i = 0; i < wVecs.size(); i++){
			pw.print("caf"); // add tag (for format consistency with CNV)
			ArrayList<Integer> basin = basins.get(i);
			int k = basin.size();
			for(int j = 0; j < k; j++){
				if(j == 0){
					pw.print("\t" + basin.get(j));
				}else{
					pw.print("," + basin.get(j));
				}
			}
			double[] fs = wVecs.get(i);
			for(int j = 0; j < m; j++){
				pw.print("\t" + fs[j]);
			}
			pw.println();
		}
		pw.close();
	}
	
	private double[] getWeightVector(double[] mi, int m, double a){
		double sum = 0;
		double[] w = new double[m];
		
		for(int i = 0; i < m; i++){
			if(mi[i] > 0){
				w[i] = Math.exp(a*Math.log(mi[i]));
				sum += w[i];
			}else{
				w[i] = 0;
			}
		}
		for(int i = 0; i < m; i++){
			w[i] /= sum;
		}
		return w;
	}
	
	private float[] getMetaGene(float[][] data, double[] w, int m, int n){
		float[] outf = new float[n];
		
		for(int j =0 ; j < n; j++){
			double o = 0;
			for(int i = 0; i < m; i++){
				o += w[i] * data[i][j];
			}
			outf[j] = (float) o;
		}
		
		return outf;
	}
	
	private double getScore(double[] mi,int m, int pos){
		double[] garbage = new double[m];
		System.arraycopy(mi, 0, garbage, 0, m);
		Arrays.sort(garbage);
		return garbage[m-pos];
	}
	
	public double[] AttractorScanning(DataFile ma, float[] vec, double powerStart, double powerEnd, double[] outPower) throws Exception{
		float[][] data = ma.getData();
		int m = data.length;
		int n = data[0].length;
		//ArrayList<String> genes = ma.getProbes();
		
		double[] mi = itc.getAllDoubleMIWith(vec, data);
		double preTenthMI = 0;
		int c = 0;
		boolean dominance = true;
		double[] winner = new double[m];
		winner[0] = Double.NaN;
		double bestScore = -1;
		double bestPow = powerStart;

		double a = powerStart;
		
	// Stage 1: identify best power
		while(a >= powerEnd){
			double preScore = 0;
			double[] deltaWindow = {0, 0, 0};
			
			while(c < maxIter){
				preTenthMI = getScore(mi, m, 10);
				double[] w = getWeightVector(mi, m, a);
				float[] metaGene = getMetaGene(data, w, m, n);
				mi = itc.getAllDoubleMIWith(metaGene, data);
				
				double tenthMI = getScore(mi, m, 10);
				double delta = Math.abs(tenthMI - preTenthMI);
				int k = c % 3;
				deltaWindow[k] = delta;
				//System.out.println("Iteration " + (c+1) + "\tDelta = " + delta);
				if(c >= 3){
					boolean scoreConverge = true;
					for(int i = 0; i < 3; i++){
						if(deltaWindow[i] >= 1E-4){
							scoreConverge = false;
							break;
						}
					}
					if(scoreConverge){
						if(dominance){
							w = getWeightVector(mi, m, a);
							double wMax = -1;
							for(int i = 0; i < m; i++){
								if(w[i] > wMax){
									wMax = w[i];
								}
							}
							if(wMax < 0.8) dominance = false;
							else{
								a--;
								System.out.println("Dominance. a = " + a);
								
								c = 0;
								continue;
							}
							
						}
						if(Double.isNaN(winner[0])) {
							System.arraycopy(mi, 0, winner, 0, m);
						}
						// check if the top 10 genes intersect
						//System.out.println("winner[0]: " + winner[0]);
						int[] orderMi = StatOps.order(mi, m);
						int[] orderWinner = StatOps.order(winner, m);
						HashSet<Integer> hs = new HashSet<Integer>();
						boolean intersect = false;
						for(int i = 0; i < 10; i++){
							hs.add(orderMi[m-1-i]);
						}
						for(int i = 0; i < 10; i++){
							if(hs.contains(orderWinner[m-1-i])) {
								intersect = true;
								break;
							}
						}
						if(!intersect) break; //no intersection, discontinuous
						if(tenthMI > bestScore){
							bestScore = tenthMI;
							bestPow = a;
							System.out.println("Best Score: " + bestScore + 
									"\tBest Power: " + bestPow);
							System.arraycopy(mi, 0, winner, 0, m);
						}
						double diffScore = tenthMI - preScore;
						preScore = tenthMI;
						System.out.println("Score: " + tenthMI + "\tDiff Score: " + diffScore);
						if(diffScore >= 0){
							a -= 0.1;
						}else{
							a -= 0.3;
						}
						if(a < 0) break;
						System.out.println("a = " + a);
						c = 0;
					}
					
					
				}
				
				
				c++;
			}
			if(c >= maxIter){
				a -= 0.1;
				c = 0;
				System.out.println("Not converged.");
				System.out.println("a = " + a);
				
			}else break;
		}
		
		if(Double.isNaN(winner[0])) return winner;
		
	// Stage 2: more sophisticated convergence
		System.out.println("Best Score: " + bestScore + "\tBest Pow:" + bestPow);
		double convergeTh = epsilon;
		a = bestPow;
		double[] premi = new double[m];
		System.arraycopy(winner, 0, premi, 0, m);
		double[] w = getWeightVector(premi, m, a);
		float[] metaGene = getMetaGene(data, w, m, n);
		mi = itc.getAllDoubleMIWith(metaGene, data);
		c = 0;
		double delta = 1;
		while(delta >= convergeTh && c <= 200){
			System.arraycopy(mi, 0, premi, 0, m);
			w = getWeightVector(mi, m, a);
			metaGene = getMetaGene(data, w, m, n);
			mi = itc.getAllDoubleMIWith(metaGene, data);
			delta = calcMSE(mi, premi, m);
			//System.out.println(delta);
			c++;
		}
		if(c > 200){
			mi[0] = Double.NaN;
		}
		outPower[0] = bestPow;
		return mi;
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
	public void setPrecision(double precision){
		Converger.precision = precision;
	}
	public void setEpsilon(double epsilon){
		Converger.epsilon = epsilon;
	}
	public void miNormalization(boolean miNorm){
		Converger.miNorm = miNorm;
	}
	public void linkITComputer(ITComputer itc){
		Converger.itc = itc;
	}

	public void setMaxIter(int maxIter){
		this.maxIter = maxIter;
	}
	
	
	
}
