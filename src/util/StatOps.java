package util;

import java.util.Arrays;
import java.util.HashMap;

/**
 * @author Wei-Yi Cheng
 * @version 0.2
 * @date 05/05/2011
 */

public class StatOps {
	
	/**
     * Holds for a gene a value and the corresponding phenotype.
     */
    private static class ValdIdx implements Comparable<ValdIdx> {
        double val;
        int idx;
        
        public ValdIdx(int idx, double d){
        	this.val = d;
        	this.idx = idx;
        }
        
        public int compareTo(ValdIdx o) {
            return Double.compare(val, o.val);
        }
    }
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
	public static float mean(float[] x){
		int n = x.length;
		float mean = 0;
		for(int i = 0; i < n; i++){
			if(!Float.isNaN(x[i])){
				mean += x[i];
			}else{
				n--;
			}
		}
		mean /= n;
		return mean;
	}
	
	public static float std(float[] x) {
        int n = x.length;
        float m = mean(x);
        float std = 0;

        for (int i = 0; i < n; i++) {
        	if(!Float.isNaN(x[i])){
        		std += (x[i] - m) * (x[i] - m);
        	}else{
        		n--;
        	}
        }
        int nn = n==0? 1:n-1;
        return (float) Math.sqrt(1 / (float) nn * std);
    }
	public static float rss(float[] x){
		int n = x.length;
		float m = mean(x);
        float ss = 0;

        for (int i = 0; i < n; i++) {
        	if(!Float.isNaN(x[i])){
        		ss += (x[i] - m) * (x[i] - m);
        	}else{
        		n--;
        	}
        }
        
        return ss;
	}
	public static float median(float[] in){
    	int numIn = in.length;
    	Arrays.sort(in);
    	if(numIn%2==0){
			return (in[numIn/2-1] + in[numIn/2])/2;
		}else{
			return in[numIn/2];
		}
    }
	public static float quantile(float[] in, float pct){
		int numIn = in.length;
		Arrays.sort(in);
		return in[ (int) (numIn*pct)];
	}
	public static float mse(float[] x, float[] y){
		int n = x.length;
		int q = 0;
		for(int i = 0; i < n; i++){
			q += (x[i] - y[i]) * (x[i] - y[i]);
		}
		q /= n;
		return (float) Math.sqrt(q);
	}
	
	public static float pearsonCorr(float[] x, float[] y){
		float xMean = 0;
		float yMean = 0;
		double xSd = 0;
		double ySd = 0;
		float rho = 0;
		int N = x.length;
		int n = N;
		for(int i = 0; i < N; i++){
			if(! (Float.isNaN(x[i])||Float.isNaN(y[i]) )){
				xMean += x[i];
				yMean += y[i];
				xSd += x[i] * x[i];
				ySd += y[i] * y[i];
				rho += x[i] * y[i];
			}else{
        		n--;
        	}
		}
		xMean /= n;
		xSd = Math.sqrt((xSd - n*xMean*xMean));
		yMean /= n;
		ySd = Math.sqrt((ySd - n*yMean*yMean));
		
		if(xSd==0 || ySd==0){
			return 0;
		}
		
		rho = (float) ((rho - n * xMean * yMean)/xSd/ySd);
		return rho;
	}
	public static float cov(float[] x, float[] y, int n){
		float xMean = 0, yMean = 0;
		float q = 0;
		for(int i = 0; i < n; i++){
			xMean += x[i];
			yMean += y[i];
			q += x[i] * y[i];
		}
		xMean /= n;
		yMean /= n;
		q = (q - n*xMean*yMean)/(n-1);
		return q;
	}
	public static float[] cov(float[] x, float[][] y, int m, int n){
		float xMean = 0;
		float[] yMean = new float[m];
		float[] rho = new float[m];
		for(int i = 0; i < n; i++){
			xMean += x[i];
			for(int j = 0; j < m; j++){
				yMean[j] += y[j][i];
				rho[j] += x[i] * y[j][i];
			}
		}
		xMean /= n;

		for(int i = 0; i < m; i++){
			yMean[i] /= n;
			rho[i] = (float) ((rho[i] - n * xMean * yMean[i])/(n-1));
		}
		//System.out.println();
		return rho;
	}
	public static float var(float[] x, int n){
		float xMean = 0;
		float v = 0;
		for(int i = 0; i < n; i++){
			xMean += x[i];
			v += x[i] * x[i];
		}
		xMean /= n;
		v = (v - n * xMean * xMean)/(n-1);
		return v;
	}
	
	
	public static float[] pearsonCorr(float[] x, float[][] y, int m, int n){
		float xMean = 0;
		float[] yMean = new float[m];
		double xSd = 0;
		double[] ySd = new double[m];
		float[] rho = new float[m];
		for(int i = 0; i < n; i++){
			xMean += x[i];
			xSd += x[i] * x[i];	
			for(int j = 0; j < m; j++){
				yMean[j] += y[j][i];
				ySd[j] += y[j][i] * y[j][i];
				rho[j] += x[i] * y[j][i];
			}
		}
		xMean /= n;
		xSd = Math.sqrt((xSd - n*xMean*xMean));

		if(xSd == 0){
			return rho;
		}
		
		for(int i = 0; i < m; i++){
			yMean[i] /= n;
			ySd[i] = Math.sqrt((ySd[i] - n*yMean[i]*yMean[i]));
			if(ySd[i]==0){
				rho[i] = 0;
			}else{
				rho[i] = (float) ((rho[i] - n * xMean * yMean[i])/xSd/ySd[i]);
			}
		}
		//System.out.println();
		return rho;
	}
	public static float beta1(float[] y, float[] x){
		float xMean = 0;
		float yMean = 0;
		float xSq = 0;
		float xy = 0;
		int N = x.length;
		int n = N;
		for(int i = 0; i < N; i++){
			if(! (Float.isNaN(x[i])||Float.isNaN(y[i]) )){
				xMean += x[i];
				yMean += y[i];
				xSq += x[i] * x[i];
				xy += x[i] * y[i];
			}else{
        		n--;
        	}
		}
		xMean /= n;
		float denom = xSq - n * xMean * xMean;
		yMean /= n;
		float num = (xy - n * xMean * yMean);
		return (num/denom);
	}
	
	public static float[] mergeSort(float[] x, int[] idx){
		int n = x.length;
		if(n <= 1) return x;
		
		float left[], right[],result[];
		int idxL[], idxR[];
		int m = n/2;
		left = new float[m];
		idxL = new int[m];
		right = new float[n - m];
		idxR = new int[n-m];
		
		for(int i = 0; i < m; i++){
			left[i] = x[i];
			idxL[i] = idx[i];
		}
		for(int i = m; i < n; i++){
			right[i-m] = x[i];
			idxR[i-m] = idx[i];
		}
		left = mergeSort(left, idxL);
		right = mergeSort(right, idxR);
		
		result = merge(left, right, idxL, idxR, idx);
		
		return result;
	}
	
	private static float[] merge(float[] left, float[] right, int[] idxL, int[] idxR, int[] mergedIdx){
		int nL = left.length;
		int nR = right.length;
		float[] result = new float[nL + nR];
		int iL = 0, iR = 0, iOut = 0;;
		while(iL < nL || iR < nR){
			if(iL < nL && iR < nR){
				if(left[iL] <= right[iR]){
					result[iOut] = left[iL];
					mergedIdx[iOut] = idxL[iL];
					iL++; 
				}else{
					result[iOut] = right[iR];
					mergedIdx[iOut] = idxR[iR];
					iR++;
				}
			}else if(iL < nL){
				result[iOut] = left[iL];
				mergedIdx[iOut] = idxL[iL];
				iL++; 
			}else if(iR < nR){
				result[iOut] = right[iR];
				mergedIdx[iOut] = idxR[iR];
				iR++;
			}
			iOut++;
		}
		return result;
	}
	
	public static float[] rank(float[] x){
		int n = x.length;
		float[] xt = new float[n];
		System.arraycopy(x, 0, xt, 0, n);
		Arrays.sort(xt);
		HashMap<Float, Float> map = new HashMap<Float, Float>();
		
		int cumRank = 0;
		int cnt = 0;
		for(int i = 0; i < n; i++){
			if(map.containsKey(xt[i])){
				cumRank += (i);
				cnt++;
			}else{
				if(cumRank != 0){
					map.put(xt[i-1], (float)(cumRank + i)/(cnt+1));
					cumRank = 0;
					cnt = 0;
				}
				map.put(xt[i], (float)i+1);
			}
		}
		
		float[] idx = new float[n];
		for(int i = 0; i < n; i++){
			idx[i] = map.get(x[i]);
		}
		
		return idx;
	}
	public static int[] rankNoTie(float[] x){
		int n = x.length;
		int order[] = new int[n];
		for(int i = 0; i < n; i++){
			order[i] = i;
		}
		mergeSort(x, order);
		int rank[] = new int[n];
		for(int i = 0; i < n; i++){
			rank[order[i]] = i+1;
		}
		return rank;
	}
	public static double[] pAdjust(double[] p, int m){
		if(m == 1) return p;
		ValdIdx[] pVals = new ValdIdx[m];
		double[] padj = new double[m];
		for(int i = 0; i < m; i++){
			pVals[i] = new ValdIdx(i, p[i]);
			padj[i] = 1;
		}
		Arrays.sort(pVals);
		/*for(IdxVal iv : pVals){
			System.out.print(iv.val + "( " + iv.idx + ") \t");
		}System.out.println();*/
		double cummin = 1;
		for(int i = (m-1) ; i >= 0; i--){
			int o = pVals[i].idx;
			padj[o] = Math.min(cummin, pVals[i].val * m / (i+1)) ;
			cummin = padj[o];
		}
		return padj;
	}
	public static double[] pAdjustBonf(double[] p, int m){
		if(m == 1) return p;
		double[] padj = new double[m];
		for(int i = 0 ; i < m; i++){
			padj[i] = m * p[i];
		}
		return padj;
	}
	/*
	 * Transforming the spearman correlation to z-score using Fisher's transformation
	 * 
	 */
	public static float rsToZ(float rs, int n){ 
		return (float)(Math.sqrt((n-3)/1.06)*1/2 * Math.log((1+rs)/(1-rs)));
	}
	
	public static float rpToZ(float rp, int n){ 
		return (float)(Math.sqrt(n-3)*1/2 * Math.log((1+rp)/(1-rp)));
	}
	
	public static float[] xToZ(float[] x, int n){
		float[] z = new float[n];
		int k = n;
		float xMean = 0;
		float xSq = 0;
		for(int i = 0; i < n; i++){
			if(Float.isNaN(x[i])){
				k--;
			}else{
				xMean += x[i];
				xSq += x[i]*x[i];
			}
		}
		xMean /= k;
		float xSd = (float) Math.sqrt((xSq - k*xMean*xMean)/(k-1));
		
		for(int i = 0; i < n; i++){
			z[i] = Float.isNaN(x[i])? Float.NaN : (x[i] - xMean) / xSd;
		}
		
		return z;
	}
}
