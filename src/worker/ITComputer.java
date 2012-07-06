package worker;

/**
 * ITComputer.java
 * @author Wei-Yi Cheng
 * @version 0.23
 * @date 06/22/2011
 */

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import obj.DataFile;

import util.SplineMI;
import util.StatOps;

public class ITComputer extends DistributedWorker{
	private static int bins;
	private static int splineOrder;
	private static int[] knots;
	private static double[] dknots;
	private static boolean output2File = false;
	private static boolean rankBased = false;
	private static boolean negateMI = false;
	private static boolean normalizeMI = true;
	
	public ITComputer(int bins, int splineOrder, int thisSeg, int totalSegs){
		super(thisSeg, totalSegs);
		ITComputer.bins = bins;
		ITComputer.splineOrder = splineOrder;
		int[] knots = new int[bins + splineOrder];
		SplineMI.splineKnots(knots, bins, splineOrder);
		ITComputer.knots = knots;
		//System.out.println("ITComputer " + (id+1) + " of " + totalComputers + " created.");
	}
	public ITComputer(int bins, int splineOrder, int thisSeg, int totalSegs, boolean miNorm){
		super(thisSeg, totalSegs);
		ITComputer.bins = bins;
		ITComputer.splineOrder = splineOrder;
		double[] knots = new double[bins + splineOrder];
		SplineMI.knotVector(knots, bins, splineOrder);
		ITComputer.dknots = knots;
		ITComputer.normalizeMI = miNorm;
		//System.out.println("ITComputer " + (id+1) + " of " + totalComputers + " created.");
	}
	public float[][][] getWeights(DataFile mset){
		int m = mset.getNumRows();
		int n = mset.getNumCols();
		float[][] data = mset.getData();
		float[][][] weights = new float[m][bins][n];
		for (int i = 0; i < m; i++) {
			float[] in = data[i];	
			if(rankBased){
				in = StatOps.rank(in);
			}
            SplineMI.findWeights(in, knots, weights[i], n, splineOrder, bins);
        }
		return weights;
	}
	
	public boolean[][] getValid(DataFile mset){
		int m = mset.getNumRows();
		int n = mset.getNumCols();
		float[][] data = mset.getData();
		boolean[][] valid = new boolean[m][n];
		for (int i = 0; i < m; i++) {
			for(int j = 0; j < n; j++){
				valid[i][j] = Float.isNaN(data[i][j])? false : true;
			}
        }
		return valid;
	}
	public float[] getEntropy1D(float[][][] weights){
		int m = weights.length;
		int n = weights[0][0].length;
		float[] e1 = new float[m];
		for(int i = 0; i < m; i++){
			e1[i] = (float)SplineMI.entropy1f(weights[i], n, bins);
		}
		return e1;
	}
	public float[] getEntropy1D(ArrayList<float[][]> weights){
		int m = weights.size();
		float[] e1 = new float[m];
		
		for(int i = 0; i < m; i++){
			float[][] w = weights.get(i);
			int b = w.length;
			e1[i] = (float)SplineMI.entropy1f(w, w[0].length, b);
		}
		return e1;
	}
	
	public float[] getEntropy1D(DataFile mset){
		return getEntropy1D(getWeights(mset));
	}
	
	public float[][] getEntropy2D(float[][][] weights) throws Exception{
		int m = weights.length;
		int n = weights[0][0].length;
		float[][] e2 = new float[m][m];
		for(int i = 0; i < m; i++){
			for(int j = i; j < m; j++){
				//System.out.println("Processing probe " + (i+1) + "," + (j+1));
				e2[i][j] = (float)SplineMI.entropy2f(weights[i], weights[j], n, bins);
				e2[j][i] = e2[i][j];
			}
		}
		if(output2File){
			PrintWriter pw = new PrintWriter(new FileWriter("output/" + jobID + "/entropy2.txt"));
			for(int i = 0; i < m; i++){
				for(int j = 0; j < m; j++){
					if(j==0) pw.print(e2[i][j]);
					else pw.print("\t" + e2[i][j]);
				}pw.println();
			}
			pw.close();
		}
		return e2;
	}
	
	public float[][] getEntropy2D(DataFile mset) throws Exception{
		return getEntropy2D(getWeights(mset));
	}
	public void getDistEntropy2D(float[][][] weights, boolean[][] valid) throws Exception{
		int m = weights.length;
		int n = weights[0][0].length;
		int start = m * id / totalComputers;
		int end = m * (id + 1)/totalComputers;
		System.out.println("Processing 2D entropies probe " + (start+1) + " to probe " + end);
		prepare("ent");
		for(int i = start; i < end; i++){
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("tmp/" + jobID + "/ent/entropy2-" + String.format("%05d", i) + ".bin")));
			for(int j = 0; j < m; j++){
				float e2 = (float)SplineMI.entropy2f(weights[i], weights[j], n, bins);
				out.writeFloat(e2);
			}
			out.close();
		}
		System.out.println("Done.");
	}
	public void getDistEntropy2D(DataFile mset) throws Exception{
		getDistEntropy2D(getWeights(mset), getValid(mset));
	}
	public void outputDistEntropy2D(int m) throws Exception{
		File dir = new File("tmp/" + jobID + "/ent");
        String[] allFiles = dir.list();
        int n = allFiles.length;
        if(m != n){
        	throw new RuntimeException("Number of entropy files is different from number of probes.");
        }
        PrintWriter pw = new PrintWriter(new FileWriter("output/" + jobID + "/entropy2.txt"));
        for(int i = 0; i < n; i++){
        	DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream("tmp/" + jobID + "/ent/entropy2-" + String.format("%05d", i) + ".bin")));
        	int l = dis.available()*Byte.SIZE/Float.SIZE; // float: 4 bytes
        	for(int j = 0; j < l; j++){
        		if(j==0) pw.print(dis.readFloat());
				else pw.print("\t" + dis.readFloat());
        	}pw.println();
        }
		pw.close();
        
	}
	public float[] getAllMIWith(int idx, float[][][] weights) throws Exception{
		int n = weights[0][0].length;
		int m = weights.length;
		
		float[] mi = new float[m];
		float miMax = 0;
		for(int i = 0; i < m; i++){
			float e1tf = 0, e1tg = 0; 
			float[] histValtf = new float[bins];
			float[] histValtg = new float[bins];
			int numSamples = n;
			
			for(int curSample = 0; curSample < n; curSample++){
				
				if(!Float.isNaN(weights[idx][0][curSample]) && !Float.isNaN(weights[i][0][curSample])){
					for (int curBin = 0; curBin < bins; curBin++) {
						histValtf[curBin] += weights[idx][curBin][curSample];
						histValtg[curBin] += weights[i][curBin][curSample];
		            }
	        	}else{
	        		numSamples--;
	        	}
	        }
			
			for (int curBin = 0; curBin < bins; curBin++){
				histValtg[curBin] /= numSamples;
				if (histValtg[curBin] > 0) {
	        		e1tg -= histValtg[curBin] * SplineMI.log2d(histValtg[curBin]);
	        	}
				histValtf[curBin] /= numSamples;
				if (histValtf[curBin] > 0) {
	        		e1tf -= histValtf[curBin] * SplineMI.log2d(histValtf[curBin]);
	        	}
			}
			
			float e2 = (float)SplineMI.entropy2f(weights[idx], weights[i], n, bins);
			if(!normalizeMI){
				mi[i] = (e1tf + e1tg - e2);
			}else{
				float e2tf = (float)SplineMI.entropy2f(weights[idx], weights[idx], n, bins);
				float e2tg = (float)SplineMI.entropy2f(weights[i], weights[i], n, bins);
				float mitf = 2*e1tf - e2tf;
				float mitg = 2*e1tg - e2tg;
				
				mi[i] = (e1tf + e1tg - e2) / Math.max(mitf, mitg);
				if(miMax < mi[i]) miMax = mi[i];
			}
		}
		if(normalizeMI){
			for(int i = 0; i < m; i++){
				mi[i] = mi[i]/miMax;
			}
		}
		return mi;
	}
	private static float getMomentSign(final float[] x, final float[] y, final int n){
		float productMoment = 0;
		float xmean = 0, ymean = 0;
		int numSamples = n;
		for(int i = 0; i < n; i++){
			if(!Float.isNaN(x[i]) && !Float.isNaN(y[i])){
				xmean += x[i];
				ymean += y[i];
				productMoment += x[i] * y[i];
			}else{
				numSamples--;
			}
		}
		
		return Math.signum(productMoment - xmean*ymean/numSamples);
		
	}
	private static double getMomentSign(final double[] x, final double[] y, final int n){
		double productMoment = 0;
		double xmean = 0, ymean = 0;
		for(int i = 0; i < n; i++){
			xmean += x[i];
			ymean += y[i];
			productMoment += x[i] * y[i];
		}
		
		return Math.signum(productMoment - xmean*ymean/n);
		
	}
	public double[] getAllDoubleMIWith(double[] fixVec, double[][] data) throws Exception{
		int n = data[0].length;
		int m = data.length;
		
	// calculate the MI of fixVec to itself (for normalization)
		
		double[][] weightFix = new double[bins][n];
		SplineMI.findWeights(fixVec, dknots, weightFix, n, splineOrder, bins);
		double e1tf = 0;
		for (int curBin = 0; curBin < bins; curBin++){
			double h = 0;
			for(int curSample = 0; curSample < n; curSample++){
				h += weightFix[curBin][curSample];
	        }
			h /= n;
			if (h > 0) {
        		e1tf -= h * SplineMI.log2d(h);
        	}
		}
		
		//float miMax = 2 * e1fix - e2fix;
		
		double[] mi = new double[m];
		for(int i = 0; i < m; i++){
			
			double[][] weightTg = new double[bins][n];
			SplineMI.findWeights(data[i], dknots, weightTg, n, splineOrder, bins);
			
			double e1tg = 0;
			
			for (int curBin = 0; curBin < bins; curBin++){
				double h = 0;
				for(int curSample = 0; curSample < n; curSample++){
					h += weightTg[curBin][curSample];
				}
				h /= n;
				if (h > 0) {
	        		e1tg -= h * SplineMI.log2d(h);
	        	}
			}
			
			double e2 = SplineMI.entropy2d(weightFix, weightTg, n, bins);
			if(!normalizeMI){
				mi[i] = (e1tf + e1tg - e2);
				if(negateMI) mi[i] = mi[i] * getMomentSign(fixVec, data[i], n);
			}else{
				double e2tf = SplineMI.entropy2d(weightFix, weightFix, n, bins);
				double e2tg = SplineMI.entropy2d(weightTg, weightTg, n, bins);
				double mitf = 2*e1tf - e2tf;
				double mitg = 2*e1tg - e2tg;
				double largerMI = mitf > mitg? mitf:mitg;
				if(largerMI == 0) largerMI = 1;
				mi[i] = (e1tf + e1tg - e2) / largerMI;
				if(negateMI) mi[i] = mi[i] * getMomentSign(fixVec, data[i], n);
			}
			
		}
		return mi;
	}
	
	public double[] getAllDoubleMIWith(float[] fixVec, float[][] data) throws Exception{
		int n = data[0].length;
		int m = data.length;
		
	// calculate the MI of fixVec to itself (for normalization)
		
		double[][] weightFix = new double[bins][n];
		SplineMI.findWeights(fixVec, dknots, weightFix, n, splineOrder, bins);
		double[] histValtf = new double[bins];
		int numSamples = n;
		
		//float miMax = 2 * e1fix - e2fix;
		
		double[] mi = new double[m];
		for(int i = 0; i < m; i++){
			
			double[][] weightTg = new double[bins][n];
			SplineMI.findWeights(data[i], dknots, weightTg, n, splineOrder, bins);
			
			double e1tf = 0, e1tg = 0;
			histValtf = new double[bins];
			double[] histValtg = new double[bins];
			numSamples = n;
			
			for(int curSample = 0; curSample < n; curSample++){
				
				if(!Double.isNaN(weightFix[0][curSample]) && !Double.isNaN(weightTg[0][curSample])){
					for (int curBin = 0; curBin < bins; curBin++) {
						histValtf[curBin] += weightFix[curBin][curSample];
						histValtg[curBin] += weightTg[curBin][curSample];
		            }
	        	}else{
	        		numSamples--;
	        	}
	        }
			
			for (int curBin = 0; curBin < bins; curBin++){
				histValtg[curBin] /= numSamples;
				if (histValtg[curBin] > 0) {
	        		e1tg -= histValtg[curBin] * SplineMI.log2d(histValtg[curBin]);
	        	}
				histValtf[curBin] /= numSamples;
				if (histValtf[curBin] > 0) {
	        		e1tf -= histValtf[curBin] * SplineMI.log2d(histValtf[curBin]);
	        	}
			}
			
			double e2 = SplineMI.entropy2d(weightFix, weightTg, n, bins);
			if(!normalizeMI){
				mi[i] = (e1tf + e1tg - e2);
				if(negateMI) mi[i] = mi[i] * getMomentSign(fixVec, data[i], n);
			}else{
				double e2tf = SplineMI.entropy2d(weightFix, weightFix, n, bins);
				double e2tg = SplineMI.entropy2d(weightTg, weightTg, n, bins);
				double mitf = 2*e1tf - e2tf;
				double mitg = 2*e1tg - e2tg;
				double largerMI = mitf > mitg? mitf:mitg;
				if(largerMI == 0) largerMI = 1;
				mi[i] = (e1tf + e1tg - e2) / largerMI;
				if(negateMI) mi[i] = mi[i] * getMomentSign(fixVec, data[i], n);
			}
			
		}
		return mi;
	}
	
	
	public float[][] getMI2D(float[][][] weights, int[] tfs, int[] targets) throws Exception{
		int n = weights[0][0].length;
		int numTfs = tfs.length;
		int numTargets = targets.length;
		float[] e1 = getEntropy1D(weights);
		float[][] mi = new float[numTfs][numTargets];
		for(int i = 0; i < numTfs; i++){
			int tf = tfs[i];
			for(int j = 0; j < numTargets; j++){
				int target = targets[j];
				float e2 = (float)SplineMI.entropy2f(weights[tf], weights[target], n, bins);
				mi[i][j] = e1[tf] + e1[target] - e2;
				
			}
		}
		if(output2File){
			PrintWriter pw = new PrintWriter(new FileWriter("output/" + jobID + "/mi2.txt"));
			for(int i = 0; i < numTfs; i++){
				for(int j = 0; j < numTargets; j++){
					if(j==0) pw.print(mi[i][j]);
					else pw.print("\t" + mi[i][j]);
				}pw.println();
			}
			pw.close();
		}
		return mi;
	}
	public float[][] getMI2D(DataFile mset, int[] tfs, int[] targets) throws Exception{
		return getMI2D(getWeights(mset), tfs, targets);
	}
	public void getDistMI2D(float[][][] weights) throws Exception{
		int n = weights[0][0].length;
		int m = weights.length;
		float[] e1 = getEntropy1D(weights);
		
		
		int start = m * id / totalComputers;
		int end = m * (id + 1)/totalComputers;
		System.out.println("Processing MI of probe " + (start+1) + " to probe " + end);
		prepare("mi");
		for(int i = start; i < end; i++){
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("tmp/" + jobID + "/mi/probe" + String.format("%05d", i) + ".bin")));
			for(int j = 0; j < m; j++){
				float e2 = (float)SplineMI.entropy2f(weights[i], weights[j], n, bins);
				out.writeFloat(e1[i] + e1[j] - e2);
			}
			out.close();
		}
	}
	public void getDistMI2D(DataFile mset) throws Exception{
		getDistMI2D(getWeights(mset));
	}
	public void getDistMI2D(DataFile mset, DataFile se, int[] tfs, HashMap<Integer, Integer> tfMap)throws Exception{
		getDistMI2D(getWeights(mset), getWeights(se), tfs, tfMap);
	}
	
	public void getDistMI2D(float[][][] weightsTg, float[][][] weightsTf, int[] tfs, HashMap<Integer, Integer> tfMap) throws Exception{
		int k = weightsTg[0][0].length; // number of samples
		int n = weightsTg.length;	// number of genes
		int m = weightsTf.length;	// number of tfs
		
		float[] e1Tf = getEntropy1D(weightsTf);
		float[] e1Tg = getEntropy1D(weightsTg);
		
		int start = n * id / totalComputers;
		int end = n * (id + 1)/totalComputers;
		System.out.println("Processing MI of probe " + (start+1) + " to probe " + end);
		prepare("mi");
		for(int i = start; i < end; i++){
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("tmp/" + jobID + "/mi/mi2-probe" + String.format("%05d", i) + ".bin")));
			for(int j = 0; j < m; j++){
				float e2 = 0;
				int tf = tfs[j];
				//if(tfMap.containsKey(j)){
					//int tfIdx = tfMap.get(j);
					int tfIdx = tfMap.get(tf); 
					e2 = (float)SplineMI.entropy2f(weightsTg[i], weightsTf[tfIdx], k, bins);
					out.writeFloat(e1Tg[i] + e1Tf[tfIdx] - e2);
					tfIdx++;
				//}
				/*
				else{
					e2 = (float)SplineMI.entropy2d(weightsTg[i], weightsTg[j], k, bins);
					out.writeFloat(e1Tg[i] + e1Tg[j] - e2);
				}
				*/
			}
			out.close();
		}
		
	}
	
	public void getDistMI2D(float[][][] weightsGe, ArrayList<float[][]> weightsClnc) throws Exception{
		int n = weightsGe[0][0].length; // number of samples
		int mg = weightsGe.length;	// number of genes
		int mf = weightsClnc.size();	// number of features
		
		int start = mg * id / totalComputers;
		int end = mg * (id + 1)/totalComputers;
		System.out.println("Processing MI of probe " + (start+1) + " to probe " + end);
		prepare("mi");
		for(int i = start; i < end; i++){
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("tmp/" + jobID + "/mi/mi2-probe" + String.format("%05d", i) + ".bin")));
			for(int j = 0; j < mf; j++){
				float[][] w = weightsClnc.get(j);
				int b = w.length;
				float e1g = 0, e1f = 0; 
				float wf = 0;
				float[] histValF = new float[b];
				float[] histValG = new float[bins];
				int numSamples = n;
				for(int curSample = 0; curSample < n; curSample++){
					
					if(!Float.isNaN(w[0][curSample])){
						for (int curBin = 0; curBin < b; curBin++) {
							wf = w[curBin][curSample];
							histValF[curBin] += wf;
			            }
						for (int curBin = 0; curBin < bins; curBin++){
							histValG[curBin] += weightsGe[i][curBin][curSample];
						}
	            	}else{
	            		numSamples--;
	            	}
					
		            
		        }
				for (int curBin = 0; curBin < b; curBin++) {
					histValF[curBin] /= numSamples;
					if (histValF[curBin] > 0) {
		        		e1f -= histValF[curBin] * SplineMI.log2d(histValF[curBin]);
		        	}
				}
				for (int curBin = 0; curBin < bins; curBin++){
					histValG[curBin] /= numSamples;
					if (histValG[curBin] > 0) {
		        		e1g -= histValG[curBin] * SplineMI.log2d(histValG[curBin]);
		        	}
				}
	            
	            
	            
				float e2 = (float)SplineMI.entropy2f(weightsGe[i], w, n, bins, b);
				out.writeFloat(e1g + e1f - e2);
			}
			out.close();
		}
		
	}
	public void getDistMI2D(DataFile ma1, DataFile ma2)throws Exception{
		getDistMI2D(getWeights(ma1), getWeights(ma2));
	}
	public void getDistMI2D(float[][][] wRegulator, float[][][] wRegulon) throws Exception{
		int n = wRegulator[0][0].length;
		if(n != wRegulon[0][0].length){
			throw new RuntimeException("Two files do not have same number of samples!");
		}
		int ntf = wRegulator.length;
		int ntg = wRegulon.length;
		
		int start =  ntf * id / totalComputers;
		int end = ntf * (id + 1)/totalComputers;
		System.out.println("Processing regulator " + (start+1) + " to " + end);
		prepare("mi");
		
		for(int i = start; i < end; i++){
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("tmp/" + jobID + "/mi/probe" + String.format("%05d", i) + ".bin")));
			for(int j = 0; j < ntg; j++){
				float e1tf = 0, e1tg = 0; 
				float[] histValtf = new float[bins];
				float[] histValtg = new float[bins];
				int numSamples = n;
				for(int curSample = 0; curSample < n; curSample++){
					
					if(!Float.isNaN(wRegulon[j][0][curSample]) && !Float.isNaN(wRegulator[i][0][curSample])){
						for (int curBin = 0; curBin < bins; curBin++) {
							histValtf[curBin] += wRegulator[i][curBin][curSample];
							histValtg[curBin] += wRegulon[j][curBin][curSample];
			            }
	            	}else{
	            		numSamples--;
	            	}
		        }
				for (int curBin = 0; curBin < bins; curBin++){
					histValtg[curBin] /= numSamples;
					if (histValtg[curBin] > 0) {
		        		e1tg -= histValtg[curBin] * SplineMI.log2d(histValtg[curBin]);
		        	}
					histValtf[curBin] /= numSamples;
					if (histValtf[curBin] > 0) {
		        		e1tf -= histValtf[curBin] * SplineMI.log2d(histValtf[curBin]);
		        	}
				}
	            
	            
	            
				float e2 = (float)SplineMI.entropy2f(wRegulon[j], wRegulator[i], n, bins);
				out.writeFloat(e1tg + e1tf - e2);
			}
			out.close();
		}
		
	}
	
	
	public void getDistMI2D(float[][][] weights, int[] tfs, int[] targets) throws Exception{
		int numTfs = tfs.length;
		int numTargets = targets.length;
		int n = weights[0][0].length;
		float[] e1 = getEntropy1D(weights);
		
		
		int start = numTfs * id / totalComputers;
		int end = numTfs * (id + 1)/totalComputers;
		System.out.println("Processing MI of tf " + (start+1) + " to tf " + end);
		prepare("mi");
		for(int i = start; i < end; i++){
			int tf = tfs[i];
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("tmp/" + jobID + "/mi/probe" + String.format("%05d", i) + ".bin")));
			for(int j = 0; j < numTargets; j++){
				int target = targets[j];
				float e2 = (float)SplineMI.entropy2f(weights[tf], weights[target], n, bins);
				out.writeFloat(e1[tf] + e1[target] - e2);
			}
			out.close();
		}
	}
	public void getDistMI2D(DataFile mset, int[] tfs, int[] targets) throws Exception{
		getDistMI2D(getWeights(mset), tfs, targets);
	}
	public void getMplusS(DataFile mset, int[] tfs, int[] tgs, int[] ptns) throws Exception{
		getMplusS(getWeights(mset), tfs, tgs, ptns);
	}
	public void getMplusS(float[][][] weights, int[] tfs, int[] tgs, int[] ptns) throws Exception{
		int numProbes = weights.length;
		int numTfs = tfs.length;
		int numTargets = tgs.length;
		int numPartners = ptns.length;

		int n = weights[0][0].length;
		float[] e1 = getEntropy1D(weights);
		int start = numTargets * id / totalComputers;
		int end = numTargets * (id + 1)/totalComputers;
		System.out.println("Processing synergies of target " + (start+1) + " to target " + end);
		prepare("mPs");
		prepare("partner");
		float[] e2Tg = new float[numProbes];
		float[] e2Tf = new float[numProbes];
		DataInputStream dis;
		DataOutputStream out, outP;
		for(int i = start; i < end; i++){
			System.out.print("\rProcessing target " + (i+1) + "...");
			int tg = tgs[i];
			out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("tmp/" + jobID + "/mPs/mPlusS-" + String.format("%05d", i) + ".bin")));
			outP = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("tmp/" + jobID + "/partner/partner-" + String.format("%05d", i) + ".bin")));
			dis = new DataInputStream(new BufferedInputStream(new FileInputStream("tmp/" + jobID + "/ent/entropy2-" + String.format("%05d", tg) + ".bin")));
			for(int k = 0; k < numProbes; k++){
				e2Tg[k] = dis.readFloat();
			}
			dis.close();
			for(int j = 0; j < numTfs; j++){
				float syn = 0;
				int p = -1;
				int tf = tfs[j];
				float mi2 = e1[tf] + e1[tg] - e2Tg[tf];
					
				if(tf != tg){
					dis = new DataInputStream(new BufferedInputStream(new FileInputStream("tmp/" + jobID + "/ent/entropy2-" + String.format("%05d", tf) + ".bin")));
					for(int k = 0; k < numProbes; k++){
						e2Tf[k] = dis.readFloat();
					}
					dis.close();
					for(int k = 0; k < numPartners; k++){
						int ptn = ptns[k];
						if(ptn != tf && ptn != tg){
							float miTfP = e1[tf] + e1[ptn] - e2Tf[ptn];
							float miTgP = e1[tg] + e1[ptn] - e2Tg[ptn];
							
							if((mi2 >= miTfP) && (miTgP >= miTfP)){
							
								float e3 = (float) SplineMI.entropy3f(weights[tf], weights[tg], weights[ptn], n, bins);
								float mi3 = e1[tf] + e1[tg] + e1[ptn] - e2Tf[ptn] - e2Tg[ptn] - e2Tg[tf] + e3;
								if(-mi3 > syn){
									syn = -mi3;
									p = ptn;
								}
							}
						}
					}
				}
				out.writeFloat(mi2 + syn);
				outP.writeInt(p);
			}
			out.close();
			outP.close();
		}
		System.out.println("Done.");
	}
	public void getMplusS(DataFile mset, int[] ptns) throws Exception{
		getMplusS(getWeights(mset), ptns);
	}
	public void getMplusS(float[][][] weights, int[] ptns) throws Exception{
		int numProbes = weights.length;
		int numPartners = ptns.length;

		int n = weights[0][0].length;
		float[] e1 = getEntropy1D(weights);
		int start = numProbes * id / totalComputers;
		int end = numProbes * (id + 1)/totalComputers;
		System.out.println("Processing synergies of probe " + (start+1) + " to probe " + end);
		prepare("mPs");
		prepare("partner");
		float[] e2Tg = new float[numProbes];
		float[] e2Tf = new float[numProbes];
		DataInputStream dis;
		DataOutputStream out, outP;
		for(int i = start; i < end; i++){
			System.out.print("\rProcessing probe " + (i+1) + "...");
			out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("tmp/" + jobID + "/mPs/mPlusS-" + String.format("%05d", i) + ".bin")));
			outP = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("tmp/" + jobID + "/partner/partner-" + String.format("%05d", i) + ".bin")));
			dis = new DataInputStream(new BufferedInputStream(new FileInputStream("tmp/" + jobID + "/ent/entropy2-" + String.format("%05d", i) + ".bin")));
			for(int k = 0; k < numProbes; k++){
				e2Tg[k] = dis.readFloat();
			}
			dis.close();
			for(int j = 0; j < numProbes; j++){
				float syn = 0;
				int p = -1;
				float mi2 = e1[i] + e1[j] - e2Tg[j];
					
				if(i != j){
					dis = new DataInputStream(new BufferedInputStream(new FileInputStream("tmp/" + jobID + "/ent/entropy2-" + String.format("%05d", j) + ".bin")));
					for(int k = 0; k < numProbes; k++){
						e2Tf[k] = dis.readFloat();
					}
					dis.close();
					for(int k = 0; k < numPartners; k++){
						int ptn = ptns[k];
						if(i != ptn && j != ptn){
							float miTfP = e1[j] + e1[ptn] - e2Tf[ptn];
							float miTgP = e1[i] + e1[ptn] - e2Tg[ptn];
							
							if((mi2 >= miTfP) && (miTgP >= miTfP)){
							
								float e3 = (float) SplineMI.entropy3f(weights[j], weights[i], weights[ptn], n, bins);
								float mi3 = e1[j] + e1[i] + e1[ptn] - e2Tf[ptn] - e2Tg[ptn] - e2Tg[j] + e3;
								
								if(-mi3 > syn){
									syn = -mi3;
									p = ptn;
								}
							}
						}
					}
				}
				
				out.writeFloat(mi2 + syn);
				outP.writeInt(p);
			}
			out.close();
			outP.close();
		}
		System.out.println("Done.");
	}
	public void outputDistMI2D(int m) throws Exception{
		File dir = new File("tmp/" + jobID + "/mi");
        String[] allFiles = dir.list();
        int n = allFiles.length;
        if(m != n){
        	throw new RuntimeException("Number of mi files is different from number of targets.");
        }
        PrintWriter pw = new PrintWriter(new FileWriter("output/" + jobID + "/mi2.txt"));
        for(int i = 0; i < n; i++){
        	DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream("tmp/" + jobID + "/mi/mi2-target" + String.format("%05d", i) + ".bin")));
        	int l = dis.available()*Byte.SIZE/Float.SIZE;; // float: 4 bytes
        	for(int j = 0; j < l; j++){
        		if(j==0) pw.print(dis.readFloat());
				else pw.print("\t" + dis.readFloat());
        	}pw.println();
        	if(clean){
            		new File("tmp/" + jobID + "/mi/mi2-target" + String.format("%05d", i) + ".bin").delete();
            	}
		dis.close();
        }
		pw.close();
        
	}
	public void generateBGMI(DataFile mset , int[] tfs, int[] targets, int numPerm) throws Exception{
		int numTfs = tfs.length;
		int numTargets = targets.length;
		prepare("permMi");
		DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("tmp/" + jobID + "/permMi/mi2-permWorker" + String.format("%05d", id) + ".bin")));
		int start = numTargets * id / totalComputers;
		int end = numTargets * (id + 1)/totalComputers;
		System.out.println("Processing target " + (start+1) + " to target " + end);
		
		for(int p = 0; p < numPerm; p++){
			System.out.print("Permutation " + p + "...");
			Random randGen = new Random(numPerm*id + p);
			System.out.print("Seed: " + (numPerm*id + p) + "...");
	        for (int i = 0; i < numTargets; i++) {
	            mset.permute(randGen.nextLong(), targets[i]);
	        }
			float[][][] weights = getWeights(mset);
			int n = weights[0][0].length;
			float[] e1 = getEntropy1D(weights);
			if(totalComputers > numTargets){
				throw new RuntimeException("Number of segments is greater than the number of targets.");
			}
			for(int i = start; i < end; i++){
				int target = targets[i];
				for(int j = 0; j < numTfs; j++){
					int tf = tfs[j];
					float e2 = (float)SplineMI.entropy2f(weights[tf], weights[target], n, bins);
					out.writeFloat(e1[tf] + e1[target] - e2);
				}
			}
			System.out.println("Done.");
		}
		out.close();
	}
	
	public void setBinsSplineOrders(int bins, int splineOrder) {
		ITComputer.bins = bins; 
		ITComputer.splineOrder = splineOrder;
		int[] knots = new int[bins + splineOrder];
		SplineMI.splineKnots(knots, bins, splineOrder);
		ITComputer.knots = knots;
	}
	public void rankBasedMI(boolean b){
		ITComputer.rankBased = b;
	}
	public void negateMI(boolean b){
		ITComputer.negateMI = b;
	}
	public void normalizeMI(boolean b){
		ITComputer.normalizeMI = b;
	}
}
