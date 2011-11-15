package caf;

import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Properties;

import org.apache.commons.math.distribution.NormalDistributionImpl;

import obj.DataFile;
import util.StatOps;
import worker.Converger;
import worker.Scheduler;

public class CorrAttractorFinder {
	private static Properties config;
	private static String configFile;
	
	private static int segment = 0;
	private static int numSegments = 1;
	private static long jobID;
	private static int minSize = 0;
	private static int maxSize = Integer.MAX_VALUE;
	private static boolean debugging = false;
	private static boolean rankBased = true;
	private static double fdrThreshold = 0.05;
	private static String breakPoint = null;
	private static int maxIter = 100;
	private static boolean rowNorm = true;
	private static float corrThreshold = 0.7f;
	
	private static DataFile ma;
	
	
	private static void SystemConfiguration() throws Exception{
		String numSegmentsProperty = System.getProperty("SGE_TASK_LAST");
		if (numSegmentsProperty != null && numSegmentsProperty.length() > 0) {
            try {
                numSegments = Integer.parseInt(numSegmentsProperty);
            } catch (NumberFormatException nfe) {
                System.out.println("WARNING: Couldn't parse number of segments property SGE_TASK_LAST=" + numSegmentsProperty + ", use 1");
            }
        }
		String jobIDProperty = System.getProperty("JOB_ID");
        jobID = System.currentTimeMillis();
        if (jobIDProperty != null && jobIDProperty.length() > 0) {
            try {
                jobID = Integer.parseInt(jobIDProperty);
            } catch (NumberFormatException nfe) {
                System.out.println("WARNING: Couldn't parse job id JOB_ID=" + jobIDProperty + ", use current time");
                jobID = System.currentTimeMillis();
            }
        }else if(numSegments != 1){
        	throw new RuntimeException("No job ID is assigned!!");
        }
        System.out.printf("%-25s%s\n", "Job ID: ",jobID);
		String segmentProperty = System.getProperty("SGE_TASK_ID");
        if (segmentProperty != null && segmentProperty.length() > 0) {
            try {
                segment = Integer.parseInt(segmentProperty) - 1;
            } catch (NumberFormatException nfe) {
                System.out.println("WARNING: Couldn't parse segment property SGE_TASK_ID=" + segmentProperty + ", use 0");
            }
        }
        System.out.println("Total Segments: " + numSegments + "\tThis Segment: " + segment);
        
		
	}
	private static void fileConfiguration()throws Exception{
		//========Gene Expression==========
			String confLine = config.getProperty("ge");
			ma = null;
			if(confLine == null){
				throw new RuntimeException("No data file specified!");
			}
			//String[] files = confLine.split(",|;");
			
			System.out.printf("%-25s%s\n", "Expression Data:", confLine);
			try {
				ma = DataFile.parse(confLine);
				ma.sortProbes();
		    } catch (Exception e) {
		    	throw new RuntimeException("ERROR:Problem parsing data file.\n" + e);
		    }
			
			confLine = config.getProperty("rank_based");
	    	if (confLine != null && confLine.length() > 0) {
	            try {
	               rankBased = Boolean.parseBoolean(confLine);
	            } catch (NumberFormatException nfe) {
	            	System.out.println("WARNING: Couldn't parse whether using rank-based MI: " + confLine + ", defaultly using it.");
	            }
	        }
	    	System.out.printf("%-25s%s\n", "Rank-based correlation:", rankBased);
	    	
	    	confLine = config.getProperty("fdr");
	    	if (confLine != null && confLine.length() > 0) {
	            try {
	               fdrThreshold = Double.parseDouble(confLine);
	            } catch (NumberFormatException nfe) {
	            	System.out.println("WARNING: Couldn't parse FDR threshold: " + confLine + ", using default = 0.05.");
	            }
	        }
	    	System.out.printf("%-25s%s\n", "FDR threshold:", fdrThreshold);
	    	
	    	confLine = config.getProperty("max_iter");
	    	if (confLine != null && confLine.length() > 0) {
	            try {
	               maxIter = Integer.parseInt(confLine);
	            } catch (NumberFormatException nfe) {
	            	System.out.println("WARNING: Couldn't parse Max Iteration: " + confLine + ", using default = 100");
	            }
	        }
	    	System.out.printf("%-25s%s\n", "Max Iteration:", maxIter);
	    	
	    	confLine = config.getProperty("row_normalization");
	    	if (confLine != null && confLine.length() > 0) {
	            try {
	               rowNorm = Boolean.parseBoolean(confLine);
	            } catch (NumberFormatException nfe) {
	            	System.out.println("WARNING: Couldn't parse Row Normalization: " + confLine + ", using default = true");
	            }
	        }
	    	System.out.printf("%-25s%s\n", "Row Normalization:", rowNorm);
	    	
	    	confLine = config.getProperty("corr_threshold");
	    	if (confLine != null && confLine.length() > 0) {
	            try {
	               corrThreshold = Float.parseFloat(confLine);
	            } catch (NumberFormatException nfe) {
	            	System.out.println("WARNING: Couldn't parse Correlation Threshold: " + confLine + ", using default = " + corrThreshold);
	            }
	        }
	    	System.out.printf("%-25s%s\n", "Correlation Threshold:", corrThreshold);
	    	
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
	private static float[] getMetaGene(float[][] data, float[] r, HashSet<Integer> idx, int n){
		float w = 0;
		float[] out = new float[n];
		for(Integer i : idx){
			for(int j = 0; j < n; j++){
				out[j] += r[i] * data[i][j];
			}
			w += r[i];
		}
		for(int j = 0; j < n; j++){
			out[j] /= w;
		}
		return out;
	}
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		long tOrigin = System.currentTimeMillis();
		
		System.out.println("==Correlation Attractor Finder============Wei-Yi Cheng==wc2302@columbia.edu===========\n");
			SystemConfiguration();
		System.out.println("\n===================================================================================\n");
			
		if(args.length < 1)	throw new RuntimeException("No configuration file specified!!!");
		
		if(args.length > 1){
			System.out.println("**** Debugging mode");
			debugging = true;
			jobID = Long.parseLong(args[1]);
			breakPoint = args[2];
			
			System.out.printf("%-25s%s\n", "ID:", jobID);
			System.out.printf("%-25s%s\n", "BreakPoint:", breakPoint);
		}
		
			configFile = args[0];
			config = new Properties();
			config.load(new FileInputStream(configFile));	
			fileConfiguration();
		System.out.println("\n===================================================================================\n");
		
		if(rowNorm){
			ma.normalizeRows();
		}
		float[][] data = ma.getData();
		
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		
		// transform the first data matrix into ranks
		float[][] val = new float[m][n];
		if(rankBased){
			for(int i = 0; i < m; i++){
				System.arraycopy(StatOps.rank(data[i]), 0, val[i], 0, n);
			}
		}else{
			val = data;
		}
		Scheduler scdr = new Scheduler(segment, numSegments, jobID);
		Converger cvg = new Converger(segment, numSegments, jobID, fdrThreshold, maxIter, corrThreshold, rankBased);
		cvg.findAttractor(val, data);
		
		
		
	}

}
