package caf;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Properties;

import obj.Annotations;
import obj.DataFile;
import obj.Genome;

public class Sharpening {
	private static Properties config;
	private static String configFile;
	
	private static int segment = 0;
	private static int numSegments = 1;
	private static long jobID;
	private static boolean normMI = true;
	private static boolean negateMI = true;
	private static String breakPoint = "";
	private static int maxIter = 100;
	private static float precision = (float) 1E-4;
	
	// general attractor parameter
	private static double weightExp = 5.0;
		
	// MI parameter
	private static int bins = 6;
	private static int splineOrder = 3;
	
	private static DataFile ma;
	private static Genome gn;
	
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
			
		//========gene location file========================
		    gn = null;
		    confLine = config.getProperty("genome");
		    if(confLine != null){
		    	System.out.printf("%-25s%s\n", "Genome File:", confLine);
		    	try {
		    		gn = Genome.parseGeneLocation(confLine);
	            } catch (Exception e) {
	                throw new RuntimeException("ERROR: Couldn't parse gene location file '" + confLine+ "\n" + e);
	            }
		    }
		    
	    	confLine = config.getProperty("max_iter");
	    	if (confLine != null && confLine.length() > 0) {
	            try {
	               maxIter = Integer.parseInt(confLine);
	            } catch (NumberFormatException nfe) {
	            	System.out.println("WARNING: Couldn't parse Max Iteration: " + confLine + ", using default " + maxIter);
	            }
	        }
	    	System.out.printf("%-25s%s\n", "Max Iteration:", maxIter);
	    	
	    	confLine = config.getProperty("MI_normalization");
	    	if (confLine != null && confLine.length() > 0) {
	            try {
	               normMI = Boolean.parseBoolean(confLine);
	            } catch (NumberFormatException nfe) {
	            	System.out.println("WARNING: Couldn't parse MI Normalization: " + confLine + ", using default = true");
	            }
	        }
	    	System.out.printf("%-25s%s\n", "MI Normalization:", normMI);
	    	
	    	confLine = config.getProperty("negate_MI");
	    	if (confLine != null && confLine.length() > 0) {
	            try {
	               negateMI = Boolean.parseBoolean(confLine);
	            } catch (NumberFormatException nfe) {
	            	System.out.println("WARNING: Couldn't parse MI negation: " + confLine + ", using default = true");
	            }
	        }
	    	System.out.printf("%-25s%s\n", "Negate MI:", negateMI);
	    	
    		confLine = config.getProperty("precision");
	    	if (confLine != null && confLine.length() > 0) {
	            try {
	               precision = Float.parseFloat(confLine);
	            } catch (NumberFormatException nfe) {
	            	System.out.println("WARNING: Couldn't parse delta: " + confLine + ", using default " + precision);
	            }
	        }
    		System.out.printf("%-25s%s\n", "Precision:", precision);

	    	confLine = config.getProperty("bins");
	    	if (confLine != null && confLine.length() > 0) {
	           try {
	               bins = Integer.parseInt(confLine);
	           } catch (NumberFormatException nfe) {
	               System.out.println("WARNING: Couldn't parse number of bins property: " + confLine + ", use 7");
	           }
	        }
	    	
	    	confLine = config.getProperty("spline_order");
	    	if (confLine != null && confLine.length() > 0) {
	            try {
	               splineOrder = Integer.parseInt(confLine);
	            } catch (NumberFormatException nfe) {
	               System.out.println("WARNING: Couldn't parse spline order property: " + confLine + ", use 3");
	            }
	        }
	    	System.out.printf("%-25s%s\n", "Bins:", bins);
	    	System.out.printf("%-25s%s\n", "Spline Order:", splineOrder);
	}
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		long tOrigin = System.currentTimeMillis();
		
		System.out.println("==Correlation Attractor Sharpener============Wei-Yi Cheng==wc2302@columbia.edu=========\n");
			SystemConfiguration();
		System.out.println("\n===================================================================================\n");
		
		String entryPoint = "/home/weiyi/workspace/javaworks/caf/output/weighted.ncbi/brca.gse2034.jetset.ncbi";
		if(!entryPoint.endsWith("/")) entryPoint += "/";
		FilenameFilter filter = new FilenameFilter() {
			public boolean accept(File dir, String name) {
				return name.startsWith("confCAF");
			}
		    
		};	
		String[] conffile = new File(entryPoint).list(filter);
		
		configFile = entryPoint + conffile[0];
		config = new Properties();
		config.load(new FileInputStream(configFile));	
		
		System.out.println("configFile:" + configFile);
		
		fileConfiguration();
		
		
		
	}

}
