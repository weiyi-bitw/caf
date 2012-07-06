package caf;

import java.io.FileInputStream;
import java.util.Properties;

import obj.Annotations;
import obj.DataFileD;
import obj.GeneSet;
import obj.Genome;
import worker.Converger;
import worker.GeneSetMerger;
import worker.ITComputer;
import worker.Scheduler;

public class CorrAttractorFinder {
	private static Properties config;
	private static String configFile;
	private static String command = "CAF";
	
	private static int segment = 0;
	private static int numSegments = 1;
	private static long jobID;
	private static int winsize = 10;
	private static boolean debugging = false;
	private static boolean normMI = true;
	private static boolean negateMI = true;
	private static String breakPoint = "";
	private static int maxIter = 100;
	private static double precision =  1E-4;
	
	// general attractor parameter
	private static double weightExp = 5.0;
	
	// MI parameter
	private static int bins = 6;
	private static int splineOrder = 3;
	
	// CNV parameter
	private static int wstart = 11;
	private static int wend = 51;
	private static int delw = 10;
	private static double pstart = 1;
	private static double pend = 6;
	private static double delp = 0.5f;
	private static int quantile = 5;
	
	private static DataFileD ma;
	private static Annotations annot;
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
				ma = DataFileD.parse(confLine);
				ma.sortProbes();
		    } catch (Exception e) {
		    	throw new RuntimeException("ERROR:Problem parsing data file.\n" + e);
		    }
		//=======Annotations for dataset 1====================
		    annot = null;
		    confLine = config.getProperty("annot");
		    if(confLine != null){
		    	System.out.printf("%-25s%s\n", "GeneExp Annotations:", confLine);
		    	try {
	                annot = Annotations.parseAnnotations(confLine);
	            } catch (Exception e) {
	                throw new RuntimeException("ERROR: Couldn't parse annotations file '" + confLine+ "\n" + e);
	            }
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
	    	
	    	confLine = config.getProperty("min_size");
	    	if (confLine != null && confLine.length() > 0) {
	            try {
	               winsize = Integer.parseInt(confLine);
	            } catch (NumberFormatException nfe) {
	            	System.out.println("WARNING: Couldn't parse minimum attractor size : " + confLine + ", using default = " + winsize);
	            }
	        }
	    	System.out.printf("%-25s%s\n", "Min Size:", winsize);
	    		
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
	
	private static void CAFConfiguration(){
		
		String confLine = config.getProperty("exp");
    	if (confLine != null && confLine.length() > 0) {
            try {
               weightExp = Double.parseDouble(confLine);
            } catch (NumberFormatException nfe) {
            	System.out.println("WARNING: Couldn't parse weight exponent: " + confLine + ", using default " + weightExp);
            }
        }
		System.out.printf("%-25s%s\n", "Exponent:", weightExp);
		
	}
	
	private static void CNVConfiguration(){
		
		String confLine = config.getProperty("win_size_start");
    	if (confLine != null && confLine.length() > 0) {
            try {
               wstart = Integer.parseInt(confLine);
            } catch (NumberFormatException nfe) {
            	System.out.println("WARNING: Couldn't parse starting window size: " + confLine + ", using default. ");
            }
        }
    	if(wstart < 0){
    		System.out.println("WARNING: window size cannot be negative! Using default." );
    		wstart = 11;
    	}
    	
		System.out.printf("%-25s%s\n", "Starting window size:", wstart);
		
		confLine = config.getProperty("win_size_end");
    	if (confLine != null && confLine.length() > 0) {
            try {
               wend = Integer.parseInt(confLine);
            } catch (NumberFormatException nfe) {
            	System.out.println("WARNING: Couldn't parse max window size: " + confLine + ", using default. ");
            }
        }
    	if(wend < 0){
    		System.out.println("WARNING: Max window size cannot be negative or smaller than starting size! Using default." );
    		wend = 51;
    	}
		System.out.printf("%-25s%s\n", "Max window size:", wend);
		
		confLine = config.getProperty("win_size_incre");
    	if (confLine != null && confLine.length() > 0) {
            try {
               delw = Integer.parseInt(confLine);
            } catch (NumberFormatException nfe) {
            	System.out.println("WARNING: Couldn't parse window size increment: " + confLine + ", using default.");
            }
        }
    	if(delw < 0){
    		System.out.println("WARNING: Window size increment cannot be negative! Using default." );
    		delw = 10;
    	}
		System.out.printf("%-25s%s\n", "Window size increment:", delw);
		
		confLine = config.getProperty("pow_start");
    	if (confLine != null && confLine.length() > 0) {
            try {
               pstart = Float.parseFloat(confLine);
            } catch (NumberFormatException nfe) {
            	System.out.println("WARNING: Couldn't parse starting weight power: " + confLine + ", using default.");
            }
        }
    	System.out.printf("%-25s%s\n", "Starting weight power:", pstart);
		
    	confLine = config.getProperty("pow_end");
    	if (confLine != null && confLine.length() > 0) {
            try {
               pend = Float.parseFloat(confLine);
            } catch (NumberFormatException nfe) {
            	System.out.println("WARNING: Couldn't parse max weight power: " + confLine + ", using default.");
            }
        }
    	System.out.printf("%-25s%s\n", "Starting weight power:", pend);
		
    	confLine = config.getProperty("pow_incre");
    	if (confLine != null && confLine.length() > 0) {
            try {
               delp = Float.parseFloat(confLine);
            } catch (NumberFormatException nfe) {
            	System.out.println("WARNING: Couldn't parse weight power increment: " + confLine + ", using default.");
            }
        }
    	System.out.printf("%-25s%s\n", "Weight power increment:", delp);
		
    	confLine = config.getProperty("max_rank");
    	if (confLine != null && confLine.length() > 0) {
            try {
               quantile = Integer.parseInt(confLine);
            } catch (NumberFormatException nfe) {
            	System.out.println("WARNING: Couldn't parse maximized rank: " + confLine + ", using default.");
            }
        }
    	System.out.printf("%-25s%s\n", "Maximized rank:", quantile);
    	
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
		if(args.length < 2) throw new RuntimeException("No command specified!!");
		
		if(args.length >= 2){
			command = args[1];
			if(command.equalsIgnoreCase("DEBUG")){
				System.out.println("**** Debugging mode");
				debugging = true;
				jobID = Long.parseLong(args[2]);
				breakPoint = args[3];
				
				System.out.printf("%-25s%s\n", "ID:", jobID);
				System.out.printf("%-25s%s\n", "BreakPoint:", breakPoint);
			}else if(command.equalsIgnoreCase("CAF")){
				System.out.println("-- Attractor finding");
				
			}else if(command.equalsIgnoreCase("CNV")){
				System.out.println("-- CNV Finding");
				
			}else{
				throw new RuntimeException("ERROR: Cannot understand command: " + command);
			}
		}
			
			configFile = args[0];
			config = new Properties();
			config.load(new FileInputStream(configFile));	
			fileConfiguration();
		System.out.println("\n===================================================================================\n");
		if(command.equalsIgnoreCase("CAF")){
			CAFConfiguration();
		}else if(command.equalsIgnoreCase("CNV")){
			CNVConfiguration();
		}
		System.out.println("\n===================================================================================\n");
		
		Scheduler scdr = new Scheduler(segment, numSegments, jobID);
		Converger cvg = new Converger(segment, numSegments, jobID, maxIter, false);
		ITComputer itc = new ITComputer(bins, splineOrder, segment, numSegments, normMI);
		itc.negateMI(negateMI);
		cvg.linkITComputer(itc);
		cvg.setPrecision(precision);
		int fold = (int) Math.round(Math.sqrt(numSegments));
		if(!debugging)
		{
			// transform the first data matrix into ranks
			if(command.equalsIgnoreCase("CAF")){
				cvg.findWeightedAttractor(ma, weightExp);
			}else if(command.equalsIgnoreCase("CNV")){
				cvg.findWeightedCNVCoef(ma, gn, wstart, wend, delw, pstart, pend, delp, quantile);
				//cvg.findWeightedCNV(ma, gn, pstart, pend, delp, quantile);
			}
			
			if(command.equalsIgnoreCase("CAF")){
				scdr.waitTillFinished(0, fold);
			}
		}
		
		if(!command.equals("CNV")){
		
				ma = null;
			
				if(!debugging || breakPoint.equalsIgnoreCase("merge"))
				{
					if(segment < fold){
						GeneSetMerger mg = new GeneSetMerger(segment, fold, jobID);
						mg.mergeWeightedGeneSets("tmp/" + jobID + "/geneset/", numSegments, precision, false);
					}else{
						System.out.println("Job finished. Exit.");
						System.exit(0);
					}
					
				}
				if(!debugging  || breakPoint.equalsIgnoreCase("output") || breakPoint.equalsIgnoreCase("merge"))
				{
					if(annot != null)GeneSet.setAnnotations(annot);
					if((scdr.allFinished(fold)|| breakPoint.equalsIgnoreCase("output")) && segment==0){
						GeneSetMerger mg = new GeneSetMerger(segment, 1, jobID);
						if(breakPoint.equalsIgnoreCase("output")){
							GeneSetMerger.addMergeCount();
						}
						mg.mergeWeightedGeneSets("tmp/" + jobID + "/merge" + (GeneSetMerger.mergeCount-1), fold, precision, true);
					}
				}
			
		
		}
		
		System.out.println("Done in " + (System.currentTimeMillis() - tOrigin) + " msecs.");
		System.out.println("\n====Thank you!!==================================@ Columbia University 2012=======\n");
	}

}
