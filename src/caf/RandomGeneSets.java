package caf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import caf.GroupWeightedAttractor2.DistPair;
import caf.GroupWeightedAttractor2.WtdAttractor;
import caf.GroupWeightedAttractor2.WtdAttractorSet;

public class RandomGeneSets {
	private static int segment = 0;
	private static int numSegments = 1;
	private static long jobID = 0;
	static int[] attrSizes = {17, 54, 23, 52, 45, 19};
	static String[] dummies = {"0","1","2","3","4","5"};
	static int geneSize = 11395;
	
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
	public static void main(String[] args) throws Exception {
		SystemConfiguration();
		
		int numPerm = 10000000;
		
		int start = numPerm * segment / numSegments;
		int end = numPerm * (segment+1) / numSegments;
		
		System.out.println("Performing permutation " + start + " to " + end );
		
		String outDir = "output/" + jobID;
		new File(outDir).mkdirs();
		PrintWriter pw = new PrintWriter(new FileWriter(outDir + "/" + "numCommonGenes." + String.format("%05d", segment) + ".txt"));
		
	//====== permutation loop ==============================
		for(int p = start; p < end; p++){
	//======================================================
		System.out.print("Permutation " + p + "...");
			
		WtdAttractor.setSeed(p + jobID);
		ArrayList<WtdAttractor> allWtdAttractors = new ArrayList<WtdAttractor>();
		
		//System.out.println("Generate 6 sets of attractors...");
		
		for(int s = 0; s < 6; s++){
			ArrayList<WtdAttractor> waInThisFile = new ArrayList<WtdAttractor>();
			while(waInThisFile.size() < attrSizes[s]){
				WtdAttractor wa = WtdAttractor.generateRandomAttractor(waInThisFile.size(), s, geneSize);
				waInThisFile.add(wa);
			}
			allWtdAttractors.addAll(waInThisFile);
			//System.out.println(" (" + waInThisFile.size() + ") ");
		}
		
		int n = allWtdAttractors.size();
		//System.out.println(n + " attractors were loaded.");
		DistPair.setTotalIdx(n);
		Collections.sort(allWtdAttractors);
		
		ArrayList<DistPair> allDistPairs = new ArrayList<DistPair>();
		
		WtdAttractorSet.setNames(dummies);
		ArrayList<WtdAttractorSet> out = new ArrayList<WtdAttractorSet>();
		//System.out.println("Calculating distance...");
		for(int i = 0; i < n; i++){
			WtdAttractorSet a = new WtdAttractorSet(allWtdAttractors.get(i));
			out.add(a);
		}
		
		for(int i = 0; i < n; i++){
			WtdAttractorSet a = out.get(i);
			
			for(int j = i+1; j < n; j++){
				WtdAttractorSet b = out.get(j);
				double ovlp = a.ovlapWith(b);
				if(ovlp > 0){
					allDistPairs.add(new DistPair(a, b , ovlp));
				}
			}
		}
		Collections.sort(allDistPairs);
		//System.out.println(allDistPairs.size() + " pairs of distances have been added.");
		
		while(allDistPairs.size() > 0){
			DistPair dp = allDistPairs.get(0);
			allDistPairs.remove(dp);
		
			//System.out.println(dp);
			WtdAttractorSet x = dp.x;
			boolean mergeable = x.merge(dp.y);
			if(mergeable){
				out.remove(dp.y);
				for(int i = allDistPairs.size() - 1; i >=0; i--){
					DistPair dpp = allDistPairs.get(i);
					if(dpp.contains(dp.x) || dpp.contains(dp.y)){
						allDistPairs.remove(i);
					}
				}
				int xi = out.indexOf(dp.x);
				for(int i = 0; i < out.size(); i++){
					if(i != xi){
						WtdAttractorSet b = out.get(i);
						double ovlp = x.ovlapWith(b);
						if(ovlp > 0){
							allDistPairs.add(new DistPair(x, b , ovlp));
						}
					}
				}
				Collections.sort(allDistPairs);
			}
		}
		int maxNumCommonGenes = 0;
		for(WtdAttractorSet was: out){
			if(was.numCommonGenes > maxNumCommonGenes){
				maxNumCommonGenes = was.numCommonGenes;
			}
		}
		
		System.out.println(maxNumCommonGenes);
		pw.println(maxNumCommonGenes);
		
		/*
		Collections.sort(out);
		System.out.println("Output to file...");
		PrintWriter pw2 = new PrintWriter(new FileWriter("test.txt"));
		int ii = 1;
		
		for(WtdAttractorSet was : out){
			pw2.print(ii + "\t");
			pw2.println(was);
			ii++;
		}
		pw2.close();
		*/
	// ===== end permutation loop ===========================
		}
	// ======================================================
		
		pw.close();
		
	}
}
