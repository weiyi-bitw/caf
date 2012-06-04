package caf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import obj.DataFile;
import obj.Genome;
import obj.ValIdx;
import worker.Converger;
import worker.ITComputer;

public class FineTuningCNV {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		int[] ws = {51};
		//int wstart = 21;
		//int wend = 201;
		float estart = 1;
		float eend = 6f;
		
		int quantile = 5;
		
		System.out.println("Loading files...");
		
		//String outPath = "/home/weiyi/workspace/javaworks/caf/output/656/";
		
		
		final String[] dataFiles={
				//"/home/weiyi/workspace/data/brca/gse3143/ge.8443x158.jetset.mean.txt",
				//"/home/weiyi/workspace/data/brca/gse3494/ge.12160x251.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/brca/gse32646/ge.19190x115.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/brca/gse36771/ge.19190x107.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/brca/gse31448/ge.19190x353.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/brca/gse2034/ge.12160x286.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/brca/tcga/ge/ge.17475x536.ncbi.txt",
				//"/home/weiyi/workspace/data/coad/gse14333/ge.19190x290.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/coad/tcga/ge/ge.17814x154.knn.txt",
				"/home/weiyi/workspace/data/ov/gse9891/ge.19190x285.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/ov/tcga/ge/ge.11963x582.ncbi.txt",
				"/home/weiyi/workspace/data/ov/gse26193/ge.19190x107.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/prad/gse17951/ge.19189x154.jetset.mean.txt",
				//"/home/weiyi/workspace/data/prad/gse8218/ge.12160x148.jetset.mean.txt"
				//"/home/weiyi/workspace/data/nbl/gse3446/ge.12160x117.jetset.mean.txt",
				//"/home/weiyi/workspace/data/nbl/gse3960/L3.9127x102.mean.txt",
				//"/home/weiyi/workspace/data/sarcoma/gse6481/ge.12160x104.jetset.mean.txt",
				//"/home/weiyi/workspace/data/sarcoma/gse21122/ge.12160x149.jetset.mean.txt",
				//"/home/weiyi/workspace/data/lihc/gse14520/ge.12160x247.jetset.mean.txt",
				//"/home/weiyi/workspace/data/lihc/gse36376/ge.34696x240.mean.txt",
				//"/home/weiyi/workspace/data/gbm/rembrandt/ge.19189x450.jetset.mean.txt",
				//"/home/weiyi/workspace/data/gbm/tcga/ge/ge.12042x545.txt"
		};
		
		final String[] outputDirs={
				//"brca.gse3494.jetset.mean",
				"brca.gse32646.jetset.ncbi",
				"brca.gse36771.jetset.ncbi",
				"brca.gse31448.jetset.ncbi",
				"brca.gse2034.jetset.ncbi",
				"brca.tcga.ncbi",
				//"coad.gse14333.jetset.ncbi",
				//"coad.tcga.ncbi",
				"ov.gse9891.jetset.ncbi",
				"ov.tcga.ncbi",
				"ov.gse26193.jetset.ncbi",
				//"prad.gse17951",
				//"prad.gse8218",
				//"nbl.gse3446",
				//"nbl.gse3960",
				//"sarcoma.gse6481",
				//"sarcoma.gse21122",
				//"lihc.gse14520",
				//"lihc.gse36376",
				//"gbm.rembrandt",
				//"gbm.tcga"
		};
		
		String[] bestSeeds = new String[outputDirs.length];
		
		String outPath = "/home/weiyi/workspace/javaworks/caf/output/window51/erbb2";
		if(!outPath.endsWith("/")){
			outPath = outPath + "/";
		}
		//new File("output").mkdir();
		new File(outPath).mkdir();
		
		final String geneLocFile = "/home/weiyi/workspace/data/annot/ncbi/gene.location.ncbi";
		//final String geneLocFile = "output/window/gene.location3";
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		
		
		for(int qq = 0; qq < dataFiles.length; qq++)
		{
		
		System.out.println("Data: " + dataFiles[qq]);
		new File(outPath + outputDirs[qq]).mkdir();
		
		
		DataFile ma = DataFile.parse(dataFiles[qq]);
		
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		
		//gn.linkToDataFile(ma);
		//String[] testList = gn.getAllGenesInChrArm(targetArm);
		String[] testList = gn.getNeighbors("ERBB2", 101);
				
		long jobID = System.currentTimeMillis();
		Converger cvg = new Converger(0, 1, jobID);
		ITComputer itc = new ITComputer(6, 3, 0, 1, true);
		itc.negateMI(true);
		cvg.linkITComputer(itc);
		HashMap<String, Integer> geneMap = ma.getRows();
		
		
		System.out.println(testList.length + " genes were selected.");
		
		String outSummaryFile = outPath + outputDirs[qq] + "/ScoreSummary.txt";
		PrintWriter pw2 = new PrintWriter(new FileWriter(outSummaryFile));
		pw2.println("Seed\tScore");
		ArrayList<String> geneNames = ma.getProbes();
		String undisputedSeed = "";
		float undisputedScore = -1;
		
		for(String gtest : testList){
			if(!geneNames.contains(gtest)){
				continue;
			}
			float power = estart;
			float bestScore = -1;
			System.out.print("Testing " + gtest + "...");
			ValIdx[] bestVec = new ValIdx[ws[0]];
			int bestWSize = 0;
			float bestPow = 0;
			while(true){
				if(power > eend){
					break;
				}
				//int wsize = wstart;
				for(int wsize:ws){
				/*while(true){
					if(wsize > wend){
						break;
					}*/
					//System.out.println("Window size " + wsize + "\tPower " + power + "...");
								
					String[] neighbors = gn.getNeighbors(gtest, wsize);
					if(neighbors == null){
						//pw.println("No neighbors");
						//pw.close();
						continue;
					}
					DataFile ma2 = ma.getSubProbes(neighbors);
					geneNames = ma2.getProbes();
					m = ma2.getNumRows();
					if(m < quantile){
						continue;
					}
					int idx = ma2.getRows().get(gtest);
					data = ma2.getData();
					float[] vec = data[idx];
					double[] out = cvg.findWeightedAttractorDouble(ma2, vec, power);
					
					
				
					if(out[0] == -1){
						//pw.println("Not converged.");
						//pw.close();
						wsize += 10;
						continue;
					}
								
					ValIdx[] vis = new ValIdx[m];
					for(int i = 0; i < m; i++){
						vis[i] = new ValIdx(i, (float)out[i]);
					}
					Arrays.sort(vis);
							
					float score = vis[quantile-1].val;
					if(score > bestScore){
						bestScore = score;
						bestVec = new ValIdx[m];
						System.arraycopy(vis, 0, bestVec, 0, m);
						bestWSize = wsize;
						bestPow = power;
					}
					//System.out.println("\tScore: " + score);
					//wsize += 50;
				}
				power += 0.5;	
			}
			if(bestScore > undisputedScore){
				undisputedScore = bestScore;
				undisputedSeed = gtest;
			}
			System.out.println("\tScore: " + bestScore + "\tPower: " + bestPow);
			pw2.println(gtest + "\t" + bestScore);
			String outFile = outPath + outputDirs[qq] + "/Attractor_" + gtest + ".txt";
			String[] neighbors = gn.getNeighbors(gtest, bestWSize);
			DataFile ma2 = ma.getSubProbes(neighbors);
			geneNames = ma2.getProbes();
			int m2 = ma2.getNumRows();
			PrintWriter pw = new PrintWriter(new FileWriter(outFile));
			pw.println("Power" + "\t" + bestPow);
			for(int i = 0; i < m2; i++){
				String gg = geneNames.get(bestVec[i].idx);
				pw.println(gg + "\t"  + bestVec[i].val);
			}
			pw.close();
			
			
		}
		pw2.close();
		System.out.println("Final best seed: " + undisputedSeed);
		bestSeeds[qq] = undisputedSeed;
		
		}
		System.out.println("Best seed:");
		for(int qq = 0; qq < dataFiles.length; qq++){
			System.out.println(outputDirs[qq] + ": " + bestSeeds[qq]);
		}
		
	}

}
