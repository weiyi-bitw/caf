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
		int wstart = 11;
		int wend = 51;
		float estart = 1;
		float eend = 5f;
		
		int quantile = 5;
		
		String targetArm = "chr8q24.3";
		
		System.out.println("Loading files...");
		
		//String outPath = "/home/weiyi/workspace/javaworks/caf/output/656/";
		new File("output").mkdir();
		String outPath = "/home/weiyi/workspace/javaworks/caf/output/cnv.finetune/8q24/";
		if(!outPath.endsWith("/")){
			outPath = outPath + "/";
		}
		//final String dataFile = "/home/weiyi/workspace/data/brca/gse2034/ge.12160x286.jetset.txt";
		//final String dataFile = "/home/weiyi/workspace/data/brca/tcga/ge/ge.17814x536.knn.txt";
		//final String dataFile = "/home/weiyi/workspace/data/coad/gse14333/ge.19189x290.jetset.txt";
		//final String dataFile = "/home/weiyi/workspace/data/coad/tcga/ge/ge.17814x154.knn.txt";
		final String dataFile = "/home/weiyi/workspace/data/ov/gse9891/ge.19189x285.jetset.txt";
		//final String dataFile = "/home/weiyi/workspace/data/ov/tcga/ge/ge.17814x584.knn.txt";
		//final String dataFile = "/home/weiyi/workspace/data/ov/tcga/ge/ge.12042x582.txt";
		final String geneLocFile = "/home/weiyi/workspace/data/annot/affy/u133p2/gene.location3";
		
		DataFile ma = DataFile.parse(dataFile);
		
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		gn.linkToDataFile(ma);
		
		ArrayList<String> gs = new ArrayList<String>();
		
		long jobID = System.currentTimeMillis();
		Converger cvg = new Converger(0, 1, jobID);
		ITComputer itc = new ITComputer(6, 3, 0, 1, true);
		//itc.negateMI(true);
		cvg.linkITComputer(itc);
		HashMap<String, Integer> geneMap = ma.getRows();
		
		//gs.add("ERBB2");
		
		String[] testList = gn.getAllGenesInChrArm(targetArm);
		//String[] testList = gn.getNeighbors("ERBB2", 51);
		
		System.out.println(testList.length + " genes were selected.");
		
		String outSummaryFile = outPath + "ScoreSummary.txt";
		PrintWriter pw2 = new PrintWriter(new FileWriter(outSummaryFile));
		pw2.println("Seed\tScore");
		ArrayList<String> geneNames = ma.getProbes();
		String undisputedSeed = "";
		float undisputedScore = -1;
		
		for(String gtest : testList){
			float power = estart;
			float bestScore = -1;
			System.out.println("Testing " + gtest + "...");
			ValIdx[] bestVec = new ValIdx[wstart];
			int bestWSize = 0;
			while(true){
				if(power > eend){
					break;
				}
				int wsize = wstart;
				while(true){
					if(wsize > wend){
						break;
					}
					System.out.println("Window size " + wsize + "\tPower " + power + "...");
								
					String[] neighbors = gn.getNeighbors(gtest, wsize);
					if(neighbors == null){
						//pw.println("No neighbors");
						//pw.close();
						continue;
					}
					DataFile ma2 = ma.getSubProbes(neighbors);
					geneNames = ma2.getProbes();
					m = ma2.getNumRows();
					int idx = ma2.getRows().get(gtest);
					data = ma2.getData();
					float[] vec = data[idx];
					float[] out = cvg.findWeightedCNV(ma2, gtest, gn, vec, wsize, power,  false, false);
								
					if(out[0] == -1){
						//pw.println("Not converged.");
						//pw.close();
						wsize += 10;
						continue;
					}
								
					ValIdx[] vis = new ValIdx[m];
					for(int i = 0; i < m; i++){
						vis[i] = new ValIdx(i, out[i]);
					}
					Arrays.sort(vis);
							
					float score = vis[quantile-1].val;
					if(score > bestScore){
						bestScore = score;
						bestVec = new ValIdx[m];
						System.arraycopy(vis, 0, bestVec, 0, m);
						bestWSize = wsize;
					}
					System.out.println("\tScore: " + score);
					wsize += 10;
				}
				power += 0.5;	
			}
			if(bestScore > undisputedScore){
				undisputedScore = bestScore;
				undisputedSeed = gtest;
			}
			pw2.println(gtest + "\t" + bestScore);
			String outFile = outPath + targetArm + "_Attractor_" + gtest + ".txt";
			String[] neighbors = gn.getNeighbors(gtest, bestWSize);
			DataFile ma2 = ma.getSubProbes(neighbors);
			geneNames = ma2.getProbes();
			PrintWriter pw = new PrintWriter(new FileWriter(outFile));
			for(int i = 0; i < bestWSize; i++){
				String gg = geneNames.get(bestVec[i].idx);
				pw.println(gg + "\t"  + bestVec[i].val);
			}
			pw.close();
			
			
		}
		pw2.close();
		System.out.println("Final best seed: " + undisputedSeed);
	}

}
