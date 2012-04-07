package caf;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import obj.Annotations;
import obj.DataFile;
import obj.Genome;
import obj.InverseAnnotations;
import obj.ValIdx;
import worker.Converger;
import worker.ITComputer;

public class FineTuningCNVProbe {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		int wstart = 11;
		int wend = 51;
		float estart = 1;
		float eend = 6f;
		
		int quantile = 5;
		
		String targetArm = "chr17q12";
		
		System.out.println("Loading files...");
		
		//String outPath = "/home/weiyi/workspace/javaworks/caf/output/656/";
		
		//final String dataFile = "/home/weiyi/workspace/data/brca/gse2034/ge.22283x286.txt";
		final String dataFile = "/home/weiyi/workspace/data/brca/tcga/ge/ge.64847x536.knn.txt";
		//final String dataFile = "/home/weiyi/workspace/data/coad/gse14333/ge.19189x290.jetset.txt";
		//final String dataFile = "/home/weiyi/workspace/data/coad/tcga/ge/ge.17814x154.knn.txt";
		//final String dataFile = "/home/weiyi/workspace/data/ov/gse9891/ge.19189x285.jetset.txt";
		//final String dataFile = "/home/weiyi/workspace/data/ov/tcga/ge/ge.17814x584.knn.txt";
		//final String dataFile = "/home/weiyi/workspace/data/ov/tcga/ge/ge.12042x582.txt";
		final String geneLocFile = "/home/weiyi/workspace/data/annot/affy/u133p2/gene.location3";
		
		String[] dataFiles = {
				"/home/weiyi/workspace/data/brca/tcga/ge/ge.64847x536.knn.txt",
				"/home/weiyi/workspace/data/coad/gse14333/ge.54675x290.txt"
				//"/home/weiyi/workspace/data/coad/tcga/ge/ge.64847x147.knn.txt",
				//"/home/weiyi/workspace/data/ov/gse9891/ge.54675x285.txt",
				//"/home/weiyi/workspace/data/ov/tcga/ge/ge.22277x586.txt"
		};
		
		String[] outdirs = {
				"brca.tcga",
				"coad.gse14333"
				//"coad.tcga",
				//"ov.gse9891",
				//"ov.tcga"
		};
		String[] annots = {
				"/home/weiyi/workspace/data/annot/tcga/4502a073/annot.csv",
				"/home/weiyi/workspace/data/annot/affy/u133p2/annot.csv"
				//"/home/weiyi/workspace/data/annot/tcga/4502a073/annot.csv",
				//"/home/weiyi/workspace/data/annot/affy/u133p2/annot.csv",
				//"/home/weiyi/workspace/data/annot/affy/u133a/annot.csv"
		};
		
		for(int qq = 0; qq < 4; qq++){
		
		DataFile ma = DataFile.parse(dataFiles[qq]);
		
		//String annotPath = "/home/weiyi/workspace/data/annot/affy/u133a/annot.csv";
		//String annotPath = "/home/weiyi/workspace/data/annot/tcga/4502a073/annot.csv";
		String annotPath = annots[qq];
		Annotations annot = Annotations.parseAnnotations(annotPath);
		InverseAnnotations invannot = annot.getInvAnnot();
		String[] allgenes = annot.getAllGenes();
		
		new File("output").mkdir();
		String outPath = "/home/weiyi/workspace/javaworks/caf/output/cnv.finetune/erbb2.pbs/" + outdirs[qq];
		new File(outPath).mkdir();
		if(!outPath.endsWith("/")){
			outPath = outPath + "/";
		}
		
		int m = ma.getNumRows();
		float[][] data;
		
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		gn.linkToGeneSet(allgenes);
		
		ArrayList<String> gs = new ArrayList<String>();
		
		long jobID = System.currentTimeMillis();
		Converger cvg = new Converger(0, 1, jobID);
		ITComputer itc = new ITComputer(6, 3, 0, 1, true);
		//itc.negateMI(true);
		cvg.linkITComputer(itc);
		
		gs.add("ERBB2");
		
		//String[] testList = gn.getAllGenesInChrArm(targetArm);
		String[] testList = gn.getNeighbors("ERBB2", 51);
		
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
								
					String[] neighborsG = gn.getNeighbors(gtest, wsize);
					int mg = neighborsG.length;
					String[] neighbors = invannot.getProbes(neighborsG);
					if(neighbors == null){
						System.out.println("No neighbors :(");
						continue;
					}
					DataFile ma2 = ma.getSubProbes(neighbors);
					geneNames = ma2.getProbes();
					m = ma2.getNumRows();
					String[] testPbs = invannot.getProbes(gtest);
					for(String tb : testPbs){
						System.out.println("Probe " + tb);
						int idx = ma2.getRows().get(tb);
						data = ma2.getData();
						float[] vec = data[idx];
						ValIdx[] out = cvg.findWeightedCNV(ma2, vec, gn, neighborsG, invannot, wsize, power, annot);
						if(out[0] == null){
							continue;
						}			
						Arrays.sort(out);
								
						float score = out[quantile-1].val;
						if(score > bestScore){
							bestScore = score;
							bestVec = new ValIdx[mg];
							System.arraycopy(out, 0, bestVec, 0, mg);
							bestWSize = wsize;
						}
						//System.out.println("\tScore: " + score);
					
					} // END for testPbs
					
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
			
			String[] neighborsG = gn.getNeighbors(gtest, bestWSize);
			String[] neighbors = invannot.getProbes(neighborsG);
			DataFile ma2 = ma.getSubProbes(neighbors);
			geneNames = ma2.getProbes();
			PrintWriter pw = new PrintWriter(new FileWriter(outFile));
			for(int i = 0; i < bestWSize; i++){
				String gg = geneNames.get(bestVec[i].idx);
				pw.println(annot.getGene(gg) + "\t"  + bestVec[i].val + "\t" + gg);
			}
			pw.close();
			
			
		}
		pw2.close();
		System.out.println("Final best seed: " + undisputedSeed);

	} // END for qq = 1:4
		
	} // END main

}
