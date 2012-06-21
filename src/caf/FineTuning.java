package caf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import obj.DataFile;
import obj.ValIdx;
import worker.Converger;
import worker.ITComputer;

public class FineTuning {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		double start = 3.5;
		int quantile = 10;
		
		System.out.println("Loading files...");
		
		//String outPath = "/home/weiyi/workspace/javaworks/caf/output/656/";
		new File("output").mkdir();
		String outPath = "/home/weiyi/workspace/javaworks/caf/output/caf.finetune/test/";
		if(!outPath.endsWith("/")){
			outPath = outPath + "/";
		}
		final String[] dataFiles={
				//"/home/weiyi/workspace/data/gbm/tcga/ge/ge.12042x545.txt"
				//"/home/weiyi/workspace/data/brca/gse3143/ge.8443x158.jetset.mean.txt",
				//"/home/weiyi/workspace/data/brca/gse3494/ge.12160x251.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/brca/gse32646/ge.19190x115.jetset.ncbi.txt",
				"/home/weiyi/workspace/data/brca/gse36771/ge.19190x107.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/brca/gse31448/ge.19190x353.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/brca/gse2034/ge.12160x286.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/brca/tcga/ge/ge.17475x536.ncbi.txt",
				//"/home/weiyi/workspace/data/coad/gse14333/ge.19190x290.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/coad/tcga/ge/ge.17475x154.ncbi.txt",
				//"/home/weiyi/workspace/data/ov/gse9891/ge.19190x285.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/ov/tcga/ge/ge.11963x582.ncbi.txt",
				//"/home/weiyi/workspace/data/ov/gse26193/ge.19189x107.jetset.mean.txt",
				//"/home/weiyi/workspace/data/prad/gse17951/ge.19189x154.jetset.mean.txt",
				//"/home/weiyi/workspace/data/prad/gse8218/ge.12160x148.jetset.mean.txt",
				//"/home/weiyi/workspace/data/ov/tcga/super.35696x511.knn.txt",
				//"/home/weiyi/workspace/data/gbm/tcga/super.40092x274.txt"
		};
			
		final String[] outputDirs={
				//"gbm.tcga"
				//"brca.gse3494.jetset.mean",
				//"brca.gse32646.jetset.mean",
				"brca.gse36771.jetset.mean",
				//"brca.gse31448.jetset.mean",
				//"brca.gse2034.jetset.ncbi",
				//"brca.tcga.ncbi",
				//"coad.gse14333.jetset.ncbi",
				//"coad.tcga.ncbi",
				//"ov.gse9891.jetset.ncbi",
				//"ov.tcga.ncbi",
				//"ov.gse26193",
				//"prad.gse17951",
				//"prad.gse8218",
				//"ov.super",
				//"gbm.super"
		};
			
		
		for(int qq = 0; qq < dataFiles.length; qq++)
		{
		System.out.println("Data: " + dataFiles[qq]);
		
		DataFile ma = DataFile.parse(dataFiles[qq]);
		
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		
		ArrayList<String> gs = new ArrayList<String>();
		
		long jobID = System.currentTimeMillis();
		Converger cvg = new Converger(0, 1, jobID);
		ITComputer itc = new ITComputer(6, 3, 0, 1, true);
		itc.negateMI(true);
		cvg.linkITComputer(itc);
		
		HashMap<String, Integer> geneMap = ma.getRows();
		
		gs.add("ESR1");
		//gs.add("FABP4");
	
		for(String g : gs){
			ArrayList<String> geneNames = ma.getProbes();
			
			if(geneNames.contains(g)){
				System.out.println("Processing " + g + " (" + (gs.indexOf(g)+1) + "/" + gs.size() + ")" + "...");

				float[] vec = new float[n];
				int idx = geneMap.get(g);
				vec = data[idx];
				double[] bestPower = {start};
				double[] tmp = cvg.AttractorScanning(ma, vec, start, bestPower);
					
				if(tmp[0] == -1){
					continue;
				}
						
				ValIdx[] vis = new ValIdx[m];
				for(int i = 0; i < m; i++){
					vis[i] = new ValIdx(i, (float)tmp[i]);
				}
				Arrays.sort(vis);
						
				new File(outPath + "mergeroom").mkdir();
				String outFile = outPath + "mergeroom/" + outputDirs[qq] + "_" + g + "_" + bestPower[0] + ".txt";
				PrintWriter pw = new PrintWriter(new FileWriter(outFile));
				pw.println("Power:\t" + bestPower[0]);
				for(int i = 0; i < m; i++){
					String gg = geneNames.get(vis[i].idx);
					pw.println(gg + "\t"  + vis[i].val);
				}
				pw.close();
						
				float score = vis[quantile-1].val;
				System.out.println("\tScore: " + score);
					
			}else{
				System.out.println("Does not contain gene " + g + "!!");
			}
		}
		
		}
	}

}
