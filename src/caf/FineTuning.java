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
		float start = 5;
		float end = 5;
		int quantile = 20;
		
		System.out.println("Loading files...");
		
		//String outPath = "/home/weiyi/workspace/javaworks/caf/output/656/";
		new File("output").mkdir();
		String outPath = "/home/weiyi/workspace/javaworks/caf/output/caf.finetune/cenpa/";
		if(!outPath.endsWith("/")){
			outPath = outPath + "/";
		}
		final String dataFile = "/home/weiyi/workspace/data/brca/gse2034/ge.12764x286.mean.txt";
		//final String dataFile = "/home/weiyi/workspace/data/brca/tcga/ge/ge.17814x536.knn.txt";
		//final String dataFile = "/home/weiyi/workspace/data/coad/gse14333/ge.19189x290.geo.jetset.mean.txt";
		//final String dataFile = "/home/weiyi/workspace/data/coad/tcga/ge/ge.17814x154.knn.txt";
		//final String dataFile = "/home/weiyi/workspace/data/ov/gse9891/ge.19177x285.geo.jetset.mean.txt";
		//final String dataFile = "/home/weiyi/workspace/data/ov/tcga/ge/ge.17814x584.knn.txt";
		//final String dataFile = "/home/weiyi/workspace/data/ov/tcga/ge/ge.12042x582.txt";
		DataFile ma = DataFile.parse(dataFile);
		
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		
		ArrayList<String> gs = new ArrayList<String>();
		
		long jobID = System.currentTimeMillis();
		Converger cvg = new Converger(0, 1, jobID);
		ITComputer itc = new ITComputer(6, 3, 0, 1, true);
		//itc.negateMI(true);
		cvg.linkITComputer(itc);
		HashMap<String, Integer> geneMap = ma.getRows();
		
		gs.add("CENPA");
		
		for(String g : gs){
			ArrayList<String> geneNames = ma.getProbes();
			
			if(geneNames.contains(g)){
				System.out.println("Processing " + g + " (" + (gs.indexOf(g)+1) + "/" + gs.size() + ")" + "...");
				ValIdx[] out = new ValIdx[m];
				float bestPower = -1;
				float bestScore = -1;
				float power = start;
				while(true){
					if(power < end){
						break;
					}
					System.out.print("Power " + power + "...");
					try{
						
						
						float[] vec = new float[n];
						int idx = geneMap.get(g);
						vec = data[idx];
						float[] tmp = cvg.findWeightedAttractor(ma, g, vec, power, false);
						
						if(tmp[0] == -1){
							/*pw.println("Not converged.");
							pw.close();*/
							power -= 0.5;
							continue;
						}
						ValIdx[] vis = new ValIdx[m];
						for(int i = 0; i < m; i++){
							vis[i] = new ValIdx(i, tmp[i]);
						}
						Arrays.sort(vis);
						
						String outFile = outPath + g + "_Attractor" + power + ".txt";
						PrintWriter pw = new PrintWriter(new FileWriter(outFile));
						for(int i = 0; i < m; i++){
							String gg = geneNames.get(vis[i].idx);
							pw.println(gg + "\t"  + vis[i].val);
						}
						pw.close();
						
						float score = vis[quantile-1].val;
						/*float score = 0;
						for(int i = 0; i < quantile; i++){
							score += vis[i].val;
						}*/
						System.out.println("\tScore: " + score);
						power -= 0.5;
						/*if(score > bestScore){
							bestScore = score;
							bestPower = power;
							System.arraycopy(vis, 0, out, 0, m);
							power -= 0.5;
						}else{
							power -= 0.5;
							if(power == 10){
								power = 15;
								continue;
							}
							break;
						}*/
						
					}catch (FileNotFoundException exc){
						System.out.println("Exception: " + exc);
						continue;
					}
				
				}
				
				
				
				
				
			}else{
				System.out.println("Does not contain gene " + g + "!!");
			}
		}
		
	}

}
