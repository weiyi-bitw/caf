package caf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import obj.DataFile;
import obj.Genome;
import obj.ValIdx;
import worker.Converger;
import worker.ITComputer;

public class GroupCNVWindow {

	private static float[] getMetaGene(float[][] data, ArrayList<Integer> idx, int n){
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
	private static ArrayList<Integer> getAttractorIdx(float[][] data, float[][] val, ArrayList<Integer> inputIdx, int m, int n, ITComputer itc, int attractorSize, int maxIter) throws Exception{
		ArrayList<Integer> metaIdx = inputIdx;
		ArrayList<Integer> preMetaIdx = new ArrayList<Integer>();
		ArrayList<Integer> cycMetaIdx = new ArrayList<Integer>();
		
		preMetaIdx.addAll(inputIdx);
		cycMetaIdx.addAll(inputIdx);
		
		float[] metagene = getMetaGene(data, metaIdx, n);
		ValIdx[] vec = new ValIdx[m];
		int cnt = 0;
		
		
		while(cnt < maxIter){
			float[] mi = itc.getAllMIWith(metagene, val);
			for(int i = 0; i < m; i++){
				vec[i] = new ValIdx(i, mi[i]);
			}
			Arrays.sort(vec);
			
			System.out.println(cnt + "\t" + vec[0].idx() + "\t" + vec[attractorSize/2].idx() + "\t" 
					+ "\t" + vec[attractorSize-1].idx());
			
			metaIdx.clear();
			for(int i = 0; i < attractorSize; i++){
				metaIdx.add(vec[i].idx());
			}
			if(metaIdx.equals(preMetaIdx) || metaIdx.equals(cycMetaIdx)){
				break;
			}
			cycMetaIdx.clear();
			cycMetaIdx.addAll(preMetaIdx);
			preMetaIdx.clear();
			preMetaIdx.addAll(metaIdx);
			cnt++;
			metagene = getMetaGene(data, metaIdx, n);
			
		}
		
		return metaIdx;
	}
	private static ArrayList<String> slidingWindowSelector(String inFileName, int winSize) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(inFileName));
		ArrayList<String> genesAl = new ArrayList<String>();
		ArrayList<Float> scoresAl = new ArrayList<Float>();
		String line = br.readLine();
		while(line != null){
			String[] tokens = line.split("\t");
			genesAl.add(tokens[0]);
			scoresAl.add(Float.parseFloat(tokens[1]));
			line = br.readLine();
		}
		
		int k = genesAl.size();
		ArrayList<String> outG = new ArrayList<String>();
		for(int i = 0; i <= k-winSize; i+=20){
			int maxIdx = -1;
			float maxF = -1;
			for(int j = 0; j < winSize; j++){
				if(scoresAl.get(i+j) > maxF){
					maxF = scoresAl.get(i+j);
					maxIdx = (i+j);
				}
			}
			if(maxIdx == (i + winSize-1) || maxIdx == i){
				continue;
			}
			String g = genesAl.get(maxIdx);
			if(!outG.contains(g)){
				outG.add(g);
			}
		}
		br.close();
		return outG;
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String path = "/home/weiyi/workspace/data/brca/gse2034/";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		
		System.out.println("Loading files...");
		
		//String outPath = "/home/weiyi/workspace/javaworks/caf/output/window41.brca.gse2034/";
		String outPath = "/home/weiyi/workspace/javaworks/caf/tmp/";
		if(!outPath.endsWith("/")){
			outPath = outPath + "/";
		}
		DataFile ma = DataFile.parse(path + "ge.13271x286.var.txt");
		
		final String geneLocFile = "/home/weiyi/workspace/data/annot/affy/u133p2/gene.location3";
		//final String geneLocFile = "/home/weiyi/workspace/javaworks/caf/output/639/gene.location3";
		
		float power = 2f;
		boolean excludeTop = false;
		boolean miDecay = false;
		int winSize = 41;
		
		//ma.normalizeRows();
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		
		ArrayList<String> gs = new ArrayList<String>();
		
		long jobID = System.currentTimeMillis();
		
		Converger cvg = new Converger(0, 1, jobID);
		
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		gn.linkToDataFile(ma);
		
		//gs.addAll(slidingWindowSelector(outPath + "basinScores.txt", 200));
		String[] nbs = gn.getNeighbors("PUF60", 20);
		for(String s : nbs){
			gs.add(s);
		}
		
		ITComputer itc = new ITComputer(6, 3, 0, 1, true);
		cvg.linkITComputer(itc);
		HashMap<String, Integer> geneMap = ma.getRows();
		
		int cnt = 0;
		new File("tmp").mkdir();
		PrintWriter pw = new PrintWriter(new FileWriter(outPath + "windowedCNVs.txt"));
		for(String g : gs){
			ArrayList<String> geneNames = ma.getProbes();
			if(geneNames.contains(g)){
				System.out.println("Processing " + g + " (" + (gs.indexOf(g)+1) + "/" + gs.size() + ")" + "...");
					float[] out = new float[m];
					float[] vec = new float[m];
					int idx = geneMap.get(g);
					String[] neighbors = gn.getNeighbors(g, winSize);
					if(neighbors == null){
						/*pw.println("No neighbors");
						pw.close();*/
						continue;
					}
					DataFile ma2 = ma.getSubProbes(neighbors);
					geneNames = ma2.getProbes();
					m = ma2.getNumRows();
					idx = ma2.getRows().get(g);
					data = ma2.getData();
					vec = data[idx];
					out = cvg.findWeightedCNV(ma2, g, gn, vec, power, excludeTop, miDecay);
					
					
					if(out[0] == -1){
						/*pw.println("Not converged.");
						pw.close();*/
						continue;
					}
					
					ArrayList<ValIdx> vi = new ArrayList<ValIdx>();
					for(int i = 0; i < m; i++){
						//out[i] = (float) Math.round(out[i] * 10000) / 10000;
						vi.add(new ValIdx(i, out[i]));
					}
					Collections.sort(vi);
					
					pw.print("Attractor_" + g + "\t" + gn.getChr(geneNames.get(vi.get(0).idx)));
					for(int i = 0; i< m; i++){
						String gg = geneNames.get(vi.get(i).idx);
						pw.print("\t" + gg + "(" + gn.getIdx(gg) + "):" + vi.get(i).val);
					}pw.println();
					
					//scores[cnt] = vi.get(9).val;
					cnt++;
					
			}else{
				System.out.println("Does not contain gene " + g + "!!");
			}
		}
		pw.close();
		System.out.println("Done.");
	}

}
