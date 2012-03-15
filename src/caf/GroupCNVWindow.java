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
	private static ArrayList<String> slidingWindowSelector(String inFileName, int winSize, int topAttractors) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(inFileName));
		ArrayList<String> genesAl = new ArrayList<String>();
		ArrayList<ValIdx> scoresAl = new ArrayList<ValIdx>();
		String line = br.readLine();
		int cnt = 0;
		
		while(line != null){
			String[] tokens = line.split("\t");
			genesAl.add(tokens[0]);
			scoresAl.add(new ValIdx(cnt,Float.parseFloat(tokens[1])));
			line = br.readLine();
			cnt ++;
		}
		br.close();
		
		int k = genesAl.size();
		Collections.sort(scoresAl);
		ArrayList<String> outG = new ArrayList<String>();
		boolean[] eliminate = new boolean[k];
		
		for(ValIdx vi : scoresAl){
			int idx = vi.idx;
			if(!eliminate[idx]){
				for(int i = (idx - winSize + 1); i <= (idx + winSize -1); i++){
					if(i >= 0 && i < k){
						eliminate[i] = true;
					}
				}
				outG.add(genesAl.get(idx));
				if(topAttractors > 0){
					if(outG.size() >= topAttractors){
						break;
					}
				}
			}
		}
		
		
		return outG;
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		//final String dataFile = "/home/weiyi/workspace/data/brca/gse2034/ge.13271x286.var.txt";
		//final String dataFile = "/home/weiyi/workspace/data/brca/tcga/ge/ge.17814x536.knn.txt";
		//final String dataFile = "/home/weiyi/workspace/data/coad/gse14333/ge.20765x290.var.txt";
		final String dataFile = "/home/weiyi/workspace/data/coad/tcga/ge/ge.17814x154.knn.txt";
		//final String dataFile = "/home/weiyi/workspace/data/ov/gse9891/ge.20765x285.var.txt";
		//final String dataFile = "/home/weiyi/workspace/data/ov/tcga/ge/ge.17814x584.knn.txt";
		
		System.out.println("Loading files...");
		DataFile ma = DataFile.parse(dataFile);
		
		String outPath = "/home/weiyi/workspace/javaworks/caf/output/window51/coad.tcga/";
		//String outPath = "/home/weiyi/workspace/javaworks/caf/tmp/";
		if(outPath.endsWith("/")){
			outPath = outPath.substring(0, outPath.length()-1);
		}
		
		final String geneLocFile = "/home/weiyi/workspace/data/annot/affy/u133p2/gene.location3";
		//final String geneLocFile = "/home/weiyi/workspace/javaworks/caf/output/639/gene.location3";
		
		float power = 2f;
		boolean excludeTop = false;
		boolean miDecay = false;
		int winSize = 51;
		
		//ma.normalizeRows();
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		
		ArrayList<String> gs = new ArrayList<String>();
		
		long jobID = System.currentTimeMillis();
		
		Converger cvg = new Converger(0, 1, jobID);
		
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		gn.linkToDataFile(ma);
		
		gs.addAll(slidingWindowSelector(outPath + "/basinScores.txt", 101, -1));
		/*String[] nbs = gn.getNeighbors("PUF60", 20);
		for(String s : nbs){
			gs.add(s);
		}*/
		
		ITComputer itc = new ITComputer(6, 3, 0, 1, true);
		cvg.linkITComputer(itc);
		HashMap<String, Integer> geneMap = ma.getRows();
		
		int cnt = 0;
		//new File("tmp").mkdir();
		String outFileName = outPath.substring(outPath.lastIndexOf("/"));
		PrintWriter pw = new PrintWriter(new FileWriter(outPath + "/../mergeroom/" + outFileName));
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
					out = cvg.findWeightedCNV(ma2, g, gn, vec, winSize, power, excludeTop, miDecay);
					
					
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
					for(int i = 0; i< 20; i++){
						String gg = geneNames.get(vi.get(i).idx);
						pw.print("\t" + gg);
					}
					pw.print("\t" +gn.getChrArm(g) + "\t" +  vi.get(9).val);
					pw.println();
					
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
