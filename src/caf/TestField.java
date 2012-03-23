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
import java.util.HashSet;

import org.apache.commons.math.distribution.NormalDistributionImpl;

import obj.Annotations;
import obj.DataFile;
import obj.Genome;
import obj.ValIdx;

import util.StatOps;
import worker.Converger;
import worker.ITComputer;

public class TestField {
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
			if(maxIdx == (i + winSize-1)){
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
	
	
	
/*	static class ValIdx implements Comparable<ValIdx>{
		float val;
		int idx;
		ValIdx(int i, float v){
			this.idx = i;
			this.val = v;
		}
		
		public int compareTo(ValIdx other) {
			return -Double.compare(this.val, other.val);
		}
	}*/
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String path = "/home/weiyi/workspace/data/coad/gse14333/";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		
		System.out.println("Loading files...");
		
		//String outPath = "/home/weiyi/workspace/javaworks/caf/output/656/";
		String outPath = "/home/weiyi/workspace/javaworks/caf/tmp/";
		if(!outPath.endsWith("/")){
			outPath = outPath + "/";
		}
		DataFile ma = DataFile.parse(path + "ge.19964x290.zmean.txt");
		
		
		final String geneLocFile = "/home/weiyi/workspace/data/annot/affy/u133p2/gene.location3";
		//final String geneLocFile = "/home/weiyi/workspace/javaworks/caf/output/639/gene.location3";
		
		String command = "CAF";
		float power = 5f;
		boolean excludeTop = false;
		boolean miDecay = false;
		int winSize = 51;
		
		//ma.normalizeRows();
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		
		ArrayList<String> gs = new ArrayList<String>();
		/*BufferedReader br = new BufferedReader(new FileReader("COL11A1_50"));
		br.readLine();
		String line = br.readLine();
		while(line != null){
			String[] tokens = line.split("\t");
			gs.add(tokens[0]);
			line = br.readLine();
		}
		br.close();*/
		
		long jobID = System.currentTimeMillis();
		
		//String annotPath = "/home/weiyi/workspace/data/annot/affy/u133p2/annot.csv";
		//Annotations annot = Annotations.parseAnnotations(annotPath);
		Annotations annot = null;
		
		Converger cvg = new Converger(0, 1, jobID);
		
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		if(command.equals("CNV")) gn.linkToDataFile(ma);
		
		
		gs.add("CENPA");
		/*gs.add("EXOSC10");
		gs.add("ASH2L");
		gs.add("EXOSC4");*/
		//gs.addAll(slidingWindowSelector(outPath + "basinScores.txt", 200));
		
		/*String[] chr8 = gn.getNeighbors("ASH2L", -1);
		for(String s : chr8){
			gs.add(s);
		}
		
		int k = chr8.length;
		float scores[] = new float[k];*/
		
		ITComputer itc = new ITComputer(6, 3, 0, 1, true);
		//itc.negateMI(true);
		cvg.linkITComputer(itc);
		HashMap<String, Integer> geneMap = ma.getRows();
		
		int cnt = 0;
		new File("tmp").mkdir();
		for(String g : gs){
			ArrayList<String> geneNames = ma.getProbes();
			if(geneNames.contains(g)){
				System.out.println("Processing " + g + " (" + (gs.indexOf(g)+1) + "/" + gs.size() + ")" + "...");
				String outFile = command.equals("CAF")? outPath + g + "_Attractor.txt" : outPath + g + "_CNV.txt";
				try{
					PrintWriter pw = new PrintWriter(new FileWriter(outFile));
					float[] out = new float[m];
					float[] vec = new float[m];
					int idx = geneMap.get(g);
					if(command.equals("CAF")){
						vec = data[idx];
						out = cvg.findWeightedAttractor(ma, g, vec, power, excludeTop);
					}else if(command.equals("CNV")){
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
						out = cvg.findWeightedCNV(ma2, g, gn, vec, winSize, power,  excludeTop, miDecay);
					}else{
						System.out.println("Wrong command!");
						System.exit(1);
					}
					
					
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
					//scores[cnt] = vi.get(9).val;
					cnt++;
					for(int i = 0; i < m; i++){
						String gg = geneNames.get(vi.get(i).idx);
						if(command.equals("CNV"))pw.println(gg + "\t" + gn.getIdx(gg) +  "\t" + vi.get(i).val);
						else pw.println(gg + "\t"  + vi.get(i).val);
					}
					pw.close();
					
				}catch (FileNotFoundException e){
					System.out.println("Exception: " + e);
					continue;
				}
			}else{
				System.out.println("Does not contain gene " + g + "!!");
			}
		}
		
		/*PrintWriter pw = new PrintWriter(new FileWriter("tmp/FixedWindowScores.txt"));
		for(int i = 0; i < k; i++){
			pw.println(gs.get(i) + "\t" + scores[i]);
		}pw.close();*/
		
		
		/*cvg.setZThreshold(8f);
		cvg.setConvergeMethos("ZSCORE");
		cvg.setMIParameter(6, 3);*/
		
		/*HashMap<String, Integer> geneMap = ma.getRows();
		ArrayList<String> attractees = new ArrayList<String>();
		ArrayList<ArrayList<ValIdx>> attractors = new ArrayList<ArrayList<ValIdx>>();
		ArrayList<String> geneNames = ma.getProbes();
		int cnt = 0;
		for(String g : gs){
			System.out.println(g + "...");
			int idx = geneMap.get(g);
			ArrayList<ValIdx> out = cvg.findAttractor(data, idx);
			if(attractors.contains(out)){
				int j = attractors.indexOf(out);
				String s = attractees.get(j);
				s = s + "," + g;
				attractees.set(j, s);
			}else{
				attractees.add(g);
				attractors.add(out);
			}
			cnt = cnt + 1;
			if(cnt == 10){
				break;
			}
		}
		
		PrintWriter pw = new PrintWriter(new FileWriter("lala.txt"));
		int kk = attractors.size();
		for(int i = 0; i < kk; i++){
			pw.print(attractees.get(i));
			ArrayList<ValIdx> vv = attractors.get(i);
			//Collections.sort(vv);
			for(ValIdx vj : vv){
				pw.print("\t" + geneNames.get(vj.idx()));
				if(annot != null){
					pw.print(":" + annot.getGene(geneNames.get(vj.idx())));
				}
				pw.print(":" + vj.val());
			}
			pw.println();
		}
		pw.close();*/
		
		/*String path = "/home/weiyi/workspace/javaworks/caf/output/207/";
		
		BufferedReader br = new BufferedReader(new FileReader(path + "leadingAttractors.txt"));
		ArrayList<String> leadingAttractors = new ArrayList<String>();
		String line = br.readLine();
		while(line != null){
			leadingAttractors.add(line.replaceAll(".txt", ""));
			line = br.readLine();
		}
		br.close();
		
		br = new BufferedReader(new FileReader(path + "attractors.gwt"));
		PrintWriter pw = new PrintWriter(new FileWriter(path + "leadingAttractors.detail.txt"));
		line = br.readLine();
		while(line != null){
			String name = line.split("\t")[0];
			if(leadingAttractors.contains(name)){
				pw.println(leadingAttractors.indexOf(name) + "\t" + line);
			}
			line = br.readLine();
		}
		pw.close();
		br.close();*/
		
		System.out.println("Done.");
	}

}
