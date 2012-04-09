package caf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

import obj.DataFile;
import obj.GeneSet;
import obj.ValIdx;
 
public class CreateAttractorDataset {

	static ArrayList<GeneSet> parseAttractorInOneFile(String file, HashMap<String, Integer> rowmap) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = br.readLine();
		ArrayList<GeneSet> allAttractors = new ArrayList<GeneSet>();
		while(line != null){ // each line is an attractor
			String[] tokens = line.split("\t");
			String name = tokens[0];
			int nt = tokens.length;
			ArrayList<ValIdx> viList = new ArrayList<ValIdx>();
			for(int j = 2; j < nt; j++){
				String[] t2 = tokens[j].split(":");
				int i = rowmap.get(t2[0]);
				float v = Float.parseFloat(t2[1]);
				viList.add(new ValIdx(i, v));
				
			}
			// no number of child information in ClusterLeaders file, set 1
			allAttractors.add(new GeneSet(name, viList.toArray(new ValIdx[0]), 1));
			line = br.readLine();
		}
		br.close();
		return allAttractors;
	}
	private static float[] getMetaGene(float[][] data, ValIdx[] idx, int n){
		int m = idx.length;
		float[] out = new float[n];
		for(int j = 0; j < n; j++){
			for(ValIdx vi : idx){
				out[j] += data[vi.idx()][j];
			}
			out[j] /= m;
		}
		return out;
	}
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String path = "/home/weiyi/workspace/javaworks/caf/output/ov.tcga.affyL3";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		System.out.println("Loading files...");
		DataFile ma = DataFile.parse("/home/weiyi/workspace/data/ov/tcga/ge/ge.12042x582.txt");
		//ma.normalizeRows();
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		
		ArrayList<GeneSet> attractors = parseAttractorInOneFile(path + "ClusterLeaders.txt", ma.getRows());
		int k = attractors.size();
		
		float[][] newData = new float[k][n];
		ArrayList<String> attractorNames = new ArrayList<String>();
		HashMap<String, Integer> newRowMap = new HashMap<String, Integer>();
		
		int cnt = 0;
		System.out.println("Calculating metagenes...");
		for(GeneSet gs: attractors){
			attractorNames.add(gs.getName());
			newRowMap.put(gs.getName(),cnt);
			System.arraycopy(getMetaGene(data, gs.getGeneIdx(), n), 0, newData[cnt], 0, n);
			cnt++;
		}
		
		DataFile out = new DataFile(newData, newRowMap, ma.getCols(), attractorNames, ma.getChipID());
		System.out.println("Output to files...");
		out.output2Gct(path + "attractorSpace." + out.getNumRows() + "x" + out.getNumCols() + ".txt");
		
		System.out.println("Done.");
	}

}
