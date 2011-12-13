package caf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.math.distribution.HypergeometricDistributionImpl;

public class MatchMaker {
	static class GeneSet implements Comparable<GeneSet>{
		static ArrayList<String> genes;
		String id;
		float[] vec;
		int sz;
		float density;
		
		GeneSet(String id, float[] vec, int sz, float density){
			this.id = id;
			this.vec = vec;
			this.sz = sz;
			this.density = density;
		}
		
		double phyperWith(GeneSet gs2){
			HypergeometricDistributionImpl hyper = new HypergeometricDistributionImpl(genes.size(), sz, gs2.sz);
			float[] vec2 = gs2.vec;
			int cnt = 0;
			int n = vec.length;
			for(int i = 0; i < n; i++){
				if(vec[i] > 0 && vec2[i] > 0){
					cnt++;
				}
			}
			return hyper.upperCumulativeProbability(cnt);
			//return cnt;
		}
		
		static float cosineSimilarity(GeneSet gs1, GeneSet gs2) throws Exception{
			float xsq=0;
			float ysq=0;
			float xy = 0;
			
			int n = gs1.vec.length;
			if(n != gs2.vec.length){
				throw new RuntimeException("ERROR: Two gene sets have different vector size!!");
			}
			float[] x = gs1.vec;
			float[] y = gs2.vec;
			for(int i = 0; i < n; i++){
				xsq += x[i] * x[i];
				ysq += y[i] * y[i];
				xy += x[i] * y[i];
			}
			
			return (float) (xy / Math.sqrt(xsq) / Math.sqrt(ysq));
		}
		float cosineSimilarityWith(GeneSet gs2) throws Exception{
			float xsq=0;
			float ysq=0;
			float xy = 0;
			int n = vec.length;
			
			if(n != gs2.vec.length){
				throw new RuntimeException("ERROR: Two gene sets have different vector size!!");
			}
			float[] x = this.vec;
			float[] y = gs2.vec;
			for(int i = 0; i < n; i++){
				xsq += x[i] * x[i];
				ysq += y[i] * y[i];
				xy += x[i] * y[i];
			}
			
			return (float) (xy / Math.sqrt(xsq) / Math.sqrt(ysq));
		}
		
		static void setGeneSpace(ArrayList<String> genes){
			GeneSet.genes = genes;
		}
		
		public int compareTo(GeneSet other) {
			return -Double.compare(this.density, other.density);
		}
	}
	
	static ArrayList<GeneSet> parse(String file, HashMap<String, Integer> gMap, int dim) throws Exception{
		ArrayList<GeneSet> gs = new ArrayList<GeneSet>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = br.readLine();
		
		while(line!= null){
			String[] tokens = line.split("\t");
			String id = tokens[0];
			String[] t2 = tokens[1].split(":");
			int sz = 0;
			int nt = tokens.length;
			float[] w = new float[dim];
			for(int i = 2; i < nt; i++){
				t2 = tokens[i].split(":");
				int idx;
				if(gMap.get(t2[0])!=null){
					sz++;
					idx = gMap.get(t2[0]);
					if(t2.length > 1){
						w[idx] = Float.parseFloat(t2[1]);
					}else{
						w[idx] = 1f;
					}
				}
				
			}
			float density = Float.parseFloat(t2[1])/sz;
			
			gs.add(new GeneSet(id, w, sz, density));
			line = br.readLine();
		}
		Collections.sort(gs);
		return gs;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{
		long tOrigin = System.currentTimeMillis();
		
		String inputFile1 = args[0];
		String inputFile2 = args[1];
		String geneSpaceFile = args[2]; // usually annotation files
		
		System.out.printf("%-25s%s\n", "File 1:" , inputFile1);
		System.out.printf("%-25s%s\n", "File 2:" , inputFile2);
		System.out.printf("%-25s%s\n", "Gene Space:" , geneSpaceFile);
		
		ArrayList<String> genes = new ArrayList<String>();
		HashMap<String, Integer> gMap = new HashMap<String, Integer>();
		BufferedReader br = new BufferedReader(new FileReader(geneSpaceFile));
		String line = br.readLine();
		int cnt = 0;
		while(line != null){
			String[] tokens = line.split(",");
			if(gMap.get(tokens[1]) == null){
				genes.add(tokens[1]);
				gMap.put(tokens[1], cnt);
				cnt++;
			}
			line = br.readLine();
		}
		GeneSet.setGeneSpace(genes);
		int spaceSize = genes.size();
		System.out.printf("%-25s%s\n", "Gene Space Size:" , spaceSize);
		
		System.out.println("Parsing attractors...");
		ArrayList<GeneSet> gss1 = parse(inputFile1, gMap, spaceSize);
		ArrayList<GeneSet> gss2 = parse(inputFile2, gMap, spaceSize);
		
		System.out.println("Calculating Distances...");
		
		if(gss2.size() > gss1.size()){
			ArrayList<GeneSet> tmp = gss2;
			gss2 = gss1;
			gss1 = tmp;
			String tmp2 = inputFile2;
			inputFile2 = inputFile1;
			inputFile1 = tmp2;
		}
		
		PrintWriter pw = new PrintWriter(new FileWriter("match.txt"));
		pw.println(inputFile1 + "\t" + inputFile2 + "\tDistance(Cosine)\tP-value");
		
		
		for(GeneSet gs : gss1){
			pw.print(gs.id);
			GeneSet closest = null;
			float closestSim = 0;
			
			for(GeneSet gs2:gss2){
				float d = gs.cosineSimilarityWith(gs2);
				if(d > closestSim){
					closest = gs2;
					closestSim = d;
				}
			}
			
			if(closest != null){
				pw.print("\t" + closest.id + "\t" + (1-closestSim) + "\t" + gs.phyperWith(closest));
				//gss2.remove(closest);
			}else{
				pw.print("\t" + "NA");
			}
			pw.println();
		}
		
		/*for(GeneSet gs : gss2){
			pw.println("NA\t" + gs.id);
		}*/
		
		pw.close();
		
		System.out.println("Done in " + (System.currentTimeMillis() - tOrigin) + " msecs. <3");
	}

}
