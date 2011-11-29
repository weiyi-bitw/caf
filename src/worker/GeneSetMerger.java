package worker;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import obj.GeneSet;

public class GeneSetMerger extends DistributedWorker{
	ArrayList<GeneSet> allGeneSets;
	public static int mergeCount = 0;
	public static int minSize = 0;
	
	public GeneSetMerger(int id, int totalComputers, long jobID){
		super(id, totalComputers, jobID);
		allGeneSets = new ArrayList<GeneSet>();
	}
	
	public void mergeGeneSets(String path, int numFiles, boolean finalOutput) throws IOException{
		if(!path.endsWith("/")) path = path + "/";
		int start = id * numFiles/totalComputers;
		int end = (id+1) * numFiles/totalComputers;
		BufferedReader br;
		System.out.println("Processing file " + start + " to file " + end );
		for(int i = start; i < end; i++){
			System.out.println(i);
			br = new BufferedReader(new FileReader(path + "caf."+ String.format("%05d", i)+".txt"));
			String line = br.readLine();
			// Greedily merge gene set
			while(line != null){
				String[] tokens = line.split("\t");
				if(!tokens[2].equals("NA")){
					//first token: attractees separated by ","
					HashSet<Integer> attr = new HashSet<Integer>();
					String[] t2 = tokens[0].split(",");
					for(String s: t2){
						attr.add(Integer.parseInt(s));
					}
					int nt = tokens.length;
					int[] gIdx = new int[nt-2];
					float[] wts = new float[nt-2];
					int numChild = Integer.parseInt(tokens[1]);
					for(int j = 2; j < nt; j++){
						t2 = tokens[j].split(",");
						gIdx[j-2] = Integer.parseInt(t2[0]);
						wts[j-2] = Float.parseFloat(t2[1]);
					}
					GeneSet rookie = new GeneSet(attr,gIdx, wts, numChild); 
					int origSize = allGeneSets.size();
					if(origSize == 0){
						allGeneSets.add(rookie);
					}else{
						boolean mergeable = false;
						for(int j = 0; j < origSize; j++){
							GeneSet gs = allGeneSets.get(j);
							if(gs.merge(rookie)){
								mergeable = true;
								break;
								// gene set merged
							}
						}
						if(!mergeable){
							allGeneSets.add(rookie);
						}
					}
				}
				line = br.readLine();
			}
			br.close();
		}
		if(finalOutput){
			new File("output").mkdir();
			new File("output/" + jobID).mkdir();
			new File("output/" + jobID + "/lists").mkdir();
			PrintWriter pw = new PrintWriter(new FileWriter("output/" + jobID + "/attractors.gct"));
			PrintWriter pw2 = new PrintWriter(new FileWriter("output/" + jobID + "/attractees.gct"));
			//PrintWriter pw3 = new PrintWriter(new FileWriter("output/" + jobID + "/weights.txt"));
			
			int cnt = 0;
			for(GeneSet gs : allGeneSets){
				if(gs.size() > minSize){
					String name = "Attractor" + String.format("%03d", cnt);
					
					gs.sort();
					gs.calcWeight();
					/*pw3.print(name + "\t" + gs.size() + ":" + gs.getAttracteeSize() + "\t");
					pw3.println(gs.getWeight());*/
					
					pw2.print(name + "\t" + gs.size() + ":" + gs.getAttracteeSize() + "\t");
					pw2.println(gs.getAttractees());
					
					pw.print(name + "\t" + gs.size() + ":" + gs.getAttracteeSize() + "\t");
					if(GeneSet.hasAnnot()){
						pw.println(gs.toGenes());
					}else{
						pw.println(gs.toProbes());
					}
					
					PrintWriter pw4 = new PrintWriter(new FileWriter("output/" + jobID + "/lists/" + name + ".txt"));
					if(GeneSet.hasAnnot()){
						pw4.println("Probe\tGene\rWeight");
						int[] indices = gs.getGeneIdx();
						for(Integer i : indices){
							pw4.println(gs.getOnePair(i) + "\t" + gs.getOneWeight(i));
						}
					}else{
						pw4.println("Gene\tWeight");
						ArrayList<String> geneNames = gs.getGeneNames();
						HashMap<String, Float> geneWeightMap = gs.getGeneWeightMap();
						for(String s: geneNames){
							pw4.println(s + "\t" + geneWeightMap.get(s));
						}
					}
					
					pw4.close();
					cnt++;
				}
			}
			//pw3.close();
			pw2.close();
			pw.close();
		}else{
			prepare("merge" + mergeCount);
			PrintWriter pw = new PrintWriter(new FileWriter("tmp/" + jobID + "/merge" + mergeCount + "/caf."+ String.format("%05d", id)+".txt"));
			for(GeneSet gs : allGeneSets){
				pw.println(gs.toString());
			}
			pw.close();
			mergeCount++;
		}
	}
	public void setMinSize(int minSize){
		GeneSetMerger.minSize = minSize;
	}
	
}
