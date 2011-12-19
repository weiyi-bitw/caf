package caf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import obj.Annotations;
import obj.DataFile;
import obj.GeneSet;
import worker.Converger.ValIdx;

public class GroupAttractors {
	
	static GeneSet parseAttractor(String file, HashMap<String, Integer> rowmap) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(file));
		br.readLine();
		String line = br.readLine();
		ArrayList<ValIdx> viList = new ArrayList<ValIdx>();
		while(line != null){
			String[] tokens = line.split("\t");
			int nt = tokens.length;
			int i = rowmap.get(tokens[0]);
			float v = Float.parseFloat(tokens[nt-1]);
			viList.add(new ValIdx(i, v));
			line = br.readLine();
		}
		br.close();
		return new GeneSet(viList.toArray(new ValIdx[0]));
	}
	
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String path = "/home/weiyi/workspace/javaworks/caf/output/207/";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		boolean annotation = true;
		
		System.out.println("Loading files...");
		DataFile ma = DataFile.parse("/home/weiyi/workspace/data/ov/gse9891/ge.54675x285.txt");
		//ma.normalizeRows();
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		HashMap<String, Integer> rowmap = ma.getRows();
		
		String annotPath = "/home/weiyi/workspace/data/annot/affy/u133p2/annot.csv";
		Annotations annot = annotation? Annotations.parseAnnotations(annotPath) : null;
		
		ArrayList<String> probeNames = ma.getProbes();
		if(annotation){
			GeneSet.setAnnotations(annot);
		}
		GeneSet.setProbeNames(probeNames);
		
		
		System.out.print("Loading gene sets...");
		ArrayList<GeneSet> allGeneSet = new ArrayList<GeneSet>();
		ArrayList<Integer> deleteIdx = new ArrayList<Integer>();
		
		String[] attractorFiles = new File(path + "lists/").list();
		int na = attractorFiles.length;
		int progress = 0;
		for(String f : attractorFiles){
			GeneSet gs = parseAttractor(path + "lists/" + f, rowmap);
			gs.setName(f);
			allGeneSet.add(gs);
		}
		int N = allGeneSet.size();
		System.out.println(N + " gene sets are loaded.");
		
		System.out.println("Filtering gene sets...");
		for(int i = 0; i < N; i++){
			progress++;
			boolean delete = false;
			GeneSet gs = allGeneSet.get(i);
			int sz = gs.size();
			float lastMI = gs.getGeneIdx()[gs.size()-1].val();
			for(int j = 0; j < N; j++){
				if(j == i){
					continue;
				}
				GeneSet gs2 = allGeneSet.get(j);
				if(gs.overlapWith(gs2)){
					int sz2 = gs2.size();
					if(sz2 > sz){
						delete=true;
						break;
					}else if(sz2 == sz){
						float lastMI2 = gs2.getGeneIdx()[gs2.size()-1].val();
						if(lastMI2 > lastMI){
							delete = true;
							break;
						}
					}
				}
			}
			if(delete){
				deleteIdx.add(i);
			}
			
			if(progress % 100 == 0) System.out.println(progress + " / " + na);
		}
		Collections.sort(deleteIdx);
		Collections.reverse(deleteIdx);
		for(Integer i : deleteIdx){
			allGeneSet.remove(i.intValue());
		}
		N = allGeneSet.size();
		System.out.println(N + " gene sets are left.");
		System.out.println("Output to files...");
		
		PrintWriter pw = new PrintWriter(new FileWriter(path + "ClusterLeaders.txt"));
		PrintWriter pw4 = new PrintWriter(new FileWriter(path + "forProf.txt"));
		
		
		int numOut = 5;
		
		for(int i = 0; i < N; i++){
			GeneSet gs = allGeneSet.get(i);
			pw.print(gs.getName() + "\t" + gs.size() + "\t");
			if(annotation){
				pw.println(gs.toGenes());
			}else{
				pw.println(gs.toProbes());
			}
			
			ValIdx[] lala = gs.getGeneIdx();
			pw4.print(gs.getName() + "\t" + gs.size());
			for(int j = 0; j < numOut; j++){
				if(annotation){
					pw4.print("\t" + probeNames.get(lala[j].idx()) + " : " + annot.getGene(probeNames.get(lala[j].idx())));
				}else{
					pw4.print("\t" + probeNames.get(lala[j].idx()));
				}
			}
			pw4.println();
			
		}
		pw.close();
		pw4.close();
		/*for(ArrayList<GeneSet> gslist : cluster){
			String name = "Cluster" + String.format("%03d", cnt) + ".txt";
			PrintWriter pw = new PrintWriter(new FileWriter(path + "clusters/" + name));
			pw3.print(name);
			Integer[] member = members.get(cnt);
			for(int i : member){
				if(annotation){
					pw3.print("\t" + probeNames.get(i) + ":" + annot.getGene(probeNames.get(i)));
				}else{
					pw3.print("\t" + probeNames.get(i));
				}
			}
			
			int sz = gslist.size();
			for(int i = 0; i < sz; i++){
				GeneSet gs = gslist.get(i);
				pw.print(gs.getName());
				pw.print("\t" + gs.size() + "\t");
				if(!annotation){
					pw.println(gs.toProbes());
				}else{
					pw.println(gs.toGenes());
				}
			}
			cnt++;
			pw.close();
		}*/
		
		System.out.println("Done.");
	}

}
