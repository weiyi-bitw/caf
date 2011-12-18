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
		String path = "/home/weiyi/workspace/javaworks/caf/output/206/";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		boolean annotation = false;
		
		System.out.println("Loading files...");
		DataFile ma = DataFile.parse("/home/weiyi/workspace/data/ov/tcga/mergeroom/ge.17814x514.common.txt");
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
		
		System.out.println("Filtering gene sets...");
		ArrayList<GeneSet> leaders = new ArrayList<GeneSet>();
		ArrayList<String> clustMembers = new ArrayList<String>();
		ArrayList<HashSet<Integer>> clustGeneIdx = new ArrayList<HashSet<Integer>>();
		
		String[] attractorFiles = new File(path + "lists/").list();
		int na = attractorFiles.length;
		int progress = 0;
		for(String f : attractorFiles){
			progress++;
			GeneSet gs = parseAttractor(path + "lists/" + f, rowmap);
			gs.setName(f);
			
			int nc = leaders.size();
			boolean mergeable = false;
			int replaceIdx = -1;
			int convergeIdx = -1;
			int convergeSize = -1;
			float convergeMI = -1;
			ArrayList<Integer> deleteIdx = new ArrayList<Integer>();
			boolean willBeDeleted = false;
			for(int i = 0; i < nc; i++){
				GeneSet gs2 = leaders.get(i);
				if(gs.overlapWith(gs2)){
					mergeable = true;
					System.out.println("Find mergeable: " + i);
					if(gs.size() > gs2.size()){
						System.out.println("Larger size.");
						if(replaceIdx >= 0){// already going to replace, delete this position
							deleteIdx.add(i);
						}else if(willBeDeleted){ // this gs is already going to be deleted, delete this position also
							deleteIdx.add(i);
						}else{
							replaceIdx = i;
							convergeIdx = replaceIdx;
							convergeSize = gs.size();
							convergeMI = gs.getGeneIdx()[gs.size()-1].val();
						}
					}else if(gs.size() == gs2.size()){ // tie; compare the last MI
						System.out.println("Equal size.");
						float v1 = gs.getGeneIdx()[gs.size()-1].val();
						float v2 = gs2.getGeneIdx()[gs2.size()-1].val();
						if(v1 > v2){
							System.out.println("higher MI.");
							if(replaceIdx >= 0){// already going to replace, delete this position
								deleteIdx.add(i);
								convergeIdx = replaceIdx;
							}else if(willBeDeleted){ // this gs is already going to be deleted, delete this position also
								deleteIdx.add(i);
							}else{
								replaceIdx = i;
								convergeIdx = replaceIdx;
								convergeSize = gs.size();
								convergeMI = gs.getGeneIdx()[gs.size()-1].val();
							}
						}else{
							System.out.println("Lower MI");
							willBeDeleted = true;
							if(convergeIdx >= 0){ // already set to converge somewhere
								if(gs2.size() > convergeSize){
									System.out.println("New convergence.");
									deleteIdx.add(convergeIdx);
									convergeIdx = i;
									convergeSize = gs2.size();
									convergeMI = gs2.getGeneIdx()[gs.size()-1].val();
								}else if(gs2.size() == convergeSize){
									float gs2MI = gs2.getGeneIdx()[gs.size()-1].val();
									if(gs2MI > convergeMI){
										System.out.println("New convergence.");
										deleteIdx.add(convergeIdx);
										convergeIdx = i;
										convergeMI = gs2MI;
									}else{
										System.out.println("Original convergence.");
										deleteIdx.add(i);
									}
								}else{
									System.out.println("Original convergence.");
									deleteIdx.add(i);
								}
							}else{ // haven't set to converge anywhere
								System.out.println("Set original convergence");
								convergeIdx = i;
								convergeSize = gs2.size();
								convergeMI = gs2.getGeneIdx()[gs.size()-1].val();
							}
							
							if(replaceIdx >= 0){ // set to replace someone, send this someone to be deleted
								deleteIdx.add(replaceIdx);
								replaceIdx = -1;
							}
						}
					}else{ // cannot replace, going to be deleted
						System.out.println("Smaller size.");
						willBeDeleted = true;
						if(convergeIdx >= 0){ // already set to converge somewhere
							if(gs2.size() > convergeSize){
								System.out.println("New convergence.");
								deleteIdx.add(convergeIdx);
								convergeIdx = i;
								convergeSize = gs2.size();
								convergeMI = gs2.getGeneIdx()[gs.size()-1].val();
							}else if(gs2.size() == convergeSize){
								float gs2MI = gs2.getGeneIdx()[gs.size()-1].val();
								if(gs2MI > convergeMI){
									System.out.println("New convergence.");
									deleteIdx.add(convergeIdx);
									convergeIdx = i;
									convergeSize = gs2.size();
									convergeMI = gs2MI;
								}else{
									System.out.println("Original convergence.");
									deleteIdx.add(i);
								}
							}else{
								System.out.println("Original convergence.");
								deleteIdx.add(i);
							}
						}else{ // haven't set to converge anywhere
							System.out.println("Set original convergence");
							convergeIdx = i;
							convergeSize = gs2.size();
							convergeMI = gs2.getGeneIdx()[gs.size()-1].val();
						}
						if(replaceIdx >= 0){ // set to replace someone, send this someone to be deleted
							deleteIdx.add(replaceIdx);
							replaceIdx = -1;
						}
					}
					
				}
			}
			System.out.println("replaceIdx=" + replaceIdx + "\tconvergeIdx=" + convergeIdx + "\tdeleteIdx size=" + deleteIdx.size() + "\tWill be deleted=" + willBeDeleted);
			if(mergeable){ // if mergeable, there will be a converge Idx
				ValIdx[] vis;
				HashSet<Integer> ggg = clustGeneIdx.get(convergeIdx);
				String origMembers = clustMembers.get(convergeIdx);

				if(replaceIdx >= 0){
					leaders.set(replaceIdx, gs);
					// if there is a replace index, the gene idx is already in the first load
				}
				if(deleteIdx.size() > 0){
					Collections.sort(deleteIdx);
					Collections.reverse(deleteIdx);
					// delete the gene set reversely
					for(Integer ii : deleteIdx){
						// add all genes into the converge gene set
						vis = leaders.get(ii).getGeneIdx();
						origMembers = origMembers + "\t" + clustMembers.get(ii);
						for(ValIdx vi : vis){
							ggg.add(vi.idx());
						}
						// remove this gene set
						leaders.remove(ii);
						clustMembers.remove(ii);
						clustGeneIdx.remove(ii);
					}
				}
				// add the gene set loaded from file into the converge gene set
				origMembers = origMembers + "\t" + gs.getName();
				clustMembers.set(convergeIdx, origMembers);
				vis = gs.getGeneIdx();
				for(ValIdx vi : vis){
					ggg.add(vi.idx());
				}
			}
			// cannot find any overlap, create new list
			else{
				System.out.println("Cannot merge.");
				leaders.add(gs);
				clustMembers.add(gs.getName());
				HashSet<Integer> genes = new HashSet<Integer>();
				for(ValIdx vi : gs.getGeneIdx()){
					genes.add(vi.idx());
				}
				clustGeneIdx.add(genes);
			}
			if(progress % 100 == 0) System.out.println(progress + " / " + na);
		}
		
		
		System.out.println("Output to files...");
		
		int nc = leaders.size();
		boolean err = false;
		if(nc != clustGeneIdx.size()){
			System.out.println("Error: clustGeneIdx size is different from leaders!");
		}
		if(nc != clustMembers.size()){
			System.out.println("Error: clustMembers size is different from leaders!");
		}
		if(err){
			throw new RuntimeException("Error occurs. Exit.");
		}
		
		PrintWriter pw = new PrintWriter(new FileWriter(path + "ClusterLeaders.txt"));
		PrintWriter pw2 = new PrintWriter(new FileWriter(path + "ClusterGenes.txt"));
		PrintWriter pw3 = new PrintWriter(new FileWriter(path + "ClusterMembers.txt"));
		PrintWriter pw4 = new PrintWriter(new FileWriter(path + "forProf.txt"));
		
		int numOut = 5;
		
		for(int i = 0; i < nc; i++){
			GeneSet gs = leaders.get(i);
			pw.print(gs.getName() + "\t" + gs.size() + "\t");
			if(annotation){
				pw.println(gs.toGenes());
			}else{
				pw.println(gs.toProbes());
			}
			HashSet<Integer> genes = clustGeneIdx.get(i);
			String name = String.format("Cluster%-3s", i);
			pw2.print(name + "\t" + genes.size());
			for(Integer ii : genes){
				if(annotation){
					pw2.print("\t" + probeNames.get(ii) + ":" + annot.getGene(probeNames.get(ii)));
				}else{
					pw2.print("\t" + probeNames.get(ii));
				}
			}
			pw2.println();
			pw3.print(name);
			pw3.println("\t" + clustMembers.get(i));
			
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
		pw2.close();
		pw3.close();
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
