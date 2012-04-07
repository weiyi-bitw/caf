package caf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import obj.Annotations;
import obj.DataFile;
import obj.Genome;
import obj.InverseAnnotations;
import obj.ValIdx;
import worker.Converger;
import worker.ITComputer;

public class TestFieldProbe {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String annotPath = "/home/weiyi/workspace/data/annot/affy/u133a/annot.csv";
		//String annotPath = "/home/weiyi/workspace/data/annot/tcga/4502a073/annot.64847.csv";
		Annotations annot = Annotations.parseAnnotations(annotPath);
		InverseAnnotations invannot = annot.getInvAnnot();
		String[] allgenes = annot.getAllGenes();
		
		String path = "/home/weiyi/workspace/data/brca/gse2034/";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		System.out.println("Loading files...");
		
		//String outPath = "/home/weiyi/workspace/javaworks/caf/output/656/";
		String outPath = "/home/weiyi/workspace/javaworks/caf/tmp/";
		if(!outPath.endsWith("/")){
			outPath = outPath + "/";
		}
		DataFile ma = DataFile.parse(path + "ge.22283x286.txt");
		
		final String geneLocFile = "/home/weiyi/workspace/data/annot/affy/u133p2/gene.location3";
		//final String geneLocFile = "/home/weiyi/workspace/javaworks/caf/output/639/gene.location3";
		
		String command = "CNV";
		float power = 2f;
		int winSize = 51;
		
		//ma.normalizeRows();
		int m = ma.getNumRows();
		int mg = allgenes.length;
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
		
		Converger cvg = new Converger(0, 1, jobID);
		
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		if(command.equals("CNV")) gn.linkToGeneSet(allgenes);
		
		gs.add("HSF1");
		
		ITComputer itc = new ITComputer(6, 3, 0, 1, true);
		//itc.negateMI(true);
		cvg.linkITComputer(itc);
		HashMap<String, Integer> geneMap = ma.getRows();
		ArrayList<String> geneNames = ma.getProbes();
		new File("tmp").mkdir();
		for(String gtmp : gs){
			String[] pbs = invannot.getProbes(gtmp);
			
			if(command.equals("CNV")){
				String[] neighborsG = gn.getNeighbors(gtmp, winSize);
				String[] neighbors = invannot.getProbes(neighborsG);
				if(neighbors == null){
					System.out.println("No neighbors :(");
					continue;
				}
				DataFile ma2 = ma.getSubProbes(neighbors);
				geneNames = ma2.getProbes();
				m = ma2.getNumRows();
				
				for(String p : pbs){
					String genename = annot.getGene(p);
					System.out.println("Processing " + genename + "(" + p + ") (" + (gs.indexOf(p)+1) + "/" + gs.size() + ")" + "...");
					String outFile = outPath + genename + "_" + p + "_CNV.txt";
					try{
						PrintWriter pw = new PrintWriter(new FileWriter(outFile));
						int idx = ma2.getRows().get(p);
						data = ma2.getData();
						float[] vec = data[idx];
						ValIdx[] out = cvg.findWeightedCNV(ma2, vec, gn, neighborsG, invannot, winSize, power, annot);
						if(out[0] == null){
							continue;
						}
						
						Arrays.sort(out);
						for(int i = 0; i < out.length; i++){
							if(out[i].idx > 0){
								String gg = geneNames.get(out[i].idx);
								pw.println(annot.getGene(gg) + "\t"  + out[i].val + "\t" + gg);
							}
						}
						pw.close();
						
					}catch (FileNotFoundException e){
						System.out.println("Exception: " + e);
						continue;
					}
					
					
				} //  END CNV probe iteration
			}else{
				for(String g : pbs){
				
				if(geneNames.contains(g)){
					String genename = annot.getGene(g);
					System.out.println("Processing " + genename + "(" + g + ") (" + (gs.indexOf(g)+1) + "/" + gs.size() + ")" + "...");
					String outFile = outPath + genename + "_"+ g + "_CA.txt";
					try{
						PrintWriter pw = new PrintWriter(new FileWriter(outFile));
						ValIdx[] out = new ValIdx[mg];
						float[] vec = new float[n];
						int idx = geneMap.get(g);
						vec = data[idx];
						out = cvg.findWeightedAttractor(ma, vec, allgenes, invannot, power, annot);
						
						if(out[0] == null){
							continue;
						}
						
						Arrays.sort(out);
						for(int i = 0; i < mg; i++){
							if(out[i].idx > 0){
								String gg = geneNames.get(out[i].idx);
								pw.println(annot.getGene(gg) + "\t"  + out[i].val  + "\t" + gg);
							}
						}
						pw.close();
						
					}catch (FileNotFoundException e){
						System.out.println("Exception: " + e);
						continue;
					}
				}else{
					System.out.println("Does not contain gene " + g + "!!");
				}
			
			} // END CAF Probe iteration
			
			}// END if command == CNV
		}// END Target iteration
		
	}

}
