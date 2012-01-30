package caf;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import obj.Chromosome;
import obj.DataFile;
import worker.Converger;
import worker.Converger.ValIdx;






public class CNVFinder {
	static ArrayList<Chromosome> parseChromGenes(String file, ArrayList<String> genes, HashMap<String, Integer> geneMap) throws IOException{
		ArrayList<Chromosome> chrs = new ArrayList<Chromosome>();
		Chromosome.setGeneMap(geneMap);
		Chromosome.setGeneNames(genes);
		BufferedReader br = new BufferedReader(new FileReader(file));
		br.readLine(); // first line header
		
		String line = br.readLine();
		while(line != null){
			//System.out.println(line);
			String[] tokens = line.split("\t");
			int nt = tokens.length;
			Chromosome chr = new Chromosome(tokens[nt-1]);
			float coord = Float.parseFloat(tokens[5]);
			boolean strand = tokens[2].equals("(+)");
			if(chrs.contains(chr)){
				chrs.get(chrs.indexOf(chr)).addGene(tokens[0], coord, strand);
			}else{
				chr.addGene(tokens[0], coord, strand);
				chrs.add(chr);
			}
			line = br.readLine();
		}
		return chrs;
	}
	public static void main(String args[]) throws Exception{
		String path = "/home/weiyi/workspace/data/ov/tcga/ge/";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		System.out.println("Loading files...");
		
		DataFile ma = DataFile.parse(path + "ge.12042x582.txt");
		//ma.normalizeRows();
		ArrayList<Chromosome> chrs = parseChromGenes("/home/weiyi/workspace/data/annot/affy/u133p2/gene.location3", 
				ma.getProbes(), ma.getRows());
		ArrayList<Chromosome> chr2 = new ArrayList<Chromosome>();
		chr2.add(chrs.get(10));
		
		Converger cvg = new Converger(0, 1, System.currentTimeMillis(), "ZSCORE", 100, false);
		System.out.println("Finding CNVs...");
		
		cvg.findCNV(ma.getData(),ma.getData(), chr2, 4.0f);
		
		System.out.println("Done.");
	}
}
