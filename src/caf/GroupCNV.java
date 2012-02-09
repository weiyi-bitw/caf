package caf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import obj.Genome;
import obj.Attractor;

public class GroupCNV {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{
		final String geneLocFile = "/home/weiyi/workspace/data/annot/affy/u133p2/gene.location3";
		String path = "/home/weiyi/workspace/javaworks/caf/output/cnv/coad.gse14333.rownorm.sz10.cnv";
		
		System.out.println("Loading gene location file...");
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		
		
		if(!path.endsWith("/")){
			path = path + "/";
		}
		
		Attractor.CNV(true);
		Attractor.setGenome(gn);
		
		ArrayList<Attractor> attractors = new ArrayList<Attractor>();
		
		BufferedReader br = new BufferedReader(new FileReader(path + "attractors.gwt"));
		String line = br.readLine();
		while(line != null){
			Attractor a = Attractor.parseAttractor(line);
			attractors.add(a);
			line = br.readLine();
		}
		
		int N = attractors.size();
		System.out.println( N + " attractors loaded.");
		br.close();
		Collections.sort(attractors);
		
		System.out.println("Filtering...");
		for(int i = N-1; i >= 1; i--){
			boolean kill = false;
			Attractor a = attractors.get(i);
			for(int j = 0; j < i; j++){
				if(a.getOvlp(attractors.get(j)).size() > 0){
					kill = true;
					break;
				}
			}
			if(kill){
				attractors.remove(a);
			}
		}
		
		N = attractors.size();
		System.out.println( N + " attractors left.");
		
		PrintWriter pw = new PrintWriter(new FileWriter(path + "attractors.filtered.gwt"));
		for(Attractor a : attractors){
			pw.println(a);
		}
		pw.close();
		
		pw = new PrintWriter(new FileWriter(path + "attractors.filtered.gwt.detail"));
		for(Attractor a : attractors){
			pw.println(a.toStringInDetail());
		}
		pw.close();
	}

}
