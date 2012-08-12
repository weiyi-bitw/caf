package caf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import obj.Annotations;
import obj.DataFile;
import obj.Genome;
import obj.ValIdx;
import obj.ValIdxD;

public class GroupWeightedAttractor {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		String path = "/home/weiyi/workspace/javaworks/caf/output/brca/weighted/";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		System.out.println("Loading files...");
		int IDX = 2;
		int outgenes = 500;
		int quantile = 50;
		
		String[] dataFiles = {
				//"/home/weiyi/workspace/data/m3d/sone/ge.4054x203.txt",
				//"/home/weiyi/workspace/data/m3d/ecoli/ge.4297x466.txt",
				//"/home/weiyi/workspace/data/m3d/yeast/ge.4515x407.txt",
				//"/home/weiyi/workspace/data/brca/gse3494/ge.12160x251.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/brca/gse32646/ge.19190x115.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/brca/gse36771/ge.19190x107.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/brca/gse31448/ge.19190x353.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/brca/gse2034/ge.12160x286.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/brca/tcga/ge/ge.17475x536.ncbi.txt",
				"/home/weiyi/workspace/data/dream7/preTraining/train/ge.37586x500.mean.txt",
				
				//"/home/weiyi/workspace/data/coad/gse14333/ge.19190x290.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/ov/gse9891/ge.19190x285.jetset.ncbi.txt",
				//"/home/weiyi/workspace/data/coad/tcga/ge/ge.17475x154.ncbi.txt",
				//"/home/weiyi/workspace/data/ov/tcga/ge/ge.11963x582.ncbi.txt"
				//"/home/weiyi/workspace/data/gbm/tcga/ge/ge.12042x545.txt",
				//"/home/weiyi/workspace/data/ov/tcga/super.35696x511.knn.txt",
				//"/home/weiyi/workspace/data/gbm/tcga/super.40092x274.txt",
				//"/home/weiyi/workspace/data/gbm/tcga/ge_mir_meth/meth.23094x278.knn.txt",
				//"/home/weiyi/workspace/data/ov/tcga/mergeroom/meth.23094x514.knn.txt",
				//"/home/weiyi/workspace/data/gbm/tcga/ge_mir_meth/mir.472x278.knn.txt",
				//"/home/weiyi/workspace/data/ov/tcga/mergeroom/mir.560x514.common.txt",
				
		};
		final String[] outputDirs={
				//"brca.gse3494",
				//"brca.gse32646",
				//"brca.gse36771",
				//"brca.gse31448",
				//"brca.gse2034",
				//"brca.tcga",
				"dream7",
				//"sone",
				//"ecoli",
				//"yeast",
				//"coad.gse14333",
				//"ov.gse9891",
				//"coad.tcga",
				//"ov.tcga",
				//"gbm.tcga"
				//"ov.super"
				//"gbm.super"
				//"meth.gbm",
				//"meth.ov",
				//"mir.gbm",
				//"mir.ov"
		};
		
		/*String[] annots = {
				"/home/weiyi/workspace/data/annot/affy/u133a/annot.jetset.csv",
				"/home/weiyi/workspace/data/annot/tcga/4502a073/annot.csv",
				"/home/weiyi/workspace/data/annot/affy/u133p2/annot.jetset.csv",
				"/home/weiyi/workspace/data/annot/tcga/4502a073/annot.csv",
				"/home/weiyi/workspace/data/annot/affy/u133p2/annot.jetset.csv",
				"/home/weiyi/workspace/data/annot/affy/u133a/annot.jetset.csv"
		};*/
		
		for(int qq = 0; qq < dataFiles.length; qq++)
		{
			
		System.out.println("Processing " + dataFiles[qq] + "...");
		//Annotations annot = Annotations.parseAnnotations(annots[IDX]);
		DataFile ma = DataFile.parse(dataFiles[qq]);
		ma.sortProbes();
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		ArrayList<String> genes = ma.getProbes();
		
		BufferedReader br = new BufferedReader(new FileReader(path + outputDirs[qq] + "/attractors.gwt"));
		BufferedReader br2 = new BufferedReader(new FileReader(path + outputDirs[qq] + "/attractees.gwt"));
		String line = br.readLine();
		String line2 = br2.readLine();
		/*String[] tokens = line.split("\t");
		int nt = tokens.length;
		String[] allgenes = new String[nt-2];
		int mg = allgenes.length;
		System.arraycopy(tokens, 2, allgenes, 0, nt-2);*/
		
		new File(path + "mergeroom").mkdirs();
		PrintWriter pw = new PrintWriter(new FileWriter(path + "mergeroom/" + outputDirs[qq]));
		PrintWriter pw2 = new PrintWriter(new FileWriter(path + outputDirs[qq] + "/attractees.decoded.gwt"));
		
		//line = br.readLine();
		
		while(line != null){
			ArrayList<ValIdxD> vec = new ArrayList<ValIdxD>();
			String[] tokens = line.split("\t");
			int numBasins = line2.split("\t").length - 2;
			String name = tokens[0];
			String[] t2 = line2.split("\t");
			pw.print(name + "\t" + numBasins);
			pw2.print(name + "\t" + numBasins);
			//double sum = 0;
			for(int i = 0; i < m; i++){
				//float f= Float.parseFloat(tokens[i+2]);
				//float w = (float) Math.exp(5 * Math.log(f));
				double w = Double.parseDouble(tokens[i+2]);
				vec.add(new ValIdxD(i, w));
					
				//sum += w;
			}
			for(int i = 0; i < numBasins; i++){
				pw2.print("\t" + genes.get(Integer.parseInt(t2[i+2])) );
			}
			Collections.sort(vec);
			for(int i = 0; i < outgenes; i++){
				String g = genes.get(vec.get(i).idx);
				pw.print("\t" + g + ":" + vec.get(i).val);
				//pw.print("\t" + g);
			}
			pw.print("\t" + vec.get(10-1).val);
			pw.println();
			pw2.println();
			
			line2 = br2.readLine();
			line = br.readLine();
		}
		br.close();
		pw.close();
		pw2.close();
		
		
		}
		System.out.println("Done.");
	}

}
