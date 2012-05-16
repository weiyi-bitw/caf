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

public class GroupWeightedAttractor {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		String path = "/home/weiyi/workspace/javaworks/caf/output/weighted/";
		if(!path.endsWith("/")){
			path = path + "/";
		}
		System.out.println("Loading files...");
		int IDX = 2;
		int outgenes = 50;
		int quantile = 50;
		
		String[] dataFiles = {
				//"/home/weiyi/workspace/data/m3d/sone/ge.4054x203.txt",
				//"/home/weiyi/workspace/data/m3d/ecoli/ge.4297x466.txt",
				//"/home/weiyi/workspace/data/m3d/yeast/ge.4515x407.txt",
				/*"/home/weiyi/workspace/data/brca/gse2034/ge.12160x286.jetset.mean.txt",
				"/home/weiyi/workspace/data/coad/gse14333/ge.19189x290.jetset.mean.txt",
				"/home/weiyi/workspace/data/ov/gse9891/ge.19189x285.jetset.mean.txt",
				"/home/weiyi/workspace/data/brca/tcga/ge/ge.17814x536.knn.txt",
				"/home/weiyi/workspace/data/coad/tcga/ge/ge.17814x154.knn.txt",
				"/home/weiyi/workspace/data/ov/tcga/ge/ge.12042x582.txt"*/
				//"/home/weiyi/workspace/data/gbm/tcga/ge/ge.12042x545.txt",
				//"/home/weiyi/workspace/data/ov/tcga/super.35696x511.knn.txt",
				//"/home/weiyi/workspace/data/gbm/tcga/super.40092x274.txt",
				//"/home/weiyi/workspace/data/gbm/tcga/ge_mir_meth/meth.23094x278.knn.txt",
				//"/home/weiyi/workspace/data/ov/tcga/mergeroom/meth.23094x514.knn.txt",
				"/home/weiyi/workspace/data/gbm/tcga/ge_mir_meth/mir.472x278.knn.txt",
				"/home/weiyi/workspace/data/ov/tcga/mergeroom/mir.560x514.common.txt",
				
		};
		final String[] outputDirs={
				//"sone",
				//"ecoli",
				//"yeast",
				/*"brca.gse2034.jetset.mean",
				"coad.gse14333.jetset.mean",
				"ov.gse9891.jetset.mean",
				"brca.tcga",
				"coad.tcga",
				"ov.tcga.affy"*/
				//"gbm.tcga"
				//"ov.super"
				//"gbm.super"
				//"meth.gbm",
				//"meth.ov",
				"mir.gbm",
				"mir.ov"
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
		
		new File(path + "mergeroom.300").mkdir();
		PrintWriter pw = new PrintWriter(new FileWriter(path + "mergeroom.300/" + outputDirs[qq]));
		PrintWriter pw2 = new PrintWriter(new FileWriter(path + outputDirs[qq] + "/attractees.decoded.gwt"));
		
		//line = br.readLine();
		
		while(line != null){
			ArrayList<ValIdx> vec = new ArrayList<ValIdx>();
			String[] tokens = line.split("\t");
			int numBasins = line2.split("\t").length - 2;
			String name = tokens[0];
			String[] t2 = line2.split("\t");
			pw.print(name + "\t" + numBasins);
			pw2.print(name + "\t" + numBasins);
			float sum = 0;
			for(int i = 0; i < m; i++){
				//float f= Float.parseFloat(tokens[i+2]);
				//float w = (float) Math.exp(5 * Math.log(f));
				float w = Float.parseFloat(tokens[i+2]);
				vec.add(new ValIdx(i, w));
					
				sum += w;
			}
			for(int i = 0; i < numBasins; i++){
				pw2.print("\t" + genes.get(Integer.parseInt(t2[i+2])) );
			}
			Collections.sort(vec);
			for(int i = 0; i < outgenes; i++){
				String g = genes.get(vec.get(i).idx);
				//pw.print("\t" + g + ":" + vec.get(i).val);
				pw.print("\t" + g);
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
