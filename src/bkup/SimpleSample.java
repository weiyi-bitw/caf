package bkup;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;

import obj.DataFile;
import util.StatOps;
import worker.ITComputer;

public class SimpleSample {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String file = args[0];
		
		System.out.println("Loading files...");
		DataFile ma = DataFile.parse(file);
		//ma.normalizeRows();
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		
		int nSamples = m/100;
		HashSet<Integer> sample = new HashSet<Integer>();
		Random r = new Random(System.currentTimeMillis());
		while(sample.size() < nSamples){
			int ii = r.nextInt(m);
			sample.add(ii);
		}
		
		ITComputer itc = new ITComputer(6, 3, 0, 1);
		float zth[] = {5f, 6f, 7f, 8f, 9f, 10f};
		int nz = zth.length;
		int sizeMat[][] = new int[nSamples][nz];
		
		System.out.println("Sample Size:" + nSamples);
		
		int cnt = 0;
		float[][][] w = itc.getWeights(ma);
		for(Integer i : sample){
			//float[] fixVec = data[i];
			float[] mi = itc.getAllMIWith(i, w);
			float[] z = StatOps.xToZ(mi, m);
			//Arrays.sort(z);
			for(int j = 0; j < m; j++){
				for(int k = 0; k < nz; k++){
					if(z[j] >= zth[k]){
						sizeMat[cnt][k]++;
					}
				}
			}
			if(cnt % 10 == 0){
				System.out.print(100*cnt/nSamples + "%\r");
			}
			cnt++;
		}
		
		
		float meanSize[] = new float[nz];
		float size75[] = new float[nz];
		float size25[] = new float[nz];
		float medianSize[] = new float[nz];
		for(int j = 0; j < nz; j++){
			float x[] = new float[nSamples];
			for(int i = 0; i < nSamples; i++){
				meanSize[j] += sizeMat[i][j];
				x[i] = sizeMat[i][j];
			}
			size75[j] = StatOps.quantile(x, 0.75f);
			size25[j] = StatOps.quantile(x, 0.25f);
			medianSize[j] = StatOps.median(x);
			meanSize[j] /= nSamples;
		}
		System.out.println("Z-Score\t.25 Size\tMedian Gene Set Size\tAverage Gene Set Size\t.75 Size");
		for(int i = 0; i < nz; i++){
			System.out.println(zth[i] + "\t" + size25[i] + "\t" + medianSize[i] + "\t" + meanSize[i] + "\t" + size75[i]);
		}
		
		
	}

}
