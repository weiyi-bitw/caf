package obj;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Random;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintStream;

import obj.ValIdx;

/**
 * @author John Watkinson
 * @revision Wei-Yi Cheng
 * @version 0.1
 * @date 06/07/2011
 * 
 * Limited version. Without summarize probe into gene function, and no output to workspace function.
 */
public class DataFile {

    static class Size {
        int n, m;
    }

    public static String stripQuotes(String s) {
        s = s.replace("'", "");
    	return s.replace("\"", "");
    }

    public static DataFile parse(String file) throws Exception {
        final Size size = new Size();
        processTokens(file, new TokenProcessor() {
            public void process(String[] tokens) {
                if (size.n == 0) {
                    size.n = tokens.length - 1;
                }
                size.m++;
            }
        }, true, "\t");
        int n = size.n;
        int m = size.m;
        if (n < 1) {
            throw new Exception("Must have at least 1 column of data.");
        }
        if (m < 1) {
            throw new Exception("Must have at least 1 row of data.");
        }
        final float[][] data = new float[m][n];
        final HashMap<String,Integer> rows = new HashMap<String, Integer>();
        final HashMap<String,Integer> cols = new HashMap<String, Integer>();
        final ArrayList<String> probes = new ArrayList<String>();
        final ArrayList<String> chipID = new ArrayList<String>();
        final ArrayList<Integer> nullIndices = new ArrayList<Integer>();
        processTokens(file, new TokenProcessor() {
            boolean start = true;
            int row = 0;
            int cnt = 0;
            public void process(String[] tokens) {
                if (start) {
                    for (int i = 1; i < tokens.length; i++) {
                    	if(cols.get(stripQuotes(tokens[i])) != null){
                    		throw new RuntimeException("Warning! Chip name " + tokens[i] + " duplicates!!");
                    	}
                    	cols.put(stripQuotes(tokens[i]), (i - 1));
                        chipID.add(stripQuotes(tokens[i]));
                    }
                    start = false;
                } else {
                	if(rows.get(stripQuotes(tokens[0]))!=null){
                		throw new RuntimeException("Warning! Probe name " + tokens[0] + " is duplicates!!");
                	}
                    rows.put(stripQuotes(tokens[0]), row);
                    probes.add(stripQuotes(tokens[0]));
                    for (int i = 1; i < tokens.length; i++) {
                        try {
                            if(stripQuotes(tokens[i]).equals("null") || stripQuotes(tokens[i]).equals("NA")) {
                            	data[row][i - 1] = Float.NaN;
                            	nullIndices.add(cnt);
                            }
                            else	data[row][i - 1] = Float.parseFloat(stripQuotes(tokens[i]));
                        } catch (NumberFormatException nfe) {
                            throw new RuntimeException("Couldn't parse number at row " + (row + 1) + ", column " + i + ":" + stripQuotes(tokens[i]));
                        }
                        cnt++;
                    }
                    row++;
                }
            }
        }, false, "\t");
        return new DataFile(data, rows, cols, probes, chipID, nullIndices);
    }

    public void parsePhenotype(final String file) throws Exception {
        final float[] phen = new float[n];
        processTokens(file, new TokenProcessor() {
            public void process(String[] tokens) {
                String id = stripQuotes(tokens[0]);
                Integer index = cols.get(id);
                if (index == null) {
                    throw new RuntimeException("Could not find phenotype '" + id + "' in phenotype file '" + file + "'.");
                }
                try {
                    phen[index] = Float.parseFloat(stripQuotes(tokens[1]));
                } catch (NumberFormatException nfe) {
                    throw new RuntimeException("Couldn't parse number for phenotype id '" + id + "': '" + stripQuotes(tokens[1]) + "'.");
                }
            }
        }, true, "\t");
        float[][] newData = new float[m + 1][];
        for (int i = 0; i < m; i++) {
            newData[i] = data[i];
        }
        newData[m] = phen;
        data = newData;
        probes.add("Phenotype");
    }
    
    
 // 06/22/10 adding output2File function
    
    public void output2File(final String file) throws Exception{
    	PrintStream ps = new PrintStream(new FileOutputStream(file));
    	ps.print("ExpID");
    	for(int j = 0; j < n; j++){
    		ps.print("\t" + chipID.get(j));
    	}ps.println();
    	for(int i = 0; i < m; i++){
    		ps.print(probes.get(i));
    		for(int j = 0; j < n ; j++){
    			ps.print("\t" + data[i][j]);
    		}
    		ps.println();
    	}
    	ps.close();
    }
// 06/07/11 output2Gct function
    public void output2Gct(final String file) throws Exception{
    	PrintStream ps = new PrintStream(new FileOutputStream(file));
    	ps.println("#1.2");
    	ps.println(m + "\t" + n);
    	ps.print("Feature\tDescription");
    	for(int j = 0; j < n; j++){
    		ps.print("\t" + chipID.get(j));
    	}ps.println();
    	for(int i = 0; i < m; i++){
    		ps.print(probes.get(i) + "\t" + "na");
    		for(int j = 0; j < n ; j++){
    			ps.print("\t" + data[i][j]);
    		}
    		ps.println();
    	}
    	ps.close();
    }
    
    public void output2Gct(final String file, int[] order) throws Exception{
    	PrintStream ps = new PrintStream(new FileOutputStream(file));
    	ps.println("#1.2");
    	ps.println(m + "\t" + n);
    	ps.print("Feature\tDescription");
    	int y;
    	for(int j = 0; j < n; j++){
    		y = order[j];
    		ps.print("\t" + chipID.get(y));
    	}ps.println();
    	for(int i = 0; i < m; i++){
    		ps.print(probes.get(i) + "\t" + "na");
    		for(int j = 0; j < n ; j++){
    			y = order[j];
    			ps.print("\t" + data[i][y]);
    		}
    		ps.println();
    	}
    	ps.close();
    }
    
    public void output2MatlabWorkspace(String path) throws Exception{
    	if(path.endsWith("/")){
    		path = path + "/";
    	}
    	PrintStream ps = new PrintStream(new FileOutputStream(path + "E.txt"));
    	for(int i = 0; i < m; i++){
    		ps.print(data[i][0]);
    		for(int j = 1; j < n; i++){
    			ps.print("\t" + data[i][j]);
    		}ps.println();
    	}
    	ps.close();
    	
    	ps = new PrintStream(new FileOutputStream(path + "probe.txt"));
    	for(String s : probes){
    		ps.println(s);
    	}
    	ps.close();
    	
    	ps = new PrintStream(new FileOutputStream(path + "samples.txt"));
    	for(String s : chipID){
    		ps.println(s);
    	}
    	ps.close();
    	
    }
// 06/23/10 adding transpose function
    
    public void transpose(final String file) throws Exception{
    	float[][] tData = new float[n][m];
    	for(int i = 0; i < n; i++){
    		for(int j = 0; j < m; j++){
    			tData[i][j] = data[j][i];
    		}
    	}
    	HashMap<String, Integer> temp = new HashMap<String, Integer>();
    	temp = rows;
    	rows = cols;
    	cols = temp;
    	
    	ArrayList<String> temp2 = new ArrayList<String>();
    	temp2 = probes;
    	probes = chipID;
    	chipID = temp2;
    	
    	data = tData;
    	int tempInt = n;
    	n = m;
    	m = tempInt;
    	
    }

 // 11/10/10 adding combine function
    
    public void combine(final DataFile dataFile2) throws Exception{
    	HashMap<String, Integer> newColMap = new HashMap<String, Integer>();
    	ArrayList<String> newColNames = new ArrayList<String>();
    	
    	System.out.print("Find common sample IDs...");
    	
    	String[] mySampleNames = chipID.toArray(new String[0]);
    	String[] yourSampleNames = dataFile2.getChipID().toArray(new String[0]);
    	int n2 = dataFile2.getNumCols();
    	int newSampleCount = 0;
    	ArrayList<Integer> newOrder1 = new ArrayList<Integer>();
    	ArrayList<Integer> newOrder2 = new ArrayList<Integer>();
    	for(int i = 0; i < n; i++){
    		//String name1 = mySampleNames[i].substring(0, 16);
    		String name1 = mySampleNames[i];
    		for(int j = 0; j < n2; j++){
    			//String name2 = yourSampleNames[j].substring(0, 16);
    			String name2 = yourSampleNames[j];
    			if(name1.equals(name2)){
    				newOrder1.add(i);
    				newOrder2.add(j);
    				newColNames.add(name1);
    				newColMap.put(name1, newSampleCount);
    				newSampleCount++;
    				break;
    			}
    		}
    	}
    	if(newSampleCount ==0) throw new RuntimeException("No common sample found!!");
    	else System.out.println(newSampleCount + " samples are in common.");
    	
    	int newM = m + dataFile2.getNumRows();
    	String[] probeName2 = dataFile2.getProbes().toArray(new String[0]);
    		float[][] newData = new float[newM][newSampleCount];
    		float[][] data2 = dataFile2.getData();
    		for(int i = 0; i < newM; i++){
    			for(int j = 0; j < newSampleCount; j++){
    				newData[i][j] = i < m ? data[i][newOrder1.get(j)] : data2[i-m][newOrder2.get(j)]; 
    			}
    			if(i >= m){
    				rows.put(probeName2[i-m], i);
    			}
    		}
    		probes.addAll(dataFile2.getProbes());
    		//rows.putAll(dataFile2.getRows());
    		data = newData;
    		m = newM;
    		n = newSampleCount;
    		chipID = newColNames;
    		cols = newColMap;
    		
    		System.out.println("New size: " + m + "x" + n);
    }
// 07/16/10 adding mergeSamples function
    public void mergeSamples(final String[] sampleIDs, final String methods){    	
    	int numSamples = sampleIDs.length;
    	int[] sampleIndexes = new int[numSamples];
    	ArrayList<Integer> indAL = new ArrayList<Integer>();
    	for(int i = 0; i < m; i++){
    		indAL.add(i);
    	}
    	for(int i = 0; i < numSamples; i++){
    		sampleIndexes[i] = cols.get(sampleIDs[i]);
    		chipID.remove(sampleIDs[i]);
    		cols.remove(sampleIDs[i]);
    		indAL.remove((Integer) sampleIndexes[i]);
    	}
    	Integer[] ind = indAL.toArray(new Integer[0]);
    	chipID.add(sampleIDs[0] + "_et" + numSamples + "_" + methods);
    	int newN = chipID.size();
    	HashMap<String, Integer> newCols = new HashMap<String, Integer>();
    	for(int i = 0; i < newN; i++){
    		newCols.put(chipID.get(i), i);
    	}
    	cols = newCols;
    	float[][] newData = new float[m][newN];
    	for(int i = 0; i < m; i++){
    		for(int j = 0; j < newN - 1; j++){
    			newData[i][j] = data[i][ind[j]];
    		}
    	}
    	for(int i = 0; i < m; i++){
    		float[] subset = new float[numSamples]; 
    		for(int j = 0; j < numSamples; j++){
    			subset[j] = data[i][sampleIndexes[j]];
    		}
    		if(methods.equals("mean")){
    			newData[i][newN-1] = (float)meanf(subset, numSamples);
    		}else if(methods.equals("median")){
    			newData[i][newN-1] = (float)medianf(subset, numSamples);
    		}else{
    			throw new RuntimeException("Cannot recognize the command " + methods);
    		}
    	}
    	n=newN;
    	data = newData;
    }
    // 08/14/2010 add mergeProbes
    public void mergeProbes(String[] probeIDs, String method, String newName){
    	int numProbes = probeIDs.length;
    	int[] rowIndexes = new int[numProbes];
    	ArrayList<Integer> indAL = new ArrayList<Integer>();
    	for(int i = 0; i < m; i++){
    		indAL.add(i);
    	}
    	for(int i = 0; i < numProbes; i++){
    		rowIndexes[i] = rows.get(probeIDs[i]);
    		probes.remove(probeIDs[i]);
    		rows.remove(probeIDs[i]);
    		indAL.remove((Integer) rowIndexes[i]);
    	}
    	Integer[] ind = indAL.toArray(new Integer[0]);
    	probes.add(newName);
    	int newM = probes.size();
    	HashMap<String, Integer> newRows = new HashMap<String, Integer>();
    	for(int i = 0; i < newM; i++){
    		newRows.put(probes.get(i), i);
    	}
    	rows = newRows;
    	float[][] newData = new float[newM][n];
    	for(int i = 0; i < newM-1; i++){
    		for(int j = 0; j < n; j++){
    			newData[i][j] = data[ind[i]][j];
    		}
    	}
    	for(int j = 0; j < n; j++){
    		float[] subset = new float[numProbes]; 
    		for(int i = 0; i < numProbes; i++){
    			subset[i] = data[rowIndexes[i]][j];
    		}
    		if(method.equals("mean")){
    			newData[newM - 1][j] = (float)meanf(subset, numProbes);
    		}else if(method.equals("median")){
    			newData[newM - 1][j] = (float)medianf(subset, numProbes);
    		}else if(method.equals("max")){
    			newData[newM - 1][j] = max(subset, numProbes);
    		}else{
    			throw new RuntimeException("Cannot recognize the command " + method);
    		}
    	}
    	m=newM;
    	data = newData;
    }
 // 09/23/2010 filter probes function: filter out probes according to their variance
    public void filterProbes(int numFilterIn){
    	float[] xMeans = new float[m];
		float[] xSds = new float[m];
		for(int i = 0; i < m; i++){
			for(int j = 0; j < n; j++){
				xMeans[i] += data[i][j];
				xSds[i] += data[i][j] * data[i][j];
			}
		}
		for(int i = 0; i < m; i++){
			xMeans[i] /= n;
			xSds[i] = (float) Math.sqrt((xSds[i] - n * xMeans[i] * xMeans[i])/(n - 1));
			if(Float.isNaN(xSds[i])){
				xSds[i] = 1;
			}
		}
		float[] xSdsSort = new float[m];
		System.arraycopy(xSds, 0, xSdsSort, 0, m);
		Arrays.sort(xSdsSort);
		float threshold = xSdsSort[m - numFilterIn];
		int[] newRowsIdx = new int[numFilterIn];
		int count = 0;
		for(int i = 0; i < m && count < numFilterIn; i++){
			if(xSds[i] >= threshold){
				newRowsIdx[count] = i;
				count++;
			}
		}
		ArrayList<String> newProbes = new ArrayList<String>();
		HashMap<String, Integer> newRows = new HashMap<String, Integer>();
    	float[][] newData = new float[numFilterIn][n];
    	for(int i = 0; i < numFilterIn; i++){
    		String probeName = probes.get(newRowsIdx[i]);
    		System.arraycopy(data[newRowsIdx[i]], 0, newData[i], 0, n);
    		newProbes.add(probeName);
    		newRows.put(probeName, i);
    	}
		probes = newProbes;
    	rows = newRows;
		m = numFilterIn;
		data = newData;
    }
    // 04/20/2011 Sort probes lexicographically
    public void sortProbes(){
    	Collections.sort(probes);
    	float[][] newData = new float[m][n];
    	int j;
    	for(int i = 0; i < m; i++){
    		String gene = probes.get(i);
    		j = rows.get(gene);
    		System.arraycopy(data[j], 0, newData[i], 0, n);
    		rows.put(gene, i);
    	}
    	data = newData;
    }
    
    // 08/18/2010 normalizeRows function
    public void normalizeRows(){
    	float[] means = new float[m];
    	float[] stds = new float[m];
    	
    	
    	for(int i = 0; i < m; i++){
    		float square = 0;
    		for(int j = 0; j < n; j++){
    			means[i] += data[i][j];
    			square = data[i][j] * data[i][j];
    			stds[i] += square;
    		}
    		means[i] /= n;
    		stds[i] = (float)Math.sqrt((stds[i] - n * means[i] * means[i])/(n - 1));
    		if(Float.isNaN(stds[i])){
    			stds[i] = 1;
    		}
    	}
    	
    	for(int i = 0; i < m; i++){
    		if(stds[i] == 0){
				System.out.println("Warning: row " + i + " has 0 standard deviation!");
			}else{
	    		for(int j = 0; j < n; j++){
	    			data[i][j] = (data[i][j] - means[i])/stds[i];
	    		}
			}
    	}
    	
    }
    
    public void normalizeCols(){
    	float[] means = new float[n];
    	float[] stds = new float[n];
    	
    	
    	for(int i = 0; i < n; i++){
    		float square = 0;
    		for(int j = 0; j < m; j++){
    			means[i] += data[j][i];
    			square = data[j][i] * data[j][i];
    			stds[i] += square;
    		}
    		means[i] /= m;
    		stds[i] = (float)Math.sqrt((stds[i] - m * means[i] * means[i])/(m - 1));
    		if(Float.isNaN(stds[i])){
    			stds[i] = 1;
    		}
    	}
    	
    	for(int i = 0; i < n; i++){
    		if(stds[i] == 0){
				System.out.println("Warning: col " + i + " has 0 standard deviation!");
			}else{
	    		for(int j = 0; j < m; j++){
	    			data[j][i] = (data[j][i] - means[i])/stds[i];
	    		}
			}
    	}
    	
    }
	
    //02/08/2011 log2normalization
    public void log2normalization(){
    	for(int i = 0; i < m; i++){
    		for(int j = 0; j < n; j++){
    			data[i][j] = (float) (Math.log(data[i][j]) / Math.log(2));
    		}
    	}
    }
    
    public DataFile(float[][] data, HashMap<String, Integer> rows, HashMap<String, Integer> cols, ArrayList<String> probes, ArrayList<String> chipID) {
        this.data = data;
        this.rows = rows;
        this.cols = cols;
        this.probes = probes;
        this.chipID = chipID;
        m = rows.size();
        n = cols.size();
    }

    public DataFile(float[][] data, HashMap<String, Integer> rows, HashMap<String, Integer> cols, ArrayList<String> probes, ArrayList<String> chipID, ArrayList<Integer> nullIndices) {
        this.data = data;
        this.rows = rows;
        this.cols = cols;
        this.probes = probes;
        this.chipID = chipID;
        this.nullIndices = nullIndices;
        m = rows.size();
        n = cols.size();
    }
    
    // 05/05/11 add new state: nullIndices
    private float[][] data;
    private HashMap<String, Integer> rows;
    private HashMap<String, Integer> cols;
    private ArrayList<String> probes;
    private ArrayList<String> chipID;
    private int n; // number of columns
    private int m; // number of rows : m-by-n !!!
    private ArrayList<Integer> nullIndices;
    
    private void swap(int a, int b) {
        float temp;
        for (int i = 0; i < m; i++) {
            temp = data[i][a];
            data[i][a] = data[i][b];
            data[i][b] = temp;
        }
    }
    
    private void swap(int i, int a, int b) {
        float temp;
            temp = data[i][a];
            data[i][a] = data[i][b];
            data[i][b] = temp;
    }
    private static double meanf(float[] data, int numSamples) {
        int curSample;
        double mean = 0;

        for (curSample = 0; curSample < numSamples; curSample++) {
            mean += data[curSample];
        }
        return mean / (double) numSamples;
    }
    private static double medianf(float[] data, int numElem) {
        float[] dataCopy = new float[numElem];
        System.arraycopy(data, 0, dataCopy, 0, numElem);
        int half = numElem / 2;
        double median;
        
        Arrays.sort(dataCopy);
        if(numElem % 2 ==0){
        	median = ((double) dataCopy[half-1] + (double) dataCopy[half]) / 2.0;
        }else{
        	median = (double) dataCopy[half];
        }
        return median;
    }
    private static float max(float[] data, int numElem){
    	float maximum = data[0];
    	for(int i = 1; i < numElem; i++){
    		if (data[i] > maximum) {maximum = data[i];} 
    	}
    	return maximum;
    }
    public void permute(long seed, int index) { // permute a row
        //System.out.print(seed + "\t");
	Random random = new Random(seed);
        for (int i = n-1; i > 0; i--) {
            int j = random.nextInt(i);
            //System.out.print(j + "\t");
            swap(index, i, j);
        }//System.out.println();
    }

    public void permute(long seed) { // permute the whole columns (entries within the columns remain untouched)
        Random random = new Random(seed);
        for (int i = n-1; i > 0; i--) {
            int j = random.nextInt(i);
            swap(i, j);
        }
    }

    public int getNumRows() {
        return m;
    }

    public int getNumCols() {
        return n;
    }

    public float[][] getData() {
        return data;
    }

    public DataFile getSubChips(String[] chipNames){
    	int numChips = chipNames.length;
    	int[] chipInd = new int[numChips];
    	HashMap<String, Integer> subCols = new HashMap<String, Integer>();
    	ArrayList<String> subChipNames = new ArrayList<String>();
    	int cnt = 0;
    	for(int i = 0; i < numChips; i++){
    		if(cols.get(chipNames[i]) != null){
    			chipInd[cnt] = cols.get(chipNames[i]);
    			subChipNames.add(chipNames[i]);
    			subCols.put(chipNames[i], cnt);
    			cnt++;
    		}else{
    			//throw new RuntimeException("At " + i + ": no chip name " + chipNames[i] + "!!");
    			System.out.println("At " + i + ": no chip name " + chipNames[i] + "!!");
    		}
    	}
    	float[][] subData = new float[m][cnt];
    	for(int i = 0; i < m; i++){
    		for(int j = 0; j < cnt; j++){
    			subData[i][j] = data[i][chipInd[j]];
    		}
    	}
    	return new DataFile(subData, rows, subCols, probes, subChipNames);
    }
    public DataFile getSubProbes(ArrayList<ValIdx> probeIdx){
    	int numProbes = probeIdx.size();
    	HashMap<String, Integer> subRows = new HashMap<String, Integer>();
    	ArrayList<String> subProbes = new ArrayList<String>();
    	int mm = numProbes;
    	
    	float[][] subData = new float[mm][n];
    	for(int i = 0; i < mm; i++){
    		int idx = probeIdx.get(i).idx();
    		String p = probes.get(idx);
    		subProbes.add(p);
    		subRows.put(p, i);
    		System.arraycopy(data[idx], 0, subData[i], 0, n);
    	}
    	return new DataFile(subData, subRows, cols, subProbes, chipID);
    }
    
    public DataFile getSubProbes(String[] probeNames){
    	int numProbes = probeNames.length;
    	HashMap<String, Integer> subRows = new HashMap<String, Integer>();
    	ArrayList<String> subProbes = new ArrayList<String>();
    	int mm = numProbes;
    	for(int i = 0; i < numProbes; i++){
    		if(rows.get(probeNames[i])==null){
    			System.out.println("At row " + i + ": no probe name " + probeNames[i] + "!!");
    			mm--;
    		}
    	}
    	float[][] subData = new float[mm][n];
    	int cnt = 0;
    	for(int i = 0; i < mm; i++){
    		if(rows.get(probeNames[i])!=null){
    			int probeInd = rows.get(probeNames[i]);
    			subProbes.add(probeNames[i]);
    			subRows.put(probeNames[i], i);
    			System.arraycopy(data[probeInd], 0, subData[cnt], 0, n);
    			cnt++;
    		}
    	}
    	
    	return new DataFile(subData, subRows, cols, subProbes, chipID);
    }
    public ArrayList<String> getProbes() {
        return probes;
    }
    
    public ArrayList<String> getChipID(){
    	return chipID;
    }
    
    // 05/13/2011 add setChipID
    public void setChipID(ArrayList<String> newChipID)throws Exception{
    	if(newChipID.size() != n){
    		throw new RuntimeException("New ID has different size from the sample size!");
    	}
    	HashMap<String, Integer> newCols = new HashMap<String, Integer>();
    	for(int i = 0; i < n; i++){
    		newCols.put(newChipID.get(i), i);
    	}
    	chipID = newChipID;
    	cols = newCols;
    }
    
    
    public HashMap<String, Integer> getRows() {
        return rows;
    }

    public HashMap<String, Integer> getCols() {
        return cols;
    }
    public int[] parseRowIndex(String[] in){
    	int n = in.length;
    	int[] idx = new int[n];
    	for(int i = 0; i < n; i++){
    		idx[i] = rows.get(in[i]);
    	}
    	return idx;
    }
    public int[] parseColIndex(String[] in){
    	int n = in.length;
    	int[] idx = new int[n];
    	for(int i = 0; i < n; i++){
    		idx[i] = cols.get(in[i]);
    	}
    	return idx;
    }
    public boolean containsNull(){
    	return nullIndices.size() > 0;
    }
    public int[] getNullIndices(){
    	int n = nullIndices.size();
    	if( n == 0){
    		return null;
    	}
    	int[] idx = new int[n];
    	for(int i = 0; i < n; i++){
    		idx[i] = nullIndices.get(i);
    	}
    	return idx;
    }
    static interface TokenProcessor {
        void process(String[] tokens);
    }

    public static void processTokens(String inFile, TokenProcessor processor, boolean skipHeader, String delimiter) throws IOException {
        BufferedReader in = new BufferedReader(new FileReader(inFile));
        try {
            if (skipHeader) {
                in.readLine();
            }
            String line = in.readLine();
            while (line != null) {
                String[] tokens = getDelimitedTokens(line, delimiter);
                processor.process(tokens);
                line = in.readLine();
            }
        } finally {
            in.close();
        }
    }

    private static String[] getDelimitedTokens(String line, String delim) {
        return line.split(delim);
    }
}
