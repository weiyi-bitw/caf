package obj;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class DataFileD {
	static class Size {
        int n, m;
    }
	
	 private double[][] data;
	 private HashMap<String, Integer> rows;
	 private HashMap<String, Integer> cols;
	 private ArrayList<String> probes;
	 private ArrayList<String> chipID;
	 private int n; // number of columns
	 private int m; // number of rows : m-by-n !!!
	 public void log2normalization(){
	    	for(int i = 0; i < m; i++){
	    		for(int j = 0; j < n; j++){
	    			data[i][j] = (Math.log(data[i][j]) / Math.log(2));
	    		}
	    	}
    }
	    
    public DataFileD(double[][] data, HashMap<String, Integer> rows, HashMap<String, Integer> cols, ArrayList<String> probes, ArrayList<String> chipID) {
        this.data = data;
        this.rows = rows;
        this.cols = cols;
        this.probes = probes;
        this.chipID = chipID;
        m = rows.size();
        n = cols.size();
    }

    public static String stripQuotes(String s) {
        s = s.replace("'", "");
    	return s.replace("\"", "");
    }
    static interface TokenProcessor {
        void process(String[] tokens);
    }
    private static String[] getDelimitedTokens(String line, String delim) {
        return line.split(delim);
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
    public static DataFileD parse(String file) throws Exception {
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
        final double[][] data = new double[m][n];
        final HashMap<String,Integer> rows = new HashMap<String, Integer>();
        final HashMap<String,Integer> cols = new HashMap<String, Integer>();
        final ArrayList<String> probes = new ArrayList<String>();
        final ArrayList<String> chipID = new ArrayList<String>();
        processTokens(file, new TokenProcessor() {
            boolean start = true;
            int row = 0;
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
                            	data[row][i - 1] = Double.NaN;
                            }
                            else	data[row][i - 1] = Double.parseDouble(stripQuotes(tokens[i]));
                        } catch (NumberFormatException nfe) {
                            throw new RuntimeException("Couldn't parse number at row " + (row + 1) + ", column " + i + ":" + stripQuotes(tokens[i]));
                        }
                    }
                    row++;
                }
            }
        }, false, "\t");
        return new DataFileD(data, rows, cols, probes, chipID);
    }
    public int getNumRows() {
        return m;
    }

    public int getNumCols() {
        return n;
    }

    public double[][] getData() {
        return data;
    }

    public DataFileD getSubChips(String[] chipNames){
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
    	double[][] subData = new double[m][cnt];
    	for(int i = 0; i < m; i++){
    		for(int j = 0; j < cnt; j++){
    			subData[i][j] = data[i][chipInd[j]];
    		}
    	}
    	return new DataFileD(subData, rows, subCols, probes, subChipNames);
    }
    public DataFileD getSubProbes(ArrayList<ValIdx> probeIdx){
    	int numProbes = probeIdx.size();
    	HashMap<String, Integer> subRows = new HashMap<String, Integer>();
    	ArrayList<String> subProbes = new ArrayList<String>();
    	int mm = numProbes;
    	
    	double[][] subData = new double[mm][n];
    	for(int i = 0; i < mm; i++){
    		int idx = probeIdx.get(i).idx();
    		String p = probes.get(idx);
    		subProbes.add(p);
    		subRows.put(p, i);
    		System.arraycopy(data[idx], 0, subData[i], 0, n);
    	}
    	return new DataFileD(subData, subRows, cols, subProbes, chipID);
    }
    
    public DataFileD getSubProbes(String[] probeNames){
    	int numProbes = probeNames.length;
    	HashMap<String, Integer> subRows = new HashMap<String, Integer>();
    	ArrayList<String> subProbes = new ArrayList<String>();
    	int mm = numProbes;
    	for(int i = 0; i < numProbes; i++){
    		if(rows.get(probeNames[i])==null){
    			//System.out.println("At row " + i + ": no probe name " + probeNames[i] + "!!");
    			mm--;
    		}
    	}
    	double[][] subData = new double[mm][n];
    	int cnt = 0;
    	for(int i = 0; i < numProbes; i++){
    		if(rows.get(probeNames[i])!=null){
    			int probeInd = rows.get(probeNames[i]);
    			subProbes.add(probeNames[i]);
    			subRows.put(probeNames[i], cnt);
    			System.arraycopy(data[probeInd], 0, subData[cnt], 0, n);
    			cnt++;
    		}
    	}
    	
    	return new DataFileD(subData, subRows, cols, subProbes, chipID);
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
    
 // 04/20/2011 Sort probes lexicographically
    public void sortProbes(){
    	Collections.sort(probes);
    	double[][] newData = new double[m][n];
    	int j;
    	for(int i = 0; i < m; i++){
    		String gene = probes.get(i);
    		j = rows.get(gene);
    		System.arraycopy(data[j], 0, newData[i], 0, n);
    		rows.put(gene, i);
    	}
    	data = newData;
    }
}
