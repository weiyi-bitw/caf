package obj;

/**
 * InverseAnnotations.java
 * @author Wei-Yi Cheng
 * @version 0.22
 * @date 01/27/2011
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class InverseAnnotations {
	HashMap<String, String[]> invAnnot;
	
	public static InverseAnnotations parseInvAnnot(Annotations annot) throws Exception{
	// toDo
		HashMap<String, String[]> invAnnot = new HashMap<String, String[]>();
		return new InverseAnnotations(invAnnot);
	}
	
	public static InverseAnnotations parseInvAnnot(String fileName) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		HashMap<String, String[]> invAnnot = new HashMap<String, String[]>();
		br.readLine();// ignore header
		String line = br.readLine();
		while(line != null){
			String[] tokens = line.split(",");
			if(invAnnot.containsKey(tokens[0])){
				throw new RuntimeException("Error: gene name " + tokens[0] + "duplicated in inverse annotation file!");
			}
			int nt = tokens.length;
			String[] probes = new String[nt-1];
			System.arraycopy(tokens, 1, probes, 0, nt-1);
			invAnnot.put(tokens[0],probes);
			line = br.readLine();
		}
		br.close();
		return new InverseAnnotations(invAnnot);
	}
	InverseAnnotations(HashMap<String, String[]> map) {
        this.invAnnot = map;
    }
	
	public boolean containsGene(String id){
		return invAnnot.containsKey(id);
	}
    public String[] getProbes(String id) {
        return invAnnot.get(id);
    }

    public int getSize() {
        return invAnnot.size();
    }
    public String[] getGenes(){
    	return invAnnot.keySet().toArray(new String[0]);
    }
    public String[] getProbes(String[] gs){
    	HashSet<String> outgene = new HashSet<String>();
    	for(String g : gs){
    		if(invAnnot.get(g)== null) continue;
    		for(String h : invAnnot.get(g)){
    			outgene.add(h);
    		}
    	}
    	String[] out = outgene.toArray(new String[0]);
    	Arrays.sort(out);
    	return out;
    }
}
