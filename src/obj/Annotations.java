package obj;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * @author John Watkinson
 */
public class Annotations {

    public static Annotations parseAnnotations(String file) throws Exception {
        final HashMap<String,String> geneName = new HashMap<String, String>();
        DataFile.processTokens(file, new DataFile.TokenProcessor() {
            boolean start = false;
            int geneIndex = -1;
            public void process(String[] tokens) {
                if (!start) {
                    if (DataFile.stripQuotes(tokens[0]).equals("Probe Set ID")||DataFile.stripQuotes(tokens[0]).equals("ProbeSet") ) {
                        start = true;
                        for (int i = 1; i < tokens.length; i++) {
                            if (DataFile.stripQuotes(tokens[i]).equals("Gene Symbol")||DataFile.stripQuotes(tokens[i]).equals("GeneSymbol")) {
                                geneIndex = i;
                                break;
                            }
                        }
                        if (geneIndex == -1) {
                            throw new RuntimeException("Couldn't find 'Gene Symbol' annotation.");
                        }
                    }
                } else {
                    String id = DataFile.stripQuotes(tokens[0]);
                    String gene = DataFile.stripQuotes(tokens[geneIndex]);
                    if (gene == null || gene.length() == 0 || gene.startsWith("-")) {
                        gene = id;
                    }
                    geneName.put(new String(id), new String(gene));
                }
            }
        }, false, "\",\"|,");
        return new Annotations(geneName);
    }

    private HashMap<String,String> geneName;

    private Annotations(HashMap<String, String> geneName) {
        this.geneName = geneName;
    }
    
    public String getGene(String id) {
        return geneName.containsKey(id)? geneName.get(id):"NA";
    }
    
    public int getSize() {
        return geneName.size();
    }
    public String[] getAllGenes(){
    	HashSet<String> allpbs = new HashSet<String>(geneName.keySet());
    	HashSet<String> allgenes = new HashSet<String>(); 
    	for(String s : allpbs){
    		String g = geneName.get(s);
    		if(g.contains(" /// ") || g.contains("_at")) continue;
    		allgenes.add(g);
    	}
    	String[] out = allgenes.toArray(new String[0]);
    	Arrays.sort(out);
    	return out;
    }
    public InverseAnnotations getInvAnnot(){
    	HashMap<String, String[]> invmap = new HashMap<String, String[]>();
    	HashSet<String> allpbs = new HashSet<String>(geneName.keySet());
    	for(String s : allpbs){
    		String g = geneName.get(s);
    		if(g.contains(" /// ") || g.contains("_at")) continue;
    		if(invmap.get(g)==null){
    			String[] ss = {s};
    			invmap.put(g, ss);
    		}else{
    			String[] ss = invmap.get(g);
    			String[] sss = new String[ss.length + 1];
    			System.arraycopy(ss, 0, sss, 0, ss.length);
    			sss[sss.length-1] = s;
    			invmap.put(g, sss);
    		}
    	}
    	return new InverseAnnotations(invmap);
    	
    }
}
