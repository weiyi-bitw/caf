package obj;

public class ValIdxD implements Comparable<ValIdxD>{
	public double val;
	public int idx;
	public ValIdxD(int i, double v){
		this.idx = i;
		this.val = v;
	}
	public int hashCode(){
		return idx;
	}
	
	public int compareTo(ValIdxD other) {
		return -Double.compare(this.val, other.val);
	}
	
	public int idx(){
		return idx;
	}
	public double val(){
		return val;
	}
	public boolean equals(Object other){
		boolean result = false;
        if (other instanceof ValIdx) {
        	ValIdx that = (ValIdx) other;
            result = (this.idx == that.idx);
        }
        return result;
	}
}
