package obj;

public class ValIdx implements Comparable<ValIdx>{
	float val;
	int idx;
	public ValIdx(int i, float v){
		this.idx = i;
		this.val = v;
	}
	public int hashCode(){
		return idx;
	}
	
	public int compareTo(ValIdx other) {
		return -Double.compare(this.val, other.val);
	}
	
	public int idx(){
		return idx;
	}
	public float val(){
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