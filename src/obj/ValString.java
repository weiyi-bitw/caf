package obj;

public class ValString implements Comparable<ValString>{
	public float val;
	public String s;
	public ValString(String s, float v){
		this.s = s;
		this.val = v;
	}
	public int hashCode(){
		return s.hashCode();
	}
	
	public int compareTo(ValString other) {
		return -Double.compare(this.val, other.val);
	}
	public String toString(){
		String ss = s + ":" + val;
		return ss;
	}
	public String s(){
		return s;
	}
	public float val(){
		return val;
	}
	public boolean equals(Object other){
		boolean result = false;
        if (other instanceof ValString) {
        	ValString that = (ValString) other;
            result = (this.s.equals(that.s));
        }
        return result;
	}
}
