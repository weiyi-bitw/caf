package obj;

public class IntPair {
	public int x;
	public int y;
	
	public IntPair(int x, int y){
		if(x <= y){
			this.x = x;
			this.y = y;
		}else{
			this.x = y;
			this.y = x;
		}
	}
	
	public boolean overlapWith(IntPair ip2){
		if(this.x > ip2.y && this.y > ip2.y){
			return false;
		}else if(this.x < ip2.x && this.y < ip2.x){
			return false;
		}else{
			return true;
		}
	}
	
	public String toString(){
		String s = "[ " + x + " , " + y +" ]";
		return s;
	}
	
	public int range(){
		return y-x;
	}
}
