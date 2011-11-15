package caf;

import java.util.HashSet;

import util.StatOps;

public class TestField {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		HashSet<Integer> hs = new HashSet<Integer>();
		HashSet<Integer> hs2 = new HashSet<Integer>();
		
		hs.add(1);
		hs.add(3);
		hs2.add(1);
		hs2.add(2);
		System.out.println(hs.equals(hs2));
		
	}

}
