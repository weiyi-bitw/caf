package worker;

/**
 * DistributedWorker.java
 * @author Wei-Yi Cheng
 * @version 0.22
 * @date 01/26/2011
 */

import java.io.File;

public class DistributedWorker {
	protected int id;
	protected static int totalComputers;
	protected boolean clean = true; 
	protected static long jobID;
	
	DistributedWorker(){};
	DistributedWorker(int thisSeg, int totalSegs){
		this.id = thisSeg;
		DistributedWorker.totalComputers = totalSegs;
	}
	DistributedWorker(int thisSeg, int totalSegs, long jobID){
		this.id = thisSeg;
		DistributedWorker.totalComputers = totalSegs;
		DistributedWorker.jobID = jobID;
	}
	protected void clearDir(File dir){
		File[] allFiles = dir.listFiles();
		for(File f : allFiles){
			f.delete();
		}
		dir.delete();
	}
	protected void clearDir(String dirName){
		File dir = new File(dirName);
		File[] allFiles = dir.listFiles();
		for(File f : allFiles){
			f.delete();
		}
		dir.delete();
	}
	public void cleanFile(boolean b){
		this.clean = b;
	}
	protected void prepare(String dir){
		boolean success = (new File("tmp")).mkdir();
		if (success) {
		      System.out.println("Create directory tmp");
		}
		success = (new File("tmp/" + jobID)).mkdir();
		success = (new File("tmp/" + jobID + "/" + dir)).mkdir();
		if (success) {
		      System.out.println("Create directory " + dir);
		} 
	}
}
