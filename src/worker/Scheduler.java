package worker;

/**
 * Scheduler.java
 * @author Wei-Yi Cheng
 * @version 0.25
 * @date 06/16/2011
 */


import java.io.File;
import java.io.PrintWriter;
import java.util.Random;

public class Scheduler extends DistributedWorker{
	
	public Scheduler(int thisSeg, int totalSegs, long jobID){
		super(thisSeg, totalSegs, jobID);
		new File(".notify").mkdir();
		new File(".notify/" + jobID).mkdir();
		System.out.println("Scheduler " + (id+1) + " of " + totalComputers + " created.");
	}
	public void outputLog() throws Exception{
		PrintWriter pw = new PrintWriter("spawn.log." + jobID);
		pw.println("Job ID: " + jobID);
		pw.println("Total jobs spawned: " + totalComputers);
		pw.println("Last finished jobs: " + (id+1));
		pw.close();
	}
	
	private void finishFlag() throws Exception{
		
		boolean success = (new File(".finish")).mkdir();
		if (success) {
		      System.out.println("Create directory .finish");
		}
		new File(".finish/" + jobID).mkdir();
		PrintWriter pw = new PrintWriter(".finish/" + jobID + "/flag" + String.format("%05d", id));
		System.out.println("Finish flag created.");
		pw.close();
	}
	private void finishFlag(int iteration) throws Exception{
		
                boolean success = (new File(".finish" + iteration)).mkdir();
                if (success) {
                      System.out.println("Create directory .finish");
                }
                new File(".finish" + iteration + "/" + jobID).mkdir();
                PrintWriter pw = new PrintWriter(".finish" + iteration + "/" + jobID + "/flag" + String.format("%05d", id));
                System.out.println("Finish flag created.");
                pw.close();


	}
	public void waitTillFinished()throws Exception{
		this.finishFlag();
		File dir = new File(".finish/" + jobID);
	    String[] allFiles = dir.list();
	    int numFlags = allFiles.length;
	    if(numFlags == totalComputers){ // last one, clear the .finish directory
	    	clearDir(dir);
		}else{
		    System.out.println("Waiting for finish...");
		    Random r = new Random(id);
		    while(true){
				Thread.sleep(10000 + r.nextInt(10000));
		    	try{
		    		if(!new File(".finish/" + jobID + "/flag" + String.format("%05d", id)).exists()){
		    			break;
		    		}
				}catch(NullPointerException e){
					System.out.println("Warning: Can't find folder (probably finished). Proceed...");
					break;
				}
			}
		}
	    System.out.println("Finished. Proceed...");
	    
		/*
		if(id == 0){
	    	System.out.println("Wait for 20 sec...");
	    	Thread.sleep(20000);
	    	clearDir(dir);
	    }else{
	    	while(dir.exists());
	    }
	    */
	}
	
// Latest version of waitTillFinished -- Use this!!
	public void waitTillFinished(int iteration)throws Exception{
		this.finishFlag(iteration);
		File dir = new File(".finish" + iteration + "/" + jobID);
	    String[] allFiles = dir.list();
	    int numFlags = allFiles.length;
	    System.out.println("Waiting for finish...");
	    Random r = new Random(id);
	    while(numFlags != totalComputers){
			if(id!=0){
				Thread.sleep(100000 + r.nextInt(10000));
			}else{
				Thread.sleep(10000);
				dir = new File(".finish" + iteration + "/" + jobID);
			}
			try{
				allFiles=dir.list();
				numFlags=allFiles.length;
				if(numFlags == totalComputers) break;
				if(id==0){System.out.print(numFlags + "\t");}
			}catch(NullPointerException e){
				System.out.println("Warning: Can't find folder (probably finished). Proceed...");
				break;
			}
		}
		System.out.println("Finished. Proceed...");
	    if(id == 0){
	    	System.out.println("Wait for 10 sec...");
	    	Thread.sleep(10000);
	    	clearDir(dir);
	    }else{
	    	while(dir.exists()){
			Thread.sleep(r.nextInt(10000));
		}
	    }
	}
	
	// job control after fold
	public void waitTillFinished(int iteration, int fold)throws Exception{
		this.finishFlag(iteration);
		if(id >= fold){
			System.out.println("Job finished. Exit.");
			return;
		}
		File dir = new File(".finish" + iteration + "/" + jobID);
	    String[] allFiles = dir.list();
	    int numFlags = allFiles.length;
	    System.out.println("Waiting for finish...");
	    Random r = new Random(id);
	    while(numFlags != totalComputers){
			if(id!=0){
				Thread.sleep(100000 + r.nextInt(10000));
			}else{
				Thread.sleep(10000);
				dir = new File(".finish" + iteration + "/" + jobID);
			}
			try{
				allFiles=dir.list();
				numFlags=allFiles.length;
				if(id==0){System.out.print(numFlags + " (" + fold + ")\t");}
			}catch(NullPointerException e){
				System.out.println("Warning: Can't find folder (probably finished). Proceed...");
				break;
			}
		}
		System.out.println("Finished. Proceed...");
	    if(id == 0){
	    	System.out.println("Wait for 10 sec...");
	    	Thread.sleep(10000);
	    	clearDir(dir);
	    }else{
	    	while(dir.exists()){
			Thread.sleep(r.nextInt(10000));
		}
	    }
	}
	
// 02/08/2011 waitTillNotified(), createNotify() 
	public void waitTillNotified() throws Exception{
		while(!new File(".notify/" + jobID + "/notify-" + String.format("%05d", id)).exists()){}
		System.out.println("Got notify, proceed...");
		new File(".notify/" + jobID + "/notify-" + String.format("%05d", id)).delete();
	}
	public void createNotify() throws Exception{
		for(int i = 0; i < totalComputers; i++){
			new File(".notify/" + jobID + "/notify-" + String.format("%05d", i)).createNewFile();
		}
		new File(".notify/" + jobID + "/notify-" + String.format("%05d", id)).delete();
	}
	public boolean allFinished()throws Exception{
		this.finishFlag();
		if(id != 0){
			System.out.println("Task in this segment is finished. Terminate.");
			return false;
		}
		File dir = new File(".finish/" + jobID);
	    String[] allFiles = dir.list();
	    int numFlags = allFiles.length;
	    System.out.println("Wait other segments to finish...");
	    Random r = new Random(id);
	    while(true){
	    	numFlags = 0;
	    	dir = new File(".finish/" + jobID);
	    	allFiles = dir.list();
		    numFlags = allFiles.length;
		    if(numFlags == totalComputers){
		    	break;
		    }
		    System.out.print(numFlags + "\t");
		    Thread.sleep(r.nextInt(10000));
		}
		System.out.println("Spawn finished. Proceed.");
		clearDir(dir);
		return true;
	}
	// job control after fold
	public boolean allFinished(int fold)throws Exception{
		if(id < fold){
			this.finishFlag();
		}
		if(id != 0){
			System.out.println("Task in this segment is finished. Terminate.");
			return false;
		}
		File dir = new File(".finish/" + jobID);
	    String[] allFiles = dir.list();
	    int numFlags = allFiles.length;
	    System.out.println("Wait other segments to finish...");
	    Random r = new Random(id);
	    while(true){
	    	numFlags = 0;
	    	dir = new File(".finish/" + jobID);
	    	allFiles = dir.list();
		    numFlags = allFiles.length;
		    if(numFlags == fold){
		    	break;
		    }
		    System.out.print(numFlags + " (" + fold + ")\t");
		    Thread.sleep(r.nextInt(10000));
		}
		System.out.println("Spawn finished. Proceed.");
		clearDir(dir);
		return true;
	}
}
