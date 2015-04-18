

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Random;
import java.util.Scanner;
import java.util.Vector;



public class PlayerSkeleton {
	
	// Global result for lotus swarm
	static int SWARM_SIZE = 1000;
	static int MAX_ITERATION = 500;
	static int PROBLEM_DIMENSION = 5;
	static double C1 = 2.0;
	static double C2 = 2.0;
	static double W_UPPERBOUND = 1.0;
	static double W_LOWERBOUND = 0.0;
	static double LOW = -100;
	static double HIGH = 3;
	static double VEL_LOW = -1;
	static double VEL_HIGH = 1;
	static double TEST_TIME = 2;
	
	
	static double[] result_particle = null;
	static int result_best = 0;
	
	public String printArrInt(int[] arr){
		String s="{";
		for(int i=0; i<arr.length; i++){
			s+=arr[i]+", ";
		}
		s+="},";
		return s;
	}
	public static String printArrDouble(double[] arr){
		String s="";
		for(int i=0; i<arr.length; i++){
			s+=arr[i]+" ";
		}
		return s;
	}
	public String print2DArr(int[][] arr){
		StringBuilder sb = new StringBuilder();
		for(int i=0; i<arr.length; i++){
			for(int j=0; j<arr[i].length; j++){
				sb.append(arr[i][j]+" ");
			}
			sb.append("\n");
		}
		return sb.toString();
	}
	public int[] getColumnHeights (State s){
		int h[] = new int[s.COLS];
		System.arraycopy(s.getTop(), 0, h, 0, s.COLS);
		return h;
	}
	public int[] getHeightDiff (State s){
		int hdiff[] = new int[s.COLS-1];
		int h[] = s.getTop();
		for(int i=0; i<hdiff.length; i++){
			hdiff[i]=Math.abs(h[i+1]-h[i]);
		}
		return hdiff;
	}
	public double getHeightMean (State s){
		int h[] = s.getTop();
		double sum = 0;
		for(int i=0; i<h.length; i++){
			sum += h[i];
		}
		return sum/h.length;
	}
	public double getHeightStd (State s){
		int h[] = s.getTop();
		double sum = 0;
		for(int i=0; i<h.length; i++){
			sum += h[i];
		}
		double mean = sum/h.length;
		double var = 0;
		for(int i=0; i<h.length; i++){
			var += (h[i]-mean)*(h[i]-mean);
		}
		return Math.sqrt(var/h.length);
	}
	public int getMaxHeight (State s){
		int h[] = s.getTop();
		int max = h[0];
		for(int i=1; i<h.length; i++){
			if(max<h[i]) max = h[i];
		}
		return max;
	}
	public int getHoleCount (int[][] field){

		int[][] holes = new int[field.length][field[0].length];
		
		int count =0;
		for(int i=field.length-2; i>=0; i--){
			for(int j=0; j<field[i].length; j++){
				if(field[i][j]==0){
					//check left
					if(!(j==0 || field[i][j-1]!=0 || holes[i][j-1]==1)){
						continue;
					}
					//check up (the field is up-down inverted)
					if(!(field[i+1][j]!=0 || holes[i+1][j]==1)){
						continue;
					}
					//check right
					if(j==field[i].length-1){
						holes[i][j]=1;
						count++;
					}else if(field[i][j+1]!=0){
						holes[i][j]=1;
						count++;
					}else{
						//right is zero, need to check if right is a hole
						for(int k=j+1; k<field[i].length; k++){
							//check up (the field is up-down inverted)
							if(!(field[i+1][k]!=0 || holes[i+1][k]==1)){
								//not a hole
								break;
							}
							//check right
							if(k==field[i].length-1 || field[i][k+1]!=0){
								holes[i][j]=1;
								count++;
								break;
							}
						}
						continue;
					}
				}	
			}
		}
		return count;
	}
	public double evaluate0(State s, State os, double[] weights){
		int[] heights 		= getColumnHeights(s);
		int[] heightDiff 	= getHeightDiff(s);
		int maxHeight 		= getMaxHeight(s);
		int holeIncreased 	= getHoleCount(s.getField());
		
		int numFeatures = 22;
		int[] features = new int[numFeatures];
		System.arraycopy(heights, 0, features, 0, heights.length);
		System.arraycopy(heightDiff, 0, features, heights.length, heightDiff.length);
		features[features.length-3]=maxHeight;
		features[features.length-2]=holeIncreased;
		features[features.length-1]=1;
		
		double V = 0;
		for(int i=0; i<heights.length; i++){
			V+=features[i]*weights[0];
		}
		for(int i=heights.length; i<heights.length+heightDiff.length; i++){
			V+=features[i]*weights[1];
		}
		V+=features[features.length-3]*weights[weights.length-3];
		V+=features[features.length-2]*weights[weights.length-2];
		V+=features[features.length-1]*weights[weights.length-1];
		return V;
	}
	
	public double evaluate1(State s, State os, double[] weights){
		int[] heights 		= getColumnHeights(s);
		int[] heightDiff 	= getHeightDiff(s);
		int maxHeight 		= getMaxHeight(s);
		int holeIncreased 	= getHoleCount(s.getField())-getHoleCount(os.getField());
		
		int numFeatures = 22;
		int[] features = new int[numFeatures];
		System.arraycopy(heights, 0, features, 0, heights.length);
		System.arraycopy(heightDiff, 0, features, heights.length, heightDiff.length);
		features[features.length-3]=maxHeight;
		features[features.length-2]=holeIncreased;
		features[features.length-1]=1;
		
		double V = 0;
		for(int i=0; i<heights.length; i++){
			V+=features[i]*weights[0];
		}
		for(int i=heights.length; i<heights.length+heightDiff.length; i++){
			V+=features[i]*weights[1];
		}
		V+=features[features.length-3]*weights[weights.length-3];
		V+=features[features.length-2]*weights[weights.length-2];
		V+=features[features.length-1]*weights[weights.length-1];
		return V;
	}
	public double evaluate2(State s, State os, double[] weights){

		int numFeatures 	= 5;
		double[] features 	= new double[numFeatures];
		features[0]			= getHeightMean(s);
		features[1]			= getHeightStd(s);
		features[2]			= getMaxHeight(s);
		features[3]			= getHoleCount(s.getField())-getHoleCount(os.getField());
		features[4]			= 1;	//bias term
		//System.out.println(printArrDouble(features));
		//double[] weights 	= {-0.0, -0.8, -0.1, -0.1, -0.0};
		
		double V = 0;
		for(int i=0; i<features.length; i++){
			V+=features[i]*weights[i];
		}
		return V;
	}
	public double evaluate3(State s, State os, double[] weights){

		int numFeatures 	= 5;
		double[] features 	= new double[numFeatures];
		features[0]			= getHeightMean(s)-getHeightMean(os);
		features[1]			= getHeightStd(s);
		features[2]			= getMaxHeight(s)-getMaxHeight(os);
		features[3]			= getHoleCount(s.getField())-getHoleCount(os.getField());
		features[4]			= 1;	//bias term
		//System.out.println(printArrDouble(features));
		//double[] weights 	= {-0.0, -0.8, -0.1, -0.1, -0.0};
		
		double V = 0;
		for(int i=0; i<features.length; i++){
			V+=features[i]*weights[i];
		}
		return V;
	}
	
	public int pickMove(State s, int[][] legalMoves, double[] weights) {

		double maxU = -1000000;
		double maxV = 0;
		int maxR =0;
		int maxUAction = 0;
		
		for(int i=0; i<legalMoves.length; i++){
			//simulate a move
			SimState ss = new SimState();
			ss.setField(s.getField());
			ss.setNextPiece(s.getNextPiece());
			ss.setTop(s.getTop());
			ss.makeMove(legalMoves[i]);
			
			//avoid lost move
			if(ss.hasLost()){
				//System.out.println("lost move idx = "+i);
				continue;
			}
			
			//evaluate move
			double V = evaluate1(ss, s, weights);
			int R = ss.getRowsCleared();
			double U = R+V;
			
			//update selection
			if(U>maxU){
				maxU = U;
				maxUAction = i;
				maxV = V;
				maxR = R;
			}
		}
		//System.out.println("Action = "+ printArrInt(legalMoves[maxUAction]));
		//System.out.println ("Action="+maxUAction+" U="+maxU+" R="+maxR+" V="+maxV);
		
		return maxUAction;
	}
	public int pickMoveAdvance(State s, int[][] legalMoves, double[] weights) {

		double maxU = -1000000;
		double maxV = 0;
		int maxR =0;
		int maxUAction = 0;
		
		double[] bestFiveU = new double[5];
		int[] bestFiveAction = new int[5];
		int currentNumberOfBestFive = 0;
		double minBestFiveU = Integer.MAX_VALUE;
		int minBestFiveAction = -1;
		int minBestFiveIndex = -1;
		
		
		
		for(int i=0; i<legalMoves.length; i++){
			//simulate a move
			SimState ss = new SimState();
			ss.setField(s.getField());
			ss.setNextPiece(s.getNextPiece());
			ss.setTop(s.getTop());
			ss.makeMove(legalMoves[i]);
			
			//avoid lost move
			if(ss.hasLost()){
				//System.out.println("lost move idx = "+i);
				continue;
			}
			
			//evaluate move
			double V = evaluate1(ss, s, weights);
			int R = ss.getRowsCleared();
			double U = R+V;
			
			//update best five move
			if(currentNumberOfBestFive <5){
				bestFiveU[currentNumberOfBestFive] = U;
				bestFiveAction[currentNumberOfBestFive] = i;
				if (U < minBestFiveU)
				{
					minBestFiveU = U;
					minBestFiveAction = i;
					minBestFiveIndex = currentNumberOfBestFive;
				}
				currentNumberOfBestFive++;
			}
			else{
				if (U > minBestFiveU)
				{
					bestFiveU[minBestFiveIndex] = U;
					bestFiveAction[minBestFiveIndex] = i;
					minBestFiveU = bestFiveU[0];
					minBestFiveAction = bestFiveAction[0];
					minBestFiveIndex = 0;
					for (int j=1; j<5; j++)
					{
						if (bestFiveU[j]<minBestFiveU)
						{
							minBestFiveU = bestFiveU[j];
							minBestFiveAction = bestFiveAction[j];
							minBestFiveIndex = j;
						}
					}
				}
			}

		}
		
		//for (int i=0; i<5; i++)
		//	System.out.print(bestFiveU[i]+" ");
		//System.out.println();
		

		//System.out.println("Action = "+ printArrInt(legalMoves[maxUAction]));
		//System.out.println ("Action="+maxUAction+" U="+maxU+" R="+maxR+" V="+maxV);
		
		int MaxU = -100000000;
		int MaxIndex = 0;
		for (int current=0; current<5; current++)
		{
			int sumCurrentMoveU = 0;
			SimState ss = new SimState();
			ss.setField(s.getField());
			ss.setNextPiece(s.getNextPiece());
			ss.setTop(s.getTop());
			ss.makeMove(bestFiveAction[current]);
			for (int i=0; i<7; i++)
			{
				SimState sss = new SimState();
				sss.setField(ss.getField());
				sss.setNextPiece(i);
				sss.setTop(ss.getTop());
				
				maxU = -10000000;
				maxV = 0;
				maxR =0;
				maxUAction = 0;
				
				for (int j=0; j<sss.legalMoves().length; j++)
				{
					SimState ssss = new SimState();
					ssss.setField(sss.getField());
					ssss.setNextPiece(i);
					ssss.setTop(sss.getTop());
					ssss.makeMove(sss.legalMoves()[j]);
					//System.out.println("J = "+j);
					
					//avoid lost move
					if(ssss.hasLost()){
						//System.out.println("lost move idx = "+i);
						continue;
					}
					
					//evaluate move
					double V = evaluate1(ssss, ss, weights);
					int R = ssss.getRowsCleared();
					double U = R+V;
					
					if(U>maxU){
						maxU = U;
						maxUAction = i;
						maxV = V;
						maxR = R;
					}
					
				}
				//update selection
				sumCurrentMoveU += maxU;
				
			}
			if (sumCurrentMoveU>MaxU)
			{
				MaxU = sumCurrentMoveU;
				MaxIndex = current;
			}
		}
		
		return bestFiveAction[MaxIndex];
	}
	
	
	public static int runOnce(double[] weights, boolean visualOn){
		if(visualOn){
			State s = new State();
			TFrame t = new TFrame(s);
			PlayerSkeleton p = new PlayerSkeleton();
			int move =0;
			while(!s.hasLost()) {
				move++;
				s.makeMove(p.pickMoveAdvance(s,s.legalMoves(),weights));
				
				s.draw();
				s.drawNext(0,0);
				
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			int score = s.getRowsCleared();
			//System.out.println("moves made = "+move);
			//System.out.println("You have completed "+score+" rows.");
			
			t.dispose();
			return score;
		}else{
			State s = new State();
			PlayerSkeleton p = new PlayerSkeleton();
			int move =0;
			while(!s.hasLost()) {
				move++;
				s.makeMove(p.pickMove(s,s.legalMoves(),weights));
			}
			int score = s.getRowsCleared();
			//System.out.println("moves made = "+move);
			//System.out.println("You have completed "+score+" rows.");
			return score;
		}
	}	
	
	public static double[] PSOExecute (){
		Vector<Particle> swarm = new Vector<Particle>();
		double[] pBest = new double[SWARM_SIZE];
		Vector<Location> pBestLocation = new Vector<Location>();
		double gBest = 0;
		Location gBestLocation = null;
		double[] fitnessValueList = new double[SWARM_SIZE];
		
		Random generator = new Random();
		
		Particle p;
		for(int i=0; i<SWARM_SIZE; i++) {
			p = new Particle();
			
			// randomize location inside a space defined in Problem Set
			double[] loc = new double[PROBLEM_DIMENSION];
			for (int j=0; j<PROBLEM_DIMENSION; j++)
				loc[j] = LOW + generator.nextDouble()*(HIGH-LOW);
			Location location = new Location(loc);
			
			// randomize velocity in the range defined in Problem Set
			double[] vel = new double[PROBLEM_DIMENSION];
			for (int j=0; j<PROBLEM_DIMENSION; j++)
				vel[j] = VEL_LOW + generator.nextDouble()*(VEL_HIGH-VEL_LOW);
			
			Velocity velocity = new Velocity(vel);
			
			p.setLocation(location);
			p.setVelocity(velocity);
			swarm.add(p);
		}
		
		for(int i=0; i<SWARM_SIZE; i++) {
			double ans = 0;
			for (int j=0; j<TEST_TIME; j++)
				ans += runOnce(swarm.get(i).getLocation().getLoc(), false);
			fitnessValueList[i] = ans/TEST_TIME;
		}
		
		for(int i=0; i<SWARM_SIZE; i++) {
			pBest[i] = fitnessValueList[i];
			pBestLocation.add(swarm.get(i).getLocation());
		}
		
		int t = 0;
		double w;
		
		while(t < MAX_ITERATION) {
			// step 1 - update pBest
			for(int i=0; i<SWARM_SIZE; i++) {
				if(fitnessValueList[i] > pBest[i]) {
					pBest[i] = fitnessValueList[i];
					pBestLocation.set(i, swarm.get(i).getLocation());
				}
			}
				
			// step 2 - update gBest
			int bestParticleIndex = PSOUtility.getMinPos(fitnessValueList);
			if(t == 0 || fitnessValueList[bestParticleIndex] > gBest) {
				gBest = fitnessValueList[bestParticleIndex];
				gBestLocation = swarm.get(bestParticleIndex).getLocation();
			}
			
			w = W_UPPERBOUND - (((double) t) / MAX_ITERATION) * (W_UPPERBOUND - W_LOWERBOUND);
			
			for(int i=0; i<SWARM_SIZE; i++) {
				double r1 = generator.nextDouble();
				double r2 = generator.nextDouble();
				
				p = swarm.get(i);
				
				// step 3 - update velocity
				double[] newVel = new double[PROBLEM_DIMENSION];
				for (int j=0; j<PROBLEM_DIMENSION; j++)
					newVel[j] = (w * p.getVelocity().getPos()[j]) + 
							(r1 * C1) * (pBestLocation.get(i).getLoc()[j] - p.getLocation().getLoc()[j]) +
							(r2 * C2) * (gBestLocation.getLoc()[j] - p.getLocation().getLoc()[j]);
				Velocity vel = new Velocity(newVel);
				p.setVelocity(vel);
				
				// step 4 - update location
				double[] newLoc = new double[PROBLEM_DIMENSION];
				for (int j = 0; j<PROBLEM_DIMENSION; j++)
					newLoc[j] = p.getLocation().getLoc()[j]+ newVel[j];
				Location loc = new Location(newLoc);
				p.setLocation(loc);
			}
			
			System.out.println("ITERATION " + t + ": ");
			System.out.println(printArrDouble(gBestLocation.getLoc()));
			//System.out.println(runOnce(gBestLocation.getLoc(), true));
			System.out.println(gBest);
			
			t++;
			
			for(int i=0; i<SWARM_SIZE; i++) {
				double ans = 0;
				for (int j=0; j<TEST_TIME; j++)
					ans += runOnce(swarm.get(i).getLocation().getLoc(), false);
				fitnessValueList[i] = ans/TEST_TIME;
			}
		}
		
		
		return gBestLocation.getLoc();
	}
	
	public static void main(String[] args) {	
		
		
		//so far best weights for evaluate2
		//double[] bestWeights= {-0.81, -0.35, -0.97, -0.04, -0.13};
		
		//so far best weights for evaluate1
		//double[] bestWeights= {-0.79, -0.15, -0.14, -0.1, -0.75};
		
		//so far best weights for evaluate3
		//double[] bestWeights= {-138.51506302626194,-26.13212253224357,-0.8955877227084024,-36.050852075521725,-76.97989990999787 };
		
		//so far best weights for evaluate0
		
		// Test
		double[] bestWeight = /*{-8.94571873, -4.3118907, -5.74438422, 0.341797307, 1.044628353}; */PSOExecute();
		//System.out.println(printArrDouble(bestWeight));
		
		//double max = 0;
		double ans = 0;
		for (int j=0; j<TEST_TIME; j++)
			ans += runOnce(bestWeight, true);
		ans = ans/TEST_TIME;
		System.out.println(ans);
	}
	
}
