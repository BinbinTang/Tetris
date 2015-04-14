package Utility;

public class ParameterSet implements Comparable<ParameterSet> {
	private int[] weight = new int[5];
	private int evaluateStep = 0;
	private double rewardGained = 0;
	public static final int MAXWEIGHT = (int) Math.pow(2, 16);
	
	public ParameterSet(float[] multiplier)
	{
		weight = normalize(multiplier);
	}
	
	public static int[] normalize(float[] multiplier)
	{
		float smallest = Float.MAX_VALUE;
		for (float mult : multiplier)
			if (mult<smallest)
				smallest = mult;
		int[] normalized = new int[5];
		for (int i=0; i<5; i++)
			normalized[i] = (int) (multiplier[i] / smallest);
		return normalized;
	}
	
	@Override
	public int compareTo(ParameterSet o) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	
	public double getValue()
	{
		if (evaluateStep == 0)
			return 0;
		return rewardGained/evaluateStep;
	}
	
}
