import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

import statistics.*;


public class SimulationOneStage_final {

	public static void main(String[] args) {		
		Locale.setDefault(new Locale("en", "US"));
		DecimalFormat df1 = new DecimalFormat("#.#");
		DecimalFormat df3 = new DecimalFormat("#.###");

		// Random number generator
		Random rng = new Random();
		rng.setSeed(5);

		// Demand scenarios
		double[] mR					= {9,5,6,7,8,9,10,11,12,13};
		double[] mC 				= {9,9,9,9,9,9,9,9,9,9,9};		
		double[] mD 				= {9,9,10,11,13,11,9,7,6,5};
		double[] mP 				= {9,7,7,7,7,10,12,17,7,7}; 

		// Capacity scenarios
		double[] c12			= {12,12,12,12,12,12,12,12,12,12};
		double[] c11			= {11,11,11,11,11,11,11,11,11,11};
		double[] c10			= {10,10,10,10,10,10,10,10,10,10};
		double[] cDr			= {13,13,13,13,4,4,13,13,13,13};
		double[] cUn			= {500,500,500,500,500,500,500,500,500,500};

		/*
		 * Adjustable variables
		 */

		double[] mean			= mR; 		// demand scenario used
		double[] capacity		= cUn;		// capacity scenario used
		double serviceLevel 	= 0.95; 	// service level
		double cv2				= 0.5; 		// squared coefficient of variation
		double INsim0 			= 12;		// net-inventory at time = 0
		double[] S1				= {0,21,25,29,33,37,41,45,48,31}; // base-stock levels for 10 periods

		double Xsim0 			= 0;		// X at time = 0;
		double Xsim1 			= 9;		// X at time = 1;

		/*
		 * Fixed variables
		 */
		double  h1 			= 1; 			// normalized
		int timeLimit		= 10;			// number of periods
		int runs 			= 50;			// number of runs
		int repetitions 	= 10000;		// number of repetitions

		/*
		 * Calculate costs of optimization
		 */

		double p 		= (-h1*serviceLevel)/(serviceLevel-1);	// backorder cost 

		double[] varianceD 		= new double[timeLimit];						
		double[] shapeD			= new double[timeLimit];
		double[] scaleD			= new double[timeLimit];		

		// Distributions for simulation
		GammaDistributionReturnValues[] d1 = new GammaDistributionReturnValues[timeLimit];
		for(int t=0; t<timeLimit; t++) {
			varianceD[t]	= cv2*mean[t]*mean[t];
			scaleD[t]		= varianceD[t]/mean[t];
			shapeD[t]		= mean[t]/scaleD[t];
			d1[t]			= new GammaDistributionReturnValues(shapeD[t], 1/scaleD[t], rng);		
		}

		/*
		 * Simulate with random demand, with original values for S1[t]
		 */
		double sumCostOriginal 	= 0;
		double sum2CostOriginal = 0;
		double[] sumI 			= new double[timeLimit];
		double[] sumB 			= new double[timeLimit];
		double[] sumCost		= new double[timeLimit];
		double[] sumIN 			= new double[timeLimit];
		double[] sumX			= new double[timeLimit];
		double[] sumDemand		= new double[timeLimit];
		for(int r=0; r<runs; r++) {
			double[] randomDemand 	= new double[timeLimit];
			double[] INsim			= new double[timeLimit];
			double[] Xsim			= new double[timeLimit];
			double[] costs			= new double[timeLimit];
			INsim[0] 		= INsim0;
			Xsim[0] 		= Xsim0;
			Xsim[1] 		= Xsim1;
			double sumCostRun		= 0;	
			for(int rp = 0; rp<repetitions; rp++) {
				sumX[1]			+= Xsim[1];
				for(int t=1; t<timeLimit; t++) {
					randomDemand[t] = Math.round(d1[t].nextRandom());	
					INsim[t] 		= INsim[t-1] + Xsim[t-1] - randomDemand[t];		
					costs[t] 		= Math.max(0, INsim[t])*h1 - Math.min(INsim[t], 0)*p;
					sumCostRun 		+= costs[t];		
					sumI[t]			+= Math.max(0, INsim[t]);
					sumB[t]			+= Math.max(0, -INsim[t]);
					sumCost[t]		+= costs[t];
					sumIN[t] 		+= INsim[t];
					sumDemand[t]	+= randomDemand[t];
					if(t<timeLimit-1) {
						Xsim[t+1] = Math.max(0, Math.min(capacity[t+1], S1[t+1]-INsim[t]-Xsim[t]));
						sumX[t+1] += Xsim[t+1];
					}
				}
			}
			double meanCostRun		= sumCostRun/repetitions;
			sumCostOriginal += meanCostRun;
			sum2CostOriginal += meanCostRun*meanCostRun;
		}
		double[] meanI	= new double[timeLimit];
		double[] meanB	= new double[timeLimit];
		double[] meanCost	= new double[timeLimit];
		double[] meanX	= new double[timeLimit];
		double[] meanIN	= new double[timeLimit];
		double[] meanDemand	= new double[timeLimit];
		for(int t=1; t<timeLimit; t++) {
			meanI[t]			= sumI[t]/(runs*repetitions);
			meanB[t]			= sumB[t]/(runs*repetitions);
			meanCost[t]			= sumCost[t]/(runs*repetitions);
			meanX[t]			= sumX[t]/(runs*repetitions);
			meanIN[t]			= sumIN[t]/(runs*repetitions);
			meanDemand[t]		= sumDemand[t]/(runs*repetitions);
		}
		double meanCostOriginal = sumCostOriginal/runs;
		double varCostOriginal	= sum2CostOriginal/runs - meanCostOriginal*meanCostOriginal;
		double halfwidthCostOriginal	= 1.96*Math.sqrt(varCostOriginal/runs);
		double[] CICostOriginal			= {meanCostOriginal-halfwidthCostOriginal, meanCostOriginal+halfwidthCostOriginal};
		System.out.println();
		System.out.println("Original S\t" + Arrays.toString(S1));
		System.out.println("Cost\t\t"+df3.format(meanCostOriginal));
		System.out.println("CI cost\t"+df3.format(CICostOriginal[0])+"\t"+df3.format(CICostOriginal[1]));

		System.out.print("\nrealMu:\t");
		for(int t=1; t<timeLimit; t++) {
			System.out.print(df1.format(meanDemand[t])+"\t");
		}
		System.out.print("\nrealX:\t");
		for(int t=1; t<timeLimit; t++) {
			System.out.print(df1.format(meanX[t])+"\t");
		}
		System.out.print("\nrealI:\t");
		for(int t=1; t<timeLimit; t++) {
			System.out.print(df1.format(meanI[t])+"\t");
		}
		System.out.print("\nrealB:\t");
		for(int t=1; t<timeLimit; t++) {
			System.out.print(df1.format(meanB[t])+"\t");
		}
		System.out.print("\nrealIN:\t");
		for(int t=1; t<timeLimit; t++) {
			System.out.print(df1.format(meanIN[t])+"\t");
		}
		System.out.print("\nrealCs:\t");
		for(int t=1; t<timeLimit; t++) {
			System.out.print(df1.format(meanCost[t])+"\t");
		}		 
		System.out.println("\n\n-----------------------------------------------------");


	}
}