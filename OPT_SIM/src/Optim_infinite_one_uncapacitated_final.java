


import ilog.concert.*;

import ilog.cplex.*;
import org.apache.commons.math3.distribution.*;
import java.sql.Timestamp;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;
public class Optim_infinite_one_uncapacitated_final {

	public static void main(String[] args) {

		/*
		 * Adjustable variables
		 * 
		 * Variables that typically will be changed are: service level, delta, cv2. Values can be changed inside arrays below
		 */

		double[] serviceLevel 	= {0.8,0.9,0.95}; 			// service level	
		double[] delta 			= {0.5, 1};					// delta used for discretization
		double[] cv2 			= {0.5,1,2}; 				// squared coefficent of variation
		
		double h1 				= 1; 						// holding cost echelon 1
		double mean 			= 50; 						// mean demand per period
		double[] maxCumul		= {0.9998};					// probability of occurence within discretization limit
		int combinations		= serviceLevel.length*cv2.length*delta.length*maxCumul.length;   // number of combinations 
		double[][] results	= new double[combinations][8];  // sL, p, cv2, delta, S1, cost, time

		int count			= 0;
		// nested for-loop to combine all specified variables
			for(int sL = 0; sL<serviceLevel.length; sL++) {			
				for(int c = 0; c<cv2.length; c++) {
					for(int d=0; d<delta.length; d++) {
						for(int mc=0; mc<maxCumul.length; mc++) {

						double p 		= (-h1*serviceLevel[sL])/(serviceLevel[sL]-1);	// backorder cost 
						double variance	= cv2[c]*mean*mean;
						double beta1 	= mean/variance;
						double alpha1 	= beta1*mean;
						double shape	= alpha1; 
						double scale	= 1/beta1;
						GammaDistribution d1 = new GammaDistribution(shape, scale);		// demand distribution of echelon 1		
						
						// determining discretization limit
						double cumulProb1 = 0;
						double discLimit = 1;			
						while(cumulProb1<maxCumul[mc]) {
							cumulProb1 = d1.cumulativeProbability(discLimit);
							discLimit++;
						}					

						// running model and saving results in arrays
						double[] result = model(delta[d], discLimit, d1, h1, p);	// S1, cost, time	
						double[] combiResult = {serviceLevel[sL], cv2[c], p, delta[d], result[0], result[1], result[2], discLimit, maxCumul[mc]}; // sL, cv2, p, delta, S1, cost, time
						results[count] = combiResult;
						count++;
					}			
				}
			}
		}
		// printing all results
		Locale.setDefault(new Locale("en", "US"));
		DecimalFormat df3 = new DecimalFormat("#.###");
		System.out.println("sL\tcv2\tp\td\tS1\tcost\ttime\tdLim\tmc");
		for(int i=0; i<results.length; i++) {
			for(int j=0; j<results[0].length;j++) {
				System.out.print(df3.format(results[i][j])+"\t");
			}
			System.out.print("\n");
		}
	}

	// the actual math programming program using CPLEX
	public static double[] model(double delta, double discLimit, GammaDistribution d1, double h1, double p) {
		double[] result = {0,0,0};
		try {
			IloCplex cplex = new IloCplex();

			// variables
			double h1prime = h1;								// calculate local holding cost
			double pPlush1prime = p + h1prime;					// introduce variable that is sum of backorder costs and h1prime
			int M = 10000; 										// arbitrary large number used in calculations
			int limitD = (int) Math.round(discLimit/delta);		// number of steps in discretization

			// initialization of CPLEX variables
			IloNumVar S1 = cplex.numVar(0, M);		
			IloNumVar y[] = new IloNumVar[limitD];

			for(int k = 0; k<(limitD); k++) {
				y[k] = cplex.boolVar("y"+k);
			}
			IloNumExpr sumYg = cplex.numExpr();
			for(int g=0; g<limitD; g++) {
				sumYg = cplex.sum(sumYg, y[g]);
			}

			IloNumExpr fullexpr = cplex.numExpr(); // objective function

			for(int l=0; l<limitD; l++) {		
				// variables used in expressions
				double p1l = (delta/2)*(d1.density((l+.5)*delta)+d1.density((l+1.5)*delta));	// p(1,l)
				double lDelta = (l+1)*delta;													// l*delta		
				double lDeltap1l = p1l*lDelta;													// l*delta*p(1,l)

				// expression 1
				IloNumExpr expr1 = cplex.prod(h1, cplex.sum(-lDeltap1l,cplex.prod(p1l, S1)));	

				// expression 2
				IloNumExpr expr2 = cplex.prod(cplex.prod(pPlush1prime, cplex.sum(1, cplex.prod(-1, y[l]))), cplex.sum(lDeltap1l, cplex.prod(-p1l, S1)));			

				// add expression 1 and 2 to objective function 
				fullexpr = cplex.sum(fullexpr, expr1, expr2);

				/*
				 * Constraints for consecutive 1's on y(g)
				 */
				if(l>0) {
					cplex.addGe(y[l-1], y[l]);		// y(g-1) >= yg
				}
			}

			/*
			 *  Additional constraint
			 */
			cplex.addEq(S1, cplex.prod(sumYg,delta));		// S1 = sum(yg)*delta
			
			/*
			 * Set objective (minimize) and export model
			 */
			cplex.addMinimize(fullexpr);
			cplex.exportModel ("optimization.lp");	
			cplex.setOut(null);
			long startTime = System.currentTimeMillis();

			// solve the model
			if(cplex.solve()) {
				long elapsedTime = System.currentTimeMillis() - startTime;
				double cost 	= cplex.getObjValue();
				double bestS1	= cplex.getValue(S1);
				double time		= elapsedTime/1000;
				result[0] = bestS1;
				result[1] = cost;
				result[2] = time;
				
				// possible to print individual solutions:
	/*			System.out.println("\nTotal cost:\t"+ cost);		
				System.out.println("S1: \t\t"+ bestS1);
				System.out.println("Comp. time\t"+ time);*/
			}
			else {
				System.out.println("No solution found");
			}
		}
		catch(IloException exc){
			exc.printStackTrace();
		}	

		return result;
	}

}
