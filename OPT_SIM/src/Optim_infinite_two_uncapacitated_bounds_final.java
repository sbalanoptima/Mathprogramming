

import ilog.concert.*;

import ilog.cplex.*;
import org.apache.commons.math3.distribution.*;
import java.sql.Timestamp;
import java.text.DecimalFormat;
import java.util.Locale;
public class Optim_infinite_two_uncapacitated_bounds_final {

	public static void main(String[] args) {
		//print settings
		Locale.setDefault(new Locale("en", "US"));
		DecimalFormat df3 = new DecimalFormat("#.###");

		/*
		 * Adjustable variables
		 * 
		 * Variables that typically will be changed are: service level, delta, cv2. Values can be changed inside arrays below
		 */


		double[] serviceLevel 	= {0.8,0.9,0.95}; 			// service level	
		double[] delta 			= {0.5,1};					// delta used for discretization
		double[] cv2 			= {0.5,1,2}; 				// squared coefficent of variation

		double h1 				= 1; 						// holding cost echelon 1
		double h2 				= 1; 						// holding cost echelon 2
		double mean				= 10;						// mean demand per period
		double[] maxCumul		= {0.9998};					// probability of occurence within discretization limit
		
		int combinations	= serviceLevel.length*cv2.length*delta.length*maxCumul.length;   // number of combinations 
		double[][] results	= new double[combinations][12]; // sL, p, cv2, delta, S1, cost, time

		int count			= 0;
		// nested for-loop to combine all specified variables
		for(int sL = 0; sL<serviceLevel.length; sL++) {
			for(int c = 0; c<cv2.length; c++) {
				for(int d=0; d<delta.length; d++) {	
					for(int mc=0; mc<maxCumul.length; mc++) {

						long startTime = System.currentTimeMillis();		// to measure running time				
						double p 		= (-h1*serviceLevel[sL])/(serviceLevel[sL]-1);	// backorder cost 
						double variance	= cv2[c]*mean*mean;
						double beta		= mean/variance;
						double alpha	= beta*mean;

						double alpha1	= 2*alpha; // delay of information adds extra period
						double shape1	= alpha1;
						double scale1	= 1/beta;

						GammaDistribution d1 = new GammaDistribution(shape1, scale1);	// demand distribution of echelon 1	

						double alpha2	= alpha;
						double shape2	= alpha2;
						double scale2	= 1/beta;
						GammaDistribution d2 = new GammaDistribution(shape2, scale2);	// demand distribution of echelon 2

						// determining discretization limit 1
						double cumulProb1 = 0;
						double discLimit1 = 1;
						while(cumulProb1<maxCumul[mc]) {
							cumulProb1 = d1.cumulativeProbability(discLimit1);
							discLimit1++;						
						}			

						// determining discretization limit 2
						double cumulProb2 = 0;
						double discLimit2 = 1;
						while(cumulProb2<maxCumul[mc]) {
							cumulProb2 = d2.cumulativeProbability(discLimit2);
							discLimit2++;
						}			

						// generating and printing intermediate results:
						System.out.println("discLimit1\t "+discLimit1+"\tdiscLimit2\t"+discLimit2);
						System.out.println("----------\titeration:\t"+count+"\t----------");
						System.out.println("----------\tsL\t"+serviceLevel[sL]+"\tcv2:\t"+cv2[c]+"\tdelta:\t"+delta[d]);

						double S1 		= modelS1(delta[d], discLimit1, d1, h1, h2, p);
						System.out.println("S1: \t"+S1);

						double lowerBound 	= S1; 				// lower bound on S2, can be set manually as well (replace S1 by any value)
						double upperBound 	= discLimit2; 		// upper bound on S2, can be set manually as well (replace discLimit2 by any value)

						double[] result	= model(delta[d], discLimit2, d1,d2, h1, h2, p, S1, lowerBound, upperBound);	
						System.out.println("S2: \t"+result[0]);

						// save results
						long elapsedTime = System.currentTimeMillis() - startTime;						
						double[] combiResult = {serviceLevel[sL], cv2[c], p, delta[d], S1, result[0], result[1], elapsedTime/1000, discLimit1, discLimit2, alpha, beta}; 
						results[count] = combiResult;
						count++;

						System.out.println("sL\tcv2\tp\td\tS1\tS2\tcost\ttime\tdLim1\tdLim2\talpha\tbeta");
						for(int i=0; i<results.length; i++) {
							for(int j=0; j<results[0].length;j++) {
								System.out.print(df3.format(results[i][j])+"\t");
							}
							System.out.print("\n");
						}
					}
				}
			}
		}

		// printing all results (now in comments since it will already be printed after the last iteration)
/*		System.out.println("sL\tcv2\tp\td\tS1\tS2\tcost\ttime\tdLim1\tdLim2\talpha\tbeta");
		for(int i=0; i<results.length; i++) {
			for(int j=0; j<results[0].length;j++) {
				System.out.print(df3.format(results[i][j])+"\t");
			}
			System.out.print("\n");
		}*/
	}

	public static double[] model(double delta, double discLimit, GammaDistribution d1, GammaDistribution d2, double h1, double h2, double p, double S1, double lowerBound, double upperBound) {
		double[] result = {0,0,0};
		try {
			IloCplex cplex = new IloCplex();

			// variables
			double h1prime = h1+h2;								// calculate local holding cost most downstream stage
			double pPlush1prime = p + h1prime;					// introduce variable that is sum of backorder costs and h1prime
			int M = 10000; 										// arbitrary large number used in calculations
			int limitD = (int) Math.round(discLimit/delta);		// number of steps in discretization

			double lowerboundDelta = lowerBound/delta;
			double upperboundDelta = upperBound/delta;
			double s1delta = S1/delta;

			// initialization of CPLEX variables
			IloNumVar S2 = cplex.numVar(0, M);		
			IloNumVar v[] = new IloNumVar[limitD];
			IloNumVar z[][] = new IloNumVar[limitD][limitD];
			IloNumVar q[][] = new IloNumVar[limitD][limitD];
			IloNumVar r[][] = new IloNumVar[limitD][limitD];
			for(int k = 0; k<(limitD); k++) {
				v[k] = cplex.boolVar("v"+k);
				for(int l= 0; l<(limitD); l++) {
					z[k][l]= cplex.boolVar("z"+k+l);
					q[k][l]= cplex.boolVar("q"+k+l);
					r[k][l] = cplex.intVar(0,M,"r"+k+l);
				}
			}
			IloNumExpr sumVg = cplex.numExpr();
			for(int g=0; g<limitD; g++) {
				sumVg = cplex.sum(sumVg, v[g]);
			}

			IloNumExpr fullexpr = cplex.numExpr(); // objective function

			for(int k=0; k<limitD; k++) {
				// variables used in expressions
				double kDelta = (k+1)*delta;														// k*delta
				double p2k = (delta/2)*(d2.density((k+.5)*delta)+d2.density((k+1.5)*delta));		// p(2,k)
				double kDp2k = kDelta*p2k;															// k*Delta*p(2,k)

				// expression 1
				IloNumExpr expr1 = cplex.prod(h2, cplex.sum(-kDp2k,cplex.prod(p2k, S2)));	

				// add expression 1 to objective function 
				fullexpr = cplex.sum(fullexpr, expr1);

				for(int l=0; l<limitD; l++) {		

					//variables used in expressions
					double lDelta = (l+1)*delta;													// l*delta
					double p1l = (delta/2)*(d1.density((l+.5)*delta)+d1.density((l+1.5)*delta));	// p(1,l)					
					double p1lp2k = p1l*p2k;														// p(1,l)*p(2,k)
					double S1minlDp1lp2k = (S1-lDelta)*p1lp2k;										// (S1-l*delta)*p(1,l)*p(2,k)	

					// expression 2
					if(k<=upperboundDelta-s1delta) {
						IloNumExpr expr2 = cplex.prod(h1, cplex.prod(S1minlDp1lp2k, v[k]));		
						fullexpr = cplex.sum(fullexpr, expr2);	
					}

					// variables used in expressions
					double lDkDp1lp2k = (lDelta+kDelta)*p1lp2k;										// (l*delta+k*delta)*p(1,l)*p(2,k)

					// expression 3
					if(k>=lowerboundDelta-s1delta) {
						IloNumExpr expr3_1= cplex.prod(-lDkDp1lp2k, cplex.sum(1, cplex.prod(-1, v[k])));
						IloNumExpr expr3_2 = cplex.prod(p1lp2k, cplex.prod(S2,cplex.sum(1, cplex.prod(-1, v[k]))));
						IloNumExpr expr3 = cplex.prod(h1, cplex.sum(expr3_1,expr3_2));	
						fullexpr = cplex.sum(fullexpr, expr3);	
					}

					// expression 4
					if(k>=lowerboundDelta-s1delta && l>=lowerboundDelta-k) {
						IloNumExpr expr4_1= cplex.prod(lDkDp1lp2k, z[k][l]);
						IloNumExpr expr4_2 = cplex.prod(-p1lp2k, cplex.prod(S2,z[k][l]));
						IloNumExpr expr4 = cplex.prod(pPlush1prime, cplex.sum(expr4_1,expr4_2));	
						fullexpr = cplex.sum(fullexpr, expr4);	
					}
					// variables used in expressions
					double lDminS1p1lp2k = (lDelta-S1)*p1lp2k;										// (l*delta-S1)*p(1,l)*p(2,k)

					// expression 5
					if(k<=upperboundDelta-s1delta && l>= s1delta ) { 
						IloNumExpr expr5 = cplex.prod(pPlush1prime, cplex.prod(lDminS1p1lp2k, v[k]));					
						fullexpr = cplex.sum(fullexpr, expr5);	
					}				

					/*
					 * Constraints on r(g,k,l)
					 */
					cplex.addLe(r[k][l], cplex.sum(cplex.sum(cplex.prod(1/delta, S2), M), cplex.prod(-M, q[k][l])));		
					cplex.addGe(r[k][l], cplex.sum(cplex.sum(cplex.prod(1/delta, S2), -M), cplex.prod(M, q[k][l])));			
					cplex.addLe(r[k][l], cplex.prod(M, q[k][l]));							
					cplex.addGe(r[k][l], cplex.prod(-M, q[k][l]));	

					/*
					 * Constraints on q(k,l)
					 */
					int kPlusl = (k+1)+(l+1);									// k+l
					cplex.addLe(kPlusl, cplex.sum(cplex.prod(1/delta, S2), cplex.prod(-1, r[k][l]), cplex.prod(M, q[k][l])));	// k+l <= sum(y(g)) - sum(r(g,k,l)) + q(k,l)*M
					int kPluslMinus1 = (k+1)+(l+1)-1;							// k+l-1
					cplex.addGe(kPluslMinus1, r[k][l]);							// k+l > sum(r(g,k,l))  	// here: k+l-1 >= sum(r(g,k,l))

					if(k>0) {
						cplex.addLe(q[k-1][l], q[k][l]);							// q(k-1,l) <= q(k,l)						
					}
					if(l>0) {
						cplex.addLe(q[k][l-1], q[k][l]);							// q(k,l-1) <= q(k,l)
					}

					/*
					 * Constraints on z(k,l)
					 */
					if(k>0) {
						cplex.addLe(z[k-1][l], z[k][l]);							// z(k-1,l) <= z(k,l)
					}
					if(l>0) {
						cplex.addLe(z[k][l-1], z[k][l]);							// z(k,l-1) <= z(k,l)
					}
					cplex.addLe(z[k][l], cplex.sum(1, cplex.prod(-1, v[k])));		// z(k,l) <= 1-v(k)
					cplex.addLe(z[k][l], q[k][l]);									// z(k,l) <= q(k,l)
					cplex.addGe(z[k][l], cplex.sum(q[k][l], cplex.prod(-1,v[k])));	// z(k,l) >= -v(k) + q(k,l)

				}

				/*
				 * Constraints for consecutive 1's on y(g), v(g)
				 */
				if(k>0) {
					//				cplex.addGe(y[k-1], y[k]);						// y(g-1) >= yg
					cplex.addGe(v[k-1], v[k]);										// v(g-1) >= vg
				}
			}			

			/*
			 *  Additional constraints
			 */
			cplex.addEq(sumVg, cplex.sum(cplex.prod(1/delta, S2), -s1delta));		// sum(vg) = sum(yg) - s1*/delta

			// Tell values of z, v, q outside of bounds

			for(int k=0; k<limitD; k++) {
				for(int l=0; l<limitD; l++) {
					if(k+l+2>upperboundDelta) {
						cplex.addEq(q[k][l], 1);
					}
					else if (k+l +2 < lowerboundDelta) {
						cplex.addEq(q[k][l], 0);
					}
					if(k +1< lowerboundDelta -s1delta) {
						cplex.addEq(z[k][l], 0);
					}
				}				
				if(k+1<lowerboundDelta-s1delta) {
					cplex.addEq(v[k], 1);
				}
				else if(k+1> upperboundDelta-s1delta) {
					cplex.addEq(v[k], 0);
				}
			}
			cplex.addMinimize(fullexpr);
			cplex.exportModel ("optimization.lp");		
			
			long startTime = System.currentTimeMillis();

			// solve the model
			if(cplex.solve()) {

				long elapsedTime = System.currentTimeMillis() - startTime;

				// saving results
				double cost 	= cplex.getObjValue();
				double bestS2	= cplex.getValue(S2);
				double time		= elapsedTime/1000;
				result[0] = bestS2;
				result[1] = cost;
				result[2] = time;
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

	public static double modelS1(double delta, double discLimit, GammaDistribution d1, double h1, double h2, double p) {
		double S1result = 0;
		try {
			IloCplex cplex = new IloCplex();

			// variables
			double h1prime = h1+h2;								// calculate local holding cost
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
				double lDeltap1l = p1l*lDelta;

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
					cplex.addGe(y[l-1], y[l], "114");		// y(g-1) >= yg
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

			// solve the model
			if(cplex.solve()) {

				// obtain result:
				double bestS1	= cplex.getValue(S1);
				S1result = bestS1;
			}
			else {
				System.out.println("No solution found");
			}
		}
		catch(IloException exc){
			exc.printStackTrace();
		}	

		return S1result;
	}
}
