

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

import org.apache.commons.math3.distribution.GammaDistribution;

//import ilog.concert.Exception;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

import statistics.*;

public class Optim_Simul_finite_one_capacitated_final3new {

	public static void main(String[] args) {		

		// Random number generator
		Random rng = new Random();
		//rng.setSeed(5);			// with this seed possible to reproduce exact results

		// Mean scenarios
		double[] meanRising					= {9,5,6,7,8,9,10,11,12,13};
		double[] meanDecline 				= {9,13,12,11,10,9,8,7,6,5};
		double[] meanConstant 				= {9,9,9,9,9,9,9,9,9,9,9};
		double[] meanPeak 					= {9,7,7,7,7,10,12,17,7,7};	
		
		// Capacity scenarios
		double[] capacity12				= {12,12,12,12,12,12,12,12,12,12};
		double[] capacity11				= {11,11,11,11,11,11,11,11,11,11};
		double[] capacity10				= {10,10,10,10,10,10,10,10,10,10};
		double[] capacityunlim			= {50,50,50,50,50,50,50,50,50,50};
		double[] capacitydrop			= {13,13,13,13,4,4,13,13,13,13};		

		ArrayList<double[]> allMean = new ArrayList<double[]>();
		ArrayList<double[]> capacity = new ArrayList<double[]>();
		
		/*
		 * Adjustable variables
		 */
		
		double[] serviceLevel 	= {0.95}; // or: 0.8,0.9
		double[] delta 			= {1};
		double[] cv2			= {0.5}; // or: 1,2   // squared coefficient of variation		
		
		/*
		 *  uncomment if you want to add a mean scenario and/or capacity scenario
		 */
		//allMean.add(meanRising);
		//allMean.add(meanDecline);
		//allMean.add(meanConstant);
		allMean.add(meanPeak);
		
		
		capacity.add(capacity10);
		//capacity.add(capacity11);
		//capacity.add(capacity12);
		//capacity.add(capacitydrop);
		//capacity.add(capacityunlim);
		
		/*
		 * Fixed variables
		 */

		double  h1 			= 1; // normalized
		int timeLimit		= 10;

		/*
		
		 * Run Optimization and simulation
		 */
		ArrayList<OptResult> optResults= new ArrayList<OptResult>();
		ArrayList<SimResult> simResults= new ArrayList<SimResult>();
		int ID =0;
		for(int sL = 0; sL<serviceLevel.length; sL++) {	
			for(int s=0; s<allMean.size(); s++) {		
				for(int d = 0; d<delta.length; d++) {
					for(int c = 0; c<cv2.length; c++) {
						for(int capac=0; capac<capacity.size(); capac++) {	
							long startTime = System.currentTimeMillis();		
							double p 		= (-h1*serviceLevel[sL])/(serviceLevel[sL]-1);	// backorder cost 

							double[] mean 	= allMean.get(s);
							double[] varianceD 	=  new double[timeLimit];
							double[] varianceDLT 	=  new double[timeLimit];						
							double[] shapeD		= new double[timeLimit];
							double[] shapeDLT	= new double[timeLimit];
							double[] scaleD		= new double[timeLimit];
							double[] scaleDLT	= new double[timeLimit];			

							// Distributions for simulation
							GammaDistributionReturnValues[] d1 = new GammaDistributionReturnValues[timeLimit];
							for(int t=0; t<timeLimit; t++) {
								varianceD[t]	= cv2[c]*mean[t]*mean[t];
								scaleD[t]		= varianceD[t]/mean[t];
								shapeD[t]		= mean[t]/scaleD[t];
								d1[t]			= new GammaDistributionReturnValues(shapeD[t], 1/scaleD[t], rng);		
							}

							// Distributions for optimization
							GammaDistribution[] d1LT = new GammaDistribution[timeLimit];
							double[] discLimit	= new double[timeLimit];
							for(int t=0; t<timeLimit; t++) {
								if(t==timeLimit-1) {
									varianceDLT[t]	= varianceD[t];
									scaleDLT[t]		= varianceDLT[t]/mean[t];
									shapeDLT[t]		= mean[t]/scaleDLT[t];
								}
								else {
									varianceDLT[t]	= varianceD[t] + varianceD[t+1];
									scaleDLT[t]		= varianceDLT[t]/(mean[t]+mean[t+1]);
									shapeDLT[t]		= (mean[t]+mean[t+1])/scaleDLT[t];
								}

								d1LT[t] = new GammaDistribution(shapeDLT[t], scaleDLT[t]);

								double cumulProb1 = 0;
								discLimit[t] = 1;
								while(cumulProb1<.9998) { //.9998
									cumulProb1 = d1LT[t].cumulativeProbability(discLimit[t]);
									discLimit[t]++;
								}	
							}

							/*
							 * Possible to print outputs by uncommenting
							 */
							/*for(int t=0; t<timeLimit; t++) {
								System.out.print(t+"\t");
							}
							System.out.print("\nMean\t");
							for(int t=0; t<timeLimit; t++) {
								System.out.print(d1[t].expectation()+"\t");
							}
							System.out.print("\nMeaLT\t");
							for(int t=0; t<timeLimit; t++) {
								System.out.print(d1LT[t].getNumericalMean()+"\t");
							}
							System.out.print("\nAlpha\t");
							for(int t=0; t<timeLimit; t++) {
								System.out.print(d1[t].alpha+"\t");
							}
							System.out.print("\nBeta\t");
							for(int t=0; t<timeLimit; t++) {
								System.out.print(d1[t].beta+"\t");
							}
							System.out.print("\nshap\t");
							for(int t=0; t<timeLimit; t++) {
								System.out.print(d1LT[t].getShape()+"\t");
							}
							System.out.print("\nscal\t");
							for(int t=0; t<timeLimit; t++) {
								System.out.print(d1LT[t].getScale()+"\t");
							}
							System.out.print("\ndLim\t");
							for(int t=0; t<timeLimit; t++) {
								System.out.print(discLimit[t]+"\t");
							}
*/
							OptResult optResult = optimizationModel(ID, delta[d], discLimit, d1LT, cv2[c], h1, p, mean, capacity.get(capac), timeLimit);	
							optResults.add(optResult);
							SimResult simResult = simulationModel(optResult, rng, d1);
							simResults.add(simResult);
							ID++;
							optResult.addTotalTime((System.currentTimeMillis()-startTime)/1000);
							printResults(optResult, simResult);   // comment if you do not want the results printed
							writeToFile(optResults, simResults);  // uncomment if you want to write the results to a text file. Make sure to edit file name and location in the class below.
						}
					}
				}
			}
		}
	}

	public static OptResult optimizationModel(int ID, double delta, double[] discLimit, GammaDistribution[] d1LT, double cv2, double h1, double p, double[] mean, double[] capacity, int timeLimit) 
	{

		// variables
		double h1prime = h1;
		int M = 100000; //arbitrary large number
		int[] limitD = new int[timeLimit];
		int maxLimitD = 0;
		for(int t=0; t<timeLimit; t++) {
			limitD[t] = (int) Math.round(discLimit[t]/delta);
			maxLimitD = Math.max(maxLimitD, limitD[t]);
		}

		OptResult result = new OptResult(ID, cv2, h1, p, mean, capacity, delta, discLimit, timeLimit);

		try { 
			IloCplex cplex = new IloCplex();

			IloNumVar S1[] = new IloNumVar[timeLimit];
			IloNumVar y[][] = new IloNumVar[maxLimitD][timeLimit];
			IloNumVar y1[][] = new IloNumVar[maxLimitD][timeLimit];
			IloNumVar E[][][]= new IloNumVar[maxLimitD][maxLimitD][timeLimit];
			IloNumVar X[] = new IloNumVar[timeLimit];
			IloNumVar z[] = new IloNumVar[timeLimit];
			IloNumVar u[] = new IloNumVar[timeLimit];
			IloNumVar v[] = new IloNumVar[timeLimit];
			IloNumVar IN[] = new IloNumVar[timeLimit];
			for(int t = 0; t<timeLimit; t++) {
				X[t] = cplex.numVar(0,M,"X"+t); //non-negativity constraint is covered
				IN[t] = cplex.numVar(-M, M,"IN"+t);
				S1[t] = cplex.numVar(0, M);
				z[t] = cplex.boolVar();
				u[t] = cplex.boolVar();
				v[t] = cplex.boolVar();
				for(int g=0; g<limitD[t]; g++) {
					y[g][t] = cplex.boolVar("y_"+g+"_"+t);
				}
				for(int q=0; q<limitD[t]; q++) {
					y1[q][t] = cplex.boolVar("y1_"+q+"_"+t);
				}
				for(int g=0; g<limitD[t]; g++) {
					for(int q=0; q<limitD[t]; q++) {
						E[g][q][t] = cplex.boolVar("E_"+g+"_"+q+"_"+t);
					}
				}
			}
			IloNumExpr sumYg[] = new IloNumExpr[timeLimit];
			IloNumExpr sumY1g[]= new IloNumExpr[timeLimit];
			for(int t=0; t<timeLimit; t++) {
				sumYg[t] =cplex.numExpr();
				sumY1g[t]=cplex.numExpr();
				for(int g=0; g<limitD[t]; g++) {
					sumYg[t] = cplex.sum(sumYg[t], y[g][t]);
				}
				for(int g=0; g<limitD[t]; g++) {
					sumY1g[t] = cplex.sum(sumY1g[t], y1[g][t]);
				}
			}
			

			// expressions
			IloNumExpr fullexpr = cplex.numExpr();

			for(int t=1; t<timeLimit; t++) {
				for(int l=0; l<limitD[t]; l++) {
					for (int q=0;q<limitD[t];q++) {
					// l=0 is just the index. We actually start from l=1, so whenever the l is used, don't forget to do (l+1) in code

					// creating expressions and adding to objective function
					double p1l = (delta/2)*(d1LT[t].density((l+.5)*delta)+d1LT[t].density((l+1.5)*delta));  					
					double lDp1l = (l+1)*p1l*delta;
					IloNumExpr expr1 = cplex.prod(h1, cplex.sum(-lDp1l, cplex.prod(p1l, S1[t])));

					double pplushprime = p+h1prime; 
				    IloNumExpr expr3 = cplex.prod(pplushprime, cplex.sum(cplex.prod(lDp1l,cplex.sum(1, cplex.prod(-1, y[l][t]))),cplex.prod(p1l,cplex.prod(-1, S1[t]))));
				    IloNumExpr expr4 = cplex.prod(pplushprime, cplex.prod(delta*p1l, E[l][q][t]));
					fullexpr = cplex.sum(fullexpr, expr1, expr3, expr4);
					}
				}				
				if(t>0) {
					cplex.addEq(cplex.sum(IN[t-1], cplex.sum(-mean[t], X[t-1])), IN[t], "con.4");

					// relaxing S
					cplex.addGe(cplex.sum(IN[t-1], cplex.sum(X[t], X[t-1])), S1[t], "con.3");
					cplex.addLe(cplex.sum(IN[t-1], cplex.sum(X[t], X[t-1])), cplex.sum(S1[t],cplex.prod(M, v[t])), "con.3");
				}
				cplex.addLe(X[t], capacity[t], "con1"); 

			}

			// initial conditions on IN and X
			cplex.addEq(IN[0], cplex.sum(S1[1], -9)); 
			cplex.addEq(X[0], 0);
			
			// S[t] = delta * sum of strips equation
			for(int t=0; t<timeLimit; t++) {
				// consecutive 1's on yg and sumyg = S1
				for(int g=1; g<limitD[t]; g++) {
					cplex.addGe(y[g-1][t], y[g][t],"con3");
				}
				cplex.addEq(cplex.prod(delta,sumYg[t]), S1[t]);
					if(t>1) {
					// S1 can be maximum increased with its capacity
					cplex.addLe(S1[t], cplex.sum(S1[t-1], capacity[t]));
				}
			}
			// Newly added constraint to equate the y and y1 variables in the linearization
			for(int t=0; t<timeLimit; t++) {
				for(int q=0; q<limitD[t]; q++) {
					cplex.addEq(y1[q][t], y[q][t],"new_r");
				}
			}
			
			// indicator function for relaxing S
			for(int t=0; t<timeLimit; t++) {
				cplex.addLe(cplex.sum(1, cplex.prod(-1,z[t])), X[t]);
				cplex.addGe(cplex.prod(M, cplex.sum(1, cplex.prod(-1,  z[t]))), X[t]);
				if(t>0) {
					cplex.addGe(cplex.sum(S1[t], cplex.prod(M, u[t])), cplex.sum(S1[t-1], -mean[t-1]));
					cplex.addLe(cplex.sum(0,S1[t]), cplex.sum(S1[t-1], cplex.sum(-mean[t-1],cplex.sum(M, cplex.prod(-M, u[t])))));					
				}	
				cplex.addLe(v[t], z[t]);
				cplex.addLe(v[t], u[t]);
				cplex.addGe(v[t], cplex.sum(z[t], cplex.sum(u[t], -1)));
				for(int l=0; l<limitD[t]; l++) {
					for(int q=0; q<limitD[t]; q++) {
						cplex.addLe(cplex.sum(y[l][t],y1[q][t]),cplex.sum(cplex.prod(2, E[l][q][t]),1) ,"new1");
						cplex.addGe(cplex.sum(y[l][t],cplex.prod(-1, y1[q][t])), cplex.sum(cplex.prod(1, E[l][q][t]),-1),"new2");
						cplex.addGe(cplex.sum(y1[q][t],cplex.prod(-1, y[l][t])), cplex.sum(cplex.prod(1, E[l][q][t]),-1),"new3");
						cplex.addGe(cplex.sum(y[l][t],y1[q][t]),cplex.prod(1, E[l][q][t]),"new4");
					}
					
				}
				
			}			

			cplex.addMinimize(fullexpr);
			cplex.exportModel ("revisedR.lp");		

			long startTime = System.currentTimeMillis();

			// parameters for CPLEX
			cplex.setParam(IloCplex.Param.TimeLimit, 10800); // time limit not required // Given 3 hours maximum run time
			
			// solve the model
			if(cplex.solve()) {
				long elapsedTime = System.currentTimeMillis() - startTime;

				double totalTime = elapsedTime/1000;
				double optCost = cplex.getObjValue();
				double[] S1result = new double[timeLimit];
				double[] Xresult = new double[timeLimit];
				double[] INresult = new double[timeLimit];
				for(int t = 0; t<timeLimit; t++) {
					S1result[t] = cplex.getValue(S1[t]);
					Xresult[t] = cplex.getValue(X[t]);
					INresult[t] = cplex.getValue(IN[t]);					
				}
				result.addResults(S1result, Xresult, INresult, optCost, totalTime);
				
				//return result;
				//System.out.print("optCost\t");
				DecimalFormat df4 = new DecimalFormat("#.####");
				System.out.println("optCost\t=\t"+df4.format(result.optCost));
				
				// possible to print some output values:
				// Print all y[g][t], y1[g][t] and E[g][q][t] values
//				for(int t=0;t<timeLimit;t++)
//				{
//					if (t!=7) {continue;}
//					for(int l=0; l<limitD[t]; l++) 
//					{
//						for(int q=0;q<limitD[t];q++) 
//						{
//						double s1 = cplex.getValue(y[l][t]);
//						double s2 = cplex.getValue(y1[l][t]);
//						double s3 = cplex.getValue(E[l][q][t]);
//					//	if(s1>0.99)							
//						{
//							System.out.println("y["+l+"]["+t+"]= "+ s1+"\n");
//						}
//						//if(s2>0.99)							
//						{
//							System.out.println("y1["+l+"]["+t+"]= "+ s2+"\n");
//						}
//						{
//							System.out.println("E["+l+"]["+q+"]["+t+"]= "+ s3+"\n");
//						}
//					}
//				}
//				}
				
				/*
				System.out.println("Elapsed time (s): "+elapsedTime/1000);
				System.out.println("STATUS: "+cplex.getStatus());

				System.out.println("-----------------------------------------------------------------------------------------------------------------------");

				System.out.print("t\t");
				for(int t=0; t<timeLimit; t++) {
					System.out.print(t+"\t");			
				}
				System.out.print("\n");

				System.out.print("S1\t");
				for(int t=0; t<timeLimit; t++) {
					System.out.print((int) Math.round(cplex.getValue(S1[t]))+"\t");			
				}
				System.out.print("\n");

				System.out.print("Mean\t");
				for(int t=0; t<timeLimit; t++) {
					System.out.print((int) mean[t]+"\t");			
				}
				System.out.print("\n");

				System.out.print("demLT\t");
				for(int t=0; t<timeLimit; t++) {
					System.out.print((int) demLT[t]+"\t");			
				}
				System.out.print("\n");

				System.out.print("Beta\t");
				for(int t=0; t<timeLimit; t++) {
					System.out.print((int) beta[t]+"\t");			
				}
				System.out.print("\n");

				System.out.print("Capac\t");
				for(int t=0; t<timeLimit; t++) {
					System.out.print((int) capacity[t]+"\t");			
				}
				System.out.print("\n");

				System.out.print("Xt\t");
				for(int t=0; t<timeLimit; t++) {
					int solXt = (int) Math.round(cplex.getValue(X[t]));		
					System.out.print(solXt+"\t");			
				}
				System.out.print("\n");

				System.out.print("IN\t");
				for(int t=0; t<timeLimit; t++) {
					System.out.print((int) Math.round(cplex.getValue(IN[t]))+"\t");			
				}
				System.out.print("\n");


				// printing objective value:
				System.out.println("\nTotal cost:\t"+ cost);
				System.out.println("\nAvg cost:\t"+ cost/timeLimit);
				System.out.println("Comp. time\t"+elapsedTime/1000);


				System.out.println();
				System.out.print("double[] alpha\t\t= {");
				for(int i=0; i<alpha.length; i++) {
					if(i== alpha.length-1) {
						System.out.print(alpha[i]);					
					}
					else {
						System.out.print(alpha[i]+", ");					
					}
				}
				System.out.print("};\n");
				System.out.print("double[] beta\t\t= {");
				for(int i=0; i<beta.length; i++) {
					if(i== beta.length-1) {
						System.out.print(beta[i]);					
					}
					else {
						System.out.print(beta[i]+", ");					
					}
				}
				System.out.print("};\n");
				System.out.print("double[] mean\t\t= {");
				for(int i=0; i<mean.length; i++) {
					if(i== mean.length-1) {
						System.out.print(mean[i]);					
					}
					else {
						System.out.print(mean[i]+", ");					
					}
				}
				System.out.print("};\n");
				System.out.print("double[] demLT\t\t= {");
				for(int i=0; i<demLT.length; i++) {
					if(i== demLT.length-1) {
						System.out.print(demLT[i]);					
					}
					else {
						System.out.print(demLT[i]+", ");					
					}
				}
				System.out.print("};\n");
				System.out.print("double[] capacity\t= {");
				for(int i=0; i<capacity.length; i++) {
					if(i== capacity.length-1) {
						System.out.print(capacity[i]);					
					}
					else {
						System.out.print(capacity[i]+", ");					
					}
				}
				System.out.print("};\n");
				System.out.print("double[] S1\t\t\t= {");
				for(int i=0; i<timeLimit; i++) {
					if(i== timeLimit-1) {
						System.out.print((int) Math.round(cplex.getValue(S1[i])));					
					}
					else {
						System.out.print((int) Math.round(cplex.getValue(S1[i]))+", ");					
					}
				}
				System.out.print("};\n");
				System.out.print("double[] IN\t\t\t= {");
				for(int i=0; i<timeLimit; i++) {
					if(i== timeLimit-1) {
						System.out.print((int) Math.round(cplex.getValue(IN[i])));					
					}
					else {
						System.out.print((int) Math.round(cplex.getValue(IN[i]))+", ");					
					}
				}
				System.out.print("};\n");
				System.out.print("double[] X\t\t\t= {");
				for(int i=0; i<timeLimit; i++) {
					if(i== timeLimit-1) {
						System.out.print((int) Math.round(cplex.getValue(X[i])));					
					}
					else {
						System.out.print((int) Math.round(cplex.getValue(X[i]))+", ");					
					}
				}
				System.out.print("};\n");
				*/ 
			}
			else {

				System.out.println("No solution found");
			}
		}
		catch(Exception exc){
			exc.printStackTrace();
		}
		return result;
	}

	public static SimResult simulationModel(OptResult optResult, Random rng, GammaDistributionReturnValues[] d1)  {

		int ID					= optResult.ID;
		double h1				= optResult.h1;
		double p				= optResult.p;
		double[] capacity		= optResult.capacity;
		int timeLimit			= optResult.timeLimit;
		double[] S1result		= optResult.S1result;
		double[] Xresult		= optResult.Xresult;
		double[] INresult		= optResult.INresult;

		int runs = 50;
		int repetitions = 10000;

		SimResult result = new SimResult(ID, h1, p, capacity, S1result, Xresult, INresult);

		double sumCost = 0;
		double sum2Cost = 0;
		for(int r=0; r<runs; r++) {
			double[] randomDemand 	= new double[timeLimit];
			double[] INsim			= new double[timeLimit];
			double[] Xsim			= new double[timeLimit];
			double[] costs			= new double[timeLimit];
			INsim[0] 				= INresult[0];
			Xsim[0] 				= Xresult[0];
			Xsim[1] 				= Xresult[1];

			double sumCostRun		= 0;			
			for(int rp = 0; rp<repetitions; rp++) {
				for(int t=1; t<timeLimit; t++) {
					randomDemand[t] = d1[t].nextRandom();	
					INsim[t] = INsim[t-1] + Xsim[t-1] - randomDemand[t];				
					costs[t] = Math.max(0, INsim[t])*h1 - Math.min(INsim[t], 0)*p;
					sumCostRun += costs[t];
					if(t<timeLimit-1) {
						Xsim[t+1] = Math.max(0, Math.min(capacity[t+1], S1result[t+1]-INsim[t]-Xsim[t]));
					}
				}					
			}
			double meanCostRun		= sumCostRun/repetitions;
			sumCost += meanCostRun;
			sum2Cost += meanCostRun*meanCostRun;
		}

		double meanCost = sumCost/runs;
		double varCost	= sum2Cost/runs - meanCost*meanCost;
		double halfwidthCost	= 1.96*Math.sqrt(varCost/runs);
		double[] CICost			= {meanCost-halfwidthCost, meanCost+halfwidthCost};
		result.addResults(meanCost, CICost);
		return result;
	}

	public static void printResults(OptResult optResult, SimResult simResult) {		

		Locale.setDefault(new Locale("en", "US"));
		DecimalFormat df0 = new DecimalFormat("#");
		DecimalFormat df1 = new DecimalFormat("#.#");
		DecimalFormat df3 = new DecimalFormat("#.###");

		int timeLimit = optResult.timeLimit;
		System.out.println("-----------------------------------------------------------------------------------------------------------------------");
		System.out.println("ID\t=\t"+df3.format(optResult.ID));
		System.out.println("cv2\t=\t"+df3.format(optResult.cv2));
		System.out.println("h1\t=\t"+df3.format(optResult.h1));
		System.out.println("p\t=\t"+df3.format(optResult.p));
		System.out.println("delta\t=\t"+df3.format(optResult.delta));


		System.out.print("t\t");
		for(int t=0; t<timeLimit; t++) {
			System.out.print(t+"\t");			
		}
		System.out.print("\n");

		System.out.print("S1\t");
		for(int t=0; t<timeLimit; t++) {
			System.out.print(df1.format(optResult.S1result[t])+"\t");			
		}
		System.out.print("\n");

		System.out.print("Mean\t");
		for(int t=0; t<timeLimit; t++) {
			System.out.print(df1.format(optResult.mean[t])+"\t");			
		}
		System.out.print("\n");

		System.out.print("Capac\t");
		for(int t=0; t<timeLimit; t++) {
			System.out.print(df1.format(optResult.capacity[t])+"\t");			
		}
		System.out.print("\n");

		System.out.print("Xt\t");
		for(int t=0; t<timeLimit; t++) {
			System.out.print(df1.format(optResult.Xresult[t])+"\t");			
		}
		System.out.print("\n");

		System.out.print("IN\t");
		for(int t=0; t<timeLimit; t++) {
			System.out.print(df1.format(optResult.INresult[t])+"\t");			
		}
		System.out.print("\n");

		System.out.print("dLimit\t");
		for(int t=0; t<timeLimit; t++) {
			System.out.print(df1.format(optResult.discLimit[t])+"\t");			
		}
		System.out.print("\n");

		System.out.println("Simulation results:");
		System.out.println("Costs:\t"+df3.format(simResult.CICost[0])+"\t"+df3.format(simResult.meanCost)+"\t"+df3.format(simResult.CICost[1]));

	}

	public static void writeToFile(ArrayList<OptResult> optResults, ArrayList<SimResult> simResults) {
		Locale.setDefault(new Locale("en", "US"));
		DecimalFormat df0 = new DecimalFormat("#");
		DecimalFormat df1 = new DecimalFormat("#.#");
		DecimalFormat df3 = new DecimalFormat("#.###");

		String destination="C:\\Users\\HP\\Desktop\\WSC19\\simulation_output";	 // insert your document location here!

		File file = new File(destination + "/opt_and_sim1120.txt");

		int cases = optResults.size();
		try {
			BufferedWriter pw = new BufferedWriter(new FileWriter(file)); 

			for(int c=0; c<cases; c++) {
				OptResult optResult = optResults.get(c);
				SimResult simResult = simResults.get(c);
				pw.newLine();
				pw.write("-------------------------------------------------------------------");pw.newLine();
				pw.write("ID:\t"+optResult.ID);pw.newLine();
				pw.write("cv2:\t"+df1.format(optResult.cv2));pw.newLine();
				pw.write("h1:\t"+df1.format(optResult.h1));pw.newLine();
				pw.write("p:\t"+df1.format(optResult.p));pw.newLine();
				pw.write("delta\t"+df1.format(optResult.delta));pw.newLine();
				pw.write("compT\t"+df1.format(optResult.compTime));pw.newLine();
				pw.write("totalT\t"+df1.format(optResult.totalTime));pw.newLine();
				pw.write("optCost\t"+df3.format(optResult.optCost));pw.newLine();
				pw.write("simCost\t"+df3.format(simResult.meanCost));pw.newLine();
				pw.write("simCI\t"+df3.format(simResult.CICost[0])+"\t"+df3.format(simResult.CICost[1]));pw.newLine();

				pw.write("t\t");
				for(int t=0; t<optResult.timeLimit; t++) {
					pw.write(t+"\t");
				}
				pw.newLine();	
				pw.write("S1\t");
				for(int t=0; t<optResult.timeLimit; t++) {
					pw.write(df1.format(optResult.S1result[t])+"\t");
				}
				pw.newLine();	
				pw.write("mean\t");
				for(int t=0; t<optResult.timeLimit; t++) {
					pw.write(df1.format(optResult.mean[t])+"\t");
				}
				pw.newLine();
				pw.write("capac\t");
				for(int t=0; t<optResult.timeLimit; t++) {
					pw.write(df1.format(optResult.capacity[t])+"\t");
				}
				pw.newLine();
				pw.write("X\t");
				for(int t=0; t<optResult.timeLimit; t++) {
					pw.write(df1.format(optResult.Xresult[t])+"\t");
				}
				pw.newLine();
				pw.write("IN\t");
				for(int t=0; t<optResult.timeLimit; t++) {
					pw.write(df1.format(optResult.INresult[t])+"\t");
				}
				pw.newLine();
				pw.write("dLim\t");
				for(int t=0; t<optResult.timeLimit; t++) {
					pw.write(df1.format(optResult.discLimit[t])+"\t");
				}
				pw.newLine();
			}					


			pw.flush();
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		} 
	}

}


