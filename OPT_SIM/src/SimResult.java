




public class SimResult {
	protected double ID;
	protected double h1;
	protected double p;
	protected double[] capacity;
	protected double[] S1result;
	protected double[] Xresult;
	protected double[] INresult;
	
	protected double meanCost;
	protected double[] CICost;
	public SimResult(int ID, double h1, double p, double[] capacity, double[] S1result, double[] Xresult, double[] INresult) {
		this.ID				= ID;
		this.h1 			= h1;
		this.p 				= p;
		this.capacity 		= capacity;
		this.S1result		= S1result;
		this.Xresult		= Xresult;
		this.INresult		= INresult;
	}

	public void addResults(double meanCost, double[] CICost) {
		this.meanCost	= meanCost;
		this.CICost		= CICost;
	}
}
