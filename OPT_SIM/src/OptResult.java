




public class OptResult {
	protected int ID;
	protected double cv2;
	protected double h1;
	protected double p;
	protected double[] mean;
	protected double[] capacity;
	protected double delta;
	protected double[] discLimit;
	protected int timeLimit;
	
	protected double[] S1result;
	protected double[] Xresult;
	protected double[] INresult;
	protected double optCost;
	protected double compTime;
	
	protected double totalTime;
	public OptResult(int ID, double cv2, double h1, double p, double[] mean, double[] capacity, double delta, double[] discLimit, int timeLimit) {
		this.ID				= ID;
		this.cv2			= cv2;
		this.h1 			= h1;
		this.p				= p;
		this.capacity 		= capacity;
		this.delta 			= delta;
		this.discLimit 		= discLimit;
		this.timeLimit 		= timeLimit;
		this.mean 			= mean;
	}

	public void addResults(double[] S1result, double[] Xresult, double[] INresult, double optCost, double compTime) {
		this.S1result 		= S1result;
		this.Xresult 		= Xresult;
		this.INresult 		= INresult;
		this.optCost		= optCost;
		this.compTime 		= compTime;
	}
	public void addTotalTime(double time) {
		this.totalTime		= time;
	}
}
