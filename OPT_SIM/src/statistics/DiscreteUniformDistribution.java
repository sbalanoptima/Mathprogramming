package statistics;
import java.util.*;

public class DiscreteUniformDistribution extends Distribution {
    
    protected int m; // lower bound
    protected int n; // upper bound
    
    public DiscreteUniformDistribution(int m, int n, Random random){
        this.m = m;
        this.n = n;
        this.random = random;
    }

    @Override
    public double expectation() {
        return 1.0*(m+n)/2;
    }

    @Override
    public double variance() {
        return 1.0*((n-m+1)*(n-m+1)-1)/12;
    }

    @Override
    public double nextRandom() {
        double U = random.nextDouble();
        return m + Math.floor((n-m+1)*U);
    }
    
    
    
}