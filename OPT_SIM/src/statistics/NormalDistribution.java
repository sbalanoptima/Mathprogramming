package statistics;
import java.util.*;


public class NormalDistribution extends Distribution{
    
    protected double mu;
    protected double sigma;
    
    public NormalDistribution(double mu, double sigma, Random random) {
        this.mu= mu;
        this.sigma = sigma;
        this.random = random;       
    }

    @Override
    public double expectation() {
        return mu;
    }

    @Override
    public double variance() {
        return sigma*sigma;
    }

    @Override
    public double nextRandom() {
        double U = random.nextGaussian();
        return mu + sigma*U; 
    }
    
    
    
    
}