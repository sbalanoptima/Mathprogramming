package statistics;
import java.util.*;

public class PoissonDistribution extends Distribution{

    protected double lambda;
    
    public PoissonDistribution(double lambda, Random random){
        this.lambda = lambda;
        this.random = random;
    }
    
    @Override
    public double expectation() {
        return lambda;
    }

    @Override
    public double variance() {
        return lambda;
    }

    @Override
    public double nextRandom() {
        // Discrete inverse transform method
        double U = random.nextDouble();
        int index = 0;
        double current = Math.exp(-lambda);
        double sum = current;
        while( sum < U ){
            index++;
            current = current*lambda/(1.0*index);
            sum += current;       
        }
        return index*1.0;
    }
    
    
}