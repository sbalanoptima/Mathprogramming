package statistics;
import java.util.Random;


public class GeometricDistribution extends Distribution {

    protected double p;   // Success probability
    
    public GeometricDistribution( double p, Random random ){
        this.p = p;
        this.random = random;
    }
    
        @Override
    public double expectation() {
        return 1/p;
    }
    
    @Override
    public double variance() {
        return (1-p)/(p*p);
    }

    @Override
    public double nextRandom() {
        double U = random.nextDouble();      
        return 1 + Math.floor( Math.log(U) / Math.log(1-p));
    }

    
    
    
}