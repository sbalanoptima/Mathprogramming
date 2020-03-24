package statistics;
import java.util.Random;


public class BernoulliDistribution extends Distribution {

    protected double p;  // Success probability
    
    public BernoulliDistribution( double p, Random random ){
        this.p = p;
        this.random = random;
    }
    
    
    @Override
    public double expectation() {
        return p;
    }

    @Override
    public double variance() {
        return p*(1-p);
    }

    @Override
    public double nextRandom() {
        double U = random.nextDouble();
        if ( U > p ){
            return 0;
        }else{
            return 1;
        }
    }
    
    
    
}