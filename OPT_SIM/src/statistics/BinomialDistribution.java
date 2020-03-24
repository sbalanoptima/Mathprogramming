package statistics;
import java.util.Random;


public class BinomialDistribution extends Distribution {

    protected int n;      // Number of experiments
    protected double p;   // Success probability
    
    public BinomialDistribution( int n, double p, Random random ){
        this.n = n;
        this.p = p;
        this.random = random;
    }
    
    
    @Override
    public double expectation() {
        return n*p;
    }

    @Override
    public double variance() {
        return n*p*(1-p);
    }

    @Override
    public double nextRandom() {
        // A binomial rv is the sum of n Bernoulli rv's with parameter p
        BernoulliDistribution bd = new BernoulliDistribution(p,random);
        double sum = 0;
        for( int i = 0; i < n; i++ ){
            sum += bd.nextRandom();
        }
        return sum;
    }
    
    
    
}