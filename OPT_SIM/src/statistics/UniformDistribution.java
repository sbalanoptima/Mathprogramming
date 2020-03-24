package statistics;
import java.util.*;

public class UniformDistribution extends Distribution{
    
    protected double a; // lower bound
    protected double b; // upper bound 
    
    public UniformDistribution(int a, int b, Random random){
        this.a = a;
        this.b = b;
        this.random = random;
    }

    @Override
    public double expectation() {
        return (a+b)/2;
    }

    @Override
    public double variance() {
        return (b-a)*(b-a)/12;
    }

    @Override
    public double nextRandom() {
        double U = random.nextDouble();
        return a + U*(b-a);
    }
}