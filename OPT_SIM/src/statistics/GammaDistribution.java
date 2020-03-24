/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package statistics;
import java.util.Random;
/**
 *
 * @author mboon
 */
public class GammaDistribution extends Distribution {

    protected double alpha, beta;
    
    public GammaDistribution(double alpha, double beta, Random random) {
        this.alpha = alpha;
        this.beta = beta;
        this.random = random;
    }
    
    @Override
    public double expectation() {
        return alpha/beta;
    }

    @Override
    public double variance() {
        return alpha/beta/beta;
    }
    
    public double nextRandomOld() {
        int k = (int) Math.floor(alpha);
        double a = alpha - k;
        double y = 0;
        double e = Math.E;
        if (a > 0) {
            while (y == 0) {
                double u = random.nextDouble();
                double b = (e + a) / e;
                double p = u * b;
                if (p <= 1) {
                    double x = Math.pow(p, 1/a);
                    double u1 = random.nextDouble();
                    if (u1 <= Math.exp(-x)) { // Accept
                        y = x;
                    }
                }
                else {
                    double x = -Math.log((b-p)/a);
                    double u1 = random.nextDouble();
                    if (u1 <= Math.pow(x, a-1)) { // Accept
                        y = x;
                    }   
                }    
            }
        }
        double product = 1;
        for (int i = 0; i < k; i++) {
            product = product * random.nextDouble();
        }
        return ((y - Math.log(product))/beta);
    }
    
    public double nextRandom()  {
        double shape = alpha;
        double scale = 1.0/beta;
        if (shape < 1) {
            // [1]: p. 228, Algorithm GS

            while (true) {
                // Step 1:
                final double u = random.nextDouble();
                final double bGS = 1 + shape / Math.E;
                final double p = bGS * u;

                if (p <= 1) {
                    // Step 2:

                    final double x = Math.pow(p, 1 / shape);
                    final double u2 = random.nextDouble();

                    if (u2 > Math.exp(-x)) {
                        // Reject
                        continue;
                    } else {
                        return scale * x;
                    }
                } else {
                    // Step 3:

                    final double x = -1 * Math.log((bGS - p) / shape);
                    final double u2 = random.nextDouble();

                    if (u2 > Math.pow(x, shape - 1)) {
                        // Reject
                        continue;
                    } else {
                        return scale * x;
                    }
                }
            }
        }

        // Now shape >= 1

        final double d = shape - 0.333333333333333333;
        final double c = 1 / (3 * Math.sqrt(d));

        while (true) {
            final double x = random.nextGaussian();
            final double v = (1 + c * x) * (1 + c * x) * (1 + c * x);

            if (v <= 0) {
                continue;
            }

            final double x2 = x * x;
            final double u = random.nextDouble();

            // Squeeze
            if (u < 1 - 0.0331 * x2 * x2) {
                return scale * d * v;
            }

            if (Math.log(u) < 0.5 * x2 + d * (1 - v + Math.log(v))) {
                return scale * d * v;
            }
        }
    }


    public static void main(String[] arg) {
        double sumX = 0;
        double sumX2 = 0;
        GammaDistribution dist = new GammaDistribution(0.5, 0.25, new Random());
        int n = 100000;
        for (int i = 0; i < n; i++) {
            double r = dist.nextRandom();
            sumX += r;
            sumX2 += r*r;
        }
        double EX = sumX/n;
        double EX2 = sumX2/n;
        System.out.println("E[X] = "+EX);
        System.out.println("Var[X] = "+(EX2-EX*EX));
    }
    
}
