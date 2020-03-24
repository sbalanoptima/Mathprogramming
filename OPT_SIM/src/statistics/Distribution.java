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
public abstract class Distribution {

    protected Random random = new Random();
    
    public abstract double expectation();
    
    public abstract double variance();
    
    public double standardDeviation() {
        return Math.sqrt(variance());
    }
    
    public abstract double nextRandom();
}
