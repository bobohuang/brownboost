package jboost.booster;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import jboost.booster.AdaBoost.BinaryBag;
import jboost.controller.Configuration;
import jboost.examples.Label;
import jboost.NotSupportedException;


/**
 * An implementation of BrownBoost.
 * 
 * @author Aaron Arvey
 */
public class BrownBoost extends AdaBoost {
    /** The total running time of the BrownBoost game. */
    protected double m_c;

    /** In case we ever want to make BrownBoost cost sensitive.  Also useful for Yaba. */
    protected double m_posc;
    protected double m_negc;

    /** The time remaining in the BrownBoost game. */
    protected double m_s;

    /** The potential loss of each example at the begining of the game. */
    protected double m_initialPotential;
  
    /** The last value of alpha */
    protected double m_lastAlpha;


    /**
     * This is true if and only if we are cost sensitive boosting.  If
     * this is set, then m_{pos,neg}c{1,2} and m_{pos,neg}Theta must  
     * be set in YabaBoost and m_c{pos,neg} must be set in BrownBoost.
     */
    protected boolean m_isCostSensitive;

  
    /** 
     * The time remaining in the BrownBoost game on the previous iteration. 
     * This value is useful for catching mistakes.
     */
    protected double m_oldS;
  
    /** 
     * The predictions from a hypothesis for an iteration. 
     * HACK: This should eventually be removed for efficiency concerns. 
     */
    protected double[] m_hypPredictions;
  
    /** 
     * Records the potentials.  Similar to m_margins and m_weights.
     */
    protected double m_totalPotential;
  
    /** 
     * Records the potentials.  Similar to m_margins and m_weights.
     */
    protected int[] m_examples;
  
    /**
     * Default constructor just calls AdaBoost to 
     * get everything initialized
     */
    public BrownBoost(){
	super(0.0);
    }


    public void finalizeData() {
	super.finalizeData();
	m_examples = new int[m_numExamples];
	for (int i=0; i<m_numExamples; i++) {
	    m_potentials[i] = calculatePotential(0,m_c);
	    m_examples[i] = i;
	}
    }


    public void setRuntime(double runtime){
	m_c = runtime;
	m_s = m_c;
	m_initialPotential = (1-erf(Math.sqrt(m_c)))/2;
	//System.out.println("BrownBoost:\n \t m_c: " + m_c + "\n \t intial potential: " + m_initialPotential);	  
    }


    public double getInitialPotential(){
	return m_initialPotential;
    }


    public boolean isFinished(){
	return m_s <= 0;
    }

    /**
     * Checked results with Matlab, and correctness is sufficient.
     * There may be issues when z is very close to zero.  This happens
     * in yababoost.
     *
     * Reference: http://www.cs.princeton.edu/introcs/21function/MyMath.java.html.
     *
     * @param z a value greater than 0 (should be greater than some epsilon)
     * @return the erf (1 / sqrt(pi) int e^(-x^2)
     */
    public static double erf(double z) {
	double t = 1.0 / (1.0 + 0.5 * Math.abs(z));

	// use Horner's method
	double ans = 1 - t * Math.exp( -z*z   -   1.26551223 +
				       t * ( 1.00002368 +
				       t * ( 0.37409196 + 
				       t * ( 0.09678418 + 
				       t * (-0.18628806 + 
				       t * ( 0.27886807 + 
				       t * (-1.13520398 + 
				       t * ( 1.48851587 + 
				       t * (-0.82215223 + 
				       t * ( 0.17087277))))))))));
	
	return sign(z)*Math.abs(ans);
    }

    /**
     * @param z
     * @return If z is negative return -1 else return 1 
     */
    public static double sign(double z){
	if( Double.compare(z, 0.0) == 0 ){
	    return 1.0;
	} else if( Double.compare(z, -0.0) == 0 ){
	    return -1.0;
	}
    	  
	if(z > 0){
	    return 1.0;
	} else {
	    return -1.0;
	}
    }
  
    public double getStep(short simple_label, double hyp_pred) {
	double step = getLabel(simple_label)*hyp_pred;
	double EPS = 0.000001;
	if(Math.abs(step) < EPS) return 0.0;
	return sign(step);	  
    }
  
    public double getLabel(short simple_label) {
	return sign(-simple_label+0.5);
    }

  
    /**
     * ErfVars is a container for the variables associated with the
     * BrownBoost constraint solver.  All variable names are 
     * based on the BrownBoost paper (Freund 1999, 2001).  
     */
    class ErfVars{
	/** Sum of weights (partial of potential constraint wrt t) */
	public double W;

	/** Partial derivative of correlation constraint wrt t */
	public double U;
	  
	/** 
	 * The correlation between hypothesis and labels (partial of 
	 * potential constraint wrt alpha) 
	 */
	public double B;
	  
	/** 
	 * Partial derivative of correlation constraint 
	 * wrt alpha 
	 */
	public double V;
	  
	/** Average difference in current potential and last iteration potential */
	public double E; 

	/** Average difference in current potential and initial potential */
	public double I; 
	  
	/** The summed potential */
	public double Potential;

	/** Default constructor initializes variables to 0. */
	protected ErfVars(){
	    clear();
	}
	  
	/** Initializes variables to given values.  */
	protected ErfVars(double w, double u, double b, double v, double e, double Pot){
	    W = w; U = u; B = b; V = v; E = e; I=0; Potential = Pot;
	}
	  
	/** Resets all values to 0.  */
	protected void clear(){
	    W = 0;
	    U = 0;
	    B = 0;
	    V = 0;
	    E = 0;
	    I = 0;
	    Potential = 0;
	}
    }

    /**
     * Implements the calculation of the constraints and partial
     * derivatives.  To better understand notation, see Freund 1999,
     * Freund 2001, and Freund & Arvey 2008.
     * @returns ErfVars which contains the values of the various parameters.
     */
    protected ErfVars calc_constraints(double alpha, double t)
    {
	ErfVars vars = new ErfVars();		
	

	double new_margin, new_time_remaining, new_weight, new_potential, orig_potential;
	for (int example= 0; example < m_hypPredictions.length; example++) {
	    double margin = m_margins[example];				
	    
	    // m_labels are 0 or 1 or arbitrary other, this moves it to +-1
	    double step  = getStep(m_labels[example], m_hypPredictions[example]);

	    /*
	    // variable names are taken from Freund 1999, Freund 2001
	    double aj = margin + m_s;
	    double bj = step;				
	    double dj = margin + bj*alpha + m_s - t;
	    double wj = Math.exp(-(dj*dj)/m_c);

	    vars.W += wj;
	    vars.U += wj*dj*bj;
	    vars.B += wj*bj;
	    vars.V += wj*dj*bj*bj;
	    vars.I += erf(dj / Math.sqrt(m_c)) - m_initialPotential;
	    vars.E += erf(dj / Math.sqrt(m_c)) - erf(aj / Math.sqrt(m_c));
	    vars.Potential += (1-erf( dj /Math.sqrt(m_c)))/2;
	    */

	    new_margin = margin + step*alpha;
	    new_time_remaining = m_s - t;
	    new_weight = calculateWeight(new_margin, new_time_remaining);
	    new_potential = calculatePotential(new_margin, new_time_remaining);
	    orig_potential = calculatePotential(margin, m_s);
	    vars.E += orig_potential - new_potential; // - orig_potential;
	    vars.B += new_weight * step;
	    vars.Potential += calculatePotential(new_margin, new_time_remaining);

	    //System.out.println("Example: " + j + ", Step: " + step + ", Label: " 
	    //+ label + ", Margin: " + margin + ", Pred: " + pred);
	    //System.out.println("W:" + W + ", U:" + U + ", B:" + B + ", V:" + V + ", E:" + E);
	}		
	
        vars.Potential /= m_hypPredictions.length;
	return vars;
    }



  
    protected double solve_constraints(double hyp_err, int[] examples)
    {   
	if( m_s < 0.001){
	    m_s = -1;
	    return 0;
	}
	  
	  
	double alpha = Math.min(0.1, hyp_err);
	double t = alpha*alpha/3;
	double c = m_c;

	ErfVars vars = new ErfVars();

	double new_alpha = 0.0;
	double new_t = 0.0;

	// try binary search
	// find the maximal value for t for which there exists a value
	// of alpha such that the erf potential does not decrease by much.
	    
	//System.out.println("try binary search");


        double t_step = 0.1;
        t=0.3;
        alpha=0;

        int NUM_ITERATIONS_FINISH_GAME = 10;
        int count_t_over_s = 0;
        double lastE = 0;
        double STEP_EPS = 0.0001;
	double CORR_EPS = 0.001;
        boolean first_iter = true;
        while(Math.abs(t_step) > STEP_EPS) {
            if(Math.abs(t_step) > m_s){
                t_step = m_s/2 * sign(t_step);
            }

            t+=t_step;

	    if (t >= m_s) {
		t = m_s - STEP_EPS;
		t_step = -t_step;
		count_t_over_s++;
		// if we keep going over m_s, the game is probably done
		if(count_t_over_s > NUM_ITERATIONS_FINISH_GAME){
		    t = m_s + 0.001;  // m_s is updated after loop
		    break;
		}
		continue;
	    }

            if (t < STEP_EPS) {
                t = STEP_EPS;
                t_step = -t_step;
                continue;
            }

            alpha = 0;
            //System.out.print("\nt: " + t);
            double alpha_step = 0.1;
            while(Math.abs(alpha_step) > STEP_EPS) {

                alpha += alpha_step;
                //System.out.print("   alpha: " + alpha);


                if(alpha < 0) {
                    alpha = STEP_EPS;
                    alpha_step = -alpha_step/2;
                    continue;
                }

                // calculate constraints for values of alpha and t
                vars = calc_constraints (alpha, t);
		
                // reverse alpha search direction
                if(sign(vars.B) != sign(alpha_step)) alpha_step /= -2;
		if(Math.abs(vars.B) < CORR_EPS) break;
            }

	    /*
	      System.out.println("(alpha:" + alpha + ", t:" + t + ") is, gamma + pot_diff = |" 
			       + vars.B + "| + |" + vars.E + "| = " 
			       + (Math.abs(vars.B) + Math.abs(vars.E)));
	    */

            // reverse t search direction
            if(sign(vars.E) != sign(t_step)) t_step /= -2;
        }

        // The bisection (binary search) alpha and t
        double bs_alpha = alpha;
        double bs_t = t;

	alpha = bs_alpha;
	t = bs_t;

	m_oldS = m_s;
	m_s -= t;

	System.out.format("\tBrownBoost: alpha=%.4f, t=%.4f, time left=%.4f, " 
			  + "potential=%.4f\n", alpha, t, m_s, vars.Potential);

	if(t<0) {
	    System.err.println("\nERROR: The value of t: " + t);
	    System.err.println("ERROR: Bad solution for t<0");
	    m_s = m_oldS;
	    return(0.0);
	} 
	    
	//System.out.println("\ns: " + m_s);
	return alpha;
    }

    public String getParamString() {
	String ret = String.format("BrownBoost r=%.4f s=%.4f ", m_c, m_s);
	return ret;
    }
  
    /**
     *  Update the examples, m_margins, and m_weights using the 
     *  brown boost update rule.  When m_s reaches
     *  the stopping criterion of m_s &lt; 0, this update returns 
     *  immediately and does not actually do any further updating.  
     *  @param predictions values for examples
     *  @param exampleIndex the list of examples to update
     */
  
    public void update(Prediction[] predictions, int[][] exampleIndex) {
	if (m_s < 0){
	    return;
	}
    
	for (int i= 0; i < exampleIndex.length; i++) {
	    double p = predictions[i].getClassScores()[1];
	    double[] value = new double[] { -p, p };
	    int[] indexes = exampleIndex[i];
	    for (int j= 0; j < indexes.length; j++) {
		int example   = indexes[j];
		m_oldMargins[example] = m_margins[example];
		m_margins[example] += value[m_labels[example]];
	    }
	}

	
	m_totalWeight = 0;
	m_totalPotential = 0;
        for (int i=0; i < m_hypPredictions.length; i++) {
	    m_oldWeights[i]= m_weights[i];

            m_weights[i] = calculateWeight(m_margins[i]);
            m_totalWeight += m_weights[i];

            m_potentials[i] = calculatePotential(m_margins[i]);
            m_totalPotential += m_potentials[i];
        }
    }

    /**
     * The theoretical bound on error is not defined for BrownBoost, thus
     * getTheoryBound is undefined.
     */
    public double getTheoryBound() {
	return -1.0;
    }
  
    /**
     * Calls calculateWeight(margin,m_s)
     */
    public double calculateWeight(double margin) {
  	return calculateWeight(margin, m_s);
    }

    /**
     * BrownBoost uses (1-erf(-(margin+s)/c))/2  as the potential function
     */
    public double calculatePotential(double margin) {
	return calculatePotential(margin, m_s);
    }

    /**
     * BrownBoost uses e^(-(margin+s)^2/c) as the weight calculation
     */
    public double calculateWeight(double margin, double time_remaining) {
	double s = time_remaining;
  	return Math.exp(-1 * Math.pow(margin+s,2)/m_c);
    }

    /**
     * BrownBoost uses (1-erf(-(margin+s)/c))/2  as the potential function
     */
    public double calculatePotential(double margin, double time_remaining) {
	double s = time_remaining;
  	return (1-erf((margin+s)/Math.sqrt(m_c)))/2;
    }

  
    /**
     * BrownBag is identical to BinaryBag, except for the method used to
     * derive the value of prediction (alpha in the literature).  BrownBag
     * uses the value of alpha determined by BrownBoost and its variants.
     * See comments for AdaBoost.BinaryBag.
     * @author Aaron Arvey
     */
    class BrownBag extends AdaBoost.BinaryBag{

	protected BrownBag(int[] list) {
	    m_w= new double[2];
	    reset();
	    this.addExampleList(list);
	}

	/** compute the binary prediction associated with this bag */
	public BinaryPrediction calcPrediction(double alpha) {
	    BinaryPrediction ret;
	    ret = new BinaryPrediction(m_w[1] > m_w[0] ? 1.0 : -1.0 );
	    ret.scale(alpha);
	    return ret;
	}
      
	/** compute the binary prediction associated with this bag */
	public BinaryPrediction calcPrediction(double posAlpha, double negAlpha) {
	    BinaryPrediction ret;
	    if (m_w[1] > m_w[0]) {
		ret = new BinaryPrediction(1.0);
		ret.scale(posAlpha);
	    } else {
		ret = new BinaryPrediction(-1.0);
		ret.scale(negAlpha);
	    } 
	    return ret;
	}
      
      
	/** Place holder to ensure that this function is not used in BrownBoost. */
	public BinaryPrediction calcPrediction(){
	    //System.err.println("Need to have alpha for prediction in BrownBag.calcPrediction()!");
	    return new BinaryPrediction(0);
	}
      
      
	/** default constructor */
	protected BrownBag() {
	    super();
	}
      
      
	/** constructor that copies an existing bag */
	protected BrownBag(BrownBag bag) {
	    super(bag);
	}
      
      
	/** Output the weights in the bag */
	public String toString() {
	    String s= "BrownBag.\t w0=" + m_w[0] + "\t w1=" + m_w[1] + "\n";
	    return s;
	}
      
	
    } /* End BrownBag */
    
    
    
    protected double getHypErr(Bag[] bags, int[][] exampleIndex) {
	double hyp_err = 0.0;
	double gamma = 0.0;
	double num_wrong = 0.0;
	double total_weight = 0.0;
	double potential = 0.0;
	int num_predictions = 0;
	int total = 0;

	// Keep track of which hypotheses had hypotheses associated with them.
	boolean[] examplesWithHyp = new boolean[m_margins.length];
	
	
	// Get all examples that have a hypothesis associated with them
	for (int i= 0; i < exampleIndex.length; i++) {
	    int[] indexes= exampleIndex[i];
	    BinaryPrediction pred = ((BrownBag) bags[i]).calcPrediction(1.0);
	    total += 1;
	    for (int j= 0; j < indexes.length; j++) {
		int example   = indexes[j];
		examplesWithHyp[example] = true;
		double step   = getStep(m_labels[example], m_hypPredictions[example]);

		double weight = m_weights[example];
		total_weight += weight;
		gamma += weight*step;
			
		if (step > 0){ // we got it right!
		    num_predictions += 1;
		} else { // We got it wrong
		    hyp_err += 1;
		    num_predictions += 1;
		}
	    }
	}
	    
	// Get all examples that have no hypothesis associated with them.
	// Also get current potential.
	for (int i=0; i < m_margins.length; i++) {
	    if(!examplesWithHyp[i]){
		int example   = i;
		m_hypPredictions[example] = 0;
		double weight = m_weights[example];
		total_weight += weight;
		//System.out.println("m_hypPredictions[" + i + "," + example + "]: " + 0 + " (No hyp for example " + example + ")");
	    }
	    potential += calculatePotential(m_margins[i]);
	}
	
	//updatePotential(exampleIndex);
	hyp_err /= num_predictions;
	gamma /= total_weight;
	potential /= m_margins.length;

	/*
	System.out.println("\tTotal number of examples: " + m_margins.length);
	System.out.println("\tNumber of predictions made: " + num_predictions);
	String out = "\tgamma (weighted correlation):" + gamma + ", potential (unweighted):" + potential + ", hyp error (unweighted):" + hyp_err;
	System.out.println(out);
	*/

	return gamma;
    }

    protected BinaryPrediction getZeroPred() {
	return new BinaryPrediction(0);
    }

    /*
     * Returns the predictions associated with a list of bags representing a
     * partition of the data.
     */
    public Prediction[] getPredictions(Bag[] bags, int[][] exampleIndex) {
	boolean bagsAreWeightless = true;
	for (int i=0; i < bags.length; i++) {
	    if (!bags[i].isWeightless()) {
		bagsAreWeightless = false;
	    }
	}
	
	
	
	Prediction[] p = new BinaryPrediction[bags.length];

	/* 
	 * If we have bags that are empty, then we do not process them.
	 * If we have used up all of our time, then we can't do 
	 * any more iterations.
	 */
	if (bagsAreWeightless || m_s < 0) {
	    for (int i=0; i < bags.length; i++) {
		p[i] = getZeroPred();
	    }
	    return p;		
	}
	
	
	// Create a prediction array to accompany the exampleIndex array
	m_hypPredictions = new double[m_margins.length];
	for (int i=0; i < exampleIndex.length; i++){
	    int[] index = exampleIndex[i];
	    BrownBag b = (BrownBag)bags[i];
	    for (int j=0; j < index.length; j++){
		int example = index[j];
		m_hypPredictions[example] = b.calcPrediction(1.0).getClassScores()[0];
	    }
	}

	// we solve the constraints associated with
	// the BrownBoost model and obtain alpha.  gamma is a good
	// initial guess for alpha.
	double gamma = getHypErr(bags, exampleIndex);

	if (m_isCostSensitive) {
	    System.out.println("Solving positive example constraints");
	    double posAlpha = solve_constraints(gamma, m_posExamples);
	    System.out.println("Solving negative example constraints");
	    double negAlpha = solve_constraints(gamma, m_negExamples);
	    for (int i= 0; i < bags.length; i++) {
		p[i]= ((BrownBag) bags[i]).calcPrediction(posAlpha, negAlpha);
		System.out.println("p[" + i + "] = " + p[i]);
	    }
	} else {
	    double alpha = solve_constraints(gamma, m_examples);
	    for (int i= 0; i < bags.length; i++) {
		p[i]= ((BrownBag) bags[i]).calcPrediction(alpha);
		System.out.println("p[" + i + "] = " + p[i]);
	    }
	}
	return p;
    }
    

    public Bag newBag(int[] list) {
	return new BrownBag(list);
    }

    public Bag newBag() {
	return new BrownBag();
    }

    public Bag newBag(Bag bag) {
	return new BrownBag((BrownBag) bag);
    }

    /**
     * Returns the prediction associated with a bag representing a subset of the
     * data.
     */
    protected Prediction getPrediction(Bag b) {
	return ((BrownBag) b).calcPrediction();
    }

} /* End BrownBoost Class */
