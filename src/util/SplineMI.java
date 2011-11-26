package util;

import java.util.Arrays;
import java.util.HashMap;
import java.io.*;

/**
 * InfoKit from CLR.
 * @revision Wei-Yi Cheng
 * @version 0.21
 * @date 05/13/2011
 */
public class SplineMI {

    public static final int DEPENDENCE_INDEPENDENT = 1;
    public static final int DEPENDENCE_POSITIVE = 2;

    public static final int NUM_INTEGRAL_PTS = 100;

    public static final int NEGATIVE = -1;
    public static final int POSITIVE = 1;

    public static final double EULER_MASCHERONI = 0.577215664901532;

    private static final double LOG_2 = Math.log(2);

    public static float simpsonNormCdf(float a, float b, int m, float[] params) {
        float[] evalResult;
        float point, h, s0, s1, s2, J;
        int i;

        h = (b - a) / (2 * m);

        evalResult = new float[2 * m + 1];
        for (i = 0; i <= 2 * m; i++) {
            point = a + i * h;
            evalResult[i] = normPdf(point, params);
        }
        s0 = evalResult[0] + evalResult[2 * m];
        s1 = 0;
        for (i = 1; i < 2 * m; i += 2) {
            s1 += evalResult[i];
        }
        s2 = 0;
        for (i = 2; i < 2 * m; i += 2) {
            s2 += evalResult[i];
        }

        J = h / 3 * (s0 + 4 * s1 + 2 * s2);
        return J;
    }

    public static float fdr(float[] p, int numElem, float alpha, int dependence) {
        float[] tmp = new float[numElem];
        int i, rightmost;
        float Cm;

        for (i = 0; i < numElem; i++) {
            tmp[i] = p[i];
        }
        Arrays.sort(tmp);
        // qsort(tmp, numElem, sizeof(float), (void*) compare_doubles);
        if (dependence == DEPENDENCE_INDEPENDENT || dependence == DEPENDENCE_POSITIVE) {
            Cm = 1;
        } else { /* NEGATIVE */
            Cm = (float) (Math.log(numElem) + EULER_MASCHERONI);
        }
        rightmost = -1;
        for (i = 0; i < numElem; i++) {
            if (tmp[i] < i / (numElem * Cm) * alpha)
                rightmost = i;
        }
        if (rightmost < 0) {
            return rightmost;
        } else {
            return tmp[rightmost];
        }
    }

    public static float normPdf(float x, float[] params) {
        float mu = params[0];
        float sigma = params[1];
        float ePower = -(x - mu) * (x - mu) / (2 * sigma * sigma);
        float pdf = (float) (1 / (sigma * Math.sqrt(2 * Math.PI)) * Math.exp(ePower));

        return pdf;
    }

    public static float normCdf(float val, float mu, float sigma) {
        float result;
        float params[] = new float[2];
        params[0] = mu;
        params[1] = sigma;
        result = simpsonNormCdf(-20, val, NUM_INTEGRAL_PTS / 2, params);
        return result;
    }


    public static float zToP(float z, int tail) {
        float p;

        if (tail == POSITIVE) { /* positive z is significant */
            p = 1 - normCdf(z, 0, 1);
        } else if (tail == NEGATIVE) { /* negative z is significant */
            p = normCdf(z, 0, 1);
        } else { /* absolute value is significant */
            p = normCdf(z, 0, 1);
            if (p > 0.5)
                p = 2 * (1 - p);
            else
                p = 2 * p;
        }
        return p;
    }

    public static float log2f(float x) {
        return (float) (Math.log(x) / LOG_2);
    }

    public static double log2d(double x) {
        return Math.log(x) / LOG_2;
    }

    public static double mean(double[] data, int numSamples) {
        int curSample;
        double mean = 0;

        for (curSample = 0; curSample < numSamples; curSample++) {
            mean += data[curSample];
        }
        return mean / (double) numSamples;
    }

    public static double std(double[] data, int numSamples) {
        int curSample;
        double m = mean(data, numSamples);
        double std = 0;

        for (curSample = 0; curSample < numSamples; curSample++) {
            std += (data[curSample] - m) * (data[curSample] - m);
        }
        return Math.sqrt(1 / (double) (numSamples - 1) * std);
    }

    public static float meanf(float[] data, int numSamples) {
        int curSample;
        float mean = 0;

        for (curSample = 0; curSample < numSamples; curSample++) {
            mean += data[curSample];
        }
        return mean / (float) numSamples;
    }

    public static float stdf(float[] data, int numSamples) {
        int curSample;
        float m = meanf(data, numSamples);
        float std = 0;

        for (curSample = 0; curSample < numSamples; curSample++) {
            std += (data[curSample] - m) * (data[curSample] - m);
        }
        return (float) Math.sqrt(1 / (float) (numSamples - 1) * std);
    }

    public static double meani(int[] data, int numSamples) {
        int curSample;
        double mean = 0;

        for (curSample = 0; curSample < numSamples; curSample++) {
            mean += (double) data[curSample];
        }
        return mean / (double) numSamples;
    }

    public static double mediani(int[] data, int numElem) {
        int[] dataCopy = new int[numElem];
        System.arraycopy(data, 0, dataCopy, 0, numElem);
        int half = numElem / 2;
        double median;

        Arrays.sort(dataCopy);
        median = ((double) dataCopy[half] + (double) dataCopy[half + 1]) / 2.0;
        return median;
    }

    public static double stdi(int[] data, int numSamples) {
        int curSample;
        double m = meani(data, numSamples);
        double std = 0;

        for (curSample = 0; curSample < numSamples; curSample++) {
            std += (double) (data[curSample] - m) * (data[curSample] - m);
        }
        return Math.sqrt(1 / (double) (numSamples - 1) * std);
    }

/* Accepts sorted input */

    public static double iqr(double[] data, int numSamples) {
        double q1, q3;
        int idx;
        numSamples = numSamples - 1; /* convert into indices */

        if ((numSamples + 1) % 2 == 0) {
            if ((numSamples + 1) % 4 == 0) {
                q1 = (data[numSamples / 4] + data[numSamples / 4 + 1]) / 2;
                q3 = (data[numSamples * 3 / 4] + data[numSamples * 3 / 4 + 1]) / 2;
            } else {
                q1 = data[(int) Math.ceil(numSamples / 4)];
                q3 = data[(int) Math.ceil(numSamples * 3 / 4)];
            }
        } else {
            idx = (int) Math.ceil(numSamples / 4);
            q1 = (data[idx] + data[idx + 1]) / 2;
            idx = (int) Math.ceil(numSamples * 3 / 4);
            q3 = (data[idx] + data[idx + 1]) / 2;
        }

        return q3 - q1;
    }

    public static void clrUnweightedStouffer(float[][] miMatrix, float[][] clrMatrix, int numVars) {
        float[] m = new float[numVars];
        float[] s = new float[numVars];
        float[][] curMI = new float[numVars][numVars];
        int i, j, k;
        float a, b;

        for (j = 0; j < numVars; j++) {
            for (k = 0; k < numVars; k++) {
                if (j == k) {
                    continue;
                }
                curMI[j + k * numVars] = miMatrix[j + k * numVars];
            }
        }
        for (i = 0; i < numVars; i++) {
            m[i] = meanf(curMI[i], numVars);
            s[i] = stdf(curMI[i], numVars);
        }
        for (i = 0; i < numVars - 1; i++) {
            for (j = i + 1; j < numVars; j++) {
                a = (curMI[i][j] - m[i]) / s[i];
                b = (curMI[i][j] - m[j]) / s[j];
                /* printf("a; %f, b: %f\n", a, b); */
                clrMatrix[i][j] = (float) ((a + b) / Math.sqrt(2));
                clrMatrix[j][i] = clrMatrix[i][j];
            }
        }
    }

    public static void clrGauss(float[][] miMatrix, float[][] clrMatrix, int numVars) {
        float[] m = new float[numVars];
        float[] s = new float[numVars];
        float[][] curMI = new float[numVars][numVars];
        int i, j, k;
        float a, b;

        for (j = 0; j < numVars; j++) {
            for (k = 0; k < numVars; k++) {
                if (j == k) {
                    continue;
                }
                curMI[j + k * numVars] = miMatrix[j + k * numVars];
            }
        }

        for (i = 0; i < numVars; i++) {
            m[i] = meanf(curMI[i], numVars);
            s[i] = stdf(curMI[i], numVars);
        }
        for (i = 0; i < numVars - 1; i++) {
            for (j = i + 1; j < numVars; j++) {
                a = (curMI[i][j] - m[i]) / s[i];
                b = (curMI[i][j] - m[j]) / s[j];
                /* printf("a; %f, b: %f\n", a, b); */
                if (a < 0)
                    a = 0;
                if (b < 0)
                    b = 0;
                clrMatrix[i][j] = (float) Math.sqrt(a * a + b * b);
                clrMatrix[j][i] = clrMatrix[i][j];
            }
        }
    }

/* zero-stage rule, see e.g. Wand M.P., Data-Based Choice of Histogram Bin Width */

    public static double binWidth(double[] data, int numSamples) {
        double s = std(data, numSamples);
        double i = iqr(data, numSamples) / 1.349;
        double shat = s < i ? s : i;
        return 3.49 * shat * (double) Math.pow((double) numSamples, (double) -1 / 3);
    }

/* Operates on a matrix of data */

    public static int[] calcNumBins(double[][] data, int numVars, int numSamples, double binMultiplier) {
        int[] binCount = new int[numVars];
        int curVar, curSample;
        double[][] sData = new double[numVars][numSamples];

        for (curVar = 0; curVar < numVars; curVar++)
            for (curSample = 0; curSample < numSamples; curSample++)
                sData[curVar][curSample] = data[curVar][curSample];


        for (curVar = 0; curVar < numVars; curVar++) {
            Arrays.sort(sData[curVar]);
            binCount[curVar] = (int) Math.ceil((sData[curVar][numSamples - 1] - sData[curVar][0]) / binWidth(sData[curVar], numSamples) * binMultiplier);
        }
        return binCount;
    }

/* "Lifted" from http://local.wasp.uwa.edu.au/~pbourke/curves/spline/ */
/*
  Calculate the blending value, this is done recursively.

  If the numerator and denominator are 0 the expression is 0.
  If the denominator is 0, the expression is 0
*/

    public static double SplineBlend(int k, int t, int[] u, double v, int n) {
        double value = 0.0;
        double d1, d2;

        if (t == 1) {
            if (((u[k] <= v) && (v < u[k + 1])) ||
                    (Math.abs(v - u[k + 1]) < 1e-10 && (k + 1 == n)))
                value = 1;
            else
                value = 0;
        } else {
            d1 = u[k + t - 1] - u[k];
            d2 = u[k + t] - u[k + 1];

            if ((d1 == 0) && (d2 == 0))
                value = 0;
            else if (d1 == 0)
                value = (u[k + t] - v) / (double) d2 * SplineBlend(k + 1, t - 1, u, v, n);
            else if (d2 == 0)
                value = (v - u[k]) / (double) d1 * SplineBlend(k, t - 1, u, v, n);
            else
                value = (v - u[k]) / (double) d1 * SplineBlend(k, t - 1, u, v, n) +
                        (u[k + t] - v) / (double) d2 * SplineBlend(k + 1, t - 1, u, v, n);
        }
        if (value < 0) {
            value = 0; /* rounding sometimes makes this < 0, e.g. -0.000000001 */
        }
        return (value);
    }

/*
  The positions of the subintervals of v and breakpoints, the position
  on the curve are called knots. Breakpoints can be uniformly defined
  by setting u[j] = j, a more useful series of breakpoints are defined
  by the function below. This set of breakpoints localises changes to
  the vicinity of the control point being modified.
*/

    public static void splineKnots(int[] u, int n, int t) {
    	int j;
        int d = n - 1;

        for (j = 0; j <= d + t; j++) {
            if (j < t)
                u[j] = 0;
            else if (j <= d)
                u[j] = u[j - 1] + 1;
            else if (j > d)
                u[j] = u[d] + 1;
        }
    }

    public static float max(float[] data, int numSamples) {
        int curSample;
        float curMax = Float.MIN_VALUE;

        for (curSample = 0; curSample < numSamples; curSample++) {
            if (data[curSample] > curMax) {
                curMax = data[curSample];
            }
        }
        return curMax;
    }

    public static float min(float[] data, int numSamples) {
        int curSample;
        float curMin = Float.MAX_VALUE;

        for (curSample = 0; curSample < numSamples; curSample++) {
            if (data[curSample] < curMin) {
                curMin = data[curSample];
            }
        }
        return curMin;
    }

    public static int maxi(int[] data, int numSamples) {
        int curSample;
        int curMax = data[0];

        for (curSample = 1; curSample < numSamples; curSample++) {
            if (data[curSample] > curMax) {
                curMax = data[curSample];
            }
        }
        return curMax;
    }

    public static int mini(int[] data, int numSamples) {
        int curSample;
        int curMin = data[0];

        for (curSample = 1; curSample < numSamples; curSample++) {
            if (data[curSample] < curMin) {
                curMin = data[curSample];
            }
        }
        return curMin;
    }

    public static void xToZ(float[] fromData, float[] toData, int numSamples, int splineOrder, int numBins) {
        int curSample;
        double xMax = max(fromData, numSamples);
        double xMin = min(fromData, numSamples);

        for (curSample = 0; curSample < numSamples; curSample++) {
        	toData[curSample] = (float) ((fromData[curSample] - xMin) * (numBins - splineOrder + 1) / (double) (xMax - xMin));
        }
    }

    public static void xToZFixed(float[] fromData, float[] toData, int numSamples, int splineOrder, int numBins, double xMax, double xMin) {
        int curSample;

        for (curSample = 0; curSample < numSamples; curSample++) {
            toData[curSample] = (float) ((fromData[curSample] - xMin) * (numBins - splineOrder + 1) / (double) (xMax - xMin));
        }
    }
    
    // 05/13/2011 added
    public static float[][] findWeights(String[] text, int n){
		HashMap<String, Integer> dict = new HashMap<String, Integer>();
		boolean[] valid = new boolean[n];
		int bins = 0;
		
		for(int i = 0; i < n; i++){
			if(!text[i].equalsIgnoreCase("null")){
				valid[i] = true;
				if(!dict.containsKey(text[i])){
					dict.put(text[i], bins);
					bins++;
				}
			}
		}
		
		float[][] weights;
		if(bins < 2){
			weights = null;
		}else{
			weights = new float[bins][n];
			for(int i = 0; i < n; i++){
				if(valid[i]){
					weights[dict.get(text[i])][i] = 1;
				}else{
					for(int j = 0; j < bins; j++){
						weights[j][i] = Float.NaN;
					}
				}
			}
		}
		return weights;
	}
    public static void findWeights(float[] x, int[] knots, float[][] weights, int numSamples, int splineOrder, int numBins) {
        int curSample;
        int curBin;
        float[] z = new float[numSamples];

        /*
          for (curSample = 0; curSample < numBins + splineOrder; curSample++) {
          mexPrintf("knot %d: %d\n", curSample, knots[curSample]);
          }
        */
        xToZ(x, z, numSamples, splineOrder, numBins);

        for (curSample = 0; curSample < numSamples; curSample++) {
            for (curBin = 0; curBin < numBins; curBin++) {
                weights[curBin][curSample] = Float.isNaN(z[curSample])? (float)1/numBins :(float) SplineBlend(curBin, splineOrder, knots, z[curSample], numBins);
                /* weights[curBin * numSamples + curSample] = splineFunction(z[curSample], splineOrder, curBin + 1, numBins); */
                /* mexPrintf("%d|%f(%f)\t", curBin, weights[curBin * numSamples + curSample],z[curSample]); */
            }
        }
    }
    
    
    public static void findWeights(float[] x, int[] knots, float[][] weights,	boolean[] valid, int numSamples, int splineOrder, int numBins) {
    	int curSample;
        int curBin;
        float[] z = new float[numSamples];

        /*
          for (curSample = 0; curSample < numBins + splineOrder; curSample++) {
          mexPrintf("knot %d: %d\n", curSample, knots[curSample]);
          }
        */
        xToZ(x, z, numSamples, splineOrder, numBins);

        for (curSample = 0; curSample < numSamples; curSample++) {
        	valid[curSample] = Float.isNaN(x[curSample])? false : true;
        	for (curBin = 0; curBin < numBins; curBin++) {
                weights[curBin][curSample] =Float.isNaN(z[curSample])? (float)1/numBins : (float) SplineBlend(curBin, splineOrder, knots, z[curSample], numBins);
            }
        }
	}
    public static void findWeightsFixed(float[] x, int[] knots, float[][] weights, int numSamples, int splineOrder, int numBins, double min, double max) {
        int curSample;
        int curBin;
        float[] z = new float[numSamples];

        /*
          for (curSample = 0; curSample < numBins + splineOrder; curSample++) {
          mexPrintf("knot %d: %d\n", curSample, knots[curSample]);
          }
        */
        xToZFixed(x, z, numSamples, splineOrder, numBins, min, max);

        for (curSample = 0; curSample < numSamples; curSample++) {
            for (curBin = 0; curBin < numBins; curBin++) {
                weights[curBin][curSample] = (float) SplineBlend(curBin, splineOrder, knots, z[curSample], numBins);
                /* weights[curBin * numSamples + curSample] = splineFunction(z[curSample], splineOrder, curBin + 1, numBins); */
                /* mexPrintf("%d|%f(%f)\t", curBin, weights[curBin * numSamples + curSample],z[curSample]); */
            }
        }
    }
    /*
    public static void hist1d(float[] hist, float[][] w, int numBins) {
        int curSample;
        int curBin;
        int numSamples = w[0].length;
        
        for (curBin = 0; curBin < numBins; curBin++) {
            hist[curBin] = 0;
            //int n = numSamples;
            for (curSample = 0; curSample < numSamples; curSample++) {
		//if(v[curSample]) hist[curBin] += w[curBin][curSample];
		//if(!v[curSample]) n--;
            	//else {hist[curBin] += w[curBin][curSample];}
		hist[curBin] += w[curBin][curSample];
            }
            //System.out.println("Valid samples: " + n);
            hist[curBin] /= (double) numSamples;
        }
    }
    public static void hist2d(float[][] wx, float[][] wy, float[][] hist, int numBins) {
        int curSample;
        int curBinX, curBinY;
        int numSamples = wx[0].length;
        
        for (curBinX = 0; curBinX < numBins; curBinX++) {
        	for (curBinY = 0; curBinY < numBins; curBinY++) {
                hist[curBinX][curBinY] = 0;
                //int n = numSamples;
                for (curSample = 0; curSample < numSamples; curSample++) {
                	//if(vx[curSample]&&vy[curSample]) hist[curBinX][curBinY] += wx[curBinX][curSample] * wy[curBinY][curSample];
			//if(! (vx[curSample]&&vy[curSample])) n--;
                	//else {hist[curBinX][curBinY] += wx[curBinX][curSample] * wy[curBinY][curSample];}
			hist[curBinX][curBinY] += wx[curBinX][curSample] * wy[curBinY][curSample];
                }
                //System.out.println("Valid samples: " + n);
                hist[curBinX][curBinY] /= (double) numSamples;
            }
        }

    }
	
    public static void hist3d(float[][] wx, float[][] wy, float[][] wz, float[][][] hist, float numEffectSamples, int numBins) {
        int curSample;
        int curBinX, curBinY, curBinZ;
        int numSamples = wx[0].length;
        
        for (curBinX = 0; curBinX < numBins; curBinX++) {
            for (curBinY = 0; curBinY < numBins; curBinY++) {
                for (curBinZ = 0; curBinZ < numBins; curBinZ++) {
                    hist[curBinX][curBinY][curBinZ] = 0;
                    for (curSample = 0; curSample < numSamples; curSample++) {
                        hist[curBinX][curBinY][curBinZ] += wx[curBinX][curSample] * wy[curBinY][curSample] * wz[curBinZ][curSample] / numEffectSamples;
                    }
                }
            }
        }
    }
    */
    public static double entropy1d(float[][] weights, float numSamples, int numBins) {
       double H = 0;
       float w;
        for (int curBin = 0; curBin < numBins; curBin++) {
            float histVal = 0;
            float n = numSamples;
            for(int curSample = 0; curSample < numSamples; curSample++){
            	w = weights[curBin][curSample];
            	if(!Float.isNaN(w)){
            		histVal += w;
            	}else{
            		n--;
            	}
            }
            histVal /= n;
            if (histVal > 0) {
        		H -= histVal * log2d(histVal);
        	}
        }
        return H;
    }
    public static double entropy2d(float[][] wx, float[][] wy, float numSamples, int numBins) {
        int curBinX, curBinY;
        double H = 0;
        float wxy;
        //int base = numBins * numBins;
        for (curBinX = 0; curBinX < numBins; curBinX++) {
            for (curBinY = 0; curBinY < numBins; curBinY++) {
                float histVal = 0;
                float n = numSamples;
                for (int curSample = 0; curSample < numSamples; curSample++) {
                	wxy = wx[curBinX][curSample] * wy[curBinY][curSample];
                	if(!Float.isNaN(wxy)){
                		histVal += wxy;
                	}else{
                		n--;
                	}
                }
                histVal /= n;
            	if (histVal > 0) {
                    H -= histVal * log2d(histVal);
                }
            }
        }
        return H;
    }
    
    public static double entropy2d(float[][] wx, float[][] wy, float numSamples, int nbx, int nby) {
        int curBinX, curBinY;
        double H = 0;
        float wxy;
        //int base = nbx * nby;
        for (curBinX = 0; curBinX < nbx; curBinX++) {
            for (curBinY = 0; curBinY < nby; curBinY++) {
                float histVal = 0;
                float n = numSamples;
                for (int curSample = 0; curSample < numSamples; curSample++) {
                	wxy = wx[curBinX][curSample] * wy[curBinY][curSample];
                	if(!Float.isNaN(wxy)){
                		histVal += wxy;
                	}else{
                		n--;
                	}
                }
                histVal /= n;
            	if (histVal > 0) {
                    H -= histVal * log2d(histVal);
                }
            }
        }
        return H;
    }
    
    public static double entropy3d(float[][] wx, float[][] wy, float[][] wz, float numSamples, int numBins) {
        int curBinX, curBinY, curBinZ, curSample;
        double H = 0;
        float wxyz;
        float n = numSamples;
        //int base = numBins * numBins * numBins;
        for (curBinX = 0; curBinX < numBins; curBinX++) {
            for (curBinY = 0; curBinY < numBins; curBinY++) {
                for (curBinZ = 0; curBinZ < numBins; curBinZ++) {
                    float histVal = 0;
                    for (curSample = 0; curSample < numSamples; curSample++) {
                        wxyz = wx[curBinX][curSample] * wy[curBinY][curSample] * wz[curBinZ][curSample];
                    	if(!Float.isNaN(wxyz)){
                    		histVal += wxyz;
                    	}else{
                    		n--;
                    	}
                    
                    }
                	histVal/=n;
                	if(histVal > 0){
                		H -= histVal * log2d(histVal);
                	}
                }
            }
        }
        return H;
    }

    
    /**
     * Args:
     * <ol>
     * <li> Microarray set (no labels).
     * <li> TF indices.
     * <li> Number of bins.
     * <li> Spline order.
     * <li> Total Segments (-1 for continuous).
     * <li> Segment Number.
     * <li> "all", or the indices of the targets.
     * <li> "all", or the indices of the partners.
     * <li> (Optional) Binary file with 2D MI matrix.
     * </ol>
     *
     * @param args
     */
    public static void main(String[] args) throws Exception {
    	int u[] = new int[7+3];
    	float x[] = {1f, 2f, 3f, 4f, 5f, 6f, 7f};
    	float y[] = {1f, 2f, 3f, 4f, 5f, 6f, 7f};
    	float z[] = {3f, 6f, 7f, 2f, 5f, 3f, 2f};
    	float wx[][] = new float[7][7];
    	float wy[][] = new float[7][7];
    	float wz[][] = new float[7][7];
    	
    	splineKnots(u, 7, 3);
    	findWeights(x, u, wx, 7, 3, 7);
    	findWeights(y, u, wy, 7, 3, 7);
    	findWeights(z, u, wz, 7, 3, 7);
    	
    	double e1x = entropy1d(wx, 7, 7);
    	double e1y = entropy1d(wy, 7, 7);
    	double e2xy = entropy2d(wx,wy, 7, 7);
    	System.out.println(e1x + "\t" + e1y);
    	System.out.println(e1x + e1y - e2xy);
    	System.out.println(entropy3d(wx,wy,wz, 7, 7));
    }
    /*
    public static void main(String[] args) throws Exception {
        //// Load microarray data
        float[][] data = AdPar3D.parseDataFile(args[0]);
        //// Load TF indices
        int[] tfs = AdPar3D.parseIndices(args[1]);
        //// Params
        int numGenes = data.length;
        int numSamples = data[0].length;
        int numTFs = tfs.length;
        int numBins = Integer.parseInt(args[2]);
        int splineOrder = Integer.parseInt(args[3]);
        System.out.println("==== Genes: " + numGenes + ", Arrays: " + numSamples + ", TFs: " + numTFs);
        System.out.println("==== Bins: " + numBins + ", Spline Order: " + splineOrder);
        //// Determine segment information
        int totalSegments = Integer.parseInt(args[4]);
        int segment = Integer.parseInt(args[5]) - 1;
        boolean fixedRange = Boolean.parseBoolean(args[6]);
        int[] targets;
        if ("all".equalsIgnoreCase(args[7])) {
            targets = new int[numGenes];
            for (int i = 0; i < numGenes; i++) {
                targets[i] = i;
            }
        } else {
            targets = AdPar3D.parseIndices(args[7]);
        }
        int numTargets = targets.length;
        int[] partners;
        if ("all".equalsIgnoreCase(args[8])) {
            partners = new int[numGenes];
            for (int i = 0; i < numGenes; i++) {
                partners[i] = i;
            }
        } else {
            partners = AdPar3D.parseIndices(args[8]);
        }
        int numPartners = partners.length;
        long start, stop;
        long n = ((long) numPartners) * numTargets * numTFs;
        if (totalSegments == -1) {
            start = 0;
            stop = n;
        } else {
            start = segment * n / totalSegments;
            stop = (segment + 1) * n / totalSegments;
        }
        System.out.println("Start: " + start + ", Stop: " + stop);
        DataOutputStream out;
        if (totalSegments == -1) {
            out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("spline3.synergy")));
        } else {
            out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("spline3-" + String.format("%04d", segment) + ".synergy")));
        }
        //// Construct knots
        int[] knots = new int[numBins + splineOrder];
        splineKnots(knots, numBins, splineOrder);

        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        if (fixedRange) {
            for (float[] values : data) {
                for (float v : values) {
                    if (v > max) {
                        max = v;
                    }
                    if (v < min) {
                        min = v;
                    }
                }
            }
        }

        //// Compute weights for all genes
        System.out.println("Computing weights...");
        float[][][] weights = new float[numGenes][numBins][numSamples];
        for (int i = 0; i < numGenes; i++) {
            if (fixedRange) {
                findWeightsFixed(data[i], knots, weights[i], numSamples, splineOrder, numBins, min, max);
            } else {
                findWeights(data[i], knots, weights[i], numSamples, splineOrder, numBins);
            }
        }
        //// Compute e1 for all genes
        System.out.println("Computing 1d entropies...");
        float[] e1 = new float[numGenes];
        {
            float[] hist1d = new float[numBins];
            for (int i = 0; i < numGenes; i++) {
                e1[i] = (float) entropy1d(weights[i], hist1d, numSamples, numBins);
            }
        }
        float[][] e2 = new float[numGenes][numGenes];
        //// Load e2 if specified
        if (args.length > 8) {
            System.out.println("Loading 2d entropies from '" + args[9] + "'...");
            DataInputStream miIn = new DataInputStream(new BufferedInputStream(new FileInputStream(args[9])));
            for (int y = 0; y < numGenes; y++) {
                for (int z = y + 1; z < numGenes; z++) {
                    float mi = miIn.readFloat();
                    float ent = e1[y] + e1[z] - mi;
                    e2[y][z] = ent;
                    // For convenience, store it both ways
                    e2[z][y] = ent;
                }
            }
        } else {
            //// Compute e2 for all gene pairs
            System.out.println("Computing 2d entropies...");
            {
                float[][] hist2d = new float[numBins][numBins];
                for (int i = 0; i < numGenes; i++) {
                    for (int j = i + 1; j < numGenes; j++) {
                        e2[i][j] = (float) entropy2d(weights[i], weights[j], hist2d, numSamples,numBins);
                        e2[j][i] = e2[i][j];
                    }
                }
            }
        }
        //// Compute mi3
        System.out.println("Computing 3d entropies...");
        {
            long index = 0;
            float[][][] hist3d = new float[numBins][numBins][numBins];
            for (int tf = 0; tf < numTFs; tf++) {
                int x = tfs[tf];
                for (int yi = 0; yi < numPartners; yi++) {
                    for (int zi = 0; zi < numTargets; zi++) {
                        if ((index >= start) && (index < stop)) {
                            int y = partners[yi];
                            int z = targets[zi];
                            float e3 = (float) entropy3d(weights[x], weights[y], weights[z], hist3d, numSamples, numBins);
                            //e3xyz - e2xy - e2xz - e2yz + e1x + e1y + e1z;
                            float mi3 = e3 - e2[x][y] - e2[x][z] - e2[y][z] + e1[x] + e1[y] + e1[z];
                            out.writeFloat(-mi3);
                        }
                        index++;
                    }
                }
            }
        }
        out.close();
    }
	*/
}
