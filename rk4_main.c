/*!
 *-------------------------------------------------------------------------
 *
 *  @file	rk4.c
 *
 *  @brief	Fourth Order Runge Kutta Integrator
 *
 *-------------------------------------------------------------------------
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

// #define FREQ  (1/60)          // Hz
// #define FREQ  (0.0000462962)    // Hz
#define FREQ (0.0000462962 * 8)
#define CAP   (0.02)        	// Farad 
#define Q_INIT (3.6 * 3600)
#define BUFFER (0.25)

double i;    		// unit = ampere
double i_amp;

// Define new variables for Equivalent Circuit Model
double Q = 3.6;
double R0 = 0.07;
double R1 = 0.014;
double C1 = 42000;
double Soc = 0.0;

// v1 is changed continiously (alongisde i) as time is updated. It is graphed on the y axis.
// q1 is changed continiously (alongside i) as time is updated. It is also graphed on the y-axis.
double v1 = 0.0;   // Volts
// double q1 = 12960;  // Amp per second
// double q1 = Q_INIT;
double q1 = 0.0;

// Lithium Ion Battery Voltage Chart
// Used to determine how State of Charge (SOC) relates to Voltage Released (VoC).
// Voc is plotted against Dod
const double Soc_array[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
// 3.6 is an outlier (replaced with 3.26)
const double Voc_array[] = {2.5, 3.0, 3.2, 3.22, 3.25, 3.26, 3.27, 3.3, 3.32, 3.35, 3.4};

// double Voc = Voc_array[9];
double Voc = Voc_array[0];
double v_term = Voc_array[0];


// const double Dod_array[] = {1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0};
// const double Voc_array[] = {3.4, 3.35, 3.32, 3.3, 3.27, 3.6, 3.25, 3.22, 3.2, 3.0, 2.5};


// Voc = 4.0;


// Array consisting of graph labels.
// The order corresponds to the order of values needed to be graphed.

//------------------------------------------------------------------------
//Function to approximate y' = y
//------------------------------------------------------------------------
// Integrate multiple differential equations 
    
// Updated Equations
// mode == 0: dV1(t)/dt (ROC of resistor + capacitor voltage with respect to time)
// V1 = Voc (voltage source) - resistance (current * resistance from R0)
double dydx(double x, double y, int index)
{
    if (index == 0) {
        // Returns dV1(t)/dt: value derived from first the given ODE
        // printf ("dv1/dt: {%f}   \n", (-1/(R1 * C1)) * y + R1/(R1 * C1) * i);
        return (-1/(R1 * C1)) * y + R1/(R1 * C1) * i;
    }
    else if (index == 1) {
        // printf ("dq/dt: {%f}   \n", Voc);
        // Returns dQ(t)/dt: value derived from second given ODE
        return -i;
    }  
}

double linear_interpolate(const double *x_array, const double *y_array, int size, double x) {
    // Check if x is out of bounds
    if (x < x_array[0] || x > x_array[size - 1]) {
        printf("Error: x is out of bounds.\n");
        return -1; // Return an error value
    }

    // Find the interval in which x lies
    for (int i = 0; i < size - 1; i++) {
        if (x >= x_array[i] && x <= x_array[i + 1]) {
            // Perform linear interpolation
            double x1 = x_array[i];
            double y1 = y_array[i];
            double x2 = x_array[i + 1];
            double y2 = y_array[i + 1];

            // Calculate the slope (m) and the interpolated y value
            double slope = (y2 - y1) / (x2 - x1);
            double y = y1 + slope * (x - x1);
            return y;
        }
    }

    return -1; // Return an error value if no interval is found (should not happen with valid inputs)
}

//------------------------------------------------------------------------
// Compute yn+1 given xn, yn and step size h
//------------------------------------------------------------------------

// [CHANGED]: Method now takes in an array address & size to account for the array used for y values.
// Method returns a new ptr array address to an array containing updated y values. 
double* runge_kutta_4th_order_array(
                double x, 
		        double *y_array, 
                int y_array_size,
			    double h, 
			    double (*f)(double, double, int)
            )
{
    
    // [CHANGED]: Allocate memory for y_nexts so we can return it later.
    double* y_nexts = (double*)malloc(y_array_size);
    // Memory allocation failure failsafe 
    if (y_nexts == NULL) {
        // printf("Memory allocation for y_nexts has failed. Check your code.\n");
        exit(1);
    } 

    // [CHANGED]: Iterate through all y values, applying rk4 to each. 
    // One double is 8 bytes. Dividing sizeof(double[]) by sizeof(double) gives us the number of array items.
    for (int iter = 0; iter < y_array_size / sizeof(double); iter++) {
        double k1 = h* (*f)(x, y_array[iter], iter);
        double k2 = h* (*f)(x + (h/2), y_array[iter] + (k1/2), iter);
        double k3 = h* (*f)(x + (h/2), y_array[iter] + (k2/2), iter);
        double k4 = h* (*f)(x + h, y_array[iter] + k3, iter);

        // "Append" values to y_nexts.
        y_nexts[iter] = y_array[iter] + (k1 + (2*k2) + (2*k3) + k4)/6.0;
    }

    // [CHANGED]: Return pointer object for y_nexts. 
    return y_nexts;
}

//------------------------------------------------------------------------
// Process an array consisting of all miscellaneous y values we're to graph that don't involve rk4 differentiation. 
//------------------------------------------------------------------------

double misc_methods(int index)
{
    if (index == 0) {
        // Returns v_terminal for our circuit model (Voc - v1 - total resistance)
        // printf ("VOC: {%f}   ", Voc);
        // printf ("V_Terminal: {%f}\n", Voc - v1 - R0 * i);
        if (Voc - v1 - R0 * i <= 0) {
            return Voc_array[10];
        }
        else {
            return Voc - v1 - R0 * i;
        }
    }
}

// To do: Pass in pointer rather than allocating new memory
// Rewrite malloc for misc_new_values
double* process_misc_y_array(double *y_array, int array_size, double (*f)(int)){ 
    // Iterate through given array. Changes each value according to different functions depending on index.
    for (int iter = 0; iter < array_size / sizeof(double); iter++) {
        y_array[iter] = (*f)(iter);
    }

    // Return updated array
    return y_array;
}



// Method used for updating min and max variables.
// Returns the new min/max variable. is_maximum is true for updating maximums & false for updating minimums.
double updateMinMax(double min_max, double new_value, bool is_maximum) {
    double bucket = min_max;
    if (is_maximum) {
        bucket = (new_value > min_max) ? new_value : min_max;
    }
    else {
        bucket = (new_value < min_max) ? new_value : min_max;
    }
    return bucket;
}


//------------------------------------------------------------------------
//  main program 
//------------------------------------------------------------------------
int main()
{
    // Define default variables needed for processing.
    double t, h;
    t = 0.0;
    h = 0.01;      // each step is 1 second
    i = 0.0;
    i_amp = 0.5;

    // Prepare to process array for differentiation (v_array) and miscellaneous array (misc_array)
    double v_array[] = {v1, q1};
    double misc_array[] = {v_term};

    // Allocate memory for v_new_values so we can return it later.
    double* v_new_values = (double*)malloc(sizeof(v_array));
    printf("..\n");
    // Memory allocation failure failsafe 
    if (v_new_values == NULL) {
        printf("Memory allocation for v_new_values has failed.\n");
        exit(1);
    }

    // Allocate memory for misc_new_values so we can return it later.
    
    
    double* misc_new_values = (double*)malloc(sizeof(misc_array));

    
    // Allocate memories for max_v & min_v arrays for range setting.
    // 0 = i (current)
    // 1 = v1 (voltage)
    // 2 = q1 (battery capacity)
    // 3 = v_term (terminal voltage)
    double max_v_array[4] = {-999999.0, -999999.0, -999999.0, -999999.0};
    double min_v_array[4] = {999999.0, 999999.0, 999999.0, 999999.0};
    
    // File processing
    printf("...\n");
    FILE *file = fopen("output.dat", "w");
    printf("....\n");
    fprintf(file, "Time Current Voltage \n");
    
    // ----  Integrating  -----

    // Let each iteration be 1 ms seconds. We need 63 seconds for 10 sin cycles.
    // This will populate the memory defined earlier through v_new_values. 
    printf("test-0: output.dat opened successfully.\n");

    // 4800000
    for (int iter = 0; iter < 9000000; iter++)
    {
        v_new_values = runge_kutta_4th_order_array(t, v_array, sizeof(v_array), h, dydx);

        // [2.0 CHANGE]: The code that dictates i updates is no longer a sin function, but instead a square function.
        i = (sin(2.0 * M_PI * FREQ * t) >= 0) ? i_amp * -1 : i_amp * 0;
        // i = 0.5;

        // Update max_v & min_v in relation to i.
        max_v_array[0] = updateMinMax(max_v_array[0], i, true);
        min_v_array[0] = updateMinMax(min_v_array[0], i, false);

        // [CHANGED]: Updates all values within v_array & q_array to be each cooresponding v_new_value.
        for (int iter2 = 0; iter2 < sizeof(v_array) / sizeof(double); iter2++) {
            // printf("ITER2: [%d]\n", iter2);
            v_array[iter2] = v_new_values[iter2];

            // Update v1 & q1 for next interation.
            if (iter2 == 0) {
                v1 = v_new_values[iter2];
                // printf("v1: [%f] \n", v1);
            }
            else if (iter2 == 1) {
                q1 = v_new_values[iter2];
                // printf("q1: {%f} \n", q1);

                // UPDATE VOC based on linear interpolation for SoC.
                Soc = q1 / Q_INIT;
                Voc = linear_interpolate(Soc_array, Voc_array, 11, Soc);

            }

            // Update max_v_array & min_v_array in relation to v.
            max_v_array[iter2 + 1] = updateMinMax(max_v_array[iter2 + 1], v_new_values[iter2], true);
            min_v_array[iter2 + 1] = updateMinMax(min_v_array[iter2 + 1], v_new_values[iter2], false);
        }

	    // Logging to file
        // Log time & current.
	    fprintf(file, "%g %g ", t, i);

        // Log processed ODEs. 
        // printf("%ld %ld\n", sizeof(double), sizeof(v_array));
        for (int iter2 = 0; iter2 < sizeof(v_array) / sizeof(double); iter2++) {
            fprintf(file, "%g ", v_array[iter2]);
        }
        
        // [4.O CHANGE]: Process misc_new_values using our misc_processing framework.
        misc_new_values = process_misc_y_array(misc_array, sizeof(misc_array), misc_methods);

        // Log misc values.
        for (int iter2 = 0; iter2 < sizeof(misc_new_values) / sizeof(double); iter2++) {
            fprintf(file, "%g ", misc_new_values[iter2]);
            max_v_array[iter2 + 3] = updateMinMax(max_v_array[iter2 + 3], misc_new_values[iter2], true);
            min_v_array[iter2 + 3] = updateMinMax(min_v_array[iter2 + 3], misc_new_values[iter2], false);
        }
        
        if (misc_new_values[0] <= 2.5 || Soc > 1.0 || q1 < 0 || q1 > Q_INIT) {
            break;
        }

        // for (int iter3 = 0; iter3 < sizeof(max_v_array) / sizeof(double); iter3++){
            // printf("MAX_V_ARRAY_VALUE %d: [%f] .... ", iter3, max_v_array[iter3]);
            // printf("MIN_V_ARRAY_VALUE %d: [%f]     \n", iter3, min_v_array[iter3]);
        // }
        // printf("VOC: %f\n", Voc);

        // Update time, v_array, & i for next iteration    
        t = t + h; 

        // Newline
        fprintf(file, "\n");

    }

    printf("test-1: rk4 iteration + misc processing has been processed & logged.\n");

    // Free allocated memory
    free (v_new_values);
    // free (misc_new_values);

    // Close file
    fclose(file);


    printf("test-2: output.dat has been freed & closed successfully.\n");

    // ----  Plotting -----
    // Each seperate unit is now plotted on its own graph.
    // new code

    for (int iter = 0; iter < sizeof(max_v_array) / sizeof(double); iter++) {
        // Create file
        FILE *gp = popen("gnuplot -persist", "w");
        if (!gp) {
            fprintf(stderr, "Unable to open gnuplot!\n");
            return 1;
        }

        // Sets title to determine in what units we're plotting. 
        
        switch(iter){
            case 0:
                fprintf(gp, "set multiplot layout 1,1 title 'Graph of Current (Amps)'\n");
                break;
            case 1:
                fprintf(gp, "set multiplot layout 1,1 title 'Graph of v1 (Volts)'\n");
                break;
            case 2:
                fprintf(gp, "set multiplot layout 1,1 title 'Graph of Charge Capacity (Ampere Seconds)'\n");
                break;
            case 3:
                fprintf(gp, "set multiplot layout 1,1 title 'Graph of v term (Volts)'\n");
                break; 
        }
        // x axis
        fprintf(gp, "set xlabel 'Time (s)'\n");
        fprintf(gp, "set grid\n"); 
        // The +.5 & -.5s exist to create blank space on top/bottom of graph. Can easily adjust to be scalable. 
        fprintf(gp, "set ylabel 'V'; set yrange[%lf:%lf]; set ytics \n", min_v_array[iter] - BUFFER, max_v_array[iter] + BUFFER); 
        // fprintf(gp, "set y2label 'A'; set y2range[%f:%f]; set y2tics \n", min_v_array[iter] - 0.01, max_v_array[iter] + 0.01); 
        
        // Plot the respective value (based on iter; see above key) from file
        switch(iter){
            case 0:
                fprintf(gp, "plot 'output.dat' using 1:%d with l lw 2 title 'Current'", iter + 2);
                break;
            case 1:
                fprintf(gp, "plot 'output.dat' using 1:%d with l lw 2 title 'Voltage'", iter + 2);
                break;
            case 2:
                fprintf(gp, "plot 'output.dat' using 1:%d with l lw 2 title 'Batt. Capacity'", iter + 2);
                break;
            case 3:
                fprintf(gp, "plot 'output.dat' using 1:%d with l lw 2 title 'Term. Voltage'", iter + 2);
                break; 
        }
        //fprintf(gp, ", '' using 1:%d with l lw 2 title 'Value'", iter);
        
        // Newline flush
        fprintf(gp, "\n");

        // Flush & close gp
        fflush(gp);
        pclose(gp);

        // Confirmation message
        printf("test-2.5: new graphing test complete. Iteration: %d | Y range: [%f - %f]\n", iter, min_v_array[iter] - 0.01, max_v_array[iter] + 0.01);
    }

    // Free up max & min v_arrays
    //free(max_v_array);
    //free(min_v_array);

    // end of new code

    return 0;
}
