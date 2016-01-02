#include "DCDeltaCalibrationStrategy.h"
#include "Kernel.h"
#include "Config.h"
#include "Robot.h"
#include "StreamOutputPool.h"
#include "Gcode.h"
#include "checksumm.h"
#include "ConfigValue.h"
#include "PublicDataRequest.h"
#include "EndstopsPublicAccess.h"
#include "PublicData.h"
#include "Conveyor.h"
#include "ZProbe.h"
#include "BaseSolution.h"

#include <fastmath.h>
#include <tuple>
#include <vector>

#define PIOVER180   0.01745329251994329576923690768489F
#define PI 3.14159265358979F

#define sample_count_checksum CHECKSUM("sample_count")
#define factors_checksum CHECKSUM("factors")
#define radius_checksum       CHECKSUM("radius")

// deprecated
#define probe_radius_checksum CHECKSUM("probe_radius")

bool DCDeltaCalibrationStrategy::handleConfig()
{
    // default is probably wrong
    float r = THEKERNEL->config->value(leveling_strategy_checksum, dc_delta_calibration_strategy_checksum, radius_checksum)->by_default(-1)->as_number();
    if (r == -1) {
        // deprecated config syntax]
        r = THEKERNEL->config->value(zprobe_checksum, probe_radius_checksum)->by_default(100.0F)->as_number();
    }
    this->probe_radius = r;

    this->sample_count = THEKERNEL->config->value(leveling_strategy_checksum, dc_delta_calibration_strategy_checksum, sample_count_checksum)->by_default(13)->as_number();

    this->factors = THEKERNEL->config->value(leveling_strategy_checksum, dc_delta_calibration_strategy_checksum, factors_checksum)->by_default(6)->as_number();
    
    return true;
}

bool DCDeltaCalibrationStrategy::handleGcode(Gcode *gcode)
{
    if (gcode->has_g) {
        // G code processing
        if (gcode->g == 32) { // DC auto calibration for delta
            // first wait for an empty queue i.e. no moves left
            THEKERNEL->conveyor->wait_for_empty_queue();
            if (!calibrate(gcode)) {
                gcode->stream->printf("Calibration failed to complete, probe not triggered\n");
                return true;
            }
            gcode->stream->printf("Calibration complete, save settings with M500\n");
            return true;
        }

    } else if(gcode->has_m) {
        // handle mcodes
    }

    return false;
}

/* Run the David Crocker calibration algorithm for delta calibration.

   The David Crocker algorithm calculates the derivative, or sensitivity,
   of each component at each sample point.  This is used to calculate a 
   normal matrix for each component.  This normal matrix, combined with the
   height errors from the samples, is used to calculate a solution matrix
   for each factor that is tuned.
   
   This approach simplifies delta kinematics: individual tower radius 
   adjustments are mathematically identical to a new circle with a single 
   delta radius for all towers, with slight adjustments to each tower angle.
   Furthermore, one of the 3 tower angles is assumed to be '0', with all 
   other tower angles relative to the first.  This reduces an 11-variable
   equation (DR, RL, 3 endstops, 3 delta radius adjustments, 3 tower angle
   adjustments) to a mathematically identical, bunt sipler, equation with
   7 variables (DR, RL, 3 endstops, 2 tower angles).
   
   The usage is as follows:

   G32 (Fx) (Sx) (Rx.xx) (K)
        Fx : Number of factors.  Valid values are 3, 4, 6, 7.  
             Default is 6.
             3 = calibrate endstops for towers A,B,C
             4 = calibrate endstops and delta radius
             6 = calibrate endstops, delta radius, and tower angles for A,B
             7 = calibrate endstops, delta radius, tower angles, and rod length
             NOTE: autocalibration of rod length generally should be avoided
        Sx : Number of sample points.  Should be a mutliple of 6 plus 1.  
             Default is 13.
        Rx.xx : Maximum probe radius.  Default is 100.0.
                Probes always include the center.  Probes are performed in
                a series of circles, 6 points per circle, varying the radius
                linearly from 0 to the probe radius.
     
 */ 
bool DCDeltaCalibrationStrategy::calibrate(Gcode *gcode)
{
    int factors = this->factors;
    if (gcode->has_letter('F')) {
        factors = gcode->get_value('F');
        if (factors != 3 && factors != 4 && factors != 6 && factors != 7) {
            gcode->stream->printf("Number of factors for DC calibration is incorrect--must be 3, 4, 6, or 7");
            return false;
        }
    }
    
    int sample_count = this->sample_count;
    if (gcode->has_letter('S')) {
        sample_count = gcode->get_value('S');
    }
    
    float probe_radius = this->probe_radius;
    if (gcode->has_letter('R')) {
        probe_radius = gcode->get_value('R');
    }
    
    bool keep = false;
    if (gcode->has_letter('K')) {
        keep = true;
    }
    
    return calibrate(factors, sample_count, probe_radius, keep, gcode->stream);
}

bool DCDeltaCalibrationStrategy::calibrate(int num_factors, int sample_count, float probe_radius, bool keep, StreamOutput* stream)
{
    float cartesian_mm[3];
    std::array<float, 3> actuator_mm;
    std::vector<std::array<float, 3>> all_actuator_mm(sample_count);
    BaseSolution::arm_options_t options;
    
    // Get or reset starting point for options, depending on 'keep'
    float trimx = 0.0F, trimy = 0.0F, trimz = 0.0F;
    if (!keep) {
        // zero trim values
        if (!set_trim(0, 0, 0, stream)) return false;
        stream->printf("Current Trim X: %.3f, Y: %.3f, Z: %.3f\n", trimx, trimy, trimz);

        options['A'] = 0;
        options['B'] = 0;
        options['C'] = 0;
        options['D'] = 0;
        options['E'] = 0;
        options['F'] = 0;
        THEKERNEL->robot->arm_solution->set_optional(options);
    } else {
        // get current trim, and continue from that
        if (get_trim(trimx, trimy, trimz)) {
            stream->printf("Current Trim X: %.3f, Y: %.3f, Z: %.3f\n", trimx, trimy, trimz);
        } else {
            stream->printf("Could not get current trim, are endstops enabled?\n");
            return false;
        }
    }
    if (THEKERNEL->robot->arm_solution->get_optional(options)) {
        stream->printf("    Rod length: %.3f", options['L']);
        stream->printf("  Delta radius: %.3f", options['R']);
        stream->printf("Tower A radius: %.3f", options['R'] + options['A']);
        stream->printf("Tower A  angle: %.3f", 210 + options['D']);
        stream->printf("Tower B radius: %.3f", options['R'] + options['B']);
        stream->printf("Tower B  angle: %.3f", 330 + options['E']);
        stream->printf("Tower C radius: %.3f", options['R'] + options['C']);
        stream->printf("Tower C  angle: %.3f", 90 + options['F']);
    }

    // Sample bed Z height
    std::vector<float> probe_heights(sample_count, 0);
    probe_bed(sample_count, probe_radius, probe_heights, stream);
    
    // Remember actuator heights for probe points
    for (int i = 0; i < sample_count; i++) {
        get_probe_point(i, cartesian_mm, sample_count, probe_radius);
        cartesian_mm[2] = probe_heights[i];
        THEKERNEL->robot->arm_solution->cartesian_to_actuator(cartesian_mm, actuator_mm);
        all_actuator_mm[i][0] = actuator_mm[0];
        all_actuator_mm[i][1] = actuator_mm[1];
        all_actuator_mm[i][2] = actuator_mm[2];
    }

    // The amount of Z height, at each probe point, changed by the solution 
    std::vector<float> corrections(sample_count, 0);   
    // Used to calculation deviation at end of run
    float initialSumOfSquares = 0;
    for (int i = 0; i < sample_count; i++) {
        initialSumOfSquares += (probe_heights[i] * probe_heights[i]);
    }
    
    // The derivative of height change for each point wrt each factor
    std::vector<float> tempInit(num_factors, 0);
    std::vector<std::vector<float>> derivative_matrix(sample_count, tempInit);
    
    // The normal matrix for each factor: N x N+1 
    std::vector<float> tempInit2(num_factors+1, 0);
    std::vector<std::vector<float>> normal_matrix(num_factors, tempInit2);
    
    // The error residuals after 
    std::vector<float> residuals(sample_count, 0);
    
    // The solution--the adjustments for each of the factors
    std::vector<float> solution(num_factors, 0);
    
    for (int iteration = 0; iteration < 2; iteration++) {
        // Build derivative matrix
        for (int i = 0; i < sample_count; i++) {
            get_probe_point(i, cartesian_mm, sample_count, probe_radius);
            cartesian_mm[2] = probe_heights[i];
            THEKERNEL->robot->arm_solution->cartesian_to_actuator(cartesian_mm, actuator_mm);
            for (int j = 0; j < num_factors; j++) {
                derivative_matrix[i][j] = compute_derivative(j, cartesian_mm, actuator_mm);
            }
        }
        
        // Debug
        stream->printf("Derivative Matrix:\n");
        for (int i = 0; i < sample_count; i++) {
            for (int j = 0; j < num_factors; j++) {
                if (j == 0) {
                    stream->printf("[");
                } else {
                    stream->printf(",");
                }
                stream->printf("%.3f", derivative_matrix[i][j]);
            }
            stream->printf("]\n");
        }
        
        // Build normal matrix from derivatives
        for (int i = 0; i < num_factors; i++) {
            for (int j = 0; j < num_factors; j++) {
                float temp = derivative_matrix[0][i] * derivative_matrix[0][j];
                for (int k = 1; k < sample_count; k++) {
                    temp += derivative_matrix[k][i] * derivative_matrix[k][j];
                }
                normal_matrix[i][j] = temp;
            }
            float temp = derivative_matrix[0][i] * -(probe_heights[i] + corrections[i]);
            for (int k = 1; k < sample_count; k++) {
                temp += derivative_matrix[k][i] * -(probe_heights[k] + corrections[k]);
            }
            normal_matrix[i][num_factors] = temp;
        }
        
        // Debug
        stream->printf("Normal matrix:\n");
        for (int i = 0; i < num_factors; i++) {
            for (int j = 0; j <= num_factors; j++) {
                if (j == 0) {
                    stream->printf("[");
                } else {
                    stream->printf(",");
                }
                stream->printf("%.3f", normal_matrix[i][j]);
            }
            stream->printf("]\n");
        }
        
        // Solve using Gauss-Jordan
        for (int i = 0; i < num_factors; i++) {
            float vmax = fabs(normal_matrix[i][i]);
            for (int j = i+1; j <= num_factors; j++) {
                float rmax = fabs(normal_matrix[j][i]);
                if (rmax > vmax) {
                    normal_matrix[i].swap(normal_matrix[j]);
                    vmax = rmax;
                }
            }
            
            float v = normal_matrix[i][i];
            for (int j = 0; j < i; j++) {
                float factor = normal_matrix[j][i] / v;
                normal_matrix[j][i] = 0;
                for (int k = i+1; k <= num_factors; k++) {
                    normal_matrix[j][k] -= normal_matrix[i][k] * factor;
                }
            }
            
            for (int j = i+1; j < num_factors; j++) {
                float factor = normal_matrix[j][i] / v;
                normal_matrix[j][i] = 0;
                for (int k = i+1; k <= num_factors; k++) {
                    normal_matrix[j][k] -= normal_matrix[i][k] * factor;
                }
            }
        }
        
        for (int i = 0; i < num_factors; i++) {
            solution[i] = normal_matrix[i][num_factors] / normal_matrix[i][i];
        }
        
        // Debug
        stream->printf("Solution:\n");
        for (int i = 0; i < num_factors; i++) {
            if (i == 0) {
                stream->printf("[");
            } else {
                stream->printf(",");
            }
            stream->printf("%.3f", solution[i]);
        }
        stream->printf("]\n");
        
        // compute residuals
        float sumOfSquares = 0;
        for (int i = 0; i < sample_count; i++) {
            residuals[i] = probe_heights[i];
            for (int j = 0; j < num_factors; j++) {
                residuals[i] += solution[j] * derivative_matrix[i][j];
            }
            sumOfSquares += residuals[i] * residuals[i];
        }
        stream->printf("Before RMS Error: %.3f", sqrtf(sumOfSquares));
        
        // apply solution
        for (int i = 0; i < num_factors; i++) {
            switch (i) {
                case 0:
                    trimx += solution[i];
                    break;
                case 1:
                    trimy += solution[i];
                    break;
                case 2:
                    trimz += solution[i];
                    break;
                case 3:
                    options['R'] += solution[i];
                    break;
                case 4:
                    options['D'] += solution[i] * PIOVER180;
                    break;
                case 5:
                    options['E'] += solution[i] * PIOVER180;
                    break;
                case 6:
                    options['L'] += solution[i];
                    break;
            }
            if (!set_trim(trimx, trimy, trimz, stream)) return false;
            THEKERNEL->robot->arm_solution->set_optional(options);
        }        
        
        // compute expected residuals
        sumOfSquares = 0;
        for (int i = 0; i < sample_count; i++) {
            actuator_mm[0] = all_actuator_mm[i][0] + trimx;
            actuator_mm[1] = all_actuator_mm[i][1] + trimy;
            actuator_mm[2] = all_actuator_mm[i][2] + trimz;
            THEKERNEL->robot->arm_solution->actuator_to_cartesian(actuator_mm, cartesian_mm);
            corrections[i] = cartesian_mm[2];
            residuals[i] = probe_heights[i] + cartesian_mm[2];
            sumOfSquares += residuals[i] * residuals[i];
        }
        stream->printf("Expected After RMS Error: %.3f", sqrtf(sumOfSquares));
    }

    return true;
}

bool DCDeltaCalibrationStrategy::set_trim(float x, float y, float z, StreamOutput *stream)
{
    float t[3] {x, y, z};
    bool ok = PublicData::set_value( endstops_checksum, trim_checksum, t);

    if (ok) {
        stream->printf("set trim to X:%f Y:%f Z:%f\n", x, y, z);
    } else {
        stream->printf("unable to set trim, is endstops enabled?\n");
    }

    return ok;
}

bool DCDeltaCalibrationStrategy::get_trim(float &x, float &y, float &z)
{
    void *returned_data;
    bool ok = PublicData::get_value( endstops_checksum, trim_checksum, &returned_data );

    if (ok) {
        float *trim = static_cast<float *>(returned_data);
        x = trim[0];
        y = trim[1];
        z = trim[2];
        return true;
    }
    return false;
}

bool DCDeltaCalibrationStrategy::probe_bed(int sample_count, float probe_radius, std::vector<float> probe_heights, StreamOutput* stream) {
    zprobe->home();

    // find bed, run at fast rate
    int s;
    if (!zprobe->run_probe(s, true)) return false;

    float bedht = zprobe->zsteps_to_mm(s) - zprobe->getProbeHeight(); // distance to move from home to 5mm above bed
    stream->printf("Bed ht is %f mm\n", bedht);

    // move to start position
    zprobe->home();
    zprobe->coordinated_move(NAN, NAN, -bedht, zprobe->getFastFeedrate(), true); // do a relative move from home to the point above the bed

    float cartesian_mm[3];
    for (int i = 0; i < sample_count; i++) {
        get_probe_point(i, cartesian_mm, sample_count, probe_radius);
        if (!zprobe->doProbeAt(s, cartesian_mm[0], cartesian_mm[1])) return false;
        probe_heights[i] = zprobe->zsteps_to_mm(s) - zprobe->getProbeHeight();
    }
    
    return true;
}

void DCDeltaCalibrationStrategy::get_probe_point(int sample_number, float cartesian_mm[], int sample_count, float probe_radius) {
    // The last sample is always the center
    if (sample_number == (sample_count-1)) {
        cartesian_mm[0] = 0;
        cartesian_mm[1] = 0;
    } else {
        int circle_count = floor((sample_count - 1) / 6);
        // Work from the outside in...this is the circle #
        int circle = circle_count - floor(sample_number / 6);
        // Each successive circle gets closer to the center
        float radius = (probe_radius * circle) / floor((sample_count - 1) / 6);
        // Probe points are at 0, 60, 120, 180, 240, and 300 degress for odd numbered circles
        //              and at 30, 90, 150, 210, 270, and 330 degrees for even numbered circles
        float angle = PI * (sample_number % 6) / 3.0;
        if ((circle % 2) == 0) {
            angle += (PI / 6.0);
        }
        
        cartesian_mm[0] = cosf(angle * PIOVER180) * radius;
        cartesian_mm[1] = sinf(angle * PIOVER180) * radius;
    }
}

float DCDeltaCalibrationStrategy::compute_derivative(int factor, float cartesian_mm[], std::array<float,3> actuator_mm) {
    float perturb = 0.2;
    float zLo = 0;
    float zHi = 0;
    float original;
    BaseSolution::arm_options_t options;
    
    if (factor < 0) {
        // error
        return 0;
    } else if (factor < 3) {
        // perturb endstops
        original = actuator_mm[factor];
        // DSK: TODO: check that this works--it's opposite the DeltaSim
        actuator_mm[factor] = original - perturb; // endstop offsets have inverted effects
        THEKERNEL->robot->arm_solution->actuator_to_cartesian(actuator_mm, cartesian_mm);
        zHi = cartesian_mm[2];
        
        actuator_mm[factor] = original + perturb; // endstop offsets have inverted effects
        THEKERNEL->robot->arm_solution->actuator_to_cartesian(actuator_mm, cartesian_mm);
        zLo = cartesian_mm[2];
        
        actuator_mm[factor] = original;
    } else if (factor < 7) {
        // perturb option
        char option = '-';
        switch (factor) {
            case 3:
                option = 'R'; // Delta radius adjustment
                break;
            case 4:
                option = 'D'; // Angle adjustment for tower A
                break;
            case 5:
                option = 'E'; // Angle adjustment for tower B
                break;
            case 6:
                option = 'L'; // Rod length adjustment
                break;
        }
        if (THEKERNEL->robot->arm_solution->get_optional(options)) {
            original = options[option];
            
            options[option] = original + perturb; 
            THEKERNEL->robot->arm_solution->set_optional(options);
            THEKERNEL->robot->arm_solution->actuator_to_cartesian(actuator_mm, cartesian_mm);
            zHi = cartesian_mm[2];
            
            options[option] = original - perturb; 
            THEKERNEL->robot->arm_solution->set_optional(options);
            THEKERNEL->robot->arm_solution->actuator_to_cartesian(actuator_mm, cartesian_mm);
            zLo = cartesian_mm[2];

            options[option] = original;
            THEKERNEL->robot->arm_solution->set_optional(options);
        }
    } else {
        // error
        return 0;
    }
    return (zHi - zLo) / (2 * perturb);
}
