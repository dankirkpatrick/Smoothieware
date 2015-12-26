#ifndef _DCDELTALEVELINGSTRATEGY
#define _DCDELTALEVELINGSTRATEGY

#include <tuple>
#include <vector>
#include "LevelingStrategy.h"

#define dc_delta_calibration_strategy_checksum CHECKSUM("dc-delta-calibration")

class StreamOutput;

class DCDeltaCalibrationStrategy : public LevelingStrategy
{
public:
    DCDeltaCalibrationStrategy(ZProbe* zprobe) : LevelingStrategy(zprobe) {};
    ~DCDeltaCalibrationStrategy() {};
    bool handleGcode(Gcode* gcode);
    bool handleConfig();
private:
    bool set_trim(float x, float y, float z, StreamOutput* sream);
    bool get_trim(float &x, float &y, float &z);
    float compute_derivative(int factor, float cartesian_mm[], std::array<float,3> actuator_mm);
    void get_probe_point(int sample_number, float cartesian_mm[], int sample_count, float probe_radius);
    bool probe_bed(int sample_count, float probe_radius, std::vector<float> probe_heights, StreamOutput* stream);
    bool calibrate(Gcode* gcode);
    bool calibrate(int numFactors, int sample_count, float probe_radius, bool keep, StreamOutput* stream);
    std::tuple<float, float> getProbeLocation(int probe_point_id);
    std::tuple<float, float, float> transform(float probe_x, float effector_y, float effector_z);
    std::tuple<float, float, float> inverseTranform(float motor_a, float motor_b, float motor_c);
        
    float probe_radius;
    int sample_count;
    int factors;
};

#endif
