// define a particle structure class

#include <vector>
#include <cmath>

struct particle {
    float x;
    float y;
    
    std::vector<float> position = {x, y};
    std::vector<float> velocity = {0., 0.};
};