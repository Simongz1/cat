using namespace std;
#include <Eigen/Sparse>

struct GeometricParameters{
    const std::string type = "rectangle";
    const float lx = 10;
    const float ly = 10;
    const float lz = 10;
    const float dx = 0.1;
    const float dy = 0.1;
    const float dz = 0.1;
};

//define an initialization function that executes the particles initialization

std::vector<particle> initParticles(GeometricParameters &geomParams){
    //retrieve the type of gemetry
    const std::string type = geomParams.type;
    std::vector<particle> particles;

    const float lx = geomParams.lx;
    const float dx = geomParams.dx;
    const float ly = geomParams.ly;
    const float dy = geomParams.dy;
    const float lz = geomParams.lz;
    const float dz = geomParams.dz;

    if (type == "rectangle"){
        //initialize a grid
        for (int i = 0; i < int(lx / dx); ++i){
            for (int j = 0; j < int(ly / dy); ++j){
                for (int k = 0; k < int(lz / dz); ++k){
                    particle currentParticle;
                    currentParticle.x = {i * dx, j * dy, k * dz};
                    currentParticle.x_old = currentParticle.x;


                    //initialize particles to the right with a nodzero velocity
                    if (j >= 0.95 * int(ly / dy)){
                        currentParticle.v = {0, -1, 0};
                    }

                    particles.emplace_back(currentParticle);
                }
            }
        }
    }

    return particles;
}