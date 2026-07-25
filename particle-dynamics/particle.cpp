// define a particle structure class

struct particle {
    float x = 0.0f, y = 0.0f, x_old = 0.0f, y_old = 0.0f;
    float vx = 0.0f, vy = 0.0f, vx_old = 0.0f, vy_old = 0.0f;
    float ax = 0.0f, ay = 0.0f, ax_old = 0.0f, ay_old = 0.0f;
    float x_ref = 0.0f, y_ref = 0.0f;
    float Fxx = 1.0f, Fxy = 0.0f, Fyx = 0.0f, Fyy = 1.0f;
    float J = 1.0f;
    float m = 1000e-6;
    int body_id = 0;
};
