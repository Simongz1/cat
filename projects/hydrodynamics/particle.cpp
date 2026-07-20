// define a particle structure class

struct particle {
    float x = 0.0f, y = 0.0f, x_old = 0.0f, y_old = 0.0f;
    float vx = 0.0f, vy = 0.0f, vx_old = 0.0f, vy_old = 0.0f;
    float ax = 0.0f, ay = 0.0f, ax_old = 0.0f, ay_old = 0.0f;
    float m = 1.0f;
    int body_id = 0;
};
