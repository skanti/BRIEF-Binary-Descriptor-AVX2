export void kernel1(uniform int p_x_a[], uniform int p_y_a[], uniform int p_x_b[], uniform int p_y_b[],
                    uniform int i_x_a[], uniform int i_y_a[], uniform int i_x_b[], uniform int i_y_b[],
                    uniform float cos_angle, uniform float sin_angle) {
    foreach(i = 0 ... 256) {
        i_x_a[i] = (int)(p_x_a[i] * cos_angle - p_y_a[i] * sin_angle);
        i_y_a[i] = (int)(p_x_a[i] * sin_angle + p_y_a[i] * cos_angle);
        i_x_b[i] = (int)(p_x_b[i] * cos_angle - p_y_b[i] * sin_angle);
        i_y_b[i] = (int)(p_x_b[i] * sin_angle + p_y_b[i] * cos_angle);
    }
}