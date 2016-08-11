#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/features2d.hpp>
#include <iostream>
#include <random>
#include "BRIEF.h"
#include "Timer.h"

#define N_CHANNELS 1
#define WIDTH_IMAGE 640
#define HEIGHT_IMAGE 480
#define STRIDE_IMAGE WIDTH_IMAGE*N_CHANNELS

#if N_CHANNELS == 1
#define CV_8UCX CV_8UC1
#define CV_LOAD_IMAGE_X 0
#elif N_CHANNELS == 3
#define CV_LOAD_IMAGE_X 1
#endif

typedef std::conditional<SIZE_BITS_HAMING == 32, int32_t, int64_t>::type haming_type;

void create_synthetic_data(std::vector<int> &x, std::vector<int> &y, std::vector<float> &angle, int n_features) {
    std::mt19937 mt(999);
    std::uniform_real_distribution<float> u_dist(0, 1);
    for (int i = 0; i < n_features; i++) {
        x[i] = (int) (u_dist(mt) * WIDTH_IMAGE);
        y[i] = (int) (u_dist(mt) * HEIGHT_IMAGE);
        angle[i] = u_dist(mt) * 3.14f;
    }
}

int main() {
    std::string dir = "/Users/amon/grive/development/BRIEF/picdump";
    std::vector<std::string> img_basenames = {"/apple.jpg", "/astronaut.jpg", "/bike.jpg", "/bmw.jpg",
                                              "/cake.jpg", "/cherry.jpg", "/dreamliner.jpg", "/macbook.jpg",
                                              "/orange.jpg"};
    int n = img_basenames.size();

    int n_dim_vec = N_DIM_BINARYDESCRIPTOR / SIZE_BITS_HAMING;
    int n_features = 500;
    std::vector<int> x(n_features), y(n_features);
    std::vector<float> angle(n_features);
    create_synthetic_data(x, y, angle, n_features);
    BRIEF brief;
    Matrix<haming_type> bd(n_dim_vec, n_features);
    double t_total = 0;
    for (int i = 0; i < n; i++) {
        cv::Mat image = cv::imread(dir + img_basenames[i], CV_LOAD_IMAGE_X);
        Timer::start();
        brief.rbrief(image.data, HEIGHT_IMAGE, WIDTH_IMAGE, N_CHANNELS, STRIDE_IMAGE,
                     x.data(), y.data(), angle.data(), n_features, bd.memptr(), bd.n_rows);
        Timer::stop();
        double t = Timer::get_timing_in_ms();
        t_total += t;
        std::cout << "timing (ms): " << t << std::endl;
    }
    std::cout << "timing average(ms): " << t_total / n << std::endl;
}