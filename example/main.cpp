#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <iostream>
#include <random>
#include "BRIEF.h"
#include "Timer.h"
#include <bitset>

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

struct Corners {
    int x, y;
    float angle;
};

void create_synthetic_data(Corners *corners, int n_features) {
    std::mt19937 mt(999);
    std::uniform_real_distribution<float> u_dist(0, 1);
    for (int i = 0; i < n_features; i++) {
        corners[i].x = (int) (u_dist(mt) * WIDTH_IMAGE);
        corners[i].y = (int) (u_dist(mt) * HEIGHT_IMAGE);
        corners[i].angle = u_dist(mt) * 3.14f;
    }
}

int main() {
    std::string dir = "/home/amon/grive/development/BRIEF/picdump";
    std::vector<std::string> img_basenames = {"/apple.jpg", "/astronaut.jpg", "/bike.jpg", "/bmw.jpg", "/cake.jpg", "/cherry.jpg", "/dreamliner.jpg", "/macbook.jpg", "/orange.jpg"};
    int n = img_basenames.size();

    int n_dim_vec = N_DIM_BINARYDESCRIPTOR / SIZE_BITS_HAMING;
    int n_features = 1 << 16;
    std::vector<Corners> corners(n_features);
    create_synthetic_data(corners.data(), n_features);
    Matrix<int64_t> bd(n_dim_vec, n_features);
    double t_total = 0;
    for (int i = 0; i < n; i++) {
        cv::Mat image = cv::imread(dir + img_basenames[i], CV_LOAD_IMAGE_X);
        BRIEF bd(4, n_features);
        Timer::start();
        bd.rbrief(image.data, HEIGHT_IMAGE, WIDTH_IMAGE, N_CHANNELS, STRIDE_IMAGE, corners.data(), n_features);
        Timer::stop();
        double t = Timer::get_timing_in_ms();
        t_total += t;
        std::cout << "timing (ms): " << t << std::endl;
        std::cout << std::bitset<64>(bd.bd(0, 0)).to_string() << " " << std::bitset<64>(bd.bd(1, 0)).to_string() << std::endl;
    }
    std::cout << "timing average(ms): " << t_total / n << std::endl;
}