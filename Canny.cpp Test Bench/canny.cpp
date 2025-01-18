// canny.cpp
/*******************************************************************************
 * PROGRAM: canny_edge
 * PURPOSE: This program implements a "Canny" edge detector. The processing
 * steps are as follows:
 *
 *   1) Convolve the image with a separable Gaussian filter.
 *   2) Compute the first derivatives in the x and y directions.
 *   3) Compute the magnitude of the gradient.
 *   4) Perform non-maximal suppression.
 *   5) Perform hysteresis thresholding.
 *
 * The SystemC model is structured with the following modules:
 *
 *   - Top
 *     - Stimulus
 *     - Platform
 *       - DataIn
 *       - DUT (Canny Edge Detector)
 *       - DataOut
 *     - Monitor
 *
 * NAME: Yash Ingle
 * DATE: 1/11/2024
 *******************************************************************************/

#include <systemc.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cerrno>

using namespace std;

#define VERBOSE 0

// Update constants for white background and black edges
#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Define the IMAGE struct
struct IMAGE
{
    static const int COLS = 2704;
    static const int ROWS = 1520;
    static const int SIZE = ROWS * COLS;
    unsigned char img[SIZE];

    IMAGE()
    {
        memset(img, 0, SIZE * sizeof(unsigned char));
    }

    IMAGE &operator=(const IMAGE &copy)
    {
        if (this != &copy)
        {
            memcpy(img, copy.img, SIZE * sizeof(unsigned char));
        }
        return *this;
    }

    operator unsigned char *()
    {
        return img;
    }

    unsigned char &operator[](const int index)
    {
        return img[index];
    }
};

// SystemC Modules

SC_MODULE(Stimulus)
{
    sc_fifo_out<IMAGE> ImgOut;
    SC_CTOR(Stimulus)
    {
        SC_THREAD(process);
        set_stack_size(128 * 1024 * 1024);
    }
    void process()
    {
        for (int i = 1; i <= 30; ++i)
        {
            IMAGE image;
            char infilename[256];
            sprintf(infilename, "video/Engineering%03d.pgm", i);
            if (VERBOSE)
                printf("Stimulus: Reading image %s.\n", infilename);
            if (!read_pgm_image(infilename, image.img))
            {
                fprintf(stderr, "Stimulus: Error reading the input image %s.\n", infilename);
                continue;
            }
            ImgOut.write(image);
        }
    }

private:
    int read_pgm_image(const char *infilename, unsigned char image[])
    {
        FILE *fp;
        char buf[71];
        int img_rows, img_cols;
        int maxval;

        if (infilename == NULL)
            fp = stdin;
        else
        {
            if ((fp = fopen(infilename, "rb")) == NULL)
            {
                fprintf(stderr, "Error reading the file %s: %s\n", infilename, strerror(errno));
                return 0;
            }
        }

        // Read the header
        fgets(buf, 70, fp);
        if (strncmp(buf, "P5", 2) != 0)
        {
            fprintf(stderr, "The file %s is not in PGM format.\n", infilename);
            if (fp != stdin)
                fclose(fp);
            return 0;
        }

        // Skip comments
        do
        {
            fgets(buf, 70, fp);
        } while (buf[0] == '#' || buf[0] == '\n');

        // Read image dimensions
        sscanf(buf, "%d %d", &img_cols, &img_rows);
        if (img_rows != IMAGE::ROWS || img_cols != IMAGE::COLS)
        {
            fprintf(stderr, "Image dimensions (%d x %d) do not match expected dimensions (%d x %d).\n",
                    img_cols, img_rows, IMAGE::COLS, IMAGE::ROWS);
            if (fp != stdin)
                fclose(fp);
            return 0;
        }

        // Read maxval
        do
        {
            fgets(buf, 70, fp);
        } while (buf[0] == '#' || buf[0] == '\n');
        sscanf(buf, "%d", &maxval);

        // Read image data
        if ((size_t)(IMAGE::ROWS * IMAGE::COLS) != fread(image, 1, IMAGE::ROWS * IMAGE::COLS, fp))
        {
            fprintf(stderr, "Error reading the image data from %s.\n", infilename);
            if (fp != stdin)
                fclose(fp);
            return 0;
        }

        if (fp != stdin)
            fclose(fp);
        return 1;
    }
};

SC_MODULE(DataIn)
{
    sc_fifo_in<IMAGE> ImgIn;
    sc_fifo_out<IMAGE> ImgOut;
    SC_CTOR(DataIn)
    {
        SC_THREAD(process);
        set_stack_size(128 * 1024 * 1024);
    }
    void process()
    {
        while (true)
        {
            IMAGE image;
            ImgIn.read(image);
            ImgOut.write(image);
        }
    }
};

SC_MODULE(DUT)
{
    sc_fifo_in<IMAGE> ImgIn;
    sc_fifo_out<IMAGE> ImgOut;

    SC_CTOR(DUT)
    {
        SC_THREAD(process);
        set_stack_size(256 * 1024 * 1024); // Increased due to large arrays
    }

    void process()
    {
        while (true)
        {
            IMAGE image_in, image_out;
            ImgIn.read(image_in);
            canny(image_in.img, image_out.img);
            ImgOut.write(image_out);
        }
    }

private:
    // Constants declared inside the class without initialization
    static const float sigma;
    static const float tlow;
    static const float thigh;
    static const int KERNEL_SIZE = 21;
    static const float BOOSTBLURFACTOR;

    void canny(unsigned char image[], unsigned char edge[])
    {
        int rows = IMAGE::ROWS;
        int cols = IMAGE::COLS;

        unsigned char *nms = new unsigned char[rows * cols];
        short int *smoothedim = new short int[rows * cols];
        short int *delta_x = new short int[rows * cols];
        short int *delta_y = new short int[rows * cols];
        short int *magnitude = new short int[rows * cols];

        // Perform Gaussian smoothing on the image
        if (VERBOSE)
            printf("Smoothing the image using a Gaussian kernel.\n");
        gaussian_smooth(image, smoothedim);

        // Compute the first derivatives
        if (VERBOSE)
            printf("Computing the X and Y first derivatives.\n");
        derivative_x_y(smoothedim, delta_x, delta_y);

        // Compute the magnitude of the gradient
        if (VERBOSE)
            printf("Computing the magnitude of the gradient.\n");
        magnitude_x_y(delta_x, delta_y, magnitude);

        // Perform non-maximal suppression
        if (VERBOSE)
            printf("Performing non-maximal suppression.\n");
        memset(nms, 0, rows * cols * sizeof(unsigned char));
        non_max_supp(magnitude, delta_x, delta_y, nms);

        // Use hysteresis to mark the edge pixels
        if (VERBOSE)
            printf("Applying hysteresis thresholding.\n");
        memset(edge, NOEDGE, rows * cols * sizeof(unsigned char)); // Initialize with NOEDGE
        apply_hysteresis(magnitude, nms, edge);

        // Free allocated memory
        delete[] nms;
        delete[] smoothedim;
        delete[] delta_x;
        delete[] delta_y;
        delete[] magnitude;
    }

    void gaussian_smooth(unsigned char image[], short int smoothedim[])
    {
        int rows = IMAGE::ROWS;
        int cols = IMAGE::COLS;
        int r, c, rr, cc;
        int windowsize, center;
        float *tempim = new float[rows * cols];
        float kernel[KERNEL_SIZE];
        float dot, sum;

        // Create a 1-dimensional Gaussian smoothing kernel
        if (VERBOSE)
            printf("Computing the Gaussian smoothing kernel.\n");
        make_gaussian_kernel(sigma, kernel, &windowsize);
        center = windowsize / 2;

        // Blur in the x-direction
        if (VERBOSE)
            printf("Blurring the image in the X-direction.\n");
        for (r = 0; r < rows; r++)
        {
            for (c = 0; c < cols; c++)
            {
                dot = 0.0;
                sum = 0.0;
                for (cc = (-center); cc <= center; cc++)
                {
                    if (((c + cc) >= 0) && ((c + cc) < cols))
                    {
                        dot += (float)image[r * cols + (c + cc)] * kernel[center + cc];
                        sum += kernel[center + cc];
                    }
                }
                tempim[r * cols + c] = dot / sum;
            }
        }

        // Blur in the y-direction
        if (VERBOSE)
            printf("Blurring the image in the Y-direction.\n");
        for (c = 0; c < cols; c++)
        {
            for (r = 0; r < rows; r++)
            {
                dot = 0.0;
                sum = 0.0;
                for (rr = (-center); rr <= center; rr++)
                {
                    if (((r + rr) >= 0) && ((r + rr) < rows))
                    {
                        dot += tempim[(r + rr) * cols + c] * kernel[center + rr];
                        sum += kernel[center + rr];
                    }
                }
                smoothedim[r * cols + c] = (short int)(dot * BOOSTBLURFACTOR / sum + 0.5);
            }
        }

        delete[] tempim;
    }

    void make_gaussian_kernel(float sigma, float kernel[], int *windowsize)
    {
        int i, center;
        float x, fx, sum = 0.0;

        *windowsize = 1 + 2 * ceil(2.5 * sigma);
        center = (*windowsize) / 2;

        if (*windowsize > KERNEL_SIZE)
        {
            fprintf(stderr, "Kernel size exceeds KERNEL_SIZE (%d). Increase KERNEL_SIZE.\n", KERNEL_SIZE);
            exit(1);
        }

        if (VERBOSE)
            printf("The kernel has %d elements.\n", *windowsize);

        for (i = 0; i < (*windowsize); i++)
        {
            x = (float)(i - center);
            fx = exp(-0.5 * x * x / (sigma * sigma)) / (sigma * sqrt(2.0 * M_PI));
            kernel[i] = fx;
            sum += fx;
        }

        // Normalize the kernel
        for (i = 0; i < (*windowsize); i++)
            kernel[i] /= sum;

        if (VERBOSE)
        {
            printf("The filter coefficients are:\n");
            for (i = 0; i < (*windowsize); i++)
                printf("kernel[%d] = %f\n", i, kernel[i]);
        }
    }

    void derivative_x_y(short int smoothedim[], short int delta_x[], short int delta_y[])
    {
        int rows = IMAGE::ROWS;
        int cols = IMAGE::COLS;
        int r, c, pos;

        // Compute the x-derivative
        if (VERBOSE)
            printf("Computing the X-direction derivative.\n");
        for (r = 0; r < rows; r++)
        {
            pos = r * cols;
            delta_x[pos] = smoothedim[pos + 1] - smoothedim[pos];
            pos++;
            for (c = 1; c < (cols - 1); c++, pos++)
            {
                delta_x[pos] = smoothedim[pos + 1] - smoothedim[pos - 1];
            }
            delta_x[pos] = smoothedim[pos] - smoothedim[pos - 1];
        }

        // Compute the y-derivative
        if (VERBOSE)
            printf("Computing the Y-direction derivative.\n");
        for (c = 0; c < cols; c++)
        {
            pos = c;
            delta_y[pos] = smoothedim[pos + cols] - smoothedim[pos];
            pos += cols;
            for (r = 1; r < (rows - 1); r++, pos += cols)
            {
                delta_y[pos] = smoothedim[pos + cols] - smoothedim[pos - cols];
            }
            delta_y[pos] = smoothedim[pos] - smoothedim[pos - cols];
        }
    }

    void magnitude_x_y(short int delta_x[], short int delta_y[], short int magnitude[])
    {
        int rows = IMAGE::ROWS;
        int cols = IMAGE::COLS;
        int r, c, pos;
        int sq1, sq2;

        for (r = 0, pos = 0; r < rows; r++)
        {
            for (c = 0; c < cols; c++, pos++)
            {
                sq1 = (int)delta_x[pos] * (int)delta_x[pos];
                sq2 = (int)delta_y[pos] * (int)delta_y[pos];
                magnitude[pos] = (short)(sqrt((float)sq1 + (float)sq2) + 0.5);
            }
        }
    }

    void non_max_supp(short *mag, short *gradx, short *grady, unsigned char *result)
    {
        int rows = IMAGE::ROWS;
        int cols = IMAGE::COLS;
        int rowcount, colcount;

        // Initialize result to NOEDGE
        memset(result, NOEDGE, rows * cols * sizeof(unsigned char));

        // Suppress non-maximum points
        for (rowcount = 1; rowcount < rows - 1; rowcount++)
        {
            for (colcount = 1; colcount < cols - 1; colcount++)
            {
                int idx = rowcount * cols + colcount;
                int gradx_val = gradx[idx];
                int grady_val = grady[idx];
                float angle = atan2((float)grady_val, (float)gradx_val) * (180.0 / M_PI);
                if (angle < 0)
                    angle += 180;

                int mag1, mag2;
                if ((angle >= 0 && angle < 22.5) || (angle >= 157.5 && angle <= 180))
                {
                    // 0 degrees
                    mag1 = mag[idx - 1];
                    mag2 = mag[idx + 1];
                }
                else if (angle >= 22.5 && angle < 67.5)
                {
                    // 45 degrees
                    mag1 = mag[idx - cols - 1];
                    mag2 = mag[idx + cols + 1];
                }
                else if (angle >= 67.5 && angle < 112.5)
                {
                    // 90 degrees
                    mag1 = mag[idx - cols];
                    mag2 = mag[idx + cols];
                }
                else // angle >= 112.5 && angle < 157.5
                {
                    // 135 degrees
                    mag1 = mag[idx - cols + 1];
                    mag2 = mag[idx + cols - 1];
                }

                if (mag[idx] >= mag1 && mag[idx] >= mag2)
                    result[idx] = POSSIBLE_EDGE;
                else
                    result[idx] = NOEDGE;
            }
        }
    }

    void apply_hysteresis(short int mag[], unsigned char nms[], unsigned char edge[])
    {
        int rows = IMAGE::ROWS;
        int cols = IMAGE::COLS;
        int r, c, pos, highcount, lowthreshold, highthreshold;
        int hist[32768] = {0};
        short int maximum_mag = 0;

        // Compute the histogram of the magnitude image
        for (r = 0, pos = 0; r < rows; r++)
        {
            for (c = 0; c < cols; c++, pos++)
            {
                if (nms[pos] == POSSIBLE_EDGE)
                {
                    hist[mag[pos]]++;
                    if (mag[pos] > maximum_mag)
                        maximum_mag = mag[pos];
                }
            }
        }

        // Compute the high threshold
        int sum = 0;
        for (r = 1; r <= maximum_mag; r++)
        {
            sum += hist[r];
        }
        highcount = (int)(sum * thigh + 0.5);

        sum = 0;
        for (r = 1; r <= maximum_mag; r++)
        {
            sum += hist[r];
            if (sum > highcount)
                break;
        }
        highthreshold = r;
        lowthreshold = (int)(highthreshold * tlow + 0.5);

        if (VERBOSE)
        {
            printf("Computed thresholds: low = %d, high = %d\n", lowthreshold, highthreshold);
        }

        // Perform edge tracking by hysteresis
        for (r = 0, pos = 0; r < rows; r++)
        {
            for (c = 0; c < cols; c++, pos++)
            {
                if ((nms[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold))
                {
                    edge[pos] = EDGE;
                    follow_edges(edge, mag, pos, lowthreshold);
                }
            }
        }

        // Set remaining possible edges to NOEDGE
        for (pos = 0; pos < rows * cols; pos++)
        {
            if (edge[pos] != EDGE)
                edge[pos] = NOEDGE;
        }
    }

    void follow_edges(unsigned char edge[], short mag[], int pos, short lowval)
    {
        int cols = IMAGE::COLS;
        int rows = IMAGE::ROWS;
        int x[8] = {1, 1, 0, -1, -1, -1, 0, 1};
        int y[8] = {0, 1, 1, 1, 0, -1, -1, -1};
        int i;
        int new_pos;

        for (i = 0; i < 8; i++)
        {
            int new_r = (pos / cols) + y[i];
            int new_c = (pos % cols) + x[i];
            if (new_r >= 0 && new_r < rows && new_c >= 0 && new_c < cols)
            {
                new_pos = new_r * cols + new_c;
                if ((edge[new_pos] == POSSIBLE_EDGE) && (mag[new_pos] >= lowval))
                {
                    edge[new_pos] = EDGE;
                    follow_edges(edge, mag, new_pos, lowval);
                }
            }
        }
    }
};

// Initialize static constants after the class definition
const float DUT::sigma = 1.0;       // Adjusted sigma value
const float DUT::tlow = 0.3;        // Adjusted tlow value
const float DUT::thigh = 0.8;       // Adjusted thigh value

const float DUT::BOOSTBLURFACTOR = 90.0;

SC_MODULE(DataOut)
{
    sc_fifo_in<IMAGE> ImgIn;
    sc_fifo_out<IMAGE> ImgOut;
    SC_CTOR(DataOut)
    {
        SC_THREAD(process);
        set_stack_size(128 * 1024 * 1024);
    }
    void process()
    {
        while (true)
        {
            IMAGE image;
            ImgIn.read(image);
            ImgOut.write(image);
        }
    }
};

SC_MODULE(Monitor)
{
    sc_fifo_in<IMAGE> ImgIn;
    SC_CTOR(Monitor)
    {
        SC_THREAD(process);
        set_stack_size(128 * 1024 * 1024);
    }
    void process()
    {
        for (int i = 1; i <= 30; ++i)
        {
            IMAGE edge;
            char outfilename[256];
            sprintf(outfilename, "Engineering%03d_edges.pgm", i);
            ImgIn.read(edge);
            if (VERBOSE)
                printf("Monitor: Writing edge image %s.\n", outfilename);
            if (!write_pgm_image(outfilename, edge.img))
            {
                fprintf(stderr, "Monitor: Error writing the edge image %s.\n", outfilename);
            }
        }
        sc_stop();
    }

private:
    int write_pgm_image(const char *outfilename, unsigned char image[])
    {
        FILE *fp;
        int cols = IMAGE::COLS;
        int rows = IMAGE::ROWS;

        if (outfilename == NULL)
        {
            fprintf(stderr, "write_pgm_image: No output filename provided.\n");
            return 0;
        }

        if ((fp = fopen(outfilename, "wb")) == NULL)
        {
            fprintf(stderr, "Error writing the file %s: %s\n", outfilename, strerror(errno));
            return 0;
        }

        // Write the header
        fprintf(fp, "P5\n%d %d\n255\n", cols, rows);

        // Write image data
        if ((size_t)(rows * cols) != fwrite(image, 1, rows * cols, fp))
        {
            fprintf(stderr, "Error writing the image data to %s.\n", outfilename);
            if (fp != stdout)
                fclose(fp);
            return 0;
        }

        if (fp != stdout)
            fclose(fp);
        return 1;
    }
};

SC_MODULE(Platform)
{
    sc_fifo_in<IMAGE> ImgIn;
    sc_fifo_out<IMAGE> ImgOut;
    sc_fifo<IMAGE> q1;
    sc_fifo<IMAGE> q2;
    DataIn data_in;
    DUT dut;
    DataOut data_out;
    SC_CTOR(Platform)
        : q1("q1", 1), q2("q2", 1), data_in("DataIn"), dut("DUT"), data_out("DataOut")
    {
        data_in.ImgIn.bind(ImgIn);
        data_in.ImgOut.bind(q1);
        dut.ImgIn.bind(q1);
        dut.ImgOut.bind(q2);
        data_out.ImgIn.bind(q2);
        data_out.ImgOut.bind(ImgOut);
    }
};

SC_MODULE(Top)
{
    Stimulus stimulus;
    Platform platform;
    Monitor monitor;
    sc_fifo<IMAGE> q1;
    sc_fifo<IMAGE> q2;
    SC_CTOR(Top)
        : q1("q1", 1), q2("q2", 1), stimulus("Stimulus"), platform("Platform"), monitor("Monitor")
    {
        stimulus.ImgOut.bind(q1);
        platform.ImgIn.bind(q1);
        platform.ImgOut.bind(q2);
        monitor.ImgIn.bind(q2);
    }
};

int sc_main(int argc, char *argv[])
{
    Top top("Top");
    sc_start();
    return 0;
}
