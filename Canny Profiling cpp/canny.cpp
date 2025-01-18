/*******************************************************************************
 * PROGRAM: canny_edge
 * PURPOSE: This program implements a "Canny" edge detector with a refined
 * hierarchical structure suitable for SystemC modeling and profiling.
 *
 * The SystemC model is structured with the following modules:
 *
 *   - Top
 *     - Stimulus
 *     - Platform
 *       - DataIn
 *       - DUT (Canny Edge Detector)
 *         - Gaussian_Smooth
 *         - Derivative_X_Y
 *         - Magnitude_X_Y
 *         - Non_Max_Supp
 *         - Apply_Hysteresis
 *       - DataOut
 *     - Monitor
 *
 * NAME: Yash Ingle
 * DATE:11/11/2024
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

// Constants for the image dimensions
#define COLS 2704
#define ROWS 1520
#define SIZE (COLS * ROWS)

// Define constants for edge detection (white background, black edges)
#define NOEDGE 255        // White background
#define POSSIBLE_EDGE 128 // Gray (intermediate value)
#define EDGE 0            // Black edges

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Constants for the algorithm

const float sigma = 1.0;
const float tlow = 0.05;
const float thigh = 0.05;
const int KERNEL_SIZE = 21;
const float BOOSTBLURFACTOR = 90.0;

// Define the IMAGE struct
struct IMAGE
{
    unsigned char img[SIZE];

    IMAGE()
    {
        memset(img, NOEDGE, SIZE * sizeof(unsigned char));
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

// Define the SIMAGE struct for short int images
struct SIMAGE
{
    short int img[SIZE];

    SIMAGE()
    {
        memset(img, 0, SIZE * sizeof(short int));
    }

    SIMAGE &operator=(const SIMAGE &copy)
    {
        if (this != &copy)
        {
            memcpy(img, copy.img, SIZE * sizeof(short int));
        }
        return *this;
    }

    operator short int *()
    {
        return img;
    }

    short int &operator[](const int index)
    {
        return img[index];
    }
};

// Define the TEMP_IMAGE struct
struct TEMP_IMAGE
{
    float data[SIZE];

    TEMP_IMAGE()
    {
        memset(data, 0, SIZE * sizeof(float));
    }

    TEMP_IMAGE &operator=(const TEMP_IMAGE &copy)
    {
        if (this != &copy)
        {
            memcpy(data, copy.data, SIZE * sizeof(float));
        }
        return *this;
    }

    float &operator[](const int index)
    {
        return data[index];
    }
};

// Stimulus Module
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
        if (img_rows != ROWS || img_cols != COLS)
        {
            fprintf(stderr, "Image dimensions (%d x %d) do not match expected dimensions (%d x %d).\n",
                    img_cols, img_rows, COLS, ROWS);
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
        if ((size_t)(ROWS * COLS) != fread(image, 1, ROWS * COLS, fp))
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

// Monitor Module
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
        int cols = COLS;
        int rows = ROWS;

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

// DataIn Module
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

// DataOut Module
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

// Gaussian_Smooth Module
SC_MODULE(Gaussian_Smooth)
{
    sc_fifo_in<IMAGE> ImgIn;
    sc_fifo_out<SIMAGE> ImgOut;

    SC_CTOR(Gaussian_Smooth)
    {
        SC_THREAD(process);
        set_stack_size(256 * 1024 * 1024);
    }

    void process()
    {
        while (true)
        {
            IMAGE image_in;
            SIMAGE smoothedim;
            float kernel[KERNEL_SIZE];
            int windowsize;
            TEMP_IMAGE tempim;

            ImgIn.read(image_in);

            gaussian_kernel(sigma, kernel, &windowsize);
            blur_x(image_in.img, tempim.data, kernel, windowsize);
            blur_y(tempim.data, smoothedim.img, kernel, windowsize);

            ImgOut.write(smoothedim);
        }
    }

    void gaussian_kernel(float sigma, float kernel[], int *windowsize) __attribute__((noinline));
    void blur_x(unsigned char image[], float tempim[], float kernel[], int windowsize) __attribute__((noinline));
    void blur_y(float tempim[], short int smoothedim[], float kernel[], int windowsize) __attribute__((noinline));
};

// Definition of Gaussian_Smooth member functions
void Gaussian_Smooth::gaussian_kernel(float sigma, float kernel[], int *windowsize)
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
}

void Gaussian_Smooth::blur_x(unsigned char image[], float tempim[], float kernel[], int windowsize)
{
    int rows = ROWS;
    int cols = COLS;
    int r, c, cc, center;
    float dot, sum;

    center = windowsize / 2;

    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            dot = 0.0;
            sum = 0.0;
            for (cc = -center; cc <= center; cc++)
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
}

void Gaussian_Smooth::blur_y(float tempim[], short int smoothedim[], float kernel[], int windowsize)
{
    int rows = ROWS;
    int cols = COLS;
    int r, c, rr, center;
    float dot, sum;

    center = windowsize / 2;

    for (c = 0; c < cols; c++)
    {
        for (r = 0; r < rows; r++)
        {
            dot = 0.0;
            sum = 0.0;
            for (rr = -center; rr <= center; rr++)
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
}

// Derivative_X_Y Module
SC_MODULE(Derivative_X_Y)
{
    sc_fifo_in<SIMAGE> ImgIn;
    sc_fifo_out<SIMAGE> DeltaXOut_Mag;
    sc_fifo_out<SIMAGE> DeltaYOut_Mag;
    sc_fifo_out<SIMAGE> DeltaXOut_NMS;
    sc_fifo_out<SIMAGE> DeltaYOut_NMS;

    SC_CTOR(Derivative_X_Y)
    {
        SC_THREAD(process_derivative);
        set_stack_size(256 * 1024 * 1024);
    }

    void process_derivative()
    {
        while (true)
        {
            SIMAGE smoothedim;
            SIMAGE delta_x, delta_y;
            ImgIn.read(smoothedim);
            derivative_x_y(smoothedim.img, delta_x.img, delta_y.img);

            // Write to both outputs
            DeltaXOut_Mag.write(delta_x);
            DeltaYOut_Mag.write(delta_y);
            DeltaXOut_NMS.write(delta_x);
            DeltaYOut_NMS.write(delta_y);
        }
    }

    void derivative_x_y(short int smoothedim[], short int delta_x[], short int delta_y[]) __attribute__((noinline));
};

// Definition of Derivative_X_Y member function
void Derivative_X_Y::derivative_x_y(short int smoothedim[], short int delta_x[], short int delta_y[])
{
    int rows = ROWS;
    int cols = COLS;
    int r, c, pos;

    // Compute the x-derivative
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

// Magnitude_X_Y Module
SC_MODULE(Magnitude_X_Y)
{
    sc_fifo_in<SIMAGE> DeltaXIn;
    sc_fifo_in<SIMAGE> DeltaYIn;
    sc_fifo_out<SIMAGE> MagnitudeOut_Mag;
    sc_fifo_out<SIMAGE> MagnitudeOut_Hyst;

    SC_CTOR(Magnitude_X_Y)
    {
        SC_THREAD(process_magnitude);
        set_stack_size(256 * 1024 * 1024);
    }

    void process_magnitude()
    {
        while (true)
        {
            SIMAGE delta_x, delta_y;
            SIMAGE magnitude;
            DeltaXIn.read(delta_x);
            DeltaYIn.read(delta_y);
            magnitude_x_y(delta_x.img, delta_y.img, magnitude.img);

            // Write magnitude to both outputs
            MagnitudeOut_Mag.write(magnitude);
            MagnitudeOut_Hyst.write(magnitude);
        }
    }

    void magnitude_x_y(short int delta_x[], short int delta_y[], short int magnitude[]) __attribute__((noinline));
};

// Definition of Magnitude_X_Y member function
void Magnitude_X_Y::magnitude_x_y(short int delta_x[], short int delta_y[], short int magnitude[])
{
    int rows = ROWS;
    int cols = COLS;
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

// Non_Max_Supp Module
SC_MODULE(Non_Max_Supp)
{
    sc_fifo_in<SIMAGE> MagnitudeIn;
    sc_fifo_in<SIMAGE> DeltaXIn;
    sc_fifo_in<SIMAGE> DeltaYIn;
    sc_fifo_out<IMAGE> NMSOut;

    SC_CTOR(Non_Max_Supp)
    {
        SC_THREAD(process_non_max_supp);
        set_stack_size(256 * 1024 * 1024);
    }

    void process_non_max_supp()
    {
        while (true)
        {
            SIMAGE magnitude, delta_x, delta_y;
            IMAGE nms;
            MagnitudeIn.read(magnitude);
            DeltaXIn.read(delta_x);
            DeltaYIn.read(delta_y);
            non_max_supp(magnitude.img, delta_x.img, delta_y.img, nms.img);
            NMSOut.write(nms);
        }
    }

    void non_max_supp(short *mag, short *gradx, short *grady, unsigned char *result) __attribute__((noinline));
};

// Definition of Non_Max_Supp member function
void Non_Max_Supp::non_max_supp(short *mag, short *gradx, short *grady, unsigned char *result)
{
    int rows = ROWS;
    int cols = COLS;
    int rowcount, colcount;
    int idx;

    // Initialize result to NOEDGE (white background)
    memset(result, NOEDGE, rows * cols * sizeof(unsigned char));

    // Suppress non-maximum points
    for (rowcount = 1; rowcount < rows - 1; rowcount++)
    {
        for (colcount = 1; colcount < cols - 1; colcount++)
        {
            idx = rowcount * cols + colcount;
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
                result[idx] = POSSIBLE_EDGE; // Gray value
            else
                result[idx] = NOEDGE; // White background
        }
    }
}

// Apply_Hysteresis Module
SC_MODULE(Apply_Hysteresis)
{
    sc_fifo_in<SIMAGE> MagnitudeIn;
    sc_fifo_in<IMAGE> NMSIn;
    sc_fifo_out<IMAGE> EdgeOut;

    SC_CTOR(Apply_Hysteresis)
    {
        SC_THREAD(process_apply_hysteresis);
        set_stack_size(256 * 1024 * 1024);
    }

    void process_apply_hysteresis()
    {
        while (true)
        {
            SIMAGE magnitude;
            IMAGE nms, edge;
            MagnitudeIn.read(magnitude);
            NMSIn.read(nms);
            apply_hysteresis(magnitude.img, nms.img, edge.img);
            EdgeOut.write(edge);
        }
    }

    void apply_hysteresis(short int mag[], unsigned char nms[], unsigned char edge[]) __attribute__((noinline));
    void follow_edges(unsigned char edge[], short mag[], int pos, short lowval) __attribute__((noinline));
};

// Definition of Apply_Hysteresis member functions
void Apply_Hysteresis::apply_hysteresis(short int mag[], unsigned char nms[], unsigned char edge[])
{
    int rows = ROWS;
    int cols = COLS;
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

    // Initialize edge image with NOEDGE (white background)
    memset(edge, NOEDGE, rows * cols * sizeof(unsigned char));

    // Perform edge tracking by hysteresis
    for (r = 0, pos = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++, pos++)
        {
            if ((nms[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold))
            {
                edge[pos] = EDGE; // Black edge
                follow_edges(edge, mag, pos, lowthreshold);
            }
        }
    }

    // Set remaining possible edges to NOEDGE
    for (pos = 0; pos < rows * cols; pos++)
    {
        if (edge[pos] != EDGE)
            edge[pos] = NOEDGE; // White background
    }
}

void Apply_Hysteresis::follow_edges(unsigned char edge[], short mag[], int pos, short lowval)
{
    int cols = COLS;
    int rows = ROWS;
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
                edge[new_pos] = EDGE; // Black edge
                follow_edges(edge, mag, new_pos, lowval);
            }
        }
    }
}

// DUT Module
SC_MODULE(DUT)
{
    sc_fifo_in<IMAGE> ImgIn;
    sc_fifo_out<IMAGE> ImgOut;

    // Internal FIFOs
    sc_fifo<SIMAGE> fifo_gauss_deriv;
    sc_fifo<SIMAGE> fifo_delta_x_mag;
    sc_fifo<SIMAGE> fifo_delta_y_mag;
    sc_fifo<SIMAGE> fifo_delta_x_nms;
    sc_fifo<SIMAGE> fifo_delta_y_nms;
    sc_fifo<SIMAGE> fifo_magnitude_mag;
    sc_fifo<SIMAGE> fifo_magnitude_hyst;
    sc_fifo<IMAGE> fifo_nms;

    // Module instances
    Gaussian_Smooth gaussian_smooth;
    Derivative_X_Y derivative_x_y;
    Magnitude_X_Y magnitude_x_y;
    Non_Max_Supp non_max_supp;
    Apply_Hysteresis apply_hysteresis;

    SC_CTOR(DUT)
        : fifo_gauss_deriv(1),
          fifo_delta_x_mag(1),
          fifo_delta_y_mag(1),
          fifo_delta_x_nms(1),
          fifo_delta_y_nms(1),
          fifo_magnitude_mag(1),
          fifo_magnitude_hyst(1),
          fifo_nms(1),
          gaussian_smooth("Gaussian_Smooth"),
          derivative_x_y("Derivative_X_Y"),
          magnitude_x_y("Magnitude_X_Y"),
          non_max_supp("Non_Max_Supp"),
          apply_hysteresis("Apply_Hysteresis")
    {
        // Connect modules
        gaussian_smooth.ImgIn(ImgIn);
        gaussian_smooth.ImgOut(fifo_gauss_deriv);

        derivative_x_y.ImgIn(fifo_gauss_deriv);
        derivative_x_y.DeltaXOut_Mag(fifo_delta_x_mag);
        derivative_x_y.DeltaYOut_Mag(fifo_delta_y_mag);
        derivative_x_y.DeltaXOut_NMS(fifo_delta_x_nms);
        derivative_x_y.DeltaYOut_NMS(fifo_delta_y_nms);

        magnitude_x_y.DeltaXIn(fifo_delta_x_mag);
        magnitude_x_y.DeltaYIn(fifo_delta_y_mag);
        magnitude_x_y.MagnitudeOut_Mag(fifo_magnitude_mag);
        magnitude_x_y.MagnitudeOut_Hyst(fifo_magnitude_hyst);

        non_max_supp.MagnitudeIn(fifo_magnitude_mag);
        non_max_supp.DeltaXIn(fifo_delta_x_nms);
        non_max_supp.DeltaYIn(fifo_delta_y_nms);
        non_max_supp.NMSOut(fifo_nms);

        apply_hysteresis.MagnitudeIn(fifo_magnitude_hyst);
        apply_hysteresis.NMSIn(fifo_nms);
        apply_hysteresis.EdgeOut(ImgOut);
    }
};

// Platform Module
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

// Top Module
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

// Main Function
int sc_main(int argc, char *argv[])
{
    Top top("Top");
    sc_start();
    return 0;
}
