/*
 * Buddhabrot program.
 *
 * Created by Finn Hakansson in September-November 2005.
 *
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>

using namespace std;

#define BMP_HEADER_SIZE 54
#define HISTOGRAM_BORDER_WIDTH 0

template<class INT>
class RGB
{
public:
  INT rgb[3];
  INT& r;
  INT& g;
  INT& b;
  RGB() : r(rgb[0]), g(rgb[1]), b(rgb[2]) { rgb[0] = rgb[1] = rgb[2] = 0; }
  INT sum() { return r + g + b; }
  INT& operator[] (int i) { return rgb[i]; }
};

void
CreateBitMapHeader(unsigned char* buf, int width, int height)
{
  int pad = (4 - (width * 3 % 4)) % 4;
  unsigned long size = height * (width * 3 + pad) + BMP_HEADER_SIZE;

  memset(buf, 0, BMP_HEADER_SIZE);
  
  buf[0] = 0x42;
  buf[1] = 0x4d;
  buf[2] = (unsigned char)((size << 24) >> 24);
  buf[3] = (unsigned char)((size << 16) >> 24);
  buf[4] = (unsigned char)((size << 8) >> 24);
  buf[5] = (unsigned char)((size << 0) >> 24);

  buf[10] = 14 + 40;

  buf[14 + 0] = 40;

  buf[14 + 4] = (unsigned char)((width << 24) >> 24);
  buf[14 + 5] = (unsigned char)((width << 16) >> 24);
  buf[14 + 6] = (unsigned char)((width << 8) >> 24);
  buf[14 + 7] = (unsigned char)((width << 0) >> 24);
  buf[14 + 8] = (unsigned char)((height << 24) >> 24);
  buf[14 + 9] = (unsigned char)((height << 16) >> 24);
  buf[14 + 10] = (unsigned char)((height << 8) >> 24);
  buf[14 + 11] = (unsigned char)((height << 0) >> 24);

  buf[14 + 12] = 1; // Number of planes
  buf[14 + 14] = 24;
}


template<class FP>
class MandelbrotSet
{
  FP xmin;
  FP xmax;
  FP ymin;
  FP ymax;
  int xpixels;
  int ypixels;
  FP xspace;
  FP yspace;
  long shallow_iter;
  long maxit;
  bool** conservative_mandelbrot_set;

  long diverge(FP x, FP y);
  bool conservatively_in_set(int xpos, int ypos, int border, long** ms) const;

public:
  MandelbrotSet(FP xmin, FP xmax, FP ymin, FP ymax,
                int xpix, int ypix, FP xspace, FP yspace,
		long shallow_iter, long maxit);
  ~MandelbrotSet();
  void create_conservative_mandelbrot_set();
  bool inSet(FP x, FP y) const;
  void dumpBMPImage(const char* filename);
};

template<class FP>
MandelbrotSet<FP>::MandelbrotSet(FP xmin, FP xmax, FP ymin, FP ymax,
                                 int xpix, int ypix,
                                 FP xspace, FP yspace,
                                 long shallow_iter, long maxit)
{
  this->xmin = xmin;
  this->xmax = xmax;
  this->ymin = ymin;
  this->ymax = ymax;
  this->xpixels = xpix;
  this->ypixels = ypix;
  this->xspace = xspace;
  this->yspace = yspace;
  this->shallow_iter = shallow_iter;
  this->maxit = maxit;

  conservative_mandelbrot_set = new bool* [ypixels];

  for (int y = 0; y < ypixels; y++)
  {
    conservative_mandelbrot_set[y] = new bool[xpixels];
    for (int x = 0; x < xpixels; x++)
    {
      conservative_mandelbrot_set[y][x] = false;
    }
  }

  create_conservative_mandelbrot_set();
}


template<class FP>
MandelbrotSet<FP>::~MandelbrotSet()
{
  for (int y = 0; y < ypixels; y++)
  {
    delete [] conservative_mandelbrot_set[y];
  }
  delete [] conservative_mandelbrot_set;
}


template<class FP> void
MandelbrotSet<FP>::create_conservative_mandelbrot_set()
{
  cout << "About to create conservative Mandelbrot Set." << endl;
  cout << "This may take a few days. ;-)" << endl;

  long** mandelbrot_set = new long* [ypixels];

  for (int y = 0; y < ypixels; y++)
  {
    mandelbrot_set[y] = new long[xpixels];
    memset(mandelbrot_set[y], 0, sizeof(long));
  }

  const int pixel_grid = 10;
  cout << "Starting calculation of Mandelbrot Set." << endl;
  for (int y = 0; y < ypixels / 2; y++)
  {
    FP ypos = ymin + y * yspace;
    float ypercent = y / (ypixels / 2.0) * 100.0;
    cout << "Calculating y = " << y << "  (" << ypercent << "%)" << endl;
    for (int x = 0; x < xpixels; x++)
    {
      FP xpos = xmin + x * xspace;
      for (int yi = 0; yi < pixel_grid; yi++)
      {
        for (int xi = 0; xi < pixel_grid; xi++)
        {
	  // Optimize here: If mandelbrot_set[y][x] already has been set,
	  // don't even bother to calculate if the position diverges.
          if (mandelbrot_set[y][x] /// RARARARARAR
	      || (mandelbrot_set[y][x] = diverge(xpos + xspace / pixel_grid / 2 + xspace / pixel_grid * xi,
                                              ypos + yspace / pixel_grid / 2 + yspace / pixel_grid * yi)))
          {
            xi = pixel_grid;
            yi = pixel_grid;
          }
        }
      }
    }
    // Optimization: Reflect half of all pixels.
    memcpy(mandelbrot_set[ypixels - y - 1], mandelbrot_set[y], sizeof(long) * xpixels);
  }
  cout << "Finished calculation of Mandelbrot Set." << endl;

  // Make Mandelbrot Set conservative.
  int pixels_in_set = 0;
  const int border = 1;
  for (int y = border; y < ypixels - border; y++)
  {
    for (int x = border; x < xpixels - border; x++)
    {
      if (conservatively_in_set(x, y, border, mandelbrot_set))
      {
        conservative_mandelbrot_set[y][x] = true;
        pixels_in_set++;
      }
      else
      {
        conservative_mandelbrot_set[y][x] = false;
      }
    }
  }

  double p = (((double)pixels_in_set) / (xpixels * ypixels)) * 100;
  cout << "Portion of pixels in Mandelbrot Set is " << p << "%." << endl;

  for (int y = 0; y < ypixels; y++)
  {
    delete [] mandelbrot_set[y];
  }
  delete [] mandelbrot_set;
  dumpBMPImage("mandelbrot_set_05_border=1.bmp");
}

template<class FP> bool
MandelbrotSet<FP>::conservatively_in_set(int xpos, int ypos, int border, long** ms) const
{
  bool r = true;
  for (int y = ypos - border; y < ypos + border; y++)
  {
    for (int x = xpos - border; x < xpos + border; x++)
    {
      r = r && (ms[y][x] == 0);
    }
  }
  return r;
}


template<class FP> long
MandelbrotSet<FP>::diverge(FP cr, FP ci)
{
  FP xr = 0;
  FP xi = 0;
  long iter = 0;

  while (iter <= maxit)
  {
    FP x2r = xr * xr - xi * xi + cr;
    FP x2i = 2 * xr * xi + ci;
    iter++;
    xr = x2r;
    xi = x2i;
    if ((xr * xr + xi * xi) >= 4)
    {
      return iter;
    }
  }

  return 0;
}


template<class FP> bool
MandelbrotSet<FP>::inSet(FP x, FP y) const
{
  if (x > xmin && x < xmax && y > ymin && y < ymax)
  {
    int xpix = (int) floor((x - xmin) / xspace);
    int ypix = (int) floor((y - ymin) / yspace);
    return conservative_mandelbrot_set[ypix][xpix];
  }
  else
  {
    cerr << "MandelbrotSet::inSet(" << x;
    cerr << ", " << y << ")" << endl;
    return 0;
  }
}


template<class FP> void
MandelbrotSet<FP>::dumpBMPImage(const char* filename)
{
  unsigned char bmpheader[BMP_HEADER_SIZE];
  int pad = (4 - (xpixels * 3 % 4)) % 4;
  std::ofstream bmpstream(filename, std::ios::out | std::ios::binary);
  if (bmpstream == 0)
  {
    std::cout << "Couldn't open " << filename << endl;
    return;
  }
  CreateBitMapHeader(bmpheader, xpixels, ypixels);
  bmpstream.write((const char*)bmpheader, BMP_HEADER_SIZE);

  for (int y = 0; y < ypixels; y++)
  {
    for (int x = 0; x < xpixels; x++)
    {
      if (!conservative_mandelbrot_set[y][x])
      {
        bmpstream.put((unsigned char)0x00);
        bmpstream.put((unsigned char)0x80);
        bmpstream.put((unsigned char)0xff);
      }
      else
      {
        bmpstream.put((unsigned char)0x00);
        bmpstream.put((unsigned char)0x00);
        bmpstream.put((unsigned char)0x00);
      }
    }
    // Add padding to each row of the picture.
    for (int p = 0; p < pad; p++) {
      bmpstream.put((unsigned char)0x00);
    }
  }
  bmpstream.close();
  cout << "Completed dumping conservative mandelbrot set." << endl;
}


template<class FP, class INT>
class Buddhabrot
{
private:
  FP xmin;
  FP xmax;
  FP ymin;
  FP ymax;
  int xpixels;
  int ypixels;
  FP xspace;
  FP yspace;
  unsigned long long samples;
  unsigned long long maxsamples;
  long miniterations;
  long maxiterations;
  RGB<INT>** counters;
  long min_red;
  long max_red;
  long min_green;
  long max_green;
  long min_blue;
  long max_blue;
  long longest_iteration;
  unsigned long long increments;
  FP cr;
  FP ci;
  FP xstep;
  FP ystep;
  INT histogram[256];
  bool done;
  MandelbrotSet<FP>* m;
  time_t start_time;

  static const int data_version = 1;

  void record(FP x, FP y, long iter);
  long calc_path(bool record_path);
  INT highest_counter(int color) const;
  bool get_sample();
  bool get_sample_random();
  bool in_cardioid();
  RGB<INT>** get_counter_matrix(int x, int y);
  void write_image(std::ofstream& bmpstream);

public:
  Buddhabrot(FP x0, FP x1, FP y0, FP y1, int xpix, int ypix,
             unsigned long long maxsamples,
             long minit, long maxit,
             FP xstep, FP ystep,
             long color_band[6]);
  Buddhabrot(const char* filename);
  ~Buddhabrot();
  void run();
  void make_image(const char * imagefile);
  void make_datafile(const char * datafile);
  void dump(const char * dumpfile);
  void make_histogram(const char * histfile);
};


template<class FP, class INT>
Buddhabrot<FP, INT>::Buddhabrot(FP x0, FP x1, FP y0, FP y1,
                                int xpix, int ypix,
                                unsigned long long maxsamp,
                                long minit, long maxit,
                                FP xs, FP ys, long color_band[6])
{
  xmin = x0;
  xmax = x1;
  ymin = y0;
  ymax = y1;
  xpixels = xpix;
  ypixels = ypix;
  xspace = (xmax - xmin) / xpixels;
  yspace = (ymax - ymin) / ypixels;
  samples = 0;
  maxsamples = maxsamp;
  miniterations = minit;
  maxiterations = maxit;
  cr = xmin + ((FP)rand() / RAND_MAX) * xspace * 0.1;
  ci = ymin + ((FP)rand() / RAND_MAX) * yspace * 0.021;
  xstep = xs;
  ystep = ys;
  increments = 0;
  longest_iteration = 0;
  done = false;
  time(&start_time);
  min_red = color_band[0];
  max_red = color_band[1];
  min_green = color_band[2];
  max_green = color_band[3];
  min_blue = color_band[4];
  max_blue = color_band[5];

  // Create a matrix to record the escape paths.
  counters = get_counter_matrix(xpixels, ypixels);

  for (int i = 0; i < 256; i++)
  {
    histogram[i] = 0;
  }

  m = new MandelbrotSet<FP>(xmin, xmax, ymin, ymax,
                            xpixels, ypixels, xspace, yspace, 500, 500000);

  time_t mandel_time;
  time(&mandel_time);
  long elapsed_time = (long)difftime(mandel_time, start_time);
  cout << "The Mandelbrot Set took " << elapsed_time;
  cout << " seconds to create." << endl;
}


template<class FP, class INT>
Buddhabrot<FP, INT>::Buddhabrot(const char* filename)
{
  std::cout << "About to create Buddhabrot object from file \"" << filename << "\"" << endl;

  done = false;
  int data_version = 0;
  std::ifstream datafile(filename, std::ios::in);
  datafile >> data_version;
  datafile >> xmin;
  datafile >> xmax;
  datafile >> ymin;
  datafile >> ymax;
  datafile >> xpixels;
  datafile >> ypixels;

  xspace = (xmax - xmin) / xpixels;
  yspace = (ymax - ymin) / ypixels;

  datafile >> samples;
  datafile >> maxsamples;
  datafile >> miniterations;
  datafile >> maxiterations;

  datafile >> cr;
  datafile >> ci;

  // The cr and ci could have wrapped around and we have
  // to set them straight.
  if (cr < xmin || cr > xmax)
  {
    cr = xmin + ((FP)rand() / RAND_MAX) * xspace * 0.0689;
  }
  if (ci < ymin || ci > ymax)
  {
    ci = ymin + ((FP)rand() / RAND_MAX) * yspace * 0.0219;
  }

  datafile >> xstep;
  datafile >> ystep;

  datafile >> longest_iteration;
  datafile >> increments;
  datafile >> min_red;
  datafile >> max_red;
  datafile >> min_green;
  datafile >> max_green;
  datafile >> min_blue;
  datafile >> max_blue;

  counters = get_counter_matrix(xpixels, ypixels);

  for (int y = 0; y < ypixels; y++)
  {
    for (int x = 0; x < xpixels; x++)
    {
      datafile >> counters[y][x].r;
      datafile >> counters[y][x].g;
      datafile >> counters[y][x].b;
    }
  }

  for (int i = 0; i <= 0xff; i++)
  {
    histogram[i] = 0;
  }

  datafile.close();
  time(&start_time);

  m = new MandelbrotSet<FP>(xmin, xmax, ymin, ymax,
                            xpixels, ypixels, xspace, yspace, 500000);

  time_t mandel_time;
  time(&mandel_time);
  long elapsed_time = (long)difftime(mandel_time, start_time);
  cout << "The Mandelbrot Set took " << elapsed_time;
  cout << " seconds to create." << endl;
  cout << " cr = " << cr << "  ci = " << ci << endl;
}


template<class FP, class INT>
Buddhabrot<FP, INT>::~Buddhabrot()
{
  for (int i = 0; i < ypixels; i++)
  {
    delete [] counters[i];
  }
  delete [] counters;
  delete m;

  time_t t;
  time(&t);
  long diff = (long)difftime(t, start_time);
  long hours = diff / 3600;
  long minutes = (diff - hours * 3600) / 60;
  long seconds = diff % 60;
  cout << "The entire calculation took " << hours << "h";
  cout << minutes << "m" << seconds << "s." << endl;
}

template<class FP, class INT> RGB<INT>**
Buddhabrot<FP, INT>::get_counter_matrix(int x, int y)
{
  RGB<INT>** counters = new RGB<INT>* [y];
  for (int i = 0; i < y; i++)
  {
    counters[i] = new RGB<INT>[x];
  }
  return counters;
}


template<class FP, class INT> void
Buddhabrot<FP, INT>::record(FP x, FP y, long iter)
{
  if (x > xmin && x < xmax && y > ymin && y < ymax)
  {
    int xpix = (int) floor((x - xmin) / xspace);
    int ypix = (int) floor((y - ymin) / yspace);

    if (iter >= min_red && iter <= max_red)
    {
      counters[ypix][xpix].r++;
    }
    else if (iter >= min_green && iter <= max_green)
    {
      counters[ypix][xpix].g++;
    }
    else if (iter >= min_blue && iter <= max_blue)
    {
      counters[ypix][xpix].b++;
    }

    increments++;
  }
}
    
template<class FP, class INT> long
Buddhabrot<FP, INT>::calc_path(bool record_path)
{
  FP xr = 0;
  FP xi = 0;
  long iterations = 0;
  
  while (iterations < maxiterations)
  {
    FP x2r = xr * xr - xi * xi + cr;
    FP x2i = 2 * xr * xi + ci;
    iterations++;
    xr = x2r;
    xi = x2i;
    if (record_path && iterations >= miniterations)
    {
      record(xr, xi, iterations);
    }
    if ((xr * xr + xi * xi) >= 4)
    {
      if (iterations > longest_iteration)
      {
        longest_iteration = iterations;
        cout << "longest iteration: ";
        cout << longest_iteration;
        cout << "  (" << cr << ", " << ci << ")" << endl;
      }
      return iterations; // We're done.
    }
  }

  return 0; // Doesn't escape.
}
    
template<class FP, class INT> void
Buddhabrot<FP, INT>::run()
{
  time_t t;
  time(&t);

  while (!done && get_sample())
  {
    long it = calc_path(false);
    if (it > 0 && it >= miniterations)
    {
      // Record all points in the escape path.
      calc_path(true);
      samples++;
      if (samples % 1000 == 0)
      {
        time_t d;
        time(&d);
        long diff = (long)difftime(d, t);
        t = d;
        //double progress = ((double)samples) / maxsamples * 100;
        cout << "samples: " << samples << "   ";
        //cout << progress << "%)     ";
        cout << "cr = " << cr << "   ci = " << ci;
        cout << " diff=" << diff << "   ";
        cout << " longest_iteration: " << longest_iteration << endl;
      }
    }
  }

  cout << "last sample number = " << samples << endl;
  cout << "last cr = " << cr << "    last ci = " << ci << endl;
}
    
template<class FP, class INT> INT
Buddhabrot<FP, INT>::highest_counter(int color) const
{
  INT h = 0;
  for (int y = 0; y < ypixels; y++)
  {
    for (int x = 0; x < xpixels; x++)
    {
      INT c = counters[y][x][color];
      if (c > h)
      {
        h = c;
      }
    }
  }
  return h;
}

template<class FP, class INT> bool
Buddhabrot<FP, INT>::get_sample()
{
  for (;;)
  {
    cr += xspace * xstep;
    if (cr > xmax)
    {
      cr -= (xmax - xmin);
      ci += yspace * ystep;
      if (ci > ymax)
      {
        //ci -= (ymax - ymin);
        cout << "ending iterations" << endl;
        done = true;
        return false;
      }
    }
    if (m->inSet(cr, ci))
    {
      // Point is in the Mandelbrot Set; get new point.
    }
    else
    {
      return true;
    }
  }
}


template<class FP, class INT> bool
Buddhabrot<FP, INT>::get_sample_random()
{
  cr = xmin + (xmax - xmin) * (((FP)rand()) / RAND_MAX);
  ci = ymin + (ymax - ymin) * (((FP)rand()) / RAND_MAX);
  return true;
}


template<class FP, class INT> bool
Buddhabrot<FP, INT>::in_cardioid()
{
  FP p = sqrt((cr - 0.25) * (cr - 0.25) + ci * ci);
  FP a = atan2(ci, cr - 0.25);
  FP pc = 0.5 - 0.5 * cos(a);
  return (p <= pc);
}


template<class FP, class INT> void
Buddhabrot<FP, INT>::make_datafile(const char *data_filename)
{
  std::ofstream datastream(data_filename, std::ios::out);
  if (datastream == 0)
  {
    std::cout << "Couldn't open " << data_filename << endl;
    return;
  }
  for (int y = 0; y < ypixels; y++)
  {
    for (int x = 0; x < xpixels; x++)
    {
      datastream << counters[y][x].r << " "
                 << counters[y][x].g << " "
                 << counters[y][x].b << " ";
    }
    datastream << endl;
  }
  datastream.close();
}


template<class FP, class INT> void
Buddhabrot<FP, INT>::make_image(const char *image_filename)
{
  unsigned char bmpheader[BMP_HEADER_SIZE];
  CreateBitMapHeader(bmpheader, xpixels, ypixels);

  std::ofstream bmpstream(image_filename, std::ios::out | std::ios::binary);
  if (bmpstream == 0)
  {
    std::cout << "Couldn't open " << image_filename << endl;
  }
  else
  {
    bmpstream.write((const char*)bmpheader, BMP_HEADER_SIZE);
    write_image(bmpstream);
    bmpstream.close();
  }
}

    
template<class FP, class INT> void
Buddhabrot<FP, INT>::write_image(std::ofstream& bmpstream)
{
  FP highest_red = (FP) highest_counter(0);
  FP highest_green = (FP) highest_counter(1);
  FP highest_blue = (FP) highest_counter(2);
  FP highest[3] = { highest_red / 3, highest_green / 3, highest_blue / 3 };
  int highest_intensity = 0;
  int pad = (4 - (xpixels * 3 % 4)) % 4;

  std::cout << "highest red = " << highest_red << endl;
  std::cout << "highest green = " << highest_green << endl;
  std::cout << "highest blue = " << highest_blue << endl;
  std::cout << "longest iteration: " << longest_iteration << endl;
  std::cout << "increments = " << increments << endl;
  //std::cout << "maxsamples = " << maxsamples << endl;

  for (int y = 0; y < ypixels; y++)
  {
    for (int x = 0; x < xpixels; x++)
    {
      for (int bgr = 2; bgr >= 0; bgr--)
      {
        int intensity = (int)((counters[y][x][bgr] / highest[bgr]) * 0xff);
        if (intensity > highest_intensity)
        {
          highest_intensity = intensity;
        }
        if (intensity > 0xff)
        {
          intensity = 0xff;
        }
        if (intensity < 0)
        {
          intensity = 0;
        }
        histogram[intensity]++;

        bmpstream.put((unsigned char)intensity);
      }
      for (int p = 0; p < pad; p++) {
        bmpstream.put((unsigned char)0x00);
      }
    }
  }

  std::cout << "highest_intensity = " << highest_intensity << endl;
}

template<class FP, class INT> void
Buddhabrot<FP, INT>::dump(const char* dumpfile)
{
  std::ofstream dumpstream(dumpfile, std::ios::out);
  if (dumpstream == 0)
  {
    std::cout << "Couldn't open " << dumpfile << endl;
    return;
  }
  dumpstream << data_version << endl;
  dumpstream << xmin << endl;
  dumpstream << xmax << endl;
  dumpstream << ymin << endl;
  dumpstream << ymax << endl;
  dumpstream << xpixels << endl;
  dumpstream << ypixels << endl;
  dumpstream << samples << endl;
  dumpstream << maxsamples << endl;
  dumpstream << miniterations << endl;
  dumpstream << maxiterations << endl;
  dumpstream << cr << endl;
  dumpstream << ci << endl;
  dumpstream << xstep << endl;
  dumpstream << ystep << endl;
  dumpstream << longest_iteration << endl;
  dumpstream << increments << endl;
  dumpstream << min_red << endl;
  dumpstream << max_red << endl;
  dumpstream << min_green << endl;
  dumpstream << max_green << endl;
  dumpstream << min_blue << endl;
  dumpstream << max_blue << endl;

  for (int y = 0; y < ypixels; y++)
  {
    for (int x = 0; x < xpixels; x++)
    {
      for (int rgb = 0; rgb < 3; rgb++)
      {
        dumpstream << counters[y][x][rgb] << " ";
      }
    }
    dumpstream << endl;
  }
  dumpstream.close();
}

template<class FP, class INT> void
Buddhabrot<FP, INT>::make_histogram(const char *histfile)
{
  int width = 2 * HISTOGRAM_BORDER_WIDTH + 256;
  int height = 2 * HISTOGRAM_BORDER_WIDTH + 256;
  unsigned char bmpheader[BMP_HEADER_SIZE];
  INT highest_in_histogram = 0;

  cout << "Creating histogram in file " << histfile << endl;

  // Normalize histogram
  for (int x = 0; x < 256; x++)
  {
    //cout << " " << histogram[x] << ", ";
    if (histogram[x] > highest_in_histogram)
    {
      highest_in_histogram = histogram[x];
    }
  }
  //cout << "\n";
  for (int x = 0; x <= 0xff; x++)
  {
    histogram[x] = (INT)((height * histogram[x]) / highest_in_histogram);
    //cout << " " << histogram[x] << ", ";
  }

  std::ofstream bmpstream(histfile, std::ios::out | std::ios::binary);
  if (bmpstream == 0)
  {
    std::cout << "Couldn't open " << histfile << endl;
    return;
  }
  CreateBitMapHeader(bmpheader, width, height);
  bmpstream.write((const char*)bmpheader, BMP_HEADER_SIZE);

  for (int y = 0; y <= 0xff; y++)
  {
    for (int x = 0; x <= 0xff; x++)
    {
      if (histogram[x] >= y)
      {
        bmpstream.put((unsigned char)0x00);
        bmpstream.put((unsigned char)0x00);
        bmpstream.put((unsigned char)0xff);
      }
      else
      {
        bmpstream.put((unsigned char)0x00);
        bmpstream.put((unsigned char)0x00);
        bmpstream.put((unsigned char)0x00);
      }
    }
  }
  bmpstream.close();
  //cout << endl;
}


int main(int argc, char** argv)
{
  long double xmin = -2.03;
  long double xmax = 1.12;
  long double ymin = -1.04;
  long double ymax = 1.04;
  int xpix = 3000;
  int ypix = 2000;
  long minit = 4;
  long maxit = 1000000;
  long color_band[6] = {minit, 26, 27, 250, 251, maxit};
  long double xstep = 0.08;
  long double ystep = 0.08;
  unsigned long long maxsamples = 1000;
  long filesize = xpix * ypix * 3 + BMP_HEADER_SIZE;
  const char* imagefile = "b3000x2000_06.bmp";
  const char* datafile = "b3000x2000_06.data";
  const char *histfile = "hist_06.bmp";
  const char *dumpfile = "dump_06.txt";
  char* arg;

  while ((arg = *argv++))
  {
    if (!strncmp(arg, "--maxit=", strlen("--maxit=")))
    {
      maxit = atoi(&arg[strlen("--maxit=")]);
      cout << "maxit = " << maxit << endl;
    }
    if (!strncmp(arg, "--xpix=", strlen("--xpix=")))
    {
      maxit = atoi(&arg[strlen("--xpix=")]);
      cout << "xpix = " << xpix << endl;
    }
    if (!strncmp(arg, "--ypix=", strlen("--ypix=")))
    {
      maxit = atoi(&arg[strlen("--ypix=")]);
      cout << "ypix = " << ypix << endl;
    }
  }

  //cout << "Starting at: " << localtime(0) << "\r\n";
  std::cout << "The bmp file is going to be of size " << filesize << endl;
  //Buddhabrot<long double, unsigned long long> b("dump05.txt");
  Buddhabrot<long double, unsigned long long> b(xmin, xmax, ymin, ymax, xpix, ypix, maxsamples, minit, maxit, xstep, ystep, color_band);
  b.run();

  b.make_image(imagefile);
  b.make_datafile(datafile);

  // Make histogram.
  b.make_histogram(histfile);
  // Dump all data to file.
  b.dump(dumpfile);

  return 0;
}

// Portion of pixels in Mandelbrot Set is 22.3356%.
