import java.awt.image.*;
public class ImageFFT {
   public static final int NO_WINDOW = 1;
  public static final int BARTLETT_WINDOW = 2;
  public static final int HAMMING_WINDOW = 3;
  public static final int HANNING_WINDOW = 4;
  private static final String NO_DATA = "no spectral data available";
  private static final String INVALID_PARAMS = "invalid filter parameters";
  private static final double TWO_PI = 2.0*Math.PI;
  
  /** Complex storage for results of FFT. */
  private Complex[] data;

  /** base-2 logarithm of transform width. */
  private int log2w;

  /** base-2 logarithm of transform height. */
  private int log2h;

  /** Width of transform. */
  private int width;

  /** Height of transform. */
  private int height;

  /** Windowing function applied to image data. */
  private int window;

  /** Indicates whether we have spectral or spatial data. */
  private boolean spectral = false;


  public static final double bartlettWindow(double r, double rmax) {
    return 1.0 - Math.min(r, rmax)/rmax;
  }


  public static final double hammingWindow(double r, double rmax) {
    double f = (rmax - Math.min(r, rmax)) / rmax;
    return 0.54 - 0.46*Math.cos(f*Math.PI);
  }

  public static final double hanningWindow(double r, double rmax) {
    double f = (rmax - Math.min(r, rmax)) / rmax;
    return 0.5 - 0.5*Math.cos(f*Math.PI);
  }
  
  public static final double butterworthLowPassFunction(
   int n, double radius, double r) {
    double p = Math.pow(r/radius, 2.0*n);
    return 1.0/(1.0 + p);
  }


  public static final double butterworthHighPassFunction(
   int n, double radius, double r) {
    try {
      double p = Math.pow(radius/r, 2.0*n);
      return 1.0/(1.0 + p);
    }
    catch (ArithmeticException e) {
      return 0.0;
    }
  }


  public static final double butterworthBandPassFunction(
   int n, double radius, double delta, double r) {
    return 1.0-butterworthBandStopFunction(n, radius, delta, r);
  }

  public static final double butterworthBandStopFunction(
   int n, double radius, double delta, double r) {
    try {
      double p = Math.pow(delta*radius/(r*r - radius*radius), 2.0*n);
      return 1.0/(1.0 + p);
    }
    catch (ArithmeticException e) {
      return 0.0;
    }
  }

  public ImageFFT(BufferedImage image) throws FFTException {
    this(image, NO_WINDOW);
  }
  
  public ImageFFT(BufferedImage image, int win) throws FFTException {

    if (image.getType() != BufferedImage.TYPE_BYTE_GRAY)
      throw new FFTException("image must be 8-bit greyscale");

    // Compute dimensions, allowing for zero padding

    log2w = powerOfTwo(image.getWidth());
    log2h = powerOfTwo(image.getHeight());
    width = 1 << log2w;
    height = 1 << log2h;
    window = win;

    // Allocate storage for results of FFT

    data = new Complex[width*height];
    for (int i = 0; i < data.length; ++i)
      data[i] = new Complex();

    Raster raster = image.getRaster();
    double xc = image.getWidth()/2.0, yc = image.getHeight()/2.0;
    double r, rmax = Math.min(xc, yc);
    switch (window) {

      case HAMMING_WINDOW:
        for (int y = 0; y < image.getHeight(); ++y)
          for (int x = 0; x < image.getWidth(); ++x) {
            r = Math.sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
            data[y*width+x].re =
             (float) (hammingWindow(r, rmax)*raster.getSample(x, y, 0));
          }
        break;

      case HANNING_WINDOW:
        for (int y = 0; y < image.getHeight(); ++y)
          for (int x = 0; x < image.getWidth(); ++x) {
            r = Math.sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
            data[y*width+x].re =
             (float) (hanningWindow(r, rmax)*raster.getSample(x, y, 0));
          }
        break;

      case BARTLETT_WINDOW:
        for (int y = 0; y < image.getHeight(); ++y)
          for (int x = 0; x < image.getWidth(); ++x) {
            r = Math.sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
            data[y*width+x].re =
             (float) (bartlettWindow(r, rmax)*raster.getSample(x, y, 0));
          }
        break;

      default:  // NO_WINDOW
        for (int y = 0; y < image.getHeight(); ++y)
          for (int x = 0; x < image.getWidth(); ++x)
            data[y*width+x].re = raster.getSample(x, y, 0);
        break;

    }

  }


  /**
   * @return width of FFT.
   */

  public int getWidth() {
    return width;
  }


  /**
   * @return height of FFT.
   */

  public int getHeight() {
    return height;
  }


  /**
   * @return current windowing function.
   */

  public int getWindow() {
    return window;
  }


  /**
   * @return true if data are spectral, false if data are spatial.
   */

  public boolean isSpectral() {
    return spectral;
  }


  /**
   * @return information string for an ImageFFT object.
   */

  public String toString() {
    String s = new String(getClass().getName() + ": " + width + "x" + height +
     (spectral ? ", frequency domain" : ", spatial domain"));
    return s;
  }
 
  public void transform() {

    int x, y, i;
    Complex[] row = new Complex[width];
    for (x = 0; x < width; ++x)
      row[x] = new Complex();
    Complex[] column = new Complex[height];
    for (y = 0; y < height; ++y)
      column[y] = new Complex();

    int direction;
    if (spectral)
      direction = -1;   // inverse transform
    else
      direction = 1;    // forward transform

    // Perform FFT on each row

    for (y = 0; y < height; ++y) {
      for (i = y*width, x = 0; x < width; ++x, ++i) {
        row[x].re = data[i].re;
        row[x].im = data[i].im;
      }
      reorder(row, width);
      fft(row, width, log2w, direction);
      for (i = y*width, x = 0; x < width; ++x, ++i) {
        data[i].re = row[x].re;
        data[i].im = row[x].im;
      }
    }

    // Perform FFT on each column

    for (x = 0; x < width; ++x) {
      for (i = x, y = 0; y < height; ++y, i += width) {
        column[y].re = data[i].re;
        column[y].im = data[i].im;
      }
      reorder(column, height);
      fft(column, height, log2h, direction);
      for (i = x, y = 0; y < height; ++y, i += width) {
        data[i].re = column[y].re;
        data[i].im = column[y].im;
      }
    }

    if (spectral)
      spectral = false;
    else
      spectral = true;

  }

  public BufferedImage toImage(BufferedImage image) throws FFTException {
    return toImage(image, 0);
  }
  
  public BufferedImage toImage(BufferedImage image, int bias)
   throws FFTException {

    if (spectral)
      throw new FFTException("cannot transfer spectral data to an image");

    if (image == null)
      image = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
    WritableRaster raster = image.getRaster();

    int w = Math.min(image.getWidth(), width);
    int h = Math.min(image.getHeight(), height);

    // If destination image is bigger, zero it

    if (w < image.getWidth() || h < image.getHeight())
      for (int y = 0; y < image.getHeight(); ++y)
        for (int x = 0; x < image.getWidth(); ++x)
          raster.setSample(x, y, 0, 0);

    // Copy real component of data to destination image

    int i = 0, value;
    for (int y = 0; y < height; ++y)
      for (int x = 0; x < width; ++x, ++i) {
        value = Math.max(0, Math.min(255, bias + Math.round(data[i].re)));
        raster.setSample(x, y, 0, value);
      }

    return image;

  }

  public BufferedImage getSpectrum() throws FFTException {

    if (!spectral)
      throw new FFTException(NO_DATA);

    // Collect magnitudes and find maximum

    float[] magData = new float[width*height];
    float maximum = calculateMagnitudes(magData);
    BufferedImage image =
     new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
    WritableRaster raster = image.getRaster();

    // Shift, rescale and copy to image

    double scale = 255.0 / Math.log(maximum + 1.0);
    int x2 = width/2, y2 = height/2;
    int sx, sy, value;
    for (int y = 0; y < height; ++y) {
      sy = shift(y, y2);
      for (int x = 0; x < width; ++x) {
        sx = shift(x, x2);
        value = (int) Math.round(scale*Math.log(magData[sy*width+sx]+1.0));
        raster.setSample(x, y, 0, value);
      }
    }

    return image;

  }

  public BufferedImage getUnshiftedSpectrum() throws FFTException {

    if (!spectral)
      throw new FFTException(NO_DATA);

    // Collect magnitudes and find maximum

    float[] magData = new float[width*height];
    float maximum = calculateMagnitudes(magData);
    BufferedImage image =
     new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
    WritableRaster raster = image.getRaster();

    // Rescale and copy to image

    double scale = 255.0 / Math.log(maximum + 1.0);
    int i = 0, value;
    for (int y = 0; y < height; ++y)
      for (int x = 0; x < width; ++x, ++i) {
        value = (int) Math.round(scale*Math.log(magData[i]+1.0));
        raster.setSample(x, y, 0, value);
      }

    return image;

  }

  public float getMagnitude(int u, int v) throws FFTException {
    if (!spectral)
      throw new FFTException(NO_DATA);
    if (u >= 0 && u < width && v >= 0 && v < height)
      return data[v*width+u].getMagnitude();
    else
      return 0.0f;
  }


  public float getPhase(int u, int v) throws FFTException {
    if (!spectral)
      throw new FFTException(NO_DATA);
    if (u >= 0 && u < width && v >= 0 && v < height)
      return data[v*width+u].getPhase();
    else
      return 0.0f;
  }
 

  public void setMagnitude(int u, int v, float mag) throws FFTException {
    if (!spectral)
      throw new FFTException(NO_DATA);
    if (u >= 0 && u < width && v >= 0 && v < height) {
      int i = v*width+u;
      data[i].setPolar(mag, data[i].getPhase());
    }
  }


  public void setPhase(int u, int v, float phase) throws FFTException {
    if (!spectral)
      throw new FFTException(NO_DATA);
    if (u >= 0 && u < width && v >= 0 && v < height) {
      int i = v*width+u;
      data[i].setPolar(data[i].getMagnitude(), phase);
    }
  }


  public void idealLowPassFilter(double radius) throws FFTException {

    if (!spectral)
      throw new FFTException(NO_DATA);

    if (radius < 0.0 || radius > 1.0)
      throw new FFTException("invalid filter radius");

    int u2 = width/2;
    int v2 = height/2;
    int su, sv, i = 0;
    double r, rmax = Math.min(u2, v2);

    for (int v = 0; v < height; ++v) {
      sv = shift(v, v2) - v2;
      for (int u = 0; u < width; ++u, ++i) {
        su = shift(u, u2) - u2;
        r = Math.sqrt(su*su + sv*sv) / rmax;
        if (r > radius)
          data[i].re = data[i].im = 0.0f;
      }
    }

  }


  public void idealHighPassFilter(double radius) throws FFTException {

    if (!spectral)
      throw new FFTException(NO_DATA);

    if (radius < 0.0 || radius > 1.0)
      throw new FFTException("invalid filter radius");

    int u2 = width/2;
    int v2 = height/2;
    int su, sv, i = 0;
    double r, rmax = Math.min(u2, v2);

    for (int v = 0; v < height; ++v) {
      sv = shift(v, v2) - v2;
      for (int u = 0; u < width; ++u, ++i) {
        su = shift(u, u2) - u2;
        r = Math.sqrt(su*su + sv*sv) / rmax;
        if (r < radius)
          data[i].re = data[i].im = 0.0f;
      }
    }

  }

  public void idealBandPassFilter(double radius, double delta)
   throws FFTException {

    if (!spectral)
      throw new FFTException(NO_DATA);

    double delta2 = delta/2.0;
    double r1 = radius - delta2;
    double r2 = radius + delta2;
    if (r1 < 0.0 || r2 > 1.0)
      throw new FFTException(INVALID_PARAMS);

    int u2 = width/2;
    int v2 = height/2;
    int su, sv, i = 0;
    double r, rmax = Math.min(u2, v2);

    for (int v = 0; v < height; ++v) {
      sv = shift(v, v2) - v2;
      for (int u = 0; u < width; ++u, ++i) {
        su = shift(u, u2) - u2;
        r = Math.sqrt(su*su + sv*sv) / rmax;
        if (r < r1 || r > r2)
          data[i].re = data[i].im = 0.0f;
      }
    }

  }

  public void idealBandStopFilter(double radius, double delta)
   throws FFTException {

    if (!spectral)
      throw new FFTException(NO_DATA);

    double delta2 = delta/2.0;
    double r1 = radius - delta2;
    double r2 = radius + delta2;
    if (r1 < 0.0 || r2 > 1.0)
      throw new FFTException(INVALID_PARAMS);

    int u2 = width/2;
    int v2 = height/2;
    int su, sv, i = 0;
    double r, rmax = Math.min(u2, v2);

    for (int v = 0; v < height; ++v) {
      sv = shift(v, v2) - v2;
      for (int u = 0; u < width; ++u, ++i) {
        su = shift(u, u2) - u2;
        r = Math.sqrt(su*su + sv*sv) / rmax;
        if (r >= r1 && r <= r2)
          data[i].re = data[i].im = 0.0f;
      }
    }

  }
 
  public void butterworthLowPassFilter(double radius) throws FFTException {
    butterworthLowPassFilter(1, radius);
  }
  
  public void butterworthLowPassFilter(int n, double radius)
   throws FFTException {

    if (!spectral)
      throw new FFTException(NO_DATA);

    if (n < 1 || radius <= 0.0 || radius > 1.0)
      throw new FFTException(INVALID_PARAMS);

    int u2 = width/2;
    int v2 = height/2;
    int su, sv, i = 0;
    double mag, r, rmax = Math.min(u2, v2);

    for (int v = 0; v < height; ++v) {
      sv = shift(v, v2) - v2;
      for (int u = 0; u < width; ++u, ++i) {
        su = shift(u, u2) - u2;
        r = Math.sqrt(su*su + sv*sv) / rmax;
        mag = butterworthLowPassFunction(n, radius, r)*data[i].getMagnitude();
        data[i].setPolar(mag, data[i].getPhase());
      }
    }

  }

  public void butterworthHighPassFilter(double radius) throws FFTException {
    butterworthHighPassFilter(1, radius);
  }
 
  public void butterworthHighPassFilter(int n, double radius)
   throws FFTException {

    if (!spectral)
      throw new FFTException(NO_DATA);

    if (n < 1 || radius <= 0.0 || radius > 1.0)
      throw new FFTException(INVALID_PARAMS);

    int u2 = width/2;
    int v2 = height/2;
    int su, sv, i = 0;
    double mag, r, rmax = Math.min(u2, v2);

    for (int v = 0; v < height; ++v) {
      sv = shift(v, v2) - v2;
      for (int u = 0; u < width; ++u, ++i) {
        su = shift(u, u2) - u2;
        r = Math.sqrt(su*su + sv*sv) / rmax;
        mag = butterworthHighPassFunction(n, radius, r)*data[i].getMagnitude();
        data[i].setPolar(mag, data[i].getPhase());
      }
    }

  }

  public void butterworthBandPassFilter(double radius, double delta)
   throws FFTException {
    butterworthBandPassFilter(1, radius, delta);
  }
   public void butterworthBandPassFilter(int n, double radius, double delta)
   throws FFTException {

    if (!spectral)
      throw new FFTException(NO_DATA);

    double delta2 = delta/2.0;
    if (n < 1 || radius-delta2 <= 0.0 || radius+delta2 > 1.0)
      throw new FFTException(INVALID_PARAMS);

    int u2 = width/2;
    int v2 = height/2;
    int su, sv, i = 0;
    double mag, r, rmax = Math.min(u2, v2);

    for (int v = 0; v < height; ++v) {
      sv = shift(v, v2) - v2;
      for (int u = 0; u < width; ++u, ++i) {
        su = shift(u, u2) - u2;
        r = Math.sqrt(su*su + sv*sv) / rmax;
        mag = butterworthBandPassFunction(n, radius, delta, r)
         * data[i].getMagnitude();
        data[i].setPolar(mag, data[i].getPhase());
      }
    }

  }

  public void butterworthBandStopFilter(double radius, double delta)
   throws FFTException {
    butterworthBandStopFilter(1, radius, delta);
  }

  public void butterworthBandStopFilter(int n, double radius, double delta)
   throws FFTException {

    if (!spectral)
      throw new FFTException(NO_DATA);

    double delta2 = delta/2.0;
    if (n < 1 || radius-delta2<= 0.0 || radius+delta2 > 1.0)
      throw new FFTException(INVALID_PARAMS);

    int u2 = width/2;
    int v2 = height/2;
    int su, sv, i = 0;
    double mag, r, rmax = Math.min(u2, v2);

    for (int v = 0; v < height; ++v) {
      sv = shift(v, v2) - v2;
      for (int u = 0; u < width; ++u, ++i) {
        su = shift(u, u2) - u2;
        r = Math.sqrt(su*su + sv*sv) / rmax;
        mag = butterworthBandStopFunction(n, radius, delta, r)
         * data[i].getMagnitude();
        data[i].setPolar(mag, data[i].getPhase());
      }
    }

  }


  private static int powerOfTwo(int n) {
   int i = 0, m = 1;
   while (m < n) {
     m <<= 1;
     ++i;
   }
   return i;
  }


  private static void reorder(Complex[] data, int n) {
    int j = 0, m;
    for (int i = 0; i < n; ++i) {
      if (i > j)
        data[i].swapWith(data[j]);
      m = n >> 1;
      while ((j >= m) && (m >= 2)) {
        j -= m;
        m >>= 1;
      }
      j += m;
    }
  }



  private static void fft(Complex[] data, int size, int log2n, int dir) {

    double angle, wtmp, wpr, wpi, wr, wi, tmpr, tmpi;
    int n = 1, n2;
    for (int k = 0; k < log2n; ++k) {

      n2 = n;
      n <<= 1;
      angle = (-TWO_PI/n) * dir;
      wtmp = Math.sin(0.5*angle);
      wpr = -2.0*wtmp*wtmp;
      wpi = Math.sin(angle);
      wr = 1.0;
      wi = 0.0;

      for (int m = 0; m < n2; ++m) {
        for (int i = m; i < size; i += n) {
          int j = i + n2;
          tmpr = wr*data[j].re - wi*data[j].im;
          tmpi = wi*data[j].re + wr*data[j].im;
          data[j].re = (float) (data[i].re - tmpr);
          data[i].re += (float) tmpr;
          data[j].im = (float) (data[i].im - tmpi);
          data[i].im += (float) tmpi;
        }
        wtmp = wr;
        wr = wtmp*wpr - wi*wpi + wr;
        wi = wi*wpr + wtmp*wpi + wi;
      }

    }

    // Rescale results of inverse transform

    if (dir == -1)
      for (int i = 0; i < size; ++i) {
        data[i].re /= size;
        data[i].im /= size;
      }

  }



  private static final int shift(int d, int d2) {
    return (d >= d2 ? d-d2 : d+d2);
  }


  private float calculateMagnitudes(float[] mag) {
    float maximum = 0.0f;
    for (int i = 0; i < data.length; ++i) {
      mag[i] = data[i].getMagnitude();
      if (mag[i] > maximum)
        maximum = mag[i];
    }
    return maximum;
  }


}
