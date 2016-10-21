#include "Halide.h"

#include <gflags/gflags.h>

DEFINE_int32(output_bpp, 8,  "output image bits per pixel, either 8 or 16");

using namespace std;
using namespace Halide;

template<bool Fast>
class CameraIspGen
{
 protected:
  Target target;
  Var x, y, c, yi, yo;
  //Func toneCorrected;
  //Func ispOutput;

  // Average two positive values rounding up
  Expr avg(Expr a, Expr b) {
    return (a + b) * 0.5f;
  }

  Func interleave_x(Func a, Func b) {
    Func out("ix");
    out(x, y) = select((x % 2) == 0, a(x, y), b(x, y));
    return out;
  }

  Func interleave_y(Func a, Func b) {
    Func out("iy");
    out(x, y) = select((y % 2) == 0, a(x, y), b(x, y));
    return out;
  }


  Expr redPixel() {
    return (x % 2) == 0 && (y % 2) == 1;
  }

  Expr greenRedPixel() {
    return (x % 2) == 1 && (y % 2) == 1;
  }


  Expr bluePixel() {
    return (x % 2) == 1 && (y % 2) == 0;
   }

  Expr greenBluePixel() {
    return (x % 2) == 0 && (y % 2) == 0;
  }

  Expr greenPixel() {
    return greenBluePixel() || greenRedPixel();
  }


  Func deinterleave(Func raw) {
    // Deinterleave the color channels
    Func deinterleaved("deinterleaved");

    // Fixed pattern of:
    // G B
    // R G
    deinterleaved(x, y, c) =
      cast<float>(
          select(c == 1 && greenBluePixel(), raw(x, y), // green on blue row
          select(c == 2 && bluePixel(),      raw(x, y), // blue
          select(c == 0 && redPixel(),       raw(x, y), // red
          select(c == 1 && greenRedPixel(),  raw(x, y), // green on red row
          0)))));

    return deinterleaved;
  }

  void bilinearDemosaic(
      Func rawRed,
      Func rawBlue,
      Func rawGreenB,
      Func rawGreenR,
      Func& redOnGreenBlue,
      Func& greenOnGreenBlue,
      Func& blueOnGreenBlue,
      Func& redOnGreenRed,
      Func& greenOnGreenRed,
      Func& blueOnGreenRed,
      Func& redOnBlue,
      Func& greenOnBlue,
      Func& blueOnBlue,
      Func& redOnRed,
      Func& greenOnRed,
      Func& blueOnRed,
      Func& demosaiced) {
    // Neighbor indices
    Expr x1 = x + 1;
    Expr x_1 = x - 1;

    Expr y1 = y + 1;
    Expr y_1 = y - 1;

    // Horizontally or vertically average the red and blue channels at
    // green-red and green-blue pixels.
    redOnGreenBlue(x, y) = avg(rawRed(x, y_1), rawRed(x,y1));
    greenOnGreenBlue(x, y) = rawGreenB(x, y);
    blueOnGreenBlue(x, y) = avg(rawBlue(x_1, y), rawBlue(x1,y));

    blueOnGreenRed(x, y) = avg(rawBlue(x, y_1), rawBlue(x,y1));
    greenOnGreenRed(x, y) = rawGreenR(x, y);
    redOnGreenRed(x, y) = avg(rawRed(x_1, y), rawRed(x1,y));

    // Box average pixels surrounding the red or blue channels when
    // they are alternatively on blue and red pixels respectively.
    redOnBlue(x, y) = avg(avg(rawRed(x_1, y_1), rawRed(x1, y_1)), avg(rawRed(x_1, y1), rawRed(x1, y1)));
    greenOnBlue(x, y) = avg(avg(rawGreenB(x_1, y), rawGreenB(x1, y)), avg(rawGreenR(x, y_1), rawGreenR(x, y1)));
    blueOnBlue(x, y) = rawBlue(x, y);

    // Average the green colors that are above, below, left and right
    // at red and blue pixels.
    redOnRed(x, y) = rawRed(x, y);
    greenOnRed(x, y) = avg(avg(rawGreenR(x_1, y), rawGreenR(x1, y)), avg(rawGreenB(x, y_1), rawGreenB(x, y1)));
    blueOnRed(x, y) = avg(avg(rawBlue(x_1, y_1), rawBlue(x1, y_1)), avg(rawBlue(x_1, y1), rawBlue(x1, y1)));

    // Interleave the resulting channels
    Func r, g, b;
    r = interleave_y(interleave_x(redOnGreenBlue,  redOnBlue),    interleave_x(redOnRed, redOnGreenRed));
    g = interleave_y(interleave_x(greenOnGreenBlue,greenOnBlue),  interleave_x(greenOnRed, greenOnGreenRed));
    b = interleave_y(interleave_x(blueOnGreenBlue, blueOnBlue),   interleave_x(blueOnRed, blueOnGreenRed));

    demosaiced(x, y, c) =
      select(
          c == 0, r(x, y),
          c == 1, g(x, y),
          b(x, y));
  }

  void edgeAwareDemosaic(
      Func rawRed,
      Func rawBlue,
      Func rawGreenB,
      Func rawGreenR,
      Func& redOnGreenBlue,
      Func& greenOnGreenBlue,
      Func& blueOnGreenBlue,
      Func& redOnGreenRed,
      Func& greenOnGreenRed,
      Func& blueOnGreenRed,
      Func& redOnBlue,
      Func& greenOnBlue,
      Func& blueOnBlue,
      Func& redOnRed,
      Func& greenOnRed,
      Func& blueOnRed,
      Func& demosaiced) {
    // Neighbor indices
    Expr x1 = x + 1;
    Expr x2 = x + 2;
    Expr x_1 = x - 1;
    Expr x_2 = x - 2;

    Expr y1 = y + 1;
    Expr y2 = y + 2;
    Expr y_1 = y - 1;
    Expr y_2 = y - 2;

    // Green horizonal & vertical value interpolates at green input pixels
    Func gV("gV"), gH("gH");
    gV(x, y) =
      select(greenBluePixel(), rawGreenB(x, y), // green on blue row
      select(bluePixel(),      avg(rawGreenR(x, y1), rawGreenR(x, y_1)) + (2.0f * rawBlue(x, y) - rawBlue(x, y2) - rawBlue(x, y_2)) / 4.0f, // blue
      select(redPixel(),       avg(rawGreenB(x, y1), rawGreenB(x, y_1)) + (2.0f * rawRed(x, y)  - rawRed(x, y2)  - rawRed(x, y_2)) / 4.0f, // red
      select(greenRedPixel(),  rawGreenR(x, y), // green on red row
      0))));

    gH(x, y) =
      select(greenBluePixel(), rawGreenB(x, y), // green on blue row
      select(bluePixel(),      avg(rawGreenB(x1, y), rawGreenB(x_1, y)) + (2.0f * rawBlue(x, y) - rawBlue(x2, y) - rawBlue(x_2, y)) / 4.0f, // blue
      select(redPixel(),       avg(rawGreenR(x1, y), rawGreenR(x_1, y)) + (2.0f * rawRed(x, y)  - rawRed(x2, y)  - rawRed(x_2, y)) / 4.0f, // red
      select(greenRedPixel(),  rawGreenR(x, y), // green on red row
      0))));

    // And the derivatives
    Func dV("dV"), dH("dh");
    dV(x, y) =
      select(greenBluePixel(), avg(absd(rawGreenB(x, y2), rawGreenB(x, y)), absd(rawGreenB(x, y_2), rawGreenB(x, y))), // green on blue row
      select(bluePixel(),      avg(absd(rawGreenR(x, y1), rawGreenR(x, y_1)), absd(rawBlue(x, y2) + rawBlue(x, y_2), 2.0f*rawBlue(x, y))), // blue
      select(redPixel(),       avg(absd(rawGreenB(x, y1), rawGreenB(x, y_1)), absd(rawRed(x, y2)  + rawRed(x, y_2),  2.0f*rawRed(x, y))), // red
      select(greenRedPixel(),  avg(absd(rawGreenR(x, y2), rawGreenR(x, y)), absd(rawGreenR(x, y_2), rawGreenR(x, y))), // green on red row
      0))));

    dH(x, y) =
      select(greenBluePixel(), avg(absd(rawGreenB(x2, y), rawGreenB(x,   y)), absd(rawGreenB(x_2, y), rawGreenB(x, y))), // green on blue row
      select(bluePixel(),      avg(absd(rawGreenB(x1, y), rawGreenB(x_1, y)), absd(rawBlue(x2, y) + rawBlue(x_2, y), 2.0f*rawBlue(x, y))), // blue
      select(redPixel(),       avg(absd(rawGreenR(x1, y), rawGreenR(x_1, y)), absd(rawRed(x2, y)  + rawRed(x_2, y),  2.0f*rawRed(x, y))), // red
      select(greenRedPixel(),  avg(absd(rawGreenR(x2, y), rawGreenR(x,   y)), absd(rawGreenR(x_2, y), rawGreenR(x, y))), // green on red row
      0))));

    // Homogenity calulation over a diameter size region
    const int w = 2;
    const int diameter = 2*w + 1;
    RDom d(0, diameter, 0, diameter);
    Expr xp = x + d.x - w;
    Expr yp = y + d.y - w;

    // Count the number of pixels in the region that are horizontal.
    // This is given by determining when the vertical gradient is
    // greater than the horizontal gradient. Remember that a strong
    // vertical gradiant means that the local area is horizontal since
    // the gradient is orthogonal to the direction of the local
    // feature.
    Func hCount("hCount");
    hCount(x, y) = sum(cast<int>(dH(xp,yp) < dV(xp,yp)), "hCount_sum");

    // The local green estimate is a blend the horizontal and vertical
    // green channel estimate based on how horizontal or vertical the
    // region is.
    Expr alpha = cast<float>(hCount(x, y)) / float(diameter * diameter);

    Func g("green");
    g(x, y) = lerp(gV(x, y), gH(x, y), alpha);
    
    // We take the log of green channel which helps to supress chroma
    // ringing and noise.
    Func gp("log green");
    gp(x, y) = log(1.0f + g(x, y));

    // Red and blue minus the log of the green channel
    Func rmg("rmg"), bmg("bmg");
    rmg(x, y) = rawRed(x, y) - gp(x, y);
    bmg(x, y) = rawBlue(x, y) - gp(x, y);

    // Use local box style filtering to estimate r-g and b-g
    // channels. Adding the g channel back in.
    Func r("r"), b("b");
    r(x, y) =
      select(redPixel(),       (rmg(x, y) + rmg(x, y2) + rmg(x, y_2) + rmg(x2, y) + rmg(x_2, y)) / 5.0f,
      select(greenRedPixel(),  (rmg(x_1, y_2) + rmg(x_1, y) + rmg(x_1, y2) + rmg(x1, y_2) + rmg(x1, y) + rmg(x1, y2)) / 6.0f,
      select(greenBluePixel(), (rmg(x_2, y_1) + rmg(x, y_1) + rmg(x2, y_1) + rmg(x_2, y1) + rmg(x, y1) + rmg(x2, y1)) / 6.0f,
      select(bluePixel(),      (rmg(x1, y1) + rmg(x1, y_1) + rmg(x_1, y1) + rmg(x_1, y_1)) / 4.0f,
      0))))
      + gp(x, y);

    b(x, y) =
      select(redPixel(),       (bmg(x1, y1) + bmg(x1, y_1) + bmg(x_1, y1) + bmg(x_1, y_1)) / 4.0f,
      select(greenRedPixel(),  (bmg(x_2, y_1) + bmg(x, y_1) + bmg(x2, y_1) + bmg(x_2, y1) + bmg(x, y1) + bmg(x2, y1)) / 6.0f,
      select(greenBluePixel(), (bmg(x_1, y_2) + bmg(x_1, y) + bmg(x_1, y2) + bmg(x1, y_2) + bmg(x1, y) + bmg(x1, y2)) / 6.0f,
      select(bluePixel(),      (bmg(x, y) + bmg(x, y2) + bmg(x, y_2) + bmg(x2, y) + bmg(x_2, y)) / 5.0f,
      0))))
      + gp(x, y);

    demosaiced(x, y, c) =
      select(
          c == 0, r(x, y),
          c == 1, g(x, y),
          b(x, y));
  }

  Func demosaic(
      Func deinterleaved,
      Func vignetteTableH,
      Func vignetteTableV,
      Expr blackLevelR,
      Expr blackLevelG,
      Expr blackLevelB,
      Expr whiteBalanceGainR,
      Expr whiteBalanceGainG,
      Expr whiteBalanceGainB) {
    // These are the values we already know from the input
    // x_y = the value of channel x at a site in the input of channel y
    // gb refers to green sites in the blue rows
    // gr refers to green sites in the red rows

    Expr minRawR = blackLevelR;
    Expr minRawG = blackLevelG;
    Expr minRawB = blackLevelB;
    Expr maxRaw = float((1 << 16) - 1);

    Expr invRangeR = whiteBalanceGainR / (maxRaw - minRawR);
    Expr invRangeG = whiteBalanceGainG / (maxRaw - minRawG);
    Expr invRangeB = whiteBalanceGainB / (maxRaw - minRawB);

    // Raw bayer pixels
    Func rawRed("rawRed");
    Func rawBlue("rawBlue");
    Func rawGreenB("rawGreenB");
    Func rawGreenR("rawGreenR");

    if (Fast) {
      rawRed(x, y)    = (deinterleaved(x, y, 0) - minRawR) * invRangeR;
      rawBlue(x, y)   = (deinterleaved(x, y, 2) - minRawB) * invRangeB;
      rawGreenB(x, y) = (deinterleaved(x, y, 1) - minRawG) * invRangeG;
      rawGreenR(x, y) = (deinterleaved(x, y, 1) - minRawG) * invRangeG;
    } else {
      rawRed(x, y)    = (deinterleaved(x, y, 0) - minRawR) * invRangeR * vignetteTableH(0, x) * vignetteTableV(0, y);
      rawBlue(x, y)   = (deinterleaved(x, y, 2) - minRawB) * invRangeB * vignetteTableH(2, x) * vignetteTableV(2, y);
      rawGreenB(x, y) = (deinterleaved(x, y, 1) - minRawG) * invRangeG * vignetteTableH(1, x) * vignetteTableV(1, y);
      rawGreenR(x, y) = (deinterleaved(x, y, 1) - minRawG) * invRangeG * vignetteTableH(1, x) * vignetteTableV(1, y);
    }

    // These are the ones we need to interpolate
    // Rename build up planar bayer channels
    Func redOnGreenBlue("redOnGreenBlue"); // Red on a greeb/blue pixel
    Func greenOnGreenBlue("greenOnGreenBlue"); // Green on a green/blue pixel
    Func blueOnGreenBlue("blueOnGreenBlue"); // Blue on a green/blue pixel

    Func redOnGreenRed("redOnGreenRed");  // Red on a green/red pixel
    Func greenOnGreenRed("greenOnGreenRed"); // Green on a green/red pixel
    Func blueOnGreenRed("blueOnGreenRed"); // Blue on a green/red pixel

    Func redOnBlue("redOnBlue");  // Red on a blue pixel
    Func greenOnBlue("greenOnBlue"); // Green on a blue pixel
    Func blueOnBlue("blueOnBlue");  // Blue on a blue pixel

    Func redOnRed("redOnRed"); // red on a red pixel
    Func greenOnRed("greenOnRed"); // Green on a red pixel
    Func blueOnRed("blueOnRed"); // Blue on a red pixel

    Func demosaiced("demosaiced");
    if (Fast) {
      bilinearDemosaic(rawRed, rawBlue, rawGreenB, rawGreenR,
          redOnGreenBlue, greenOnGreenBlue, blueOnGreenBlue,
          redOnGreenRed, greenOnGreenRed, blueOnGreenRed,
          redOnBlue, greenOnBlue, blueOnBlue,
          redOnRed, greenOnRed, blueOnRed,
          demosaiced);
    } else {
      edgeAwareDemosaic(rawRed, rawBlue, rawGreenB, rawGreenR,
          redOnGreenBlue, greenOnGreenBlue, blueOnGreenBlue,
          redOnGreenRed, greenOnGreenRed, blueOnGreenRed,
          redOnBlue, greenOnBlue, blueOnBlue,
          redOnRed, greenOnRed, blueOnRed,
          demosaiced);
    }

    return demosaiced;
  }

  static const int kToneCurveLutSize = 4096;

  Func applyCCM(Func input, ImageParam ccm) {
    Expr ir = input(x, y, 0);
    Expr ig = input(x, y, 1);
    Expr ib = input(x, y, 2);

    Expr r = ccm(0, 0) * ir + ccm(1, 0) * ig + ccm(2, 0) * ib;
    Expr g = ccm(0, 1) * ir + ccm(1, 1) * ig + ccm(2, 1) * ib;
    Expr b = ccm(0, 2) * ir + ccm(1, 2) * ig + ccm(2, 2) * ib;

    Func corrected("corrected");
    corrected(x, y, c) =
      cast<uint16_t>(
          clamp(
              select(
                  c == 0, r,
                  c == 1, g,
                  b), 0, kToneCurveLutSize - 1));
    return corrected;
  }


  Func gaussian_blur(Func input) {
    Func k("gaussian_params"), gaussian_x("gaussian_x"), gaussian_y("gaussian_y");
    RDom r(-2,5);

    k(x) = 0;
    k(-2) = 1;    k(-1) = 3;   k(0) = 7;   k(1) = 3;  k(2)= 1;
    gaussian_x(x,y,c) = sum(input(x+r.x, y, c) *  k(r)) / 15;
    gaussian_y(x,y,c) = sum(gaussian_x(x, y+r, c) * k(r)) / 15;

    return gaussian_y;
  }

  
  Func lp_2d_fast(
      Func input,
      Expr height) {

    Func out("gaussian_2d");
    out = gaussian_blur(input);

    return out;
  }

  Func applyUnsharpMask(
      Func input,
      Expr width,
      Expr height,
      Expr sR,
      Expr sG,
      Expr sB,
      Param<int> BGR,
      int outputBpp) {

    Func lowPass("lowPass");
    Func lowPass0("lowPass0");

    lowPass = lp_2d_fast(BoundaryConditions::mirror_image(input,0,width,0,height), height);

    Func highPass("hp");
    highPass(x, y, c) =  cast<float>(input(x, y, c)) - lowPass(x, y, c);

    Expr cp = select(BGR == 0, c, 2 - c);
    const int range = (1 << outputBpp) - 1;

    Func unsharp("unsharp");
    unsharp(x,y,c) = lowPass(x,y,c) + highPass(x,y,c) * select(cp == 0, sR,
							       cp == 1, sG,
							       sB);
    
    Func outputVal("outputVal");
    outputVal(x,y,c) = clamp( unsharp(x,y,c), 0, range );

    Func sh("sh");
    sh(x,y,c) = outputVal(x,y,c);

    Expr shr;
    shr  = sh(x,y,c);
    
    Func sharpened("sharpened");
    if (outputBpp == 8) {
      sharpened(x, y, c) = cast<uint8_t>(shr);
    } else {
      sharpened(x, y, c) = cast<uint16_t>(shr);
    }

    return sharpened;
  }

 public:
  Func generate(
      Target target,
      Func raw,
      Param<int> width,
      Param<int> height,
      Func vignetteTableH,
      Func vignetteTableV,
      Param<float> blackLevelR,
      Param<float> blackLevelG,
      Param<float> blackLevelB,
      Param<float> whiteBalanceGainR,
      Param<float> whiteBalanceGainG,
      Param<float> whiteBalanceGainB,
      Param<float> sharpenningR,
      Param<float> sharpenningG,
      Param<float> sharpenningB,
      ImageParam ccm,
      ImageParam toneTable,
      Param<int> BGR,
      int outputBpp) {

    this->target = target;
    Expr sR = 1.0f + cast<float>(sharpenningR);
    Expr sG = 1.0f + cast<float>(sharpenningG);
    Expr sB = 1.0f + cast<float>(sharpenningB);

    // This is the ISP pipeline
    Func vignettingGain("vignettingGain");
    Func deinterleaved("deinterleaved");
    Func demosaiced("demosaiced");
    Func colorCorrected("ccm");
    Func sharpened("sharpened");

    deinterleaved = deinterleave(raw);
    deinterleaved.compute_root();
    
    demosaiced = demosaic(deinterleaved, vignetteTableH, vignetteTableV, blackLevelR, blackLevelG, blackLevelB, whiteBalanceGainR, whiteBalanceGainG, whiteBalanceGainB);

    colorCorrected = applyCCM(demosaiced, ccm);
    Func toneCorrected("toneCorrected");
    toneCorrected(x, y, c) = toneTable(c, colorCorrected(x, y, c));

    demosaiced.compute_root();
    toneCorrected.compute_root();
    
    Func ispOutput("IspOutput");
    if (Fast) {
      Expr cp = select(BGR == 0, c, 2 - c);
      ispOutput(x, y, c) = toneCorrected(x , y, cp);
    } else {
      ispOutput = applyUnsharpMask(toneCorrected, width, height, sR, sG, sB, BGR, outputBpp);
    }
    
    ispOutput.output_buffer()
      .set_stride(0,3)
      .set_stride(2,1)
      .set_bounds(2,0,3);
    
    return ispOutput;
  }

  CameraIspGen() {};
};


int main(int argc, char **argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  ImageParam input(UInt(16), 2);
  Param<int> width;
  Param<int> height;
  ImageParam vignetteTableH(Float(32), 2, "vignetteH");
  ImageParam vignetteTableV(Float(32), 2, "vignetteV");
  Param<float> blackLevelR("blackLevelR");
  Param<float> blackLevelG("blackLevelG");
  Param<float> blackLevelB("blackLevelB");
  Param<float> whiteBalanceGainR("whiteBalanceGainR");
  Param<float> whiteBalanceGainG("whiteBalanceGainG");
  Param<float> whiteBalanceGainB("whiteBalanceGainB");
  Param<float> sharpenningR("sharpenningR");
  Param<float> sharpenningG("sharpenningG");
  Param<float> sharpenningB("sharpenninngB");
  ImageParam ccm(Float(32), 2, "ccm");
  ImageParam toneTable(UInt(FLAGS_output_bpp), 2, "toneTable");
  Param<int> BGR;
  

  // Pick a target architecture
  Target target = get_target_from_environment();

  target.set_feature( Target::Debug, true);
  target.set_feature( Target::Profile, true);
  target.set_feature( Target::OpenCL, true );
  
  // Build variants of the pipeline
  CameraIspGen<false> cameraIspGen;
  CameraIspGen<true> cameraIspGenFast;

  Func rawM("rawm");
  rawM = BoundaryConditions::mirror_interior(input);

  Func vignetteTableHM("vhm");
  vignetteTableHM = BoundaryConditions::mirror_image(vignetteTableH);

  Func vignetteTableVM("vvm");
  vignetteTableVM = BoundaryConditions::mirror_image(vignetteTableV);

  Func cameraIsp("cameraIsp");
  cameraIsp = cameraIspGen.generate(
      target, rawM, width, height, vignetteTableHM, vignetteTableVM, blackLevelR, blackLevelG, blackLevelB, whiteBalanceGainR,
      whiteBalanceGainG, whiteBalanceGainB, sharpenningR, sharpenningG, sharpenningB, ccm, toneTable, BGR, FLAGS_output_bpp);

  Func cameraIspFast("cameraIspFast");
  cameraIspFast =
    cameraIspGenFast.generate(
        target, rawM, width, height, vignetteTableHM, vignetteTableVM, blackLevelR, blackLevelG, blackLevelB, whiteBalanceGainR,
        whiteBalanceGainG, whiteBalanceGainB, sharpenningR, sharpenningG, sharpenningB, ccm, toneTable, BGR, FLAGS_output_bpp);

  std::vector<Argument> args = {
    input, width, height, vignetteTableH, vignetteTableV, blackLevelR, blackLevelG, blackLevelB, whiteBalanceGainR, whiteBalanceGainG,
    whiteBalanceGainB, sharpenningR, sharpenningG, sharpenningB,  ccm, toneTable,  BGR};

  // Compile the pipelines
  // Use to cameraIsp.print_loop_nest() here to debug loop unrolling
  std::cout << "Halide: " << "Generating " << FLAGS_output_bpp << " bit isp" << std::endl;
  cameraIsp.compile_to_static_library("CameraIspGen" + to_string(FLAGS_output_bpp), args, "CameraIspGen" + to_string(FLAGS_output_bpp) , target);
  cameraIsp.compile_to_assembly("CameraIspGen" + to_string(FLAGS_output_bpp) + ".s", args, target);
  cameraIsp.compile_to_lowered_stmt("CameraIspGen" + to_string(FLAGS_output_bpp) + ".html", args, HTML , target);

  
  std::cout << "Halide: " << "Generating " << FLAGS_output_bpp << " bit fastest isp" << std::endl;
  cameraIspFast.compile_to_static_library("CameraIspGenFast" + to_string(FLAGS_output_bpp), args, "CameraIspGenFast" + to_string(FLAGS_output_bpp) , target);
  cameraIspFast.compile_to_assembly("CameraIspGenFast" + to_string(FLAGS_output_bpp) + ".s", args, target);
  return EXIT_SUCCESS;
}
