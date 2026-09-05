#include <felitronics/dynamics/Compressor.h>
using namespace felitronics::dynamics;
int main(){ CompressorParams a { Detector::Rms, LinkMode::MeanPower, 7.0, Mode::UpCompress, -30.0 }; (void)a; }
