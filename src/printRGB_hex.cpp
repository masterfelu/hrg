#include "TColor.h"

void printRGB_hex()
{
	int colorIndex = kBlue;
	
	TColor *color = gROOT->GetColor(colorIndex);
	int r = color->GetRed()*255;
	int g = color->GetGreen()*255;
	int b = color->GetBlue()*255;
	std::cout << '#' << std::setfill('0') << std::setw(2) << std::hex << r << std::setw(2) << g << std::setw(2) << b << std::endl;
}