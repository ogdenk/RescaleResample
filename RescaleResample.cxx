// RescaleResample.cxx
//
//  This program will rescale a brain MRI such that the mean signal from the two hippocampi and two thalami is 512.  A label map from Freesurfer is 
//  required to locate the structures.
//
//  The label map will also be resampled to match the origin, spacing, direction, and matrix of the MRI scan.  This is required by some ITK filters
//  used to extract radiomics features from the MRI scan

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkBinaryImageToStatisticsLabelMapFilter.h"
#include "itkBinaryMedianImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include <vector>
#include <math.h>

// These are the Freesurfer label values for the structures of interest
#define LeftHippocampus 17
#define RightHippocampus 53
#define LeftThalamus 10
#define RightThalamus 49

int main(int argc, char * argv[])
{
  const unsigned int Dimension = 3;
  size_t foundchar;  // for fileNameNoExt parsing
  unsigned short int labelValue; // holds the current label value during the loop
  float sum=0, mean=0, factor=1;  // used for rescaling the MRI 
  bool interpolate = false;
 
  // both image data and labels can be opened as unsigned short
  typedef  unsigned short  int								   PixelType;
  typedef itk::Image< PixelType, Dimension >                   ImageType;
  int i = 0;  // just a counter
  
  // the file names for the brain scan, label image, and output file
  std::string inputImageFile, fileNameNoExt, inputLabelFile, Path, outputLabelFile, outputImageFile;
  
  ImageType::SizeType medianRadius; // used for binary median image filter
  medianRadius[0] = 1;
  medianRadius[1] = 1;
  medianRadius[2] = 1;

  std::vector<unsigned short int> labelValues;
  labelValues.reserve(4);  //storage for four values
  labelValues.push_back(LeftHippocampus);
  labelValues.push_back(RightHippocampus);
  labelValues.push_back(LeftThalamus);
  labelValues.push_back(RightThalamus);

   // Get the input file specified on the command line
  
  if (argc > 1)
  {
	 	inputImageFile = argv[1];
  }
  else
  {
	std::cout << std::endl;
	std::cout << std::endl;
  	std::cout <<  "Input an image file to resample the corresponding label file to match." << std::endl;
	std::cout << "The label file should be in the same directory and have the name:" << std::endl;
	std::cout << "<fileName>_label.nrrd where <fileName> is the input file." << std::endl;
	std::cout << "The output label map file will be <fileName>_label_resamp.nrrd" << std::endl;
	std::cout << "To cause the label map to be resampled, put 'resample' after the input file name." << std::endl;
	std::cout << "If there is nothing after the file name, the input file will still be rescaled." << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	exit(1);
  }

  if (argc > 2)  interpolate = true; // anything on the command line after the file name will cause labelmap interpolation

  foundchar = inputImageFile.find_last_of("/\\");
  Path = inputImageFile.substr(0, foundchar+1);

  fileNameNoExt = inputImageFile.substr(foundchar + 1);
  foundchar = fileNameNoExt.find_last_of(".");
  inputLabelFile = Path + fileNameNoExt.substr(0, foundchar) +"_label.nrrd";

  outputLabelFile = Path + fileNameNoExt.substr(0, foundchar) + "_rescaled_label.nrrd";
  outputImageFile = Path + fileNameNoExt.substr(0, foundchar) + "_rescaled.nrrd";

  // open the image data
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName(inputImageFile );
  imageReader->Update();

  // open the labels
  typedef itk::ImageFileReader< ImageType >  LabelReaderType;
  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName(inputLabelFile);
  labelReader->Update();

  // check image geometries
  ImageType::RegionType region = imageReader->GetOutput()->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();
  ImageType::RegionType labelregion = labelReader->GetOutput()->GetLargestPossibleRegion();
  ImageType::SizeType labelsize = labelregion.GetSize();
  ImageType::SpacingType spacing = imageReader->GetOutput()->GetSpacing();
  ImageType::SpacingType labelspacing = labelReader->GetOutput()->GetSpacing();
  ImageType::PointType origin = imageReader->GetOutput()->GetOrigin();
  ImageType::PointType labelorigin = labelReader->GetOutput()->GetOrigin();

  if (!interpolate)
  {
	  if (abs(int(size[0] - labelsize[0])) > 0 | abs(int(size[1] - labelsize[1])) > 0 | abs(int(size[2] - labelsize[2])) > 0 |
		  abs(spacing[0] - labelspacing[0]) > 0.01 | abs(spacing[1] - labelspacing[1]) > 0.01 | abs(spacing[2] - labelspacing[2]) > 0.01 |
		  abs(origin[0] - labelorigin[0]) > 0.01 | abs(origin[1] - labelorigin[1]) > 0.01 | abs(origin[2] - labelorigin[2]) > 0.01)
	  {
		  std::cout << "There is a problem with the data." << std::endl;
		  std::cout << "Interpolation of the label map data is turned off but the label map" << std::endl;
		  std::cout << "don't have the same geometry.  Make sure you have interpolation turned on" << std::endl;
		  std::cout << "if the label map and image data have different size, spacing, or origin." << std::endl;
		  exit(1);
	  }
	  labelReader->GetOutput()->SetSpacing(imageReader->GetOutput()->GetSpacing());
	  labelReader->GetOutput()->SetOrigin(imageReader->GetOutput()->GetOrigin());
	
  }
  
  // Need an interpolator, should use nearest neighbor for label maps
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double >  InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  

  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput(labelReader->GetOutput());
  resampler->SetInterpolator(interpolator);
  resampler->SetSize(size);
  resampler->SetOutputSpacing(spacing);
  resampler->SetOutputOrigin(origin);
  resampler->SetOutputDirection(imageReader->GetOutput()->GetDirection());
  

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputLabelFile);
  
  if(interpolate)
	  writer->SetInput(resampler->GetOutput());
  else
	  writer->SetInput(labelReader->GetOutput()); // if no interpolate then just copy the input file
  writer->Update();

  if (interpolate) // only need to overwrite the input label file if interpolation occurred
  {
	  writer->SetFileName(inputLabelFile);
	  writer->Update(); // Save two copies with different names to support the RFE program which uses the input name
  }					// to derive the label file name.  Original label file gets overwritten with resampled image

  // Now we use the resampled label image to get the stats for the four structures

  //  Median filter to smooth the edges of the labels
  typedef itk::BinaryMedianImageFilter<ImageType, ImageType> MedianFilterType;
  MedianFilterType::Pointer medianFilter = MedianFilterType::New();
  medianFilter->SetRadius(medianRadius);
  medianFilter->SetBackgroundValue(0);
  medianFilter->SetInput(resampler->GetOutput());

  //  Statistics filter to get the gray level statistics from the MRI image for a particular label location
  typedef itk::BinaryImageToStatisticsLabelMapFilter<ImageType, ImageType> BinaryImageToStatisticsLabelMapFilterType;
  BinaryImageToStatisticsLabelMapFilterType::Pointer BinaryToStatisticsFilter = BinaryImageToStatisticsLabelMapFilterType::New();
  BinaryToStatisticsFilter->SetInput1(medianFilter->GetOutput());
  BinaryToStatisticsFilter->SetComputeHistogram(TRUE);
  BinaryToStatisticsFilter->SetCoordinateTolerance(0.01);
  BinaryImageToStatisticsLabelMapFilterType::OutputImageType::LabelObjectType* StatlabelObject;

  std::vector<unsigned short int>::iterator labelValuesIterator = labelValues.begin();

  while (labelValuesIterator != labelValues.end())
  {
	  labelValue = labelValues[i];
	  medianFilter->SetForegroundValue(labelValue);

	  BinaryToStatisticsFilter->SetInput2(imageReader->GetOutput());
	  BinaryToStatisticsFilter->SetInputForegroundValue(labelValue);
	  BinaryToStatisticsFilter->Update();
	  // should only be one object
	  std::cout << "There is " << BinaryToStatisticsFilter->GetOutput()->GetNumberOfLabelObjects() << " object with statistics." << std::endl;
	  StatlabelObject = BinaryToStatisticsFilter->GetOutput()->GetNthLabelObject(0);
	  sum += StatlabelObject->GetMean();

	  i++;
	  labelValuesIterator++;
  }
  mean = sum / 4.0;
  factor = 512. / mean;

  typedef itk::ShiftScaleImageFilter<ImageType, ImageType> ScalerType;
  ScalerType::Pointer scaler = ScalerType::New();
  scaler->SetInput(imageReader->GetOutput());
  scaler->SetScale(factor);
  scaler->SetShift(0);
  scaler->Update();

  writer->SetFileName(outputImageFile);
  writer->SetInput(scaler->GetOutput());
  writer->Update();
 
}