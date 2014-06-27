/*
*
*  Copyright 2011, 2012 by the CALATK development team
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*
*
*/

#include "CALATKCommon.h"
#include "ApplicationUtils.h"

//
// Create a filename with numbering
//

namespace CALATK
{

std::string GetCALATKVersionString()
{
  std::string calatkVersionMajor;
  std::string calatkVersionMinor;

#ifdef CALATK_VERSION_MAJOR
  calatkVersionMajor = CreateIntegerString( CALATK_VERSION_MAJOR, 0);
#else
  calatkVersionMajor = "?";
#endif

#ifdef CALATK_VERSION_MINOR
  calatkVersionMinor = CreateIntegerString( CALATK_VERSION_MINOR, 0 );
#else
  calatkVersionMinor = "?";
#endif

  std::string versionString = calatkVersionMajor + "." + calatkVersionMinor;

  return versionString;
}

std::string GetCALATKJsonHeaderString()
{
  std::string headerString = "// CALATK version " + GetCALATKVersionString() + "; JSON configuration file";

  return headerString;
}

std::string CreateNumberedFileName( std::string strPrefix, unsigned int uiNr, std::string postFix )
{
  std::stringstream ss;
  ss << strPrefix << CreateIntegerString( (int)uiNr, 4 ) << postFix;
  return ss.str();
}

std::string CreateIntegerString( int iNr, unsigned int uiW )
{
  std::stringstream ss;//create a stringstream
  ss << std::setfill('0');
  ss << std::setw( uiW ) << iNr;
  return ss.str(); //return a string with the contents of the stream
}

unsigned int GetNonSingletonImageDimensionFromFile( const std::string & sourceImageOrig )
{

  // expand file name
  std::string sourceImage = ApplicationUtils::findDataFileName( sourceImageOrig );

  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO( sourceImage.c_str(),
                                   itk::ImageIOFactory::ReadMode );

  if( !imageIO )
    {
    std::string message = "No itk::ImageIO was found for image: " + sourceImage;
    throw std::runtime_error( message.c_str() );
    }

  // Now that we found the appropriate ImageIO class, ask it to
  // read the meta data from the image file.
  imageIO->SetFileName( sourceImage.c_str() );
  imageIO->ReadImageInformation();

  unsigned int uiDim = imageIO->GetNumberOfDimensions();

  unsigned int uiEff = 0;

  for ( unsigned int iI=0; iI<uiDim; ++iI )
    {
    //std::cout << "dim " << iI << " = " << imageIO->GetDimensions( iI ) << std::endl;
    if ( imageIO->GetDimensions( iI ) > 1 )
      {
      ++uiEff;
      }
    }

  return uiEff;

}

void CheckIfSameHeader(const std::string & imageOrig1, const std::string & imageOrig2, unsigned int uiDim)
{

    std::string imageName1 = ApplicationUtils::findDataFileName( imageOrig1 );
    itk::ImageIOBase::Pointer image1 = itk::ImageIOFactory::CreateImageIO( imageName1.c_str(), itk::ImageIOFactory::ReadMode );
    if( !image1 )
    {
      std::string message = "No itk::ImageIO was found for image: " + imageName1;
      throw std::runtime_error( message.c_str() );
    }
    image1->SetFileName( imageName1.c_str() );
    image1->ReadImageInformation();
    std::string imageName2 = ApplicationUtils::findDataFileName( imageOrig2 );
    itk::ImageIOBase::Pointer image2 = itk::ImageIOFactory::CreateImageIO( imageName2.c_str(), itk::ImageIOFactory::ReadMode );
    if( !image2 )
    {
      std::string message = "No itk::ImageIO was found for image: " + imageName2;
      throw std::runtime_error( message.c_str() );
    }
    image2->SetFileName( imageName2.c_str() );
    image2->ReadImageInformation();

    const double diffMax = 1.0e-6; // TODO : modify value if too strict (we have experienced diff of 1.0e-8 for spacing)

    for ( unsigned int iI=0; iI<uiDim; ++iI )
    {
        itk::SizeValueType dim1 = image1->GetDimensions( iI );
        itk::SizeValueType dim2 = image2->GetDimensions( iI );
        if (dim1 != dim2)
        {
            std::cout << imageOrig1 << " : dim(" << iI << ") = " << dim1 << std::endl;
            std::cout << imageOrig2 << " : dim(" << iI << ") = " << dim2 << std::endl;
            throw std::runtime_error("Images size do not match ");
        }

        std::vector<double> dir1 = image1->GetDirection( iI );
        std::vector<double> dir2 = image2->GetDirection( iI );
        if (dir1 != dir2)
        {
            std::cout << imageOrig1 << " : spaceDirection(" << iI << ") = " << dir1;
            std::cout << imageOrig2 << " : spaceDirection(" << iI << ") = " << dir2;
            throw std::runtime_error("Images space direction do not match ");
        }

        double spacing1 = image1->GetSpacing( iI );
        double spacing2 = image2->GetSpacing( iI );
        double diffSpacing = std::abs(1-(spacing1/spacing2));
        if (diffSpacing > diffMax)
        {
            std::cout << imageOrig1 << " : spacing(" << iI << ") = " << spacing1 << std::endl;
            std::cout << imageOrig2 << " : spacing(" << iI << ") = " << spacing2 << std::endl;
            throw std::runtime_error("Images spacing do not match ");
        }

        double origin1 = image1->GetOrigin( iI );
        double origin2 = image2->GetOrigin( iI );
        if (origin1 != origin2)
        {
            std::cout << imageOrig1 << " : origin(" << iI << ") = " << origin1 << std::endl;
            std::cout << imageOrig2 << " : origin(" << iI << ") = " << origin2 << std::endl;
            throw std::runtime_error("Images origin do not match ");
        }
    }
}

}
