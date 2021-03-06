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

/**
  * Tests the behavior of metamorphosis on ND Images.
  */

#include <iostream>
#include "CALATKCommon.h"
#include "CMetamorphosisGeodesicShootingInitialImageMomentumRegistration.h"
#include "CLDDMMSimplifiedMetamorphosisGeodesicShootingRegistration.h"
#include "CLDDMMVelocityFieldWithMomentumRegistration.h"
#include "VectorImageUtils.h"
#include "CImageManager.h"
#include "CStateInitialMomentum.h"
#include "CStateInitialImageMomentum.h"

#include "CJSONConfiguration.h"

template< class TFLOAT, unsigned int VImageDimension >
int DoIt( std::string MetamorphosisType, char* sourceImage, char* targetImage, char* resultImage, const std::string & configFileName )
{
  // define the type of state
  typedef CALATK::CStateInitialImageMomentum< TFLOAT, VImageDimension > TStateInitialImageMomentum;
  typedef CALATK::CStateInitialMomentum< TFLOAT, VImageDimension >      TStateInitialMomentum;

  // define the registration method based on this state
  typedef CALATK::CMetamorphosisGeodesicShootingInitialImageMomentumRegistration< TStateInitialMomentum > regTypeFull;
  typedef CALATK::CLDDMMSimplifiedMetamorphosisGeodesicShootingRegistration< TStateInitialMomentum > regTypeSimplified;

  typedef CALATK::CMetamorphosisGeodesicShootingInitialImageMomentumRegistration< TStateInitialImageMomentum > regTypeFullInitialImage;
  typedef CALATK::CLDDMMSimplifiedMetamorphosisGeodesicShootingRegistration< TStateInitialImageMomentum > regTypeSimplifiedInitialImage;

  // general typedefs
  typedef CALATK::VectorImageUtils< TFLOAT, VImageDimension > VectorImageUtilsType;
  typedef CALATK::CImageManager< TFLOAT, VImageDimension >    ImageManagerType;
  typedef CALATK::LDDMMUtils< TFLOAT, VImageDimension >       LDDMMUtilsType;
  typedef CALATK::VectorImage< TFLOAT, VImageDimension > VectorImageType;
  typedef CALATK::VectorField< TFLOAT, VImageDimension > VectorFieldType;

  typedef CALATK::CAlgorithmBase< TFLOAT, VImageDimension > TReg;

  typename TReg::Pointer plddmm = NULL;

  // set the registration type
  if ( MetamorphosisType.compare( "MetamorphosisFullAdjoint" ) == 0 )
    {
    plddmm = new regTypeFull;
    }
  else if ( MetamorphosisType.compare( "MetamorphosisSimplified" ) == 0 )
    {
    plddmm = new regTypeSimplified;
    }
  else if ( MetamorphosisType.compare( "MetamorphosisFullAdjointInitialImage" ) == 0 )
    {
    plddmm = new regTypeFullInitialImage;
    }
  else if ( MetamorphosisType.compare( "MetamorphosisSimplifiedInitialImage" ) == 0 )
    {
    plddmm = new regTypeSimplifiedInitialImage;
    }
  else
    {
    std::cerr << "Unknown solver type " << MetamorphosisType << ". ABORT." << std::endl;
    return EXIT_FAILURE;
    }

  ImageManagerType* ptrImageManager = dynamic_cast<ImageManagerType*>( plddmm->GetImageManagerPointer() );

  ptrImageManager->AddImage( sourceImage, 0.0, 0 );
  ptrImageManager->AddImage( targetImage, 1.0, 0 );

  CALATK::CJSONConfiguration::Pointer combinedConfiguration = new CALATK::CJSONConfiguration;
  combinedConfiguration->ReadJSONConfigurationFile( configFileName );
  CALATK::CJSONConfiguration::Pointer cleanedConfiguration = new CALATK::CJSONConfiguration;
  plddmm->SetAutoConfiguration( combinedConfiguration, cleanedConfiguration );

  plddmm->Solve();

  std::string filename(resultImage);
  filename=filename.substr( 0,filename.length()-5 )+"Map.mha";

  // generating map
  const typename VectorFieldType::Pointer ptrMap1 = new VectorFieldType( plddmm->GetMap( 1.0 ) );
  VectorImageUtilsType::writeMapITK( ptrMap1, filename, false );

  const typename VectorImageType::Pointer ptrI0W1 = new VectorImageType( plddmm->GetSourceImage( 1.0 ) );
  // generating warped image
  VectorImageUtilsType::writeImageITK( ptrI0W1, resultImage );

  return EXIT_SUCCESS;
}


int calatkMetamorphosisTest( int argc, char * argv[] )
{
  const std::string configFileName = argv[7];
  const std::string metamorphosisType = argv[1];
  const std::string dimStr = argv[2];
  std::string sFloatingPointType = argv[3];

  if ( argc !=8 )
  {
    std::cout << "Found " << argc << " arguments." << std::endl;
    std::cout << "Usage: calatkMetamorphosisTest metamorphosisType dim floatType sourceImage targetImage resultingImage configFileName" << std::endl;
    return EXIT_FAILURE;
  }

  int iImageDimension = atoi( dimStr.c_str() );

  std::for_each( sFloatingPointType.begin(), sFloatingPointType.end(), ::tolower); 
  if ( sFloatingPointType.compare( "float" )==0 )
  {
    switch ( iImageDimension )
    {
    case 1:
      return DoIt<float, 1>( metamorphosisType, argv[4], argv[5], argv[6], configFileName);
      break;
    case 2:
      return DoIt<float, 2>( metamorphosisType, argv[4], argv[5], argv[6], configFileName);
      break;
    case 3:
      return DoIt<float, 3>( metamorphosisType, argv[4], argv[5], argv[6], configFileName);
      break;
    default:
      std::cerr << "Unsupported image dimension = " << iImageDimension << std::endl;
    }
  }
  else if ( sFloatingPointType.compare( "double" )==0 )
  {
    switch ( iImageDimension )
    {
    case 1:
      return DoIt<double, 1>( metamorphosisType, argv[4], argv[5], argv[6], configFileName);
      break;
    case 2:
      return DoIt<double, 2>( metamorphosisType, argv[4], argv[5], argv[6], configFileName);
      break;
    case 3:
      return DoIt<double, 3>( metamorphosisType, argv[4], argv[5], argv[6], configFileName);
      break;
    default:
      std::cerr << "Unsupported image dimension = " << iImageDimension << std::endl;
    }

  }
  else
  {
    std::cerr << "Unsupported floating point type = " << sFloatingPointType << std::endl;
  }

  return EXIT_FAILURE;

}
