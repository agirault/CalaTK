include( ExternalProject )

set( base "${CMAKE_BINARY_DIR}" )
set_property( DIRECTORY PROPERTY EP_BASE ${base} )

set( shared OFF )
set( testing ON )
set( build_type "Debug" )
if( CMAKE_BUILD_TYPE )
  set( build_type "${CMAKE_BUILD_TYPE}" )
endif()

set( CALATK_DEPENDS "" )

set( gen "${CMAKE_GENERATOR}" )

##
## Check if sytem ITK or superbuild ITK
##
if( NOT USE_SYSTEM_ITK )

  if( NOT GIT_EXECUTABLE )
    find_package( Git REQUIRED )
  endif( NOT GIT_EXECUTABLE )

  option( GIT_PROTOCOL_HTTP
    "Use HTTP for git access (useful if behind a firewall)" OFF )
  if( GIT_PROTOCOL_HTTP )
    set( GIT_PROTOCOL "http" CACHE STRING "Git protocol for file transfer" )
  else( GIT_PROTOCOL_HTTP )
    set( GIT_PROTOCOL "git" CACHE STRING "Git protocol for file transfer" )
  endif( GIT_PROTOCOL_HTTP )
  mark_as_advanced( GIT_PROTOCOL )


set ( COMMON_PROJECT_COMPILER_FLAGS
  -DCMAKE_CXX_COMPILER:PATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
  -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
  -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
  -DCMAKE_C_COMPILER:PATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
  -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
  -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
  -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
  -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
  -DCMAKE_MODULE_LINKER_FLAGS:STRING=${CMAKE_MODULE_LINKER_FLAGS}
  -DCMAKE_GENERATOR:STRING=${CMAKE_GENERATOR}
  -DCMAKE_EXTRA_GENERATOR:STRING=${CMAKE_EXTRA_GENERATOR}
)

  ##
  ## ITK
  ##
  set( proj ITK )
  ExternalProject_Add( ${proj}
    GIT_REPOSITORY "${GIT_PROTOCOL}://itk.org/ITK.git"
    GIT_TAG "v4.4.1"
    SOURCE_DIR "${CMAKE_BINARY_DIR}/ITK"
    BINARY_DIR ITK-Build
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${COMMON_PROJECT_COMPILER_FLAGS}
      -DBUILD_SHARED_LIBS:BOOL=${shared}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF

      -DUSE_FFTWF:BOOL=ON
      -DUSE_FFTWD:BOOL=ON
      -DUSE_SYSTEM_FFTW:BOOL=OFF

      -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF
      -DITKGroup_Core:BOOL=OFF
      -DITKGroup_Filtering:BOOL=OFF
      -DITKGroup_Nonunit:BOOL=OFF
      -DITKGroup_IO:BOOL=OFF
      -DITKGroup_Video:BOOL=OFF
      -DITKGroup_Segmentation:BOOL=OFF
      -DITKGroup_Registration:BOOL=OFF
      -DITKGroup_Numerics:BOOL=OFF
      -DITKGroup_Compatibility:BOOL=OFF
      -DITKGroup_Bridge:BOOL=OFF
      -DITKGroup_ThirdParty:BOOL=OFF

      -DModule_ITKTestKernel:BOOL=ON #Core
      -DModule_ITKTransform:BOOL=ON #Core
      -DModule_ITKImageFunction:BOOL=ON #Core
      -DModule_ITKImageGrid:BOOL=ON #Filtering
      -DModule_ITKIOTransformBase:BOOL=ON #IO


    INSTALL_COMMAND ""
    )
  set( ITK_DIR "${base}/ITK-Build" )
  set( CALATK_DEPENDS ${CALATK_DEPENDS} "ITK" )
endif( NOT USE_SYSTEM_ITK )

##
## A conventient 2D/3D image viewer that can handle anisotropic spacing.
##
set( ImageViewer_DEPENDS )
if( NOT USE_SYSTEM_ITK )
  set( ImageViewer_DEPENDS ITK )
endif()
set( proj ImageViewer )
ExternalProject_Add( ImageViewer
  GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/TubeTK/ImageViewer.git"
  GIT_TAG "5e20586417e6681d4c44cc6bcc61f4b35b9888da"
  SOURCE_DIR "${CMAKE_BINARY_DIR}/ImageViewer"
  BINARY_DIR ImageViewer-Build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
      ${COMMON_PROJECT_COMPILER_FLAGS}
    "-DCMAKE_BUILD_TYPE:STRING=${build_type}"
    "-DITK_DIR:PATH=${ITK_DIR}"
  INSTALL_COMMAND ""
  DEPENDS
    ${ImageViewer_DEPENDS}
  )

##
## calatk - Normal Build
##
set( proj CalaTK )
ExternalProject_Add( ${proj}
  DOWNLOAD_COMMAND ""
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  BINARY_DIR CalaTK-Build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
      ${COMMON_PROJECT_COMPILER_FLAGS}
    -DBUILD_TESTING:BOOL=${testing}
    -DBUILD_DOCUMENTATION:BOOL=${BUILD_DOCUMENTATION}
    -DCALATK_USE_SUPERBUILD:BOOL=FALSE
    -DITK_DIR:PATH=${ITK_DIR}
    -DFFTWF_LIB:FILEPATH=${CMAKE_BINARY_DIR}/ITK-Build/fftw/lib/libfftw3f.a
    -DFFTW_LIB:FILEPATH=${CMAKE_BINARY_DIR}/ITK-Build/fftw/lib/libfftw3.a
    -DFFTW_PATH:PATH=${CMAKE_BINARY_DIR}/ITK-Build/fftw/include
  INSTALL_COMMAND ""
  DEPENDS
    ${CALATK_DEPENDS}
  )
