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

#ifndef VECTOR_FIELD_TXX
#define VECTOR_FIELD_TXX

#include "VectorField.h"

namespace CALATK
{

//////////////////////////////////
// Constructors and Destructors //
//////////////////////////////////

template <class T, unsigned int VImageDimension >
VectorField< T, VImageDimension >::VectorField() : VectorImage<T, VImageDimension >::VectorImage()
{}


template <class T, unsigned int VImageDimension >
VectorField< T, VImageDimension >::VectorField(unsigned int sizeX)
  : VectorImage<T, VImageDimension >::VectorImage(sizeX, 1) {}


template <class T, unsigned int VImageDimension >
VectorField< T, VImageDimension >::VectorField(unsigned int sizeX, unsigned int sizeY )
  : VectorImage<T, VImageDimension >::VectorImage(sizeX, sizeY, 2) {}


template <class T, unsigned int VImageDimension >
VectorField< T, VImageDimension >::VectorField(unsigned int sizeX, unsigned int sizeY, unsigned int sizeZ)
  : VectorImage<T, VImageDimension >::VectorImage(sizeX, sizeY, sizeZ, 3) {}


template <class T, unsigned int VImageDimension >
VectorField< T, VImageDimension >::VectorField( const VectorField* source )
  : VectorImage<T, VImageDimension >::VectorImage(source) {}


template <class T, unsigned int VImageDimension >
VectorField< T, VImageDimension >::VectorField( const VectorImageType* source, T dVal )
{
  this->m_SizeX = source->GetSizeX();
  this->m_SizeY = source->GetSizeY();
  this->m_SizeZ = source->GetSizeZ();
  this->m_Dimension = VImageDimension;
  this->m_Length = this->m_SizeX*this->m_SizeY*this->m_SizeZ*this->m_Dimension;
  this->m_DataPtr = NULL;
  this->Allocate();

  this->m_SpacingX = source->GetSpacingX();
  this->m_SpacingY = source->GetSpacingY();
  this->m_SpacingZ = source->GetSpacingZ();

  this->SetOrigin( source->GetOrigin() );
  this->SetDirection( source->GetDirection() );

   // initialize everthing with zero, vector image is just used to get the dimensions and spacings
  for (unsigned int uiIndx = 0; uiIndx < this->m_Length; ++uiIndx )
    {
    this->m_DataPtr[ uiIndx ] = dVal;
    }

}

////////////////////
// Public Methods //
////////////////////

//
// 1D get x
//
template <class T, unsigned int VImageDimension >
T VectorField< T, VImageDimension >::GetX(unsigned int x) const
{
  return this->GetValue(x,0);
}

//
// 2D get x
//
template <class T, unsigned int VImageDimension >
T VectorField< T, VImageDimension  >::GetX(unsigned int x, unsigned int y) const
{
  return this->GetValue(x,y,0);
}

//
// 3D get x
//
template <class T, unsigned int VImageDimension >
T VectorField< T, VImageDimension >::GetX(unsigned int x, unsigned int y, unsigned int z) const
{
  return this->GetValue(x,y,z,0);
}

//
// 2D get y
//
template <class T, unsigned int VImageDimension >
T VectorField< T, VImageDimension >::GetY(unsigned int x, unsigned int y) const
{
  return this->GetValue(x,y,1);
}

//
// 3D get y
//
template <class T, unsigned int VImageDimension >
T VectorField< T, VImageDimension >::GetY(unsigned int x, unsigned int y, unsigned int z) const
{
  return this->GetValue(x,y,z,1);
}

//
// 3D get z
//
template <class T, unsigned int VImageDimension >
T VectorField< T, VImageDimension >::GetZ(unsigned int x, unsigned int y, unsigned int z) const
{
  return this->GetValue(x,y,z,2);
}

//
// 1D set x
//
template <class T, unsigned int VImageDimension >
void VectorField< T, VImageDimension  >::SetX(unsigned int x, T value) {

  this->SetValue(x,0, value);
}

//
// 2D set x
//
template <class T, unsigned int VImageDimension >
void VectorField< T, VImageDimension >::SetX(unsigned int x, unsigned int y, T value) {

  this->SetValue(x,y,0, value);
}

//
// 3D set x
//
template <class T, unsigned int VImageDimension >
void VectorField< T, VImageDimension >::SetX(unsigned int x, unsigned int y, unsigned int z, T value) {

  this->SetValue(x,y,z,0, value);
}

//
// 2D set y
//
template <class T, unsigned int VImageDimension >
void VectorField< T, VImageDimension >::SetY(unsigned int x, unsigned int y, T value) {

  this->SetValue(x,y,1, value);
}

//
// 3D set y
//
template <class T, unsigned int VImageDimension >
void VectorField< T, VImageDimension >::SetY(unsigned int x, unsigned int y, unsigned int z, T value) {

  this->SetValue(x,y,z,1, value);
}

//
// 3D set z
//
template <class T, unsigned int VImageDimension >
void VectorField< T, VImageDimension >::SetZ(unsigned int x, unsigned int y, unsigned int z, T value) {

  this->SetValue(x,y,z,2, value);
}

//
// Get Displacement
//
template <class T, unsigned int VImageDimension >
T VectorField< T, VImageDimension  >::GetDisplacement()
{

    T mean_value = 0;

    if(this->m_Dimension == 1)
    {
        for (unsigned int uiIndx = 0; uiIndx < this->m_SizeX; ++uiIndx )
        {
            T x_value = this->GetX(uiIndx);
            mean_value += x_value;
        }
        return mean_value/this->m_SizeX;
    }
    else if(this->m_Dimension == 2)
    {
        for (unsigned int uiIndx = 0; uiIndx < this->m_SizeX; ++uiIndx )
        {
            for (unsigned int uiIndy = 0; uiIndy < this->m_SizeY; ++uiIndy )
            {
                T x_value = this->GetX(uiIndx,uiIndy);
                T y_value = this->GetY(uiIndx,uiIndy);
                T norm_value = sqrt( x_value*x_value + y_value*y_value );
                mean_value += norm_value;
            }
        }
        return mean_value/(this->m_SizeX*this->m_SizeY);
    }
    else if(this->m_Dimension == 3)
    {
        for (unsigned int uiIndx = 0; uiIndx < this->m_SizeX; ++uiIndx )
        {
            for (unsigned int uiIndy = 0; uiIndy < this->m_SizeY; ++uiIndy )
            {
                for (unsigned int uiIndz = 0; uiIndz < this->m_SizeZ; ++uiIndz )
                {
                    T x_value = this->GetX(uiIndx,uiIndy,uiIndz);
                    T y_value = this->GetY(uiIndx,uiIndy,uiIndz);
                    T z_value = this->GetZ(uiIndx,uiIndy,uiIndz);
                    T norm_value = sqrt( x_value*x_value + y_value*y_value + z_value*z_value);
                    mean_value += norm_value;
                }
            }
        }
        return mean_value/(this->m_SizeX*this->m_SizeY*this->m_SizeZ);
    }
}

} // end namespace CALATK

#endif
