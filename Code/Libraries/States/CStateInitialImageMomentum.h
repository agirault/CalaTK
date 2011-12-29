/*
*
*  Copyright 2011 by the CALATK development team
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

#ifndef C_STATE_INITIAL_IMAGE_MOMENTUM_H
#define C_STATE_INITIAL_IMAGE_MOMENTUM_H

#include "CStateImageDomain.h"
#include "VectorField.h"
#include "CResamplerLinear.h"

namespace CALATK
{

template <class T, unsigned int VImageDimension=3, class TResampler=CResamplerLinear<T, VImageDimension> >
class CStateInitialImageMomentum: public CStateImageDomain< T, VImageDimension, TResampler >
{
public:

    /* some useful typedefs */
    typedef CStateImageDomain< T, VImageDimension, TResampler > Superclass;
    typedef typename Superclass::TState SuperclassTState;
    typedef CStateInitialImageMomentum< T, VImageDimension, TResampler > TState;

    typedef typename Superclass::VectorImageType VectorImageType;
    typedef VectorImageType* VectorImagePointerType;
    typedef const VectorImageType* ConstVectorImagePointerType;

    /**
     * Empty constructor
     */
    CStateInitialImageMomentum();

    /**
     * Constructor which takes a pointer to the initial image and initial momentum;
     * Does not copy the data, just stores the pointers to it;
     * Will destroy the data in Destructor.
     */
    CStateInitialImageMomentum( VectorImagePointerType pInitialImage, VectorImagePointerType pInitialMomentum );

    /**
      * copy constructor, creation of the image and momentum for the first time, need to allocate memory
      */
    CStateInitialImageMomentum( const CStateInitialImageMomentum & c );

    /**
      * Destructor
      */
    ~CStateInitialImageMomentum();

    /*
     * Allow for upsampling of the state
     */
    SuperclassTState* CreateUpsampledStateAndAllocateMemory( const VectorImageType* pGraftImage ) const;

    // declare operators to be able to do computations with this state, which are needed in the numerical solvers
    /**
     * assignment
     */
    CStateInitialImageMomentum & operator=( const CStateInitialImageMomentum & p);

    CStateInitialImageMomentum & operator+=( const CStateInitialImageMomentum & p);

    CStateInitialImageMomentum & operator-=( const CStateInitialImageMomentum & p);

    CStateInitialImageMomentum & operator*=( const T & p);

    CStateInitialImageMomentum operator+( const CStateInitialImageMomentum & p) const;

    CStateInitialImageMomentum operator-( const CStateInitialImageMomentum & p) const;

    CStateInitialImageMomentum operator*( const T & p) const;

    VectorImagePointerType GetPointerToInitialImage() const;
    VectorImagePointerType GetPointerToInitialMomentum() const;

    void SetPointerToInitialImage( ConstVectorImagePointerType ptrImage );
    void SetPointerToInitialMomentum( ConstVectorImagePointerType ptrMomentum );

    T SquaredNorm();

protected:
    void ClearDataStructure();
    void CopyDataStructure( ConstVectorImagePointerType ptrImage, ConstVectorImagePointerType ptrMomentum );

private:
    VectorImagePointerType m_ptrInitialImage;
    VectorImagePointerType m_ptrInitialMomentum;

};

#include "CStateInitialImageMomentum.txx"

} // end namespace

#endif // CSTATEINITIALIMAGEMOMENTUM_H
