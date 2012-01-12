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

#ifndef C_LDDMM_SIMPIFIED_GEODESIC_SHOOTING_OBJECTIVE_FUNCTION_TXX
#define C_LDDMM_SIMPIFIED_GEODESIC_SHOOTING_OBJECTIVE_FUNCTION_TXX

template < class TState >
CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState >::CLDDMMSimplifiedGeodesicShootingObjectiveFunction()
{
  // storage for the map
  m_ptrMapIn = NULL;
  m_ptrMapOut = NULL;
  m_ptrMapTmp = NULL;

  m_ptrDeterminantOfJacobian = NULL;
}

template< class TState >
void CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState >::DeleteData()
{
  this->m_ptrKernel->DeallocateMemory();

  SaveDelete< VectorFieldPointerType >::Pointer( m_ptrMapIn );
  SaveDelete< VectorFieldPointerType >::Pointer( m_ptrMapOut );
  SaveDelete< VectorFieldPointerType >::Pointer( m_ptrMapTmp );

  SaveDelete< VectorImagePointerType >::Pointer( m_ptrDeterminantOfJacobian );

  m_vecMeasurementTimepoints.clear();
  m_vecTimeDiscretization.clear();
  m_vecTimeIncrements.clear();

  SaveDelete< TState* >::Pointer( this->m_pState );
  SaveDelete< TState* >::Pointer( this->m_pGradient );
}

template < class TState >
CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState >::~CLDDMMSimplifiedGeodesicShootingObjectiveFunction()
{
  DeleteData();
}

template < class TState >
void CLDDMMSimpleGeodesicShootingObjectiveFunction< TState >::CreateNewStateStructures()
{
    assert( this->m_pState == NULL );
    assert( m_vecTimeDiscretization.size() > 1 );

    // get the subject ids
    std::vector< unsigned int > vecSubjectIndices;
    this->m_ptrImageManager->GetAvailableSubjectIndices( vecSubjectIndices );

    assert( vecSubjectIndices.size()>0 );

    // obtain image from which to graft the image information for the data structures

    SImageInformation* pImInfo;
    // get information from the first image to figure out the dimensions
    this->m_ptrImageManager->GetPointerToSubjectImageInformationByIndex( pImInfo, vecSubjectIndices[0], 0 );

    VectorImageType* ptrInitialImage = new VectorImageType( pImInfo->pIm );
    VectorImageType* ptrInitialMomentum = new VectorImageType( pImInfo->pIm );
    ptrInitialMomentum->setConst(0);

    this->m_pState = new TState( ptrInitialImage, ptrInitialMomentum );
}

template< class TState >
void CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState >::ShallowCopyStateStructures( TState* pState )
{
    assert ( this->m_pState == NULL );

    VectorImageType* ptrInitialImage = pState->GetPointerToInitialImage();
    VectorImageType* ptrInitialMomentum = pState->GetPointerToInitialMomentum();

    this->m_pState = new TState( ptrInitialImage, ptrInitialMomentum );
}

template < class TState >
void CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState>::CreateGradientAndAuxiliaryStructures()
{

    // get the subject ids
    std::vector< unsigned int > vecSubjectIndices;
    this->m_ptrImageManager->GetAvailableSubjectIndices( vecSubjectIndices );

    assert( vecSubjectIndices.size()>0 );

    // obtain image from which to graft the image information for the data structures

    SImageInformation* pImInfo;
    // get information from the first image to figure out the dimensions
    this->m_ptrImageManager->GetPointerToSubjectImageInformationByIndex( pImInfo, vecSubjectIndices[0], 0 );

    // create the gradient
    VectorImageType* ptrI0Gradient = new VectorImageType( pImInfo->pIm );
    ptrI0Gradient->setConst(0);

    VectorImageType* ptrP0Gradient = new VectorImageType( pImInfo->pIm );
    ptrP0Gradient->setConst(0);

    this->m_pGradient = new TState( ptrI0Gradient, ptrP0Gradient );

    // storage for the maps

    m_ptrMapIn = new VectorFieldType( pImInfo->pIm );
    m_ptrMapOut = new VectorFieldType( pImInfo->pIm );
    m_ptrMapTmp = new VectorFieldType( pImInfo->pIm );

    // storage for the determinant of the Jacobian
    m_ptrDeterminantOfJacobian = new VectorImageType( pImInfo0>pIm );

}

template < class TState >
void CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState>::InitializeDataStructuresFromState( TState* pState )
{
    DeleteData();

    CreateTimeDiscretization();

    ShallowCopyStateStructures( pState );

    CreateGradientAndAuxiliaryStructures();

}

template < class TState >
void CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState >::InitializeDataStructures()
{
    DeleteData();

    assert (this->m_pGradient == NULL );

    CreateTimeDiscretization();

    // allocate state structures
    CreateNewStateStructures();

    // gradient and everything else
    CreateGradientAndAuxiliaryStructures();
}

template < class TState >
void CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState >::InitializeState()
{
  InitializeDataStructures();
}

template < class TState >
void CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState >::InitializeState(TState* pState)
{
  InitializeDataStructuresFromState( pState );
}

template < class TState >
void CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState >::GetMomentum( VectorImageType* ptrMomentum, T dTime )
{
  GetMap( m_ptrMapTmp, dTime );
  // now compute the momentum by interpolation
  VectorImageType* ptrInitialMomentum = this->m_pState->GetPointerToInitialMomentum();
  LDDMMUtils< T, TState::VImageDimension >::applyMap( m_ptrMapTmp, ptrInitialMomentum, ptrMomentum );
  LDDMMUtils< T, TState::VImageDimension >::computeDeterminantOfJacobian( m_ptrMapTmp, m_ptrDeterminantOfJacobian );
  ptrMomentum->multElementwise( m_ptrDeterminantOfJacobian );

}

template < class TState >
void CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState >::GetImage( VectorImageType* ptrIm, T dTime )
{
  // TODO: account for appearance changes, based on closeby images
  GetMap( m_ptrMapTmp, dTime );
  // now compute the image by interpolation
  VectorImageType* ptrInitialImage = this->m_pState->GetPointerToInitialImage();
  LDDMMUtils< T, TState::VImageDimension >::applyMap( m_ptrMapTmp, ptrInitialImage, ptrIm );

}

template < class TState >
void CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState >::GetMap( VectorFieldType* ptrMap, T dTime )
{
  T dTimeFrom = m_vecTimeDiscretization[0].dTime;
  GetMapFromTo( ptrMap, dTimeFrom, dTime );
}

template < class TState >
void CLDDMMAdjointGeodesicShootingObjectiveFunction< TState >::GetMapFromTo( VectorFieldType* ptrMap, T dTimeFrom, T dTimeTo )
{
  // TODO: implement, need to integrate to get the map
  /*CALATK::LDDMMUtils< T, TState::VImageDimension >::GetMapFromToFromSpatioTemporalVelocityField(
        ptrMap,
        dTimeFrom,
        dTimeTo,
        m_vecTimeDiscretization,
        m_ptrVelocityField,
        this->m_ptrEvolver );*/
}

template< class TState >
void CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState >::CreateTimeDiscretization()
{
    if ( this->m_ptrImageManager == NULL )
    {
      throw std::runtime_error( "ERROR: No image manager specified." );
      return;
    }

    std::vector< STimePoint > vecTimePointData;
    unsigned int uiNrOfMeasurements = CALATK::LDDMMUtils< T, TState::VImageDimension >::DetermineTimeSeriesTimePointData( this->m_ptrImageManager, 0, vecTimePointData );

    if ( uiNrOfMeasurements != 2 )
    {
      throw std::runtime_error( "CLDDMMSimplifiedGeodesicShootingObjectiveFunction only supports two measurements (one source and one target image" );
      return;
    }

    CALATK::LDDMMUtils< T, TState::VImageDimension >::CreateTimeDiscretization( vecTimePointData, m_vecTimeDiscretization, m_vecTimeIncrements, this->m_NumberOfDiscretizationVolumesPerUnitTime );

    // now add the weights, default weights are all constants here
    // TODO: make this more flexible to support local geodesic regression

    typename std::vector< STimePoint >::iterator iter;
    for ( iter = m_vecTimeDiscretization.begin(); iter != m_vecTimeDiscretization.end(); ++iter )
    {
      iter->vecWeights.clear();
      for ( unsigned int iI=0; iI<iter->vecMeasurementImages.size(); ++iI )
      {
        // TODO: make this at least a constant variable
        iter->vecWeights.push_back( 1.0/this->m_SigmaSqr );
      }
    }
    // the time discretization vector has all the N timepoint. There will be N-1 vector fields in between
}

template < class TState >
void CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState >::ComputeImageMomentumAndMapForward( )
{
  /**
    * Solves the EPDiff equation forward in time using a map-based approach
    * and computes the map, mapping to the initial condition
    *
    * \f$ I_t + \nabla I^T v = 0, \f$
    * \f$ p_t + div( p v ) = 0, \f$
    * \f$ v = -K*(p\nabla I) \f$
    *
    */

  VectorImageType* ptrInitialImage = this->m_pState->GetPointerToInitialImage();
  VectorImageType* ptrInitialMomentum = this->m_pState->GetPointerToInitialMomentum();

  LDDMMUtils< T, TState::VImageDimension>::identityMap( m_ptrMapIn );

  m_ptrCurrentI->copy( ptrInitialImage );
  m_ptrCurrentP->copy( ptrInitialMomentum );

  LDDMMUtils< T, TState::VImageDimension>::identityMap( m_ptrCurrentBackMap );
  LDDMMUtils< T, TState::VImageDimension>::identityMap( m_ptrMapIdentity );


  for ( unsigned int iI = 0; iI < m_vecTimeDiscretization.size()-1; iI++ )
  {
      ComputeVelocity( m_ptrCurrentI, m_ptrCurrentP, m_ptrCurrentVelocity, m_ptrTmpField );

      this->m_ptrEvolver->SolveForward( m_ptrCurrentVelocity, m_ptrMapIn, m_ptrMapOut, m_ptrMapTmp, this->m_vecTimeIncrements[ iI ] );

      LDDMMUtils< T, TState::VImageDimension >::applyMap( m_ptrMapOut, ptrInitialImage, m_ptrCurrentI );
      LDDMMUtils< T, TState::VImageDimension >::applyMap( m_ptrMapOut, ptrInitialMomentum, m_ptrCurrentP );

      LDDMMUtils< T, TState::VImageDimension >::computeDeterminantOfJacobian( m_ptrMapOut, m_ptrDeterminantOfJacobian );

      m_ptrCurrentP->multElementwise( m_ptrDeterminantOfJacobian );

      m_ptrMapIn->copy( m_ptrMapOut );

      // keep track of the backward map
      // negate the velocity field
      m_ptrCurrentVelocity->multConst( -1.0 );

      this->m_ptrEvolver->SolveForward( m_ptrCurrentVelocityField, m_ptrMapIdentity, m_ptrMapBackIncremental, m_ptrMapTmp, this->m_vecTimeIncrements[ iI ] );

      m_ptrMapTmp->copy( m_ptrCurrentBackMap );
      LDDMMUtils< T, TState::VImageDimension >::applyMap( m_ptrMapTmp, m_ptrMapBackIncemental, m_ptrCurrentBackMap )
  }

  // TODO: implement this for multiple time points

  // now we have the back map and can compute the adjoint in the initial frame by pulling it back

  this->m_pMetric->GetAdjointMatchingDifferenceImage( m_ptrCurrentFinalAdjoint, m_ptrCurrentI, ptrI1 );
  m_ptrCurrentAdjointIDifference->multConst( m_vecTimeDiscretization[uiNrOfTimePoints-1].vecWeights[iM] );

  LDDMMUtils< T, TState::VImageDimension >::applyMap( m_ptrCurrentBackMap, ptr_CurrentFinalAdjoint, ptr_WarpedFinalToInitialAdjoint );

}

template < class TState >
void CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState>::ComputeGradient()
{
  ComputeImageMomentumAndMapForward();

  /**
    * The gradient is now simply
    *
    \f[
      \nabla_{I_0(t_0)}E = -\lambda(0) -\frac{1}{\sigma^2}\nabla d^2(I(0),I_0)
    \f]
    *
    \f[
      \nabla_{p(t_0)}E = (\lambda(0)-p(0))
    \f]
    */

  VectorImageType* ptrInitialImage = this->m_pState->GetPointerToInitialImage();
  VectorImageType* ptrInitialMomentum = this->m_pState->GetPointerToInitialMomentum();

  VectorImageType* ptrI0Gradient = this->m_pGradient->GetPointerToInitialImage();
  VectorImageType* ptrP0Gradient = this->m_pGradient->GetPointerToInitialMomentum();

  ptrP0Gradient->copy( ptrInitialMomentum );
  ptrP0Gradient->multConst(-1);

  ptrP0Gradient->add( m_ptrWarpedFinalToIntialAdjoint );

  if ( this->m_EstimateInitialImage )
  {
    ptrI0Gradient->copy( m_ptrWarpedFinalToInitialAdjoint );
    ptrI0Gradient->multConst(-1);
  }
  else
  {
    ptrI0Gradient->setConst( 0.0 );
  }
}

template < class TState >
typename CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState >::CEnergyValues
CLDDMMSimplifiedGeodesicShootingObjectiveFunction< TState>::GetCurrentEnergy()
{
  /**
    * Computes the energy for the simplified shooting method.
    * Since everything is determined by the initial condition and the geodesics are explicitly
    * enforced we can compute the energy just as for the full adjoint formulation
    *
    \f[
      E = 0.5 \langle p(t_0)\nabla I(t_0),K*(p(t_0)\nabla I(t_0)\rangle + \frac{1}{\sigma}^2 d^2(I(t_1),Y)
    \f]
    */



}

#endif
