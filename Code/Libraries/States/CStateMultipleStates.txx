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

#ifndef C_STATE_MULTIPLE_STATES_TXX
#define C_STATE_MULTIPLE_STATES_TXX

namespace CALATK
{

//
// empty constructor
//
template < class TState >
CStateMultipleStates< TState  >::CStateMultipleStates()
{
}

//
// copy constructor
//
template < class TState >
CStateMultipleStates< TState >::CStateMultipleStates( const CStateMultipleStates & c )
{
  if ( this != &c )
    {
      assert ( this->m_IndividualStatesCollection.empty() );
      typename IndividualStatesCollectionType::const_iterator iter;
      for ( iter=c.m_IndividualStatesCollection.begin(); iter!=c.m_IndividualStatesCollection.end(); ++iter )
        {
          typename IndividualStateType::Pointer copiedState = new IndividualStateType( **iter );
          this->m_IndividualStatesCollection.push_back( copiedState );
        }
    }
}

//
// copy constructor (from individual states which have been allocated externally)
//
template < class TState >
CStateMultipleStates< TState >::CStateMultipleStates( const IndividualStatesCollectionType & individualStatesCollection )
{
  typename IndividualStatesCollectionType::const_iterator iter;
  for ( iter = individualStatesCollection.begin(); iter != individualStatesCollection.end(); ++iter )
    {
      this->m_IndividualStatesCollection.push_back( *iter );
    }
}


//
// clear data structure
//
template < class TState >
void CStateMultipleStates< TState >::ClearDataStructure()
{
  this->m_IndividualStatesCollection.clear();
}

//
// destructor
//
template < class TState>
CStateMultipleStates< TState >::~CStateMultipleStates()
{
  ClearDataStructure();
}

//
// Upsampling
//
template < class TState >
typename CStateMultipleStates< TState >::Superclass*
CStateMultipleStates< TState >::CreateUpsampledStateAndAllocateMemory( const VectorImageType* graftImage ) const
{
    IndividualStatesCollectionType individualStatesCollection;

    // upsample all the individual state components
    typename IndividualStatesCollectionType::const_iterator iter;
    for ( iter = m_IndividualStatesCollection.begin(); iter != m_IndividualStatesCollection.end(); ++iter )
    {
      typename IndividualStateType::Pointer individualState = dynamic_cast< IndividualStateType * >( (*iter)->CreateUpsampledStateAndAllocateMemory( graftImage ) );
      individualStatesCollection.push_back( individualState );
    }

    Superclass * upsampledState = new Self( individualStatesCollection );

    return upsampledState;

}

// Here come the algebraic operators and assignment

template < class TState >
CStateMultipleStates< TState > &
CStateMultipleStates< TState >::operator=(const CStateMultipleStates & p )
{
  if ( this != &p )
    {
    // now do a deep copy
    
    // check if we already have the same number of elements. If so overwrite, otherwise recreate
    if ( m_IndividualStatesCollection.size() == p.m_IndividualStatesCollection.size() )
      {
      // already memory of appropriate size allocated, so just copy
      // iterate and copy
      
      typename IndividualStatesCollectionType::const_iterator iterSource;
      typename IndividualStatesCollectionType::iterator iterTarget;
      for ( iterSource = p.m_IndividualStatesCollection.begin(), iterTarget = m_IndividualStatesCollection.begin();
            iterSource != p.m_IndividualStatesCollection.end();
            ++iterSource, ++iterTarget )
        {
        // copy the current state
        *(*iterTarget) = *(*iterSource);
        }
      }
    else
      {
      // clear and then allocate, this should typically not be necessary
      std::cerr << "WARNING: reallocating memory, should already have been assigned." << std::endl;
      ClearDataStructure();

      typename IndividualStatesCollectionType::const_iterator iter;
      for ( iter=p.m_IndividualStatesCollection.begin(); iter!=p.m_IndividualStatesCollection.end(); ++iter )
        {
          typename IndividualStateType::Pointer copiedState = new IndividualStateType( **iter );
          this->m_IndividualStatesCollection.push_back( copiedState );
        }
      }
    return *this;
    }
  else
    {
    return *this;
    }
}

template < class TState >
CStateMultipleStates< TState > &
CStateMultipleStates< TState >::operator+=(const CStateMultipleStates & p )
{
  if ( m_IndividualStatesCollection.size() != p.m_IndividualStatesCollection.size() )
    {
    throw std::runtime_error( "Size mismatch of state vectors. ABORT." );
    }

  typename IndividualStatesCollectionType::const_iterator iterSource;
  typename IndividualStatesCollectionType::iterator iterTarget;
  for ( iterSource = p.m_IndividualStatesCollection.begin(), iterTarget = m_IndividualStatesCollection.begin();
        iterSource != p.m_IndividualStatesCollection.end();
        ++iterSource, ++iterTarget )
    {
    // add the source to the target
    *(*iterTarget) += *(*iterSource );
    }

  return *this;
}

template < class TState >
CStateMultipleStates< TState > &
CStateMultipleStates< TState >::operator-=(const CStateMultipleStates & p )
{

  if ( m_IndividualStatesCollection.size() != p.m_IndividualStatesCollection.size() )
    {
    throw std::runtime_error( "Size mismatch of vector of vector fields. ABORT." );
    //return;
    }

  typename IndividualStatesCollectionType::const_iterator iterSource;
  typename IndividualStatesCollectionType::iterator iterTarget;
  for ( iterSource = p.m_IndividualStatesCollection.begin(), iterTarget = m_IndividualStatesCollection.begin();
        iterSource != p.m_IndividualStatesCollection.end(), iterTarget != m_IndividualStatesCollection.end();
        ++iterSource, ++iterTarget )
    {
     // subtract the source from the target
    *(*iterTarget) -= *( *iterSource );
    }

  return *this;
}

template < class TState >
CStateMultipleStates< TState > &
CStateMultipleStates< TState >::operator*=(const FloatType & p )
{

  typename IndividualStatesCollectionType::iterator iterTarget;
  for ( iterTarget = m_IndividualStatesCollection.begin(); iterTarget != m_IndividualStatesCollection.end(); ++iterTarget )
    {
    // multiply by the value
    *(*iterTarget)*= p;
    }

  return *this;
}

template < class TState >
CStateMultipleStates< TState >
CStateMultipleStates< TState >::operator+(const CStateMultipleStates & p ) const
{
  CStateMultipleStates r = *this;
  return r += p;
}

template < class TState >
CStateMultipleStates< TState >
CStateMultipleStates< TState >::operator-(const CStateMultipleStates & p ) const
{
  CStateMultipleStates r = *this;
  return r -= p;
}

template < class TState >
CStateMultipleStates< TState >
CStateMultipleStates< TState >::operator*(const FloatType & p ) const
{
  CStateMultipleStates r = *this;
  return r*= p;
}

//
// returns one of the individual states
//
template < class TState >
typename CStateMultipleStates< TState >::IndividualStateType*
CStateMultipleStates< TState >::GetIndividualStatePointer( unsigned int idx )
{
  const size_t numberOfStates = m_IndividualStatesCollection.size();
  if ( numberOfStates == 0 || numberOfStates <= idx )
    {
    return NULL;
    }
  else
    {
    return this->m_IndividualStatesCollection[ idx ];
    }
}

//
// computes the squared norm of the state, by adding all the individual square norm components
//
template < class TState >
typename CStateMultipleStates< TState >::FloatType
CStateMultipleStates< TState >::SquaredNorm()
{
  FloatType squaredNorm = 0;

  typename IndividualStatesCollectionType::iterator iter;
  for ( iter = m_IndividualStatesCollection.begin(); iter != m_IndividualStatesCollection.end(); ++iter )
  {
    squaredNorm += (*iter)->SquaredNorm();
  }

  return squaredNorm;

}

//
// Allows to query if the state contains the initial image
//
template < class TState >
bool CStateMultipleStates< TState >::StateContainsInitialImage()
{
  return false;
}

//
// writes the state to a file
//
template < class TState >
void CStateMultipleStates< TState >::WriteStateToFile( std::string prefix, std::string postfix, int offset )
{
  typename IndividualStatesCollectionType::iterator iter;
  int iI = 0;
  for ( iter = m_IndividualStatesCollection.begin(); iter != m_IndividualStatesCollection.end(); ++iter )
  {
    (*iter)->WriteStateToFile( prefix, postfix, offset + iI );
    ++iI;
  }
}

} // end namespace CALATK

#endif
