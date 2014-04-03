
#ifndef C_SOLVER_STEP_LENGTH_SELECTION_H
#define C_SOLVER_STEP_LENGTH_SELECTION_H

#include "CALATKCommon.h"
#include "CObjectiveFunction.h"
#include "CSolver.h"


namespace CALATK
{
/**
 * Base class for a generic line search algorithm
 */
template < class TState >
class CSolverStepLengthSelection : public CSolver< TState >
{
public:
  /** Standard class typedefs. */
  typedef CSolverStepLengthSelection      Self;
  typedef CSolver< TState >               Superclass;
  typedef itk::SmartPointer< Self >       Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  typedef typename TState::FloatType                 T;
  typedef typename Superclass::ObjectiveFunctionType ObjectiveFunctionType;
  typedef typename Superclass::CEnergyValues         CEnergyValues;

/**
 *Constructor setting the default values
 */
  CSolverStepLengthSelection();

/**
 * Destructor
 */
  virtual ~CSolverStepLengthSelection();

  SetMacro( MinGradAllowed, T );
  GetMacro( MinGradAllowed, T );
  SetMacro( MinDisplacementAllowed, T );
  GetMacro( MinDisplacementAllowed, T );
  SetMacro( DecreaseConstant, T );
  GetMacro( DecreaseConstant, T );
  SetMacro( MaxNumberOfIterations, unsigned int );
  GetMacro( MaxNumberOfIterations, unsigned int );
  SetMacro( MaxNumberOfTries, unsigned int );
  GetMacro( MaxNumberOfTries, unsigned int );

  virtual void SetAutoConfiguration( CJSONConfiguration * combined, CJSONConfiguration * cleaned );
  bool SolvePreInitialized();

protected:

  T m_MinGradAllowed;
  T m_MinDisplacementAllowed;
  T m_DecreaseConstant;
  unsigned int m_MaxNumberOfIterations;
  unsigned int m_MaxNumberOfTries;


private:

  const T DefaultMinGradAllowed;
  const T DefaultMinDisplacementAllowed;
  const T DefaultDecreaseConstant;
  const unsigned int DefaultMaxNumberOfIterations;
  const unsigned int DefaultMaxNumberOfTries;

  bool m_ExternallySetMinGradAllowed;
  bool m_ExternallySetMinDisplacementAllowed;
  bool m_ExternallySetDecreaseConstant;
  bool m_ExternallySetMaxNumberOfIterations;
  bool m_ExternallySetMaxNumberOfTries;

};

#include "CSolverStepLengthSelection.txx"

} // end namespace

#endif
