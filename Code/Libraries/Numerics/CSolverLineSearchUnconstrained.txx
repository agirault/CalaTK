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

#ifndef C_SOLVER_LINE_SEARCH_UNCONSTRAINED_TXX
#define C_SOLVER_LINE_SEARCH_UNCONSTRAINED_TXX

template < class TState >
CSolverLineSearchUnconstrained< TState >::CSolverLineSearchUnconstrained()
{
  this->m_TempState = NULL;
}

//
// minimizes the objective function
//
template < class TState >
bool CSolverLineSearchUnconstrained< TState >::SolvePreInitialized()
{

    std::cout << "[Unconstrained linesearch] Max iterations = " << this->m_MaxNumberOfIterations<< std::endl;

    unsigned int cpt_f = 0;
    unsigned int cpt_g = 0;

  ObjectiveFunctionType * objectiveFunction = this->GetObjectiveFunction();

  unsigned int uiNrOfIterationsWithImmediateDecrease = 0;
  unsigned int uiNrOfIterationsWithoutImmediateDecrease = 0;

  //* Linesearch Inputs *//
  // Energy
  CEnergyValues InitialEnergy = objectiveFunction->GetCurrentEnergy(); cpt_f++;
  std::cout << "Initial energy = " << InitialEnergy.dEnergy << std::endl;
  CEnergyValues CurrentEnergy = InitialEnergy;
  // Step Size
  T dDesiredStepSize = this->m_InitialStepSize;
  // New temp state
  this->m_TempState = new TState( *objectiveFunction->GetStatePointer() );

  // Output the initial state if desired
  //std::string sStatePrefix = "S" + CreateIntegerString( (int)this->GetExternalSolverState() ) + "-";
  //this->OutputStateInformation( 0, sStatePrefix );

  // Linesearch Outputs
  T dAlpha;
  CEnergyValues ResultingEnergy;
  unsigned int uiRequiredIterations;

  for ( unsigned int uiIter = 0; uiIter < this->m_MaxNumberOfIterations; ++uiIter )
    {

    std::cout << std::endl << ">> ITERATION " << uiIter+1 << std::endl;
    /************************ BIG CHANGEMENT ************************************/
    //bool bSufficientlyDecreasedEnergy = this->LineSearchWithBacktracking( CurrentEnergy, dDesiredStepSize, dAlpha, ResultingEnergy, uiRequiredIterations, this->m_TempState );
    bool bSufficientlyDecreasedEnergy;

    // get current energy
    CEnergyValues InitialEnergy = CurrentEnergy;
    CEnergyValues ComputedEnergy;

    T dAdjustedEnergy = std::numeric_limits< T >::infinity();

    // save the current state
    TState* pTempState = this->m_TempState;
    *pTempState = *objectiveFunction->GetStatePointer();

    // get a pointer to the state (which will be updated throughout the iterations)
    TState *pState = objectiveFunction->GetStatePointer();

    // compute the current gradient
    objectiveFunction->ComputeGradient(); cpt_g++;

    // get current gradient
    TState *pCurrentGradient = objectiveFunction->GetGradientPointer();

    // compute the norm of the gradient (required for line search with gradient descent)
    T dSquaredNorm = pCurrentGradient->SquaredNorm();

    //std::cout << "dSquaredNorm = " << dSquaredNorm << std::endl;

    // now see if we can reduce the energy by backtracking
    // FIXME: Add sufficient decrease condition: for now just see if it is decreasing

    dAlpha = dDesiredStepSize;

    bool bHitLowerStepSizeBound = false;
    bool bTerminate = false;
    uiRequiredIterations = 0;

    do
      {

      // if it has reached the smallest allowable step size
      if ( bHitLowerStepSizeBound ) bTerminate = true;

      // doing the gradient step
      //*pState = *pTempState - (*pCurrentGradient)*dAlpha;
      // here comes a more memory efficient version
      // (should need no new reallocation of memory, but simply overwrites *pState all the time)

      *pState = *pCurrentGradient;
      *pState *= -dAlpha;
      *pState += *pTempState;

      // recompute the energy

      ComputedEnergy = objectiveFunction->GetCurrentEnergy(); cpt_f++;

      /*std::cout << "initE = " << dInitialEnergy << std::endl;
      std::cout << "dc = " << m_DecreaseConstant << std::endl;
      std::cout << "alpha = " << dAlpha << std::endl;
      std::cout << "sqNorm = " << dSquaredNorm << std::endl;*/

      dAdjustedEnergy = InitialEnergy.dEnergy - this->m_DecreaseConstant*dAlpha*dSquaredNorm;

      //std::cout << "computed energy = " << dComputedEnergy << "; dAdjustedEnergy = " << dAdjustedEnergy << std::endl;
      //std::cout << "dComputedEnergy = " << dComputedEnergy << std::endl;

      if ( ComputedEnergy.dEnergy >= dAdjustedEnergy )
        {
        dAlpha *= this->m_ReductionFactor;
        if ( dAlpha < this->m_MinAllowedStepSize )
          {
          dAlpha = this->m_MinAllowedStepSize;
          bHitLowerStepSizeBound = true;
          }
        }

      uiRequiredIterations++;

      }
    while ( ( ComputedEnergy.dEnergy >= dAdjustedEnergy ) && ( uiRequiredIterations <= this->m_MaxNumberOfTries ) && ( !bTerminate ) );

    if ( ComputedEnergy.dEnergy > dAdjustedEnergy )
      {
      // could not reduce the energy, so keep the original one
      *pState = *pTempState;
      ResultingEnergy = InitialEnergy;
      bSufficientlyDecreasedEnergy = false;
      }
    else
      {
      // energy was successfully reduced, we can keep the updated state (in pState)
      ResultingEnergy = ComputedEnergy;
      bSufficientlyDecreasedEnergy = true;
      }

    /************************************************************************************************************/

    // Output the current energy information //std::setw(10)
    std::cout << "reqIter =  " << uiRequiredIterations << std::endl;
    std::cout << "Alpha   =  " << dAlpha << std::endl;
    std::cout << "E(tot)  =  " << ResultingEnergy.dEnergy << std::endl;
    std::cout << "E(I)    =  " << ResultingEnergy.dMatchingEnergy << std::endl;
    std::cout << "E(v)    =  " << ResultingEnergy.dRegularizationEnergy << std::endl;

    //Output the state if desired
    //this->OutputStateInformation( uiIter + 1, sStatePrefix );

    if ( bSufficientlyDecreasedEnergy )
      {

        CurrentEnergy = ResultingEnergy;

      if ( dDesiredStepSize == dAlpha ) // could be decreased immediately
        {
        uiNrOfIterationsWithImmediateDecrease++;
        uiNrOfIterationsWithoutImmediateDecrease = 0;
        }
      else
        {
        uiNrOfIterationsWithImmediateDecrease = 0;
        uiNrOfIterationsWithoutImmediateDecrease++;
        }

      if ( uiNrOfIterationsWithImmediateDecrease >= this->m_AdjustStepSizeUpNumber )
        {
        std::cout << "Adjusting step size up" << std::endl;
        dDesiredStepSize *= this->m_AdjustStepSizeUpFactor;
        uiNrOfIterationsWithImmediateDecrease = 0;
        uiNrOfIterationsWithoutImmediateDecrease = 0;
        }

      if ( uiNrOfIterationsWithoutImmediateDecrease >= this->m_AdjustStepSizeDownNumber )
        {
        std::cout << "Adjusting step size down" << std::endl;
        dDesiredStepSize *= this->m_AdjustStepSizeDownFactor;
        uiNrOfIterationsWithImmediateDecrease = 0;
        uiNrOfIterationsWithoutImmediateDecrease = 0;
        }

      }
    else // could not decrease energy
      {
      // terminate if smallest step size has been tried
      if ( dAlpha == this->m_MinAllowedStepSize )
        {
        std::cout << "Smallest allowable step size did not yield an energy reduction. Stopping iterations." << std::endl;
        break;
        }
      else
        {
        uiNrOfIterationsWithImmediateDecrease = 0;
        uiNrOfIterationsWithoutImmediateDecrease++;
        // set the desired step size to the last tried one
        std::cout << "Could not decrease energy. Trying again with smaller step size." << std::endl;
        dDesiredStepSize = dAlpha;
        }

      }
    }

  // write RegularizationEnergy
    std::ofstream file;
    file.open ("RegEnergy.txt");
    file <<ResultingEnergy.dRegularizationEnergy;
    file.close();

    std::cout << std::endl << "[Unconstrained linesearch] End of minimization" << std::endl; // COUT
    std::cout << "# of Energies comp  =  " << cpt_f << std::endl;
    std::cout << "# of Gradient comp  =  " << cpt_g << std::endl;


  // clean up
  if ( ResultingEnergy.dEnergy < InitialEnergy.dEnergy )
    {
    // could reduce the energy
    return true;
    }
  else
    {
    // could not reduce the energy
    return false;
    }

}


#endif
