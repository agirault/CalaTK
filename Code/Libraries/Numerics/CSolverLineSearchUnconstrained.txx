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

  ObjectiveFunctionType * objectiveFunction = this->GetObjectiveFunction();

  unsigned int uiNrOfIterationsWithImmediateDecrease = 0;
  unsigned int uiNrOfIterationsWithoutImmediateDecrease = 0;

  //* Linesearch Inputs *//
  // Energy
  CEnergyValues InitialEnergy = objectiveFunction->GetCurrentEnergy();
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
    bool bSufficientlyDecreasedEnergy = this->LineSearchWithBacktracking( CurrentEnergy, dDesiredStepSize, dAlpha, ResultingEnergy, uiRequiredIterations, this->m_TempState );

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
