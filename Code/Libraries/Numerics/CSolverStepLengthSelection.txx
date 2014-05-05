
#ifndef C_SOLVER_STEP_LENGTH_SELECTION_TXX
#define C_SOLVER_STEP_LENGTH_SELECTION_TXX

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>

template < class TState >
CSolverStepLengthSelection< TState>::CSolverStepLengthSelection()
  : DefaultMinGradAllowed( 0.1 ), //TODO analyze this value
    DefaultMinDisplacementAllowed( 0.5 ), //TODO and this value (link to resolution/sizing?)
    DefaultDecreaseConstant( 0.0001 ),
    DefaultMaxNumberOfIterations( 100 ),
    DefaultMaxNumberOfTries( 10 ),
    m_ExternallySetMinGradAllowed( false ),
    m_ExternallySetMinDisplacementAllowed( false ),
    m_ExternallySetDecreaseConstant( false ),
    m_ExternallySetMaxNumberOfIterations( false ),
    m_ExternallySetMaxNumberOfTries( false )
{
  // default setting for the parameters
  // minimum value allowed for the norm of the gradient
  m_MinGradAllowed = DefaultMinGradAllowed;

  // minimum value allowed for the displacement
  m_MinDisplacementAllowed = DefaultMinDisplacementAllowed;

  // constant for the sufficient decrease condition
  m_DecreaseConstant = DefaultDecreaseConstant;

  // maximal number of iterations
  m_MaxNumberOfIterations = DefaultMaxNumberOfIterations;

  // maximal number of tries in one line search
  m_MaxNumberOfTries = DefaultMaxNumberOfTries;

}


template < class TState >
CSolverStepLengthSelection< TState>::~CSolverStepLengthSelection()
{
}


template < class TState >
void CSolverStepLengthSelection< TState>::SetAutoConfiguration( CJSONConfiguration * combined, CJSONConfiguration * cleaned )
{
  Superclass::SetAutoConfiguration( combined, cleaned );

  Json::Value& currentConfigurationIn = this->m_CombinedJSONConfig->GetFromKey( "LineSearch", Json::nullValue );
  Json::Value& currentConfigurationOut = this->m_CleanedJSONConfig->GetFromKey( "LineSearch", Json::nullValue, CONF_NORMAL );

  SetJSONHelpForRootKey( LineSearch, "setting for the linesearch algorithm", CONF_NORMAL );

  SetJSONFromKeyDouble( currentConfigurationIn, currentConfigurationOut, MinGradAllowed, CONF_ADVANCED );
  SetJSONFromKeyDouble( currentConfigurationIn, currentConfigurationOut, MinDisplacementAllowed, CONF_ADVANCED );
  SetJSONFromKeyDouble( currentConfigurationIn, currentConfigurationOut, DecreaseConstant, CONF_EXPERT );
  SetJSONFromKeyUInt( currentConfigurationIn, currentConfigurationOut, MaxNumberOfIterations, CONF_NORMAL );
  SetJSONFromKeyUInt( currentConfigurationIn, currentConfigurationOut, MaxNumberOfTries, CONF_ADVANCED );
  SetJSONFromKeyBool( currentConfigurationIn, currentConfigurationOut, OutputStateInformation, CONF_EXPERT );
  SetJSONFromKeyUInt( currentConfigurationIn, currentConfigurationOut, OutputStateInformationFrequency, CONF_EXPERT );

  SetJSONHelpForKey( currentConfigurationIn, currentConfigurationOut, MinGradAllowed,
                     "minimum value allowed for the norm of the gradient", CONF_ADVANCED );
  SetJSONHelpForKey( currentConfigurationIn, currentConfigurationOut, MinDisplacementAllowed,
                     "minimum value allowed for the displacement", CONF_ADVANCED );
  SetJSONHelpForKey( currentConfigurationIn, currentConfigurationOut, DecreaseConstant,
                     "require sufficient decrease of energy", CONF_EXPERT );
  SetJSONHelpForKey( currentConfigurationIn, currentConfigurationOut, MaxNumberOfIterations,
                     "maximal number of iterations", CONF_NORMAL );
  SetJSONHelpForKey( currentConfigurationIn, currentConfigurationOut, MaxNumberOfTries,
                     "number of tries to find a smaller solution", CONF_ADVANCED );
  SetJSONHelpForKey( currentConfigurationIn, currentConfigurationOut, OutputStateInformation,
                     "output internal information of algorithm (for debugging)", CONF_EXPERT );
  SetJSONHelpForKey( currentConfigurationIn, currentConfigurationOut, OutputStateInformationFrequency,
                     "controls at which iterations the state information should be output", CONF_EXPERT );
}


template < class TState >
bool CSolverStepLengthSelection< TState>::SolvePreInitialized()
{
    std::cout << std::endl << "[Step Length Selection linesearch] Max iterations = " << this->m_MaxNumberOfIterations<< std::endl; // COUT

    unsigned int cpt_f = 0;
    unsigned int cpt_g = 0;

    /* Current values for x, f(x) and g(x) */
    ObjectiveFunctionType * f = this->GetObjectiveFunction();
    TState x_cur = *f->GetStatePointer();
    T f_cur = f->GetCurrentEnergy().dEnergy; cpt_f++;
    VectorFieldType * map_cur = new VectorFieldType(f->GetCurrentMap());
    f->ComputeGradient(); cpt_g++;
    TState g_cur = *f->GetGradientPointer();
    T g_norm = g_cur.SquaredNorm();

    /* variables to store previous values */
    T f_prev;
    T f_new;
    T alpha_cur;
    T alpha_prev;
    T alpha_new;
    T displacement;

    std::cout << "Energy Init   =  " << f_cur << std::endl; // COUT

    /* Iterations */
    unsigned int it_count = 0;
    do
    {
        it_count++;
        std::cout << std::endl << ">> ITERATION " << it_count << std::endl; // COUT

        /* compute alpha0 */
        if (it_count == 1) alpha_cur = 10.0/g_norm;
        else alpha_cur *=5;
        if ( alpha_cur > 0.1 ) alpha_cur = 0.1; // to avoid big steps
        else if ( alpha_cur < 0.00001 ) alpha_cur = 0.00001; // to avoid small steps (TODO)

        std::cout << "Alpha test    =  " << alpha_cur << std::endl; // COUT

        /* get the new position x and its energy f(x) */
        TState *x_new = f->GetStatePointer();
        *x_new = x_cur - g_cur * alpha_cur;
        f_new = f->GetCurrentEnergy().dEnergy; cpt_f++;

        std::cout << "Energy value  =  " << f_new << std::endl; // COUT

        /* Test Armijo condition (tries) */
        unsigned int it_armijo = 0;
        while( (f_new > f_cur - m_DecreaseConstant * alpha_cur * g_norm) && (it_armijo < m_MaxNumberOfTries) ) //TODO : is it correct for g_norm?
        {
            it_armijo++;
            std::cout << ">>> TRY " << it_armijo << std::endl; // COUT

            /* compute new alpha with step length selection */
            if(it_armijo == 1) //quadratic model
            {
                alpha_new = ( g_norm * alpha_cur * alpha_cur );
                alpha_new /= 2 * ( f_new - f_cur + g_norm * alpha_cur );
            }
            else // cubic model
            {
                vnl_matrix_fixed<T,2,2> M1;
                M1(0,0) = alpha_prev * alpha_prev;
                M1(0,1) = - alpha_cur * alpha_cur;
                M1(1,0) = - alpha_prev * alpha_prev * alpha_prev;
                M1(1,1) = alpha_cur * alpha_cur * alpha_cur;
                vnl_vector_fixed<T,2> M2( f_new - f_cur + g_norm * alpha_cur ,  f_prev - f_cur + g_norm * alpha_prev );
                vnl_vector_fixed<T,2> M = M1 * M2 / ( alpha_prev*alpha_prev * alpha_cur*alpha_cur * ( alpha_cur - alpha_prev ) );
                T a = M(0);
                T b = M(1);
                alpha_new = ( -b + sqrt(b*b + 3*a*g_norm) )/(3*a);
            }
            std::cout << "Alpha test    =  " << alpha_new << std::endl; // COUT

            /* store previous values */
            f_prev      =   f_new;
            alpha_prev  =   alpha_cur;
            alpha_cur   =   alpha_new;

            /* calculate new value to test for Armijo */
            *x_new = x_cur - g_cur * alpha_cur;
            f_new = f->GetCurrentEnergy().dEnergy; cpt_f++;

            std::cout << "Energy value  =  " << f_new << std::endl; // COUT

        }
        if( it_armijo == m_MaxNumberOfTries)
        {
            std::cout << "[!] Could not reduce Energy"<< std::endl;
            return false; // TODO : make sure this is correct.
        }

        /* Get maximum displacement */
        VectorFieldType * map_new = new VectorFieldType(f->GetCurrentMap());
        map_cur->SubtractCellwise(map_new);
        displacement = map_cur->GetDisplacement();
        delete map_cur;

        /* Update new values map, x, f(x) and g(x) */
        map_cur = map_new;
        x_cur   = *x_new;
        f_cur   = f_new;
        f->ComputeGradient(); cpt_g++;
        g_cur   =  *f->GetGradientPointer();
        g_norm  = g_cur.SquaredNorm();

        std::cout << ">> RESULTS" << std::endl; // COUT
        std::cout << "Displacement  =  " << displacement << std::endl;
        std::cout << "Gradient Norm =  " << sqrt(g_norm) << std::endl;

    }while( (displacement > m_MinDisplacementAllowed) && (g_norm > m_MinGradAllowed*m_MinGradAllowed) && (it_count < m_MaxNumberOfIterations) );

    delete map_cur;

    std::cout << std::endl << "[Step Length Selection linesearch] End of minimization" << std::endl; // COUT
    std::cout << "# of Energies comp  =  " << cpt_f << std::endl;
    std::cout << "# of Gradient comp  =  " << cpt_g << std::endl;
    return true;
}

#endif
