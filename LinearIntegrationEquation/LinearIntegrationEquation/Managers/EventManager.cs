using System;
using LinearIntegrationEquation.EventArgsExtensions;

namespace LinearIntegrationEquation.Managers
{
    class EventManager
    {
        ////SolveButtonPressed
        //public delegate void SolveButtonPressedEventHandler(object source, EventArgs eventArgs);

        //public static event SolveButtonPressedEventHandler SolveButtonPressed;

        //public static void OnSolveButtonPressed(object source, EventArgs eventArgs)// Instead EventArgs could be any type data.
        //{
            
        //    if (SolveButtonPressed != null)
        //    {
        //        SolveButtonPressed(source, new EventArgs());
        //    }
        //}
        //Matrix is formed
        public delegate void MatrixFormedEventHandler(object source, MatrixEquationEventArgs eventArgs);

        public static event MatrixFormedEventHandler MatrixFormed;

        public static void OnMatrixFormed(object source, MatrixEquationEventArgs eventArgs)
        {
            if (MatrixFormed != null)
            {
                MatrixFormed(source, eventArgs);
            }
        }
        //Matrix is solved
        public delegate void MatrixSolvedEventHandler(object source, SolutionOfMatrixEventArgs eventArgs);

        public static event MatrixSolvedEventHandler MatrixSolved;

        public static void OnMatrixSolved(object source, SolutionOfMatrixEventArgs eventArgs)
        {
            if (MatrixSolved != null)
            {
                MatrixSolved(source, eventArgs);
            }
        }
        //Solution are formed
        public delegate void SolutionFormedEventHandler(object source, EventArgs eventArgs);

        public static event SolutionFormedEventHandler SolutionFormed;

        public static void OnSolutionFormed(object source, EventArgs eventArgs)
        {
            if (SolutionFormed != null)
            {
                SolutionFormed(source, new EventArgs());
            }
        }  
    }
}
