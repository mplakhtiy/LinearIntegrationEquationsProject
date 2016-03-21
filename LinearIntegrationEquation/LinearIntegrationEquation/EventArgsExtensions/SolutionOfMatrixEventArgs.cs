using System;
namespace LinearIntegrationEquation.EventArgsExtensions
{
    class SolutionOfMatrixEventArgs:EventArgs
    {
        public double[] SolutionVector { get; set; }

        public SolutionOfMatrixEventArgs(double[] inputSolutionVector)
        {
            SolutionVector = inputSolutionVector;
        }
    }
}
