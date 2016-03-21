using System;
using LinearIntegrationEquation.BusinessLogic;

namespace LinearIntegrationEquation.EventArgsExtensions
{
    class MatrixEquationEventArgs:EventArgs
    {
        public MatrixEquation MatrixEquation { get; set; }

        public MatrixEquationEventArgs(MatrixEquation matrixEquation)
        {
            MatrixEquation = matrixEquation;
        }
    }
}
