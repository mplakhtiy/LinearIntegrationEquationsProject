using System;
using LinearIntegrationEquation.EventArgsExtensions;
using EventManager = LinearIntegrationEquation.Managers.EventManager;

namespace LinearIntegrationEquation.BusinessLogic
{
    internal class MatrixSolver
    {
        private double[] returnAnswerVector;

        private delegate void Solve(double[,] inputMatrix, double[] inputVector);

        private readonly Solve solve;

        public MatrixSolver()
        {
            solve = Gaus;
            EventManager.MatrixFormed += onMatrixFormed;
        }

        private void onMatrixFormed(object source, MatrixEquationEventArgs eventArgs)
        {
            solve.BeginInvoke(eventArgs.MatrixEquation.returnMatrix, eventArgs.MatrixEquation.returnVector,
                onMatrixSolved, null);
        }

        private void onMatrixSolved(IAsyncResult asyncResult)
        {
            EventManager.OnMatrixSolved(this, new SolutionOfMatrixEventArgs(returnAnswerVector));
        }

        private void Gaus( double[,] inputMatrix, double[] inputVector)
        {
            int count = inputVector.Length;
            returnAnswerVector = new double[count];
            for (int k = 0; k < count - 1; k++)
            {
                for (int i = k + 1; i < count; i++)
                {
                    double m = -inputMatrix[i, k] / inputMatrix[k, k];
                    inputVector[i] += m * inputVector[k];
                    for (int j = k + 1; j < count; j++)
                        inputMatrix[i, j] += m * inputMatrix[k, j];
                }

            }
            returnAnswerVector[count - 1] = inputVector[count - 1] / inputMatrix[count - 1, count - 1];
            for (int k = count - 2; k >= 0; k--)
            {
                double sum = 0;
                for (int j = k + 1; j < count; j++)
                    sum += inputMatrix[k, j] * returnAnswerVector[j];
                returnAnswerVector[k] = (inputVector[k] - sum) / inputMatrix[k, k];
            }
        }

        //private void SolveWithGauss(double[,] inputMatrix, double[] inputVector)
        //{
        //    int count = inputVector.Length;
        //    int i, j;
        //    double tempValue;
        //    returnAnswerVector = new double[count];
        //    for (int k = 0; k < count; k++)
        //    {
        //        double checkValue = Math.Abs(inputMatrix[k, k]);

        //        i = k;

        //        for (int l = k + 1; l < count; l++)
        //        {
        //            if (Math.Abs(inputMatrix[l, k]) > checkValue)
        //            {
        //                i = l;
        //                checkValue = Math.Abs(inputMatrix[l, k]);
        //            }
        //        }

        //        if (checkValue < 0.0000000000000000000000000000000000001)
        //        {
        //            throw new Exception("Insoluble Matrix Equation");
        //        }

        //        if (i != k)
        //        {
        //            for (j = k; j < count; j++)
        //            {
        //                tempValue = inputMatrix[k, j];
        //                inputMatrix[k, j] = inputMatrix[i, j];
        //                inputMatrix[i, j] = tempValue;
        //            }
        //            tempValue = inputVector[k];
        //            inputVector[k] = inputVector[i];
        //            inputVector[i] = tempValue;
        //        }

        //        checkValue = inputMatrix[k, k];

        //        for (j = k; j < count; j++)
        //        {
        //            inputMatrix[k, j] = inputMatrix[k, j]/checkValue;
        //        }

        //        inputVector[k] = inputVector[k]/checkValue;

        //        for (i = k + 1; i < count; i++)
        //        {
        //            tempValue = inputMatrix[i, k];
        //            inputMatrix[i, k] = 0;
        //            if (Math.Abs(tempValue) > Math.Pow(10, -10))
        //            {
        //                for (j = k + 1; j < count; j++)
        //                    inputMatrix[i, j] = inputMatrix[i, j] - tempValue*inputMatrix[k, j];
        //                inputVector[i] = inputVector[i] - tempValue*inputVector[k];
        //            }
        //        }

        //    }

        //    for (i = count - 1; i > -1; i--)
        //    {
        //        returnAnswerVector[i] = 0;
        //        tempValue = inputVector[i];
        //        for (j = count - 1; j > i; j--)
        //            tempValue = tempValue - inputMatrix[i, j]*returnAnswerVector[j];
        //        returnAnswerVector[i] = tempValue;
        //    }
        //}

        //private void copyMatrix(MatrixEquation matrix)
        //{

        //    int n = matrix.returnVector.Length;

        //    equation = new MatrixEquation(n/3);
        //    for (int i = 0; i <n; i++)
        //    {
        //        for (int j = 0; j < n; j++)
        //        {
        //            equation.returnMatrix[i, j] = matrix.returnMatrix[i, j];
        //        }
        //        equation.returnVector[i] = matrix.returnVector[i];
        //    }

        //}

        //private void show(MatrixEquation matrix, double[] res)
        //{

        //    using (StreamWriter writer = File.CreateText("Soft.txt"))
        //    {
        //        for (int i = 0; i < matrix.returnVector.Length; i++)
        //        {
        //            for (int j = 0; j < matrix.returnVector.Length; j++)
        //            {
        //                writer.Write("[{0},{1}] ->", i, j);
        //                writer.Write(matrix.returnMatrix[i, j]);
        //                writer.Write("                    \t ");

        //            }
        //            writer.Write(matrix.returnVector[i]);
        //            writer.WriteLine("\t");
        //        }
        //        writer.WriteLine();
        //        for (int i = 0; i < res.Length; i++)
        //        {
        //            writer.WriteLine("Result[i]  =  " + res[i]);
        //        }
        //    }

        //}
    }
}
    
