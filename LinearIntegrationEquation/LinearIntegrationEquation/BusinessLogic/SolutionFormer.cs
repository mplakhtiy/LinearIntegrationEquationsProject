using System;
using System.IO;
using EventManager = LinearIntegrationEquation.Managers.EventManager;
using LinearIntegrationEquation.EventArgsExtensions;
namespace LinearIntegrationEquation.BusinessLogic
{
    class SolutionFormer
    {
        public SolutionFormer(InputEquationData inputEquation)
        {
            equation = inputEquation;
            EventManager.MatrixSolved += onMatrixSolved;
        }

        private static InputEquationData equation;
        private static double[] solutions;

        private void onMatrixSolved(object source, SolutionOfMatrixEventArgs eventargs)
        {
            solutions = eventargs.SolutionVector;   //showSolutions();
            EventManager.OnSolutionFormed(this,new EventArgs());
        }

        public static double U1(double[] inputValue)
        {
            double resultValue = 0;
            for (int i = 0; i < solutions.Length / 3.0; i++)
            {
                resultValue += solutions[i] * HsubIJ(1, inputValue, 1, MatrixFormer.S[i]) +
                               solutions[i + solutions.Length / 3] * HsubIJ(1, inputValue, 2, MatrixFormer.S[i]);
            }
            return (3.0 / solutions.Length) * resultValue;
        }

        public static double U2(double[] inputValue)
        {
            double resultValue = 0;
            for (int i = 0; i < solutions.Length / 3.0; i++)
            {
                resultValue += solutions[i + 2 * (solutions.Length / 3)] * HsubIJ(2, inputValue, 2, MatrixFormer.S[i]);
            }
            return (3.0 / solutions.Length) * resultValue;
        }

        private const double EpsilonConst = 0.0000000000000001;

        private static double HsubIJ(int inputLimitI, double[] inputValue, int inputLimitJ, double inputSsubK)
        {
            const int numberOfFunction1 = 1, numberOfFunction2 = 2;
            return
                (-Math.Log(equation.KapaI(inputLimitI) * absOfDifference(inputValue, inputLimitJ, inputSsubK)) *
                commonSeries1(inputLimitI, inputValue, inputLimitJ, inputSsubK) +
                 commonSeries2(inputLimitI, inputValue, inputLimitJ, inputSsubK)) *
                lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitJ, numberOfFunction1, inputSsubK),
                    equation.FirstDerivativeOfXsubIJ(inputLimitJ, numberOfFunction2, inputSsubK));
        }

        private static double lengthOfVector(double inputValue1, double inputValue2)
        {
            return Math.Sqrt(Math.Pow(inputValue1, 2) + Math.Pow(inputValue2, 2));
        }

        private static double absOfDifference(double[] inputValue1, int inputLimit, double inputValue2)
        {
            return
                Math.Sqrt(Math.Pow(inputValue1[0] - equation.XsubIJ(inputLimit, 1, inputValue2), 2) +
                          Math.Pow(inputValue1[1] - equation.XsubIJ(inputLimit, 2, inputValue2), 2));
        }

        private static double commonSeries1(int inputLimit1, double[] inputValue, int inputLimit2, double inputValue2)
        {
            double returnValue = 1.0, checkValue;
            int k = 1;
            double fractionValue =
                Math.Pow((equation.KapaI(inputLimit1) / 2.0) * absOfDifference(inputValue, inputLimit2, inputValue2), 2);
            double tempValue = 1;
            do
            {
                tempValue *= fractionValue / (k * k);
                checkValue = returnValue;
                returnValue += tempValue;
                k++;

            } while (Math.Abs(returnValue - checkValue) > EpsilonConst);

            return returnValue;
        }

        private static double commonSeries2(int inputLimit1, double[] inputValue, int inputLimit2, double inputValue2)
        {
            double returnValue = MatrixFormer.Psi(1), tempValue = 1, checkValue;
            int k = 1;
            double fractionValue =
                Math.Pow((equation.KapaI(inputLimit1) / 2) * absOfDifference(inputValue, inputLimit2, inputValue2), 2);
            do
            {
                tempValue *= fractionValue / (k * k);
                checkValue = returnValue;
                returnValue += (Math.Log(2) + MatrixFormer.Psi(k + 1)) * tempValue;
                k++;
            } while (Math.Abs(returnValue - checkValue) > EpsilonConst);
            return returnValue;
        }

        public static double exactUI(int inputLimit, double[] inputValue1, double[] inputValue2)
        {
            //return (-Math.Log(equation.KapaI(inputLimit) * absOfDifference(inputValue1, inputValue2)) *
            //        commonSeries1(inputLimit, inputValue1, inputValue2) +
            //        commonSeries2(inputLimit, inputValue1, inputValue2)) / (2 * Math.PI);
            //return (0.5 - (Math.Log(0.5 * equation.KapaI(inputLimit) * absOfDifference(inputValue1, inputValue2))) *
            //         commonSeries1(inputLimit, inputValue1, inputValue2) + 2 * commonSeries2(inputLimit, inputValue1, inputValue2)) /
            //        (2 * Math.PI);
            //return( (Math.Pow(equation.KapaI(inputLimit)*Math.Abs(inputValue1[0] - inputValue2[0]), 3) -
            //        Math.Pow(equation.KapaI(inputLimit)*Math.Abs(inputValue1[1] - inputValue2[1]),3))/(2*Math.PI));
            return Math.Pow(equation.KapaI(inputLimit) * Math.Abs(inputValue1[0] - inputValue2[0]), 2) +
                   Math.Pow(equation.KapaI(inputLimit) * Math.Abs(inputValue1[1] - inputValue2[1]), 2);
        }
       
        private static double absOfDifference(double[] inputValue1, double[] inputValue2)
        {
            return Math.Sqrt(Math.Pow(inputValue1[0] - inputValue2[0], 2) + Math.Pow(inputValue1[1] - inputValue2[1], 2));
        }

        private static double commonSeries1(int inputLimit, double[] inputValue1, double[] inputValue2)
        {
            double returnValue = 1.0, checkValue;
            int k = 1;
            double fractionValue =
                Math.Pow((equation.KapaI(inputLimit) / 2.0) * absOfDifference(inputValue1, inputValue2), 2);
            double tempValue = 1;
            do
            {
                tempValue *= fractionValue / (k * k);
                checkValue = returnValue;
                returnValue += tempValue;
                k++;

            } while (Math.Abs(returnValue - checkValue) > EpsilonConst);

            return returnValue;
        }

        private  static double commonSeries2(int inputLimit, double[] inputValue1, double[] inputValue2)
        {
            double returnValue = MatrixFormer.Psi(1), tempValue = 1, checkValue;
            int k = 1;
            double fractionValue =
                Math.Pow((equation.KapaI(inputLimit) / 2) * absOfDifference(inputValue1, inputValue2), 2);
            do
            {
                tempValue *= fractionValue / (k * k);
                checkValue = returnValue;
                returnValue += (Math.Log(2) + MatrixFormer.Psi(k + 1)) * tempValue;
                k++;
            } while (Math.Abs(returnValue - checkValue) > EpsilonConst);
            return returnValue;
        }

        private void showSolutions()
        {
            using (StreamWriter writer = File.CreateText("Soft.txt"))
            {
                double[] x1 = { 1, 0.7 };
                double[] x2 = { 0, 0.1};
                double[] y = { 3, 3};
                writer.WriteLine(U1(x1));
                writer.WriteLine(U2(x2));

                writer.WriteLine(exactUI(1, x1, y));
                writer.WriteLine(exactUI(2, x2, y));
            }
        }
    }
}
