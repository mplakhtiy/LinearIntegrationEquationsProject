using System;
using LinearIntegrationEquation.EventArgsExtensions;
using EventManager = LinearIntegrationEquation.Managers.EventManager;

namespace LinearIntegrationEquation.BusinessLogic
{
    /// <summary>
    /// This class fill the matrix.
    /// </summary>
    class MatrixFormer
    {
        private static int count;
        private readonly MatrixEquation returnMatrixEquation;
        public static double[] S { get; private set; }
        private const double EpsilonConst = 0.0000000000000001, GamaConst = 0.57721;

        private readonly InputEquationData equation;

        #region Delegates

        private int countOfThreads;

        private delegate void EquationOne();

        private delegate void EquationTwo();

        private delegate void EquationThree();

        private readonly EquationOne equationOne;
        private readonly EquationTwo equationTwo;
        private readonly EquationThree equationThree;

        private void isMatrixFormed(IAsyncResult asyncResult)
        {
            countOfThreads--;
            if (countOfThreads == 0)
            {
                EventManager.OnMatrixFormed(this, new MatrixEquationEventArgs(returnMatrixEquation));
            }
        }

        #endregion

        public MatrixFormer(int inputCount, InputEquationData inputEquation)
        {
            count = inputCount;
            returnMatrixEquation = new MatrixEquation(count);
            fillS();
            equation = inputEquation;

            equationOne = fillTheMatrixWithEquationOne;
            equationTwo = fillTheMatrixWithEquationTwo;
            equationThree = fillTheMatrixWithEquationThree;

            FillMatrixEquation();
        }

        private void FillMatrixEquation()
        {
            equationOne.BeginInvoke(isMatrixFormed, null);
            countOfThreads++;
            equationTwo.BeginInvoke(isMatrixFormed, null);
            countOfThreads++;
            equationThree.BeginInvoke(isMatrixFormed, null);
            countOfThreads++;
        }

        #region EquationOne

        private void fillTheMatrixWithEquationOne()
        {

            for (int j = 0; j < count; j++) // The cycle fills matrix from 0 -> count ( Our matrix is 3*count^3*count ) 
            {
                for (int k = 0; k < count; k++)
                {
                    returnMatrixEquation.returnMatrix[j, k] = ((R(j, k))*
                                                               H2subII(inputLimitI: 1, inputSsubJ: S[j],
                                                                   inputSsubK: S[k]) +
                                                               H1subII(inputLimitI: 1, inputSsubJ: S[j],
                                                                   inputSsubK: S[k])/count);
                }

                for (int k = count; k < 2*count; k++)
                {
                    returnMatrixEquation.returnMatrix[j, k] =  (HsubIJ(inputLimitI: 1, inputSsubJ: S[j], inputLimitJ: 2,
                                                                  inputSsubK: S[k - count]))/count;
                }

                for (int k = 2*count; k < 3*count; k++)
                {
                    returnMatrixEquation.returnMatrix[j, k] = 0;
                }

                returnMatrixEquation.returnVector[j] = equation.function(S[j]);
            }
        }

        #endregion 

        #region EquationTwo

        private void fillTheMatrixWithEquationTwo()
        {
            for (int j = count,jj = 0; j < 2 * count; j++,jj++) // The cycle fills matrix from count -> 2*count ( Our matrix is 3*count^3*count )
            {

                for (int k = 0; k < count; k++)
                {
                    returnMatrixEquation.returnMatrix[j, k] =HsubIJ(inputLimitI: 2, inputSsubJ: S[jj], inputLimitJ: 1,
                                                                  inputSsubK: S[k])/count;
                }
                for (int k = count; k < 2*count; k++)
                {
                    returnMatrixEquation.returnMatrix[j, k] = (H1subII(inputLimitI: 2, inputSsubJ: S[jj], inputSsubK: S[k-count])/count +
                                                           R(jj, k-count) *
                                                           H2subII(inputLimitI: 2, inputSsubJ: S[jj], inputSsubK: S[k-count]));
                }

                for (int k = 2*count; k < 3*count; k++)
                {
                    returnMatrixEquation.returnMatrix[j, k] = -(H1subII(inputLimitI: 2, inputSsubJ: S[jj], inputSsubK: S[k-2*count]) / count +
                                                           R(jj, k-2*count) *
                                                           H2subII(inputLimitI: 2, inputSsubJ: S[jj], inputSsubK: S[k-2*count]));
                }

                returnMatrixEquation.returnVector[j] = 0;
                
            }
        }

        #endregion

        #region EquationThree

        private void fillTheMatrixWithEquationThree()
        {
            for (int j = 2*count, jj = 0; j < 3*count; j++,jj++)
                // The cycle fills matrix from 2*count -> 3*count ( Our matrix is 3*count^3*count )
            {
                for (int k = 0; k < count; k++)
                {
                    returnMatrixEquation.returnMatrix[j,k] = (equation.KapaI(1)/count)*
                                                              LsubIJ(inputLimitI: 2, inputSsubJ: S[jj], inputLimitJ: 1,
                                                                  inputSsubK: S[k]);
                }
                for (int k = count; k < 2*count; k++)
                {
                    if (j != k-count)
                    {
                        returnMatrixEquation.returnMatrix[j, k] = (equation.KapaI(1))*
                                                                      (L1subII(inputLimitI: 2, inputSsubJ: S[jj],
                                                                          inputSsubK: S[k-count])/count +
                                                                       R(jj, k-count)*
                                                                       L2subII(inputLimitI: 2, inputSsubJ: S[jj],
                                                                           inputSsubK: S[k-count]));

                    }
                    else
                    {
                        returnMatrixEquation.returnMatrix[j, k] = (equation.KapaI(1))*
                                                                      (L1subII(inputLimitI: 2, inputSsubJ: S[jj],
                                                                          inputSsubK: S[k-count])/count +
                                                                       R(jj, k-count)*
                                                                       L2subII(inputLimitI: 2, inputSsubJ: S[jj],
                                                                           inputSsubK: S[k-count]) + 0.5);
                    }
                }

                for (int k = 2*count; k < 3*count; k++)
                {
                    if (j != k-2*count)
                    {
                        returnMatrixEquation.returnMatrix[j, k] = -(equation.KapaI(2))*
                                                                      (L1subII(inputLimitI: 2, inputSsubJ: S[jj],
                                                                          inputSsubK: S[k-2*count])/count +
                                                                       R(jj, k - 2 * count) *
                                                                       L2subII(inputLimitI: 2, inputSsubJ: S[jj],
                                                                           inputSsubK: S[k - 2 * count]));

                    }
                    else
                    {
                        returnMatrixEquation.returnMatrix[j, k] = -(equation.KapaI(2))*
                                                                      (L1subII(inputLimitI: 2, inputSsubJ: S[jj],
                                                                          inputSsubK: S[k - 2 * count]) / count +
                                                                       R(jj, k - 2 * count) *
                                                                       L2subII(inputLimitI: 2, inputSsubJ: S[jj],
                                                                           inputSsubK: S[k - 2 * count]) + 0.5);
                    }
                }

                returnMatrixEquation.returnVector[j] = 0;
            }
        }

        #endregion

        private double H2subII(int inputLimitI, double inputSsubJ, double inputSsubK)
        {
            const int numberOfFunction1 = 1, numberOfFunction2 = 2;
            if (Math.Abs(inputSsubJ - inputSsubK) < EpsilonConst)
            {
                return -0.5 * lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ),
                    equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ));
            }
            return -0.5*commonSeries1(inputLimitI, inputSsubJ, inputLimitI, inputSsubK) *
                   lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubK),
                       equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubK));
        }

        private double H1subII(int inputLimitI, double inputSsubJ, double inputSsubK)
        {
            const int numberOfFunction1 = 1, numberOfFunction2 = 2;
            if (Math.Abs(inputSsubK - inputSsubJ) < EpsilonConst)
            {
                return ((-(0.5) *
                       Math.Log(Math.E * Math.Pow(equation.KapaI(inputLimitI) *
                                    lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ),
                                        equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ)), 2))) + Math.Log(2) + Psi(1)) *
                                        lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ),
                    equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ));
            }
            return (-(((0.5) *
                       Math.Log(Math.E * Math.Pow(equation.KapaI(inputLimitI) *
                                    lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ),
                                        equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ)), 2))) *
                       commonSeries1(inputLimitI, inputSsubJ, inputLimitI, inputSsubK))+
                     commonSeries2(inputLimitI, inputSsubJ, inputLimitI, inputSsubK)) *
                    lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubK),
                    equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubK));
        }

        private double HsubIJ(int inputLimitI, double inputSsubJ, int inputLimitJ, double inputSsubK)
        {
            const int numberOfFunction1 = 1, numberOfFunction2 = 2;
            return
                (-Math.Log(equation.KapaI(inputLimitI) * absOfDifference(inputLimitI, inputSsubJ, inputLimitJ, inputSsubK)) *
                commonSeries1(inputLimitI, inputSsubJ, inputLimitJ, inputSsubK) +
                 commonSeries2(inputLimitI, inputSsubJ, inputLimitJ, inputSsubK)) *
                lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitJ, numberOfFunction1, inputSsubK),
                    equation.FirstDerivativeOfXsubIJ(inputLimitJ, numberOfFunction2, inputSsubK));
        }

        private double LsubIJ(int inputLimitI, double inputSsubJ, int inputLimitJ, double inputSsubK)
        {
            const int numberOfFunction1 = 1, numberOfFunction2 = 2;

            return (1.0/(equation.KapaI(inputLimitI)*absOfDifference(inputLimitI, inputSsubJ, inputLimitJ, inputSsubK)) +
                    (Math.Log(equation.KapaI(inputLimitI)*
                              absOfDifference(inputLimitI, inputSsubJ, inputLimitJ, inputSsubK))*
                     commonSeries4(inputLimitI, inputSsubJ, inputLimitJ, inputSsubK)) -
                    commonSeries3(inputLimitI, inputSsubJ, inputLimitJ, inputSsubK))*
                   ((-equation.KapaI(inputLimitI))*
                    ((equation.XsubIJ(inputLimitI, numberOfFunction2, inputSsubJ) -
                      equation.XsubIJ(inputLimitJ, numberOfFunction2, inputSsubK))*
                     equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ) -
                     (equation.XsubIJ(inputLimitI, numberOfFunction1, inputSsubJ) -
                      equation.XsubIJ(inputLimitJ, numberOfFunction1, inputSsubK))*
                     equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ)))*
                   lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitJ, numberOfFunction1, inputSsubK),
                       equation.FirstDerivativeOfXsubIJ(inputLimitJ, numberOfFunction2, inputSsubK))/
                   ((absOfDifference(inputLimitI, inputSsubJ, inputLimitJ, inputSsubK)*
                     lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ),
                         equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ))));
        }

        private double L1subII(int inputLimitI, double inputSsubJ, double inputSsubK)
        {
            const int numberOfFunction1 = 1, numberOfFunction2 = 2;

            if (Math.Abs(inputSsubK - inputSsubJ) < EpsilonConst)
            {
                return (-(equation.SecondDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ) *
                        equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ) -
                        equation.SecondDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ) *
                        equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ)) / (2 *
                        Math.Pow(lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ),
                           equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ)), 2)));
            }

            return (1.0 / (equation.KapaI(inputLimitI) * absOfDifference(inputLimitI, inputSsubJ, inputLimitI, inputSsubK)) -
                    0.5 * (Math.Log(Math.E * Math.Pow(equation.KapaI(inputLimitI)*lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ),
                           equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ)), 2))*
                     commonSeries4(inputLimitI, inputSsubJ, inputLimitI, inputSsubK)) -
                    commonSeries3(inputLimitI, inputSsubJ, inputLimitI, inputSsubK)) *
                   ((-equation.KapaI(inputLimitI)) *
                    ((equation.XsubIJ(inputLimitI, numberOfFunction2, inputSsubJ) -
                      equation.XsubIJ(inputLimitI, numberOfFunction2, inputSsubK)) *
                     equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ) -
                     (equation.XsubIJ(inputLimitI, numberOfFunction1, inputSsubJ) -
                      equation.XsubIJ(inputLimitI, numberOfFunction1, inputSsubK)) *
                     equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ))) *
                   lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubK),
                       equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubK)) /
                   ((absOfDifference(inputLimitI, inputSsubJ, inputLimitI, inputSsubK) *
                     lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ),
                         equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ))));
        }

        private double L2subII(int inputLimitI, double inputSsubJ, double inputSsubK)
        {
            const int numberOfFunction1 = 1, numberOfFunction2 = 2;

            if (Math.Abs(inputSsubK - inputSsubJ) < EpsilonConst)
            {
                return 0;
            }

            return 0.5*commonSeries4(inputLimitI, inputSsubJ, inputLimitI, inputSsubK)*
                   ((-equation.KapaI(inputLimitI))*
                    ((equation.XsubIJ(inputLimitI, numberOfFunction2, inputSsubJ) -
                      equation.XsubIJ(inputLimitI, numberOfFunction2, inputSsubK))*
                     equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ) -
                     (equation.XsubIJ(inputLimitI, numberOfFunction1, inputSsubJ) -
                      equation.XsubIJ(inputLimitI, numberOfFunction1, inputSsubK))*
                     equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ)))*
                   lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubK),
                       equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubK))/
                   ((absOfDifference(inputLimitI, inputSsubJ, inputLimitI, inputSsubK)*
                     lengthOfVector(equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction1, inputSsubJ),
                         equation.FirstDerivativeOfXsubIJ(inputLimitI, numberOfFunction2, inputSsubJ))));
        }

        private double commonSeries1(int inputLimit1, double inputValue1, int inputLimit2, double inputValue2)
        {
            double returnValue = 1.0, checkValue;
            int k = 1;
            double fractionValue =
                Math.Pow((equation.KapaI(inputLimit1) / 2.0) * absOfDifference(inputLimit1, inputValue1, inputLimit2, inputValue2), 2);
            double tempValue = 1;
            do
            {
                tempValue *= fractionValue / (k*k);
                checkValue = returnValue;
                returnValue += tempValue ;
                k++;
                
            } while (Math.Abs(returnValue - checkValue) > EpsilonConst);

            return returnValue;
        }

        private double commonSeries2(int inputLimit1, double inputValue1, int inputLimit2, double inputValue2)
        {
            double returnValue = Psi(1),tempValue=1 , checkValue;
            int k = 1;
            double fractionValue =
                Math.Pow((equation.KapaI(inputLimit1) / 2) * absOfDifference(inputLimit1, inputValue1, inputLimit2, inputValue2), 2);
            do
            {
                tempValue *= fractionValue/(k*k);
                checkValue = returnValue;
                returnValue += (Math.Log(2) + Psi(k + 1)) * tempValue;
                k++;
            } while (Math.Abs(returnValue - checkValue) > EpsilonConst);
            return returnValue;
        }

        private double commonSeries3(int inputLimit1, double inputValue1, int inputLimit2, double inputValue2)
        {
            double checkValue;
            int k = 1;
            double fractionValue =
                ((equation.KapaI(inputLimit1) / 2) * absOfDifference(inputLimit1, inputValue1, inputLimit2, inputValue2));
            double tempValue = fractionValue;
            double returnValue = (Psi(1) + Psi(2)) * fractionValue;
            double psiTempValue1 = Psi(k + 1);
            double psiTempValue2;
            do
            {
                checkValue = returnValue;
                psiTempValue2 = Psi(k + 2);
                tempValue *= Math.Pow(fractionValue, 2) / ((k + 1) * k);
                returnValue += (Math.Log(2) + 0.5 * (psiTempValue1 + psiTempValue2)) * tempValue;
                psiTempValue1 = psiTempValue2;
                k++;

            } while (Math.Abs(returnValue - checkValue) > EpsilonConst);

            return returnValue;
        }

        private double commonSeries4(int inputLimit1, double inputValue1, int inputLimit2, double inputValue2)
        {
            double  checkValue;
            int k = 1;
            double fractionValue =(equation.KapaI(inputLimit1) / 2) * absOfDifference(inputLimit1, inputValue1, inputLimit2, inputValue2);
            double returnValue = fractionValue;
            double tempValue = fractionValue;
            do
            {
                checkValue = returnValue;
                tempValue *= Math.Pow(fractionValue, 2) / (k*(k+1));
                returnValue += tempValue;
                k++;
   
            } while (Math.Abs(returnValue - checkValue) > EpsilonConst);

            return returnValue;
        }

        public static double Psi(int inputValue)
        {
            double returnValue = -GamaConst;

            for (int k = 1; k < inputValue; k++)
            {
                returnValue += Math.Pow(k, -1);
            }
            return returnValue;
      
        }

        private double R(int inputJ, int inputK)
        {
            double returnValue = 0;
            for (int m = 1; m < count / 2.0; m++)
            {
                returnValue += (1.0 / m) * Math.Cos(m * (S[inputJ] - S[inputK]));
            }
            return ((-1.0 /count) * (1 + 2 * returnValue + (2.0/count)*Math.Cos((count/2.0)*(S[inputJ] - S[inputK])))) ;
        }

        private double absOfDifference(int inputLimit1, double inputValue1, int inputLimit2, double inputValue2)
        {
            const int numberOfFunction1 = 1, numberOfFunction2 = 2;
            return lengthOfVector(equation.XsubIJ(inputLimit1, numberOfFunction1, inputValue1) - equation.XsubIJ(inputLimit2, numberOfFunction1, inputValue2),
                                      equation.XsubIJ(inputLimit1, numberOfFunction2, inputValue1) - equation.XsubIJ(inputLimit2, numberOfFunction2, inputValue2));
        }

        private double lengthOfVector(double inputValue1, double inputValue2)
        {
            return Math.Sqrt(Math.Pow(inputValue1, 2) + Math.Pow(inputValue2, 2));
        }

        private void fillS()
        {
            S = new double[count];
            for (int i = 0; i < count; i++)
            {
                S[i] = 2 * Math.PI * i / count;
            }
        }

        ////////////////////////////////////////////////// Патерн головного мозгу виглядає так:
        //SClass s;

        //private class SClass
        //{
        //    private double[] sArray;
        //    private SClass sInstance;
        //    private SClass(int count)
        //    {
        //        sArray = new double[count];
        //        for (int i = 0; i < count; i++)
        //        {
        //            sArray[i] = 2 * Math.PI * i / count;
        //        }
        //    }

        //    public double S(int inputI)
        //    {
        //        if (sInstance != null)
        //        {
        //            return sArray[inputI];
        //        }
        //        sInstance = new SClass(count);
        //        return sArray[inputI];
        //    }
        //}
        /////////////////////////////////////////////////////////////////////////////////////////////
    }
}
