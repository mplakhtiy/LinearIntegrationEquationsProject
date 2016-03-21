using System;
namespace LinearIntegrationEquation.BusinessLogic
{
    class InputEquationData
    {
        public double XsubIJ(int inputI, int inputJ, double inputValue)
        {
            double returnValue = 0;
            switch (inputI)
            {
                case 1:
                    switch (inputJ)
                    {
                        case 1:
                            returnValue = 0.5 * Math.Cos(inputValue); 
                           // returnValue = 1.3 * Math.Cos(inputValue);
                            //returnValue = 0.5 * Math.Cos(inputValue); 
                            break;
                        case 2:
                            returnValue = 0.5 * Math.Sin(inputValue); 
                            //returnValue = Math.Sin(inputValue);
                            //returnValue = 0.4 * Math.Sin(inputValue) - 0.3 * Math.Pow(Math.Sin(inputValue), 2); 
                            break;
                    }
                    break;
                case 2:
                    switch (inputJ)
                    {
                        case 1:
                            returnValue = 1.5 * Math.Cos(inputValue); 
                            //returnValue = 0.5 * Math.Cos(inputValue);
                            //returnValue = 1.3 * Math.Cos(inputValue); 
                            break;
                        case 2:
                            returnValue = 1.5 * Math.Sin(inputValue); 
                            //returnValue = 0.4 * Math.Sin(inputValue) - 0.3 * Math.Pow(Math.Sin(inputValue), 2);
                            //returnValue = Math.Sin(inputValue); 
                            break;

                    }
                    break;

            }
            return returnValue;
        }

        public double FirstDerivativeOfXsubIJ(int inputI, int inputJ, double inputValue)
        {
            double returnValue = 0;
            switch (inputI)
            {
                case 1:
                    switch (inputJ)
                    {
                        case 1:
                            returnValue = -0.5 * Math.Sin(inputValue); 
                            //returnValue = -1.3 * Math.Sin(inputValue);
                            //returnValue = -0.5 * Math.Sin(inputValue); 
                            break;
                        case 2:
                            returnValue = 0.5 * Math.Cos(inputValue); 
                            //returnValue = Math.Cos(inputValue);
                            //returnValue = 0.4 * Math.Cos(inputValue) - 0.3 * Math.Sin(2 * inputValue); 
                            break;
                    }
                    break;
                case 2:
                    switch (inputJ)
                    {
                        case 1:
                            returnValue = -1.5 * Math.Sin(inputValue); 
                            //returnValue = -0.5 * Math.Sin(inputValue);
                            //returnValue = -1.3 * Math.Sin(inputValue); 
                            break;
                        case 2:
                            returnValue = 1.5 * Math.Cos(inputValue); 
                            //returnValue = 0.4 * Math.Cos(inputValue) - 0.3 * Math.Sin(2 * inputValue);
                            //returnValue = Math.Cos(inputValue); 
                            break;
                    }
                    break;

            }
            return returnValue;
        }

        public double SecondDerivativeOfXsubIJ(int inputI, int inputJ, double inputValue)
        {
            double returnValue = 0;
            switch (inputI)
            {
                case 1:
                    switch (inputJ)
                    {
                        case 1:
                            returnValue = -0.5 * Math.Cos(inputValue); 
                            //returnValue = -1.3 * Math.Cos(inputValue);
                            //returnValue = -0.5 * Math.Cos(inputValue); 
                            break;
                        case 2:
                            returnValue = -0.5 * Math.Sin(inputValue); 
                            //returnValue = -Math.Sin(inputValue);
                            //returnValue = -0.4*Math.Sin(inputValue) - 0.6*Math.Cos(2*inputValue); 
                            break;
                    }
                    break;
                case 2:
                    switch (inputJ)
                    {
                        case 1:
                            returnValue = -1.5 * Math.Cos(inputValue); 
                            //returnValue = -0.5 * Math.Cos(inputValue);
                            //returnValue = -1.3 * Math.Cos(inputValue); 
                            break;
                        case 2:
                            returnValue = -1.5 * Math.Sin(inputValue); 
                            //returnValue = -0.4 * Math.Sin(inputValue) - 0.6 * Math.Cos(2 * inputValue);
                            //returnValue = -Math.Sin(inputValue); 
                            break;

                    }
                    break;

            }
            return returnValue;
        }

        public double KapaI(int inputLimit)
        {
            double returnKapa = 0;
            switch (inputLimit)
            {
                case 1:
                    returnKapa = 2;
                    break;
                case 2:
                    returnKapa = 2;
                    break;
            }
            return returnKapa;
        }

        public double function(double inputSsubJ)
        {
            double[] y = { 3, 3 };
            return Math.Pow(1.5 * Math.Cos(inputSsubJ) - 3, 2) + Math.Pow(1.5 * Math.Sin(inputSsubJ) - 3, 2);
            //return (Math.Pow(KapaI(1) * (1.5 * Math.Cos(inputSsubJ) - 3), 3) -
            //   (Math.Pow(KapaI(1) * (1.5 * Math.Sin(inputSsubJ) - 3), 3))) / (2 * Math.PI);
            //return (-
            //Math.Log(KapaI(1) * absOfDifference(inputSsubJ, y)) *
            //commonSeries1(1, inputSsubJ, y) +
            //commonSeries2(1, inputSsubJ, y)) / (2 * Math.PI);
            ////return 7;

        }

        private double absOfDifference(double inputValueS, double[] inputValue)
        {
            return Math.Sqrt(Math.Pow(XsubIJ(1, 1, inputValueS) - inputValue[0], 2) + Math.Pow(XsubIJ(1, 2, inputValueS) - inputValue[1], 2));
        }

        private const double EpsilonConst = 0.0000000000000001;

        private double commonSeries1(int inputLimit1, double inputValueS, double[] inputValue)
        {
            double returnValue = 1.0, checkValue;
            int k = 1;
            double fractionValue =
                Math.Pow((KapaI(inputLimit1) / 2.0) * absOfDifference(inputValueS, inputValue), 2);
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

        private double commonSeries2(int inputLimit1, double inputValueS, double[] inputValue)
        {
            double returnValue = MatrixFormer.Psi(1), tempValue = 1, checkValue;
            int k = 1;
            double fractionValue =
                Math.Pow((KapaI(inputLimit1) / 2) * absOfDifference(inputValueS, inputValue), 2);
            do
            {
                tempValue *= fractionValue / (k * k);
                checkValue = returnValue;
                returnValue += (Math.Log(2) + MatrixFormer.Psi(k + 1)) * tempValue;
                k++;
            } while (Math.Abs(returnValue - checkValue) > EpsilonConst);
            return returnValue;
        }
    }
}
