namespace LinearIntegrationEquation.BusinessLogic
{
    /// <summary>
    /// This class represent the equation system.
    /// </summary>
    class MatrixEquation
    {
        public double[,] returnMatrix { get; set; }
        public double[] returnVector { get; set; }
        public MatrixEquation(int count)
        {
            returnMatrix = new double[3 * count, 3 * count];
            returnVector = new double[3 * count];
        }
    }
}
