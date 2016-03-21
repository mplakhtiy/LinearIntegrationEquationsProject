using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Input;
using Infragistics.Controls.Charts;
using LinearIntegrationEquation.BusinessLogic;
using LinearIntegrationEquation.Commands;
using LinearIntegrationEquation.Models;

namespace LinearIntegrationEquation.ViewModels
{
    class MainWindowViewModel:INotifyPropertyChanged
    {
        public ICommand Solve { get; private set; }

        public bool CanSolve
        {
            get
            {
                if (InputValues == null)
                {
                    return false;
                }
                return !String.IsNullOrWhiteSpace(InputValues.Count);
            }
        }

        public InputValues InputValues { get; private set; }

        public MainWindowViewModel()
        {
            Solve=new SolveCommand(this);
            InputValues=new InputValues("Write Count Here(16,32,64,128,..)");
            Data=new Data();
        }

        public void SolveAction()
        {
            InputEquationData inputEquation = new InputEquationData();
            new MatrixFormer(int.Parse(InputValues.Count),inputEquation);
            new MatrixSolver();
            new SolutionFormer(inputEquation);
            // App.Current.Shutdown();
        }

        public Data Data
        {
            get { return data; }
            private set
            {
                data = value;
                OnPropertyChanged(new PropertyChangedEventArgs("Data"));
            }
        }
        private Data data;

        public event PropertyChangedEventHandler PropertyChanged;

        private void OnPropertyChanged(PropertyChangedEventArgs ea)
        {
            if (PropertyChanged != null)
            {
                PropertyChanged(this, ea);
            }
        }

    }
}
