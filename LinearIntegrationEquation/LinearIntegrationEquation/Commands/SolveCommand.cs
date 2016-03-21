using System;
using System.Windows.Input;
using LinearIntegrationEquation.ViewModels;

namespace LinearIntegrationEquation.Commands
{
    internal class SolveCommand : ICommand
    {
        private readonly MainWindowViewModel viewModel;

        public SolveCommand(MainWindowViewModel viewModel)
        {
            this.viewModel = viewModel;
        }

        public bool CanExecute(object parameter)
        {
            return viewModel.CanSolve;
        }

        public void Execute(object parameter)
        {
            viewModel.SolveAction();
        }

        public event EventHandler CanExecuteChanged
        {
            add { CommandManager.RequerySuggested += value; }
            remove { CommandManager.RequerySuggested -= value; }
        }
    }
}
