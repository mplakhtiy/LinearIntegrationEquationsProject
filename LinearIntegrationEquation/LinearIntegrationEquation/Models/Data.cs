using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using Infragistics.Common.Util;
using LinearIntegrationEquation.BusinessLogic;
using LinearIntegrationEquation.EventArgsExtensions;
using EventManager = LinearIntegrationEquation.Managers.EventManager;
using LinearIntegrationEquation.ViewModels;
namespace LinearIntegrationEquation.Models
{
    public class Data : ObservableCollection<Data.Row>
    {
        public Data()
        {
            EventManager.SolutionFormed += onSolutionFormed;
        }
        
        private void onSolutionFormed(object source, EventArgs eventargs)
        {
            int count = 16;
            double[] S = new double[count];
            for (int i = 0; i < count; i++)
            {
                S[i] = 2 * Math.PI * i / count;
            }
            for (int k = 0; k < 5; k++)
            {
                for (int l = 0; l < count; l++)
                {
                    double[] a = new double[] { (0.1 * k) * Math.Sin(S[l]), (0.1 * k) * Math.Cos(S[l]) };
                    App.Current.Dispatcher.Invoke((Action)delegate { Add(new Row(a[0], a[1], SolutionFormer.U2(a))); });

                    if (a[0] == 0 & a[1] == 0)
                    {
                        l = count - 1;
                    }
                }
            }
            for (int i = 1; i < 10; i++)
            {
                for (int l = 0; l < count; l++)
                {
                    double[] a = { (0.1 * i + 0.5) * Math.Sin(S[l]), (0.1 * i + 0.5) * Math.Cos(S[l]) };
                    App.Current.Dispatcher.Invoke((Action)delegate { Add(new Row(a[0], a[1], SolutionFormer.U1(a))); });

                    if (a[0] == 0 & a[1] == 0)
                    {
                        l = count - 1;
                    }
                }
            }

            //for (int k = 0; k < 5; k++)
            //{
            //    for (int l = 0; l < 16; l++)
            //    {
            //        double[] a = { (0.1 * k) * Math.Sin(MatrixFormer.S[l]), (0.1 * k) * Math.Cos(MatrixFormer.S[l]) };
            //        App.Current.Dispatcher.Invoke((Action)
            //            delegate { Add(new Row(a[0], a[1], 3 * SolutionFormer.exactUI(2, a, new double[] { 3, 3 }))); });
            //        if (a[0] == 0 & a[1] == 0)
            //        {
            //            l = count - 1;
            //        }
            //    }
            //}
            //for (int i = 1; i < 10; i++)
            //{
            //    for (int l = 0; l < 16; l++)
            //    {
            //        double[] a = { (0.1 * i + 0.5) * Math.Sin(MatrixFormer.S[l]), (0.1 * i + 0.5) * Math.Cos(MatrixFormer.S[l]) };
            //        App.Current.Dispatcher.Invoke((Action)
            //            delegate { Add(new Row(a[0], a[1], 3 * SolutionFormer.exactUI(1, a, new double[] { 3, 3 }))); });
            //        if (a[0] == 0 & a[1] == 0)
            //        {
            //            l = count - 1;
            //        }
            //    }
            //}
        }

        public class Row : INotifyPropertyChanged
        {
            private static Random random = new Random();

            public Row(double x, double y, double z)
            {

                this.x = x;
                this.y = y;
                this.z = z;
            }

            public double X
            {
                get { return x; }
                set { x = value; OnPropertyChanged("X"); }
            }
            private double x;

            public double Y
            {
                get { return y; }
                set { y = value; OnPropertyChanged("Y"); }
            }
            private double y;

            public double Z
            {
                get { return z; }
                set { z = value; OnPropertyChanged("Z"); }
            }
            private double z;

            private void OnPropertyChanged(string propertyName)
            {
                if (PropertyChanged != null)
                {
                    PropertyChanged(this, new PropertyChangedEventArgs(propertyName));
                }
            }
            public event PropertyChangedEventHandler PropertyChanged;
        }
    }
}
