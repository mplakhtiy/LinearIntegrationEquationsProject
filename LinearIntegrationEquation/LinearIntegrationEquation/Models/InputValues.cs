using System.ComponentModel;
using System.Runtime.CompilerServices;
using LinearIntegrationEquation.Annotations;

namespace LinearIntegrationEquation.Models
{
    class InputValues:INotifyPropertyChanged
    {
        private string count;

        public string Count
        {
            get { return count; }
            set
            {
                count = value;
                OnPropertyChanged();
            }
        }

        public InputValues(string inputValueCount)
        {
            count = inputValueCount;
        }
        public event PropertyChangedEventHandler PropertyChanged;

        [NotifyPropertyChangedInvocator]
        protected virtual void OnPropertyChanged([CallerMemberName] string propertyName = null)
        {
            var handler = PropertyChanged;
            if (handler != null) handler(this, new PropertyChangedEventArgs(propertyName));
        }
    }
}
