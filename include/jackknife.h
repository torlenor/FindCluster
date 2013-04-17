template <class T>
void Jackknife(T *in, T &mean, T &error, int Nmeas) {
	mean=0;
        for(int j=0;j<Nmeas;j++){
                mean+=in[j];
		// cout << "j=" << j << " " << in[j] << endl;
        }   
        mean=mean/(double)Nmeas;
    
    	error=0;
        for(int j=0;j<Nmeas;j++){
                error += pow(mean-in[j],2);
        }   
    
        error=sqrt((double)(Nmeas-1)*error/(double)(Nmeas));
}

