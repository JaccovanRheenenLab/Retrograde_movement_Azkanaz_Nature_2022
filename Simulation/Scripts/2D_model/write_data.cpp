void write_data()
{
	string datname;
	

	
	int filenichtda=0;
	
	datname="surv";
	datname+= ".txt";
	ifstream testles(datname.c_str());
	if (!testles.good()) {
		filenichtda=1;
	}
	
	ofstream data(datname.c_str(),ios::out|ios::app);
	while (!data.good())
	{
		data.clear();
		data.open(datname.c_str(),ios::out|ios::app);
		cerr << "writing to measurements.txt,ios::out|ios::app faild" << endl;
		if (data.good()) cerr << " writing works" << endl;
	}
	

    // we write, for each time point, 1) the clonal survival probability of each starting position, 2) the probability of monoclonality), and 3) the average clone size
    // we write it for each of the 1000 repetition of clonal induction, which will then be averaged in "write_data.cpp".
        for(int i=0;i<T;i++)
        {
                data  << i << " " ;
                for(int p=0;p<4;p++)
                {
                    data  << surv[p][i]/count_pers[p] << " " ;
                }
                data << monoclonal[i]/(surv[0][i]+surv[1][i]+surv[2][i]+surv[3][i])*clones_num << " " << av_size[i][0]/av_size[i][1];
            data << "\n";

        }

        

  data.close();
 
}
