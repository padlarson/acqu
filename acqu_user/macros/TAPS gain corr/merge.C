#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>      // std::setprecision
#include <vector>
#include <stdlib.h>
#include <map>

int main()
{

//g++ merge.C -o merge --std=c++11


	const unsigned int start = 4952;
	const unsigned int end = 5906;
	double gain;

	std::vector<double> f1;
	std::vector<double> f2;

    typedef std::pair<int, std::vector<double>> corr_pair_1;
    typedef std::pair<int, std::vector<double>> corr_pair_2;
    std::map<int, std::vector<double>> TAPS_corr_1;
    std::map<int, std::vector<double>> TAPS_corr_2;


    // read in file0
    std::ifstream file1("/home/adlarson/data2014.07/macros/TAPSgains/taps_lg_e1.txt");
    std::string   line;
    while(std::getline(file1,line))
    {
        std::string         buffer;
        std::stringstream   ss;
        int  runnr;
        std::vector<double> vec;
        ss << line;
        std::getline(ss,buffer, '\t');
        runnr = atoi(buffer.c_str());
        while(std::getline(ss, buffer, '\t')){
            vec.push_back((std::stod(buffer)));
        }
        if(vec.size() != 438 )
        	std::cout << "For run " << runnr << " the vector length is " << vec.size() << std::endl;
        TAPS_corr_1.insert(corr_pair_1(runnr,vec));
    }
    file1.close();

    std::ifstream file2("/home/adlarson/data2014.07/macros/TAPSgains/TAPS_gain_Jun16_2015.txt");
    while(std::getline(file2,line))
    {
        std::string         buffer;
        std::stringstream   ss;
        int  runnr;
        std::vector<double> vec;
        ss << line;
        std::getline(ss,buffer, '\t');
        runnr = atoi(buffer.c_str());
        while(std::getline(ss, buffer, '\t')){
            vec.push_back((std::stod(buffer)));
        }
        TAPS_corr_2.insert(corr_pair_2(runnr,vec));
        if(vec.size() != 438 )
        	std::cout << "For run " << runnr << " the vector length is " << vec.size() << std::endl;
    }
    file2.close();


	std::ofstream out;
	out.open("/home/adlarson/data2014.07/macros/TAPSgains/taps_lg_e1_a.txt");    	

    for(auto irun = start; irun <= end; irun++ )
    {
    	out << irun << '\t';
    	if( (TAPS_corr_1.find(irun) == TAPS_corr_1.end()) && ( TAPS_corr_2.find(irun) == TAPS_corr_2.end()) ) // no gain corrections found in any of the files
    	{

    		 for(int j = 0; j < 438; j++)
    		{
    			out << std::setprecision(5) << 1 << '\t';
    		}
    		out << '\n';
		}
    	else if( (TAPS_corr_1.find(irun) != TAPS_corr_1.end()) && (TAPS_corr_2.find(irun) == TAPS_corr_2.end()) ) // no gain corrections found in file 1
    	{
    		std::cout << "irun "<< irun << " not found in TAPS_corr_2 !  " << std::endl;
    		f1= TAPS_corr_1.find(irun)->second;
    		
    		for(int j = 0; j < 438; j++)
    		{
    			out << std::setprecision(5) << f1[j] << '\t';
    		}
    		out << '\n';
    	}
    	else if( (TAPS_corr_2.find(irun) != TAPS_corr_2.end()) && (TAPS_corr_1.find(irun) == TAPS_corr_1.end()) ) // no gain corrections found in file 2
    	{
    		std::cout << "irun "<< irun << " not found in TAPS_corr_1 ! " << std::endl;
    		f2 = TAPS_corr_2.find(irun)->second;
    		
    		for(int j = 0; j < 438; j++)
    		{
    			out << std::setprecision(5) << f2[j] << '\t';
    		}
    		out << '\n';
    	}
    	else   // both gain corrections found 
    	{
    		f1 = TAPS_corr_1.find(irun)->second;
    		f2 = TAPS_corr_2.find(irun)->second;
    		
    		for(int j = 0; j < 438; j++)
    		{
    			out << std::setprecision(5) << f1[j]*f2[j] << '\t';

    			// out << std::setprecision(5) << std::stod( f1[j] )*std::stod( f2[j] ) << '\t';
    		}
    		out << '\n';
    	}

    }

out.close();


	return 0;

}