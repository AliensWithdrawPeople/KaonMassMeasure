#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <math.h>

double getrunlum_v(int run, int tr_ph_v)
{
	double buf=0;
	double lum=0;
	char pointlist[299]="";
	int sid=0;
	
	if(run>=8997 && run<=13516) sid=1;
	if(run>=14522 && run<=16597) sid=2;
	if(run>=17405 && run<=18808) sid=3;
	if(run>=18809 && run<=24395) sid=4;
	if(run>=24396 && run<=32075) sid=5;
	if(run>=36905 && run<=47683) sid=6;
	if(run>=48938 && run<=55741) sid=7;
	if(run>=55742 && run<=69268) sid=8;
	if(run>=70014 && run<=84943) sid=9;
	if(run>=85224 && run<=89844) sid=10;
	if(run>=89973 && run<=96669) sid=11;
	if(run>=98116 && run<=106704) sid=12;
	if(run>=107490 && run<=131392) sid=13;
	if(run>=131913 && run<=150158) sid=14;
	
	int version=0;
	//version=tr_ph_v;
	if(tr_ph_v==8) version=91;
	if(tr_ph_v==9) version=931;
	
	sprintf(pointlist,"/storeA/ryzhenenkov/database/luminosity/run%d/runs/v%d_ee_run_%d_tcc1.0.dat",sid,version,run);
	//sprintf(pointlist,"/store20/ryzhenenkov/tree/lum/v931_ee/runs/%d/v%d_ee_run_%d_tcc1.0.dat",sid,version,run);
	//sprintf(pointlist,"/outstage/ryzhenenkov/40107/v%d_ee_run_%d_tcc1.0.dat",version,run);
	//if(sid==7 && run>55700)std::cout<<pointlist<<endl;
	std::ifstream csf (pointlist);
	if ( csf.is_open() ){
		for (int i=0; i<=15; i++){csf>>buf;if (i==9) lum=buf;}
		csf.close();
	}
	return lum;
}

double getrunlumerr_v(int run, int tr_ph_v)
{
	double buf=0;
	double err=0;
	char pointlist[299]="";
	int sid=0;	
	if(run>=8997 && run<=13516) sid=1;
	if(run>=14522 && run<=16597) sid=2;
	if(run>=17405 && run<=18808) sid=3;
	if(run>=18809 && run<=24395) sid=4;
	if(run>=24396 && run<=32075) sid=5;
	if(run>=36905 && run<=47683) sid=6;
	if(run>=48938 && run<=55741) sid=7;
	if(run>=55742 && run<=69268) sid=8;
	if(run>=70014 && run<=84943) sid=9;
	if(run>=85224 && run<=89844) sid=10;
	if(run>=89973 && run<=96669) sid=11;
	if(run>=98116 && run<=106704) sid=12;
	if(run>=107490 && run<=131392) sid=13;
	if(run>=131913 && run<=150158) sid=14;
	
	int version=0;
	//version=tr_ph_v;
	if(tr_ph_v==8) version=91;
	if(tr_ph_v==9) version=931;

    sprintf(pointlist,"/storeA/ryzhenenkov/database/luminosity/run%d/runs/v%d_ee_run_%d_tcc1.0.dat",sid,version,run);
	std::ifstream csf (pointlist);
	if ( csf.is_open() ){
		for (int i=0; i<=15; i++){csf>>buf;if (i==12) err=buf;}
		csf.close();
	}
	return err;
}

int main()
{
    std::vector<double> points = {501, 503, 505, 508, 508.5, 509, 509.5, 510, 510.5, 511, 511.5, 514, 517, 520, 525, 530};
	std::vector<double> lumiTot = {};
    std::vector<double> lumiErrTot = {};    
	for(auto& point : points)
    {
        auto fname = "/spoolA/divanov/phi2018_runLumi/runs/runs_" + std::to_string(point) + ".txt";
        std::ifstream input(fname);
        double lumi = -1;
        double lumiErr = -1;
        int run = -1;
		double totLumi = 0;
		double totLumiErr = 0;
		std::ofstream output("/spoolA/divanov/phi2018_runLumi/lumi/lumi_" + std::to_string(point).substr(0, 5) + ".txt",
								std::ofstream::out | std::ofstream::trunc);

        std::cout << "Is input file open? " << input.is_open() << std::endl;
        while (input >> run)
        { 
            std::cout << run << std::endl;
            lumi = getrunlum_v(run, 9);
            lumiErr = getrunlumerr_v(run, 9);
			output << lumi << " " << lumiErr << std::endl; 
			totLumi += lumi;
			totLumiErr += lumiErr * lumiErr;
        }
        input.close();
		output.close();

		lumiTot.push_back(totLumi);
		lumiErrTot.push_back(sqrt(totLumiErr));
    }

	std::cout << "Lumi in energy points = {";
	for(auto& val : lumiTot)
	{ std::cout << val << ", "; }

	std::cout << "}\nLumi error in energy points = {";
	for(auto& val : lumiErrTot)
	{ std::cout << val << ", "; }
	std::cout << "}" << std::endl;
    return 0;
}