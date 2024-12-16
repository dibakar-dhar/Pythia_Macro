#include<iostream>

#include "Pythia8/Pythia.h"

int main()
{
	int nevents = 1;

	Pythia8::Pythia pythia;

	pythia.readString("Beams:idA = 2212");
	pythia.readString("Beams:idB = 2212");
	pythia.readString("Beams:eCM = 13.e3");
	pythia.readString("SoftQCD:all = on");
	pythia.readString("HardQCD:all = off");

	pythia.init();

	for(int i = 0; i < nevents; i++)
	{
		if(!pythia.next())	continue;		
					
		int entries = pythia.event.size();
		std::cout << i+1 <<"\t \t "<< entries << std::endl;
//		std::cout << "Event size: " << entries << std::endl;
//		std::cout <<"Part. Number"<<"\t Part. ID"<< "\t Part. Mass"<<"\t Part. Eta"<< std::endl;

		for(int j = 0; j < entries; j++)
		{			
			int id = pythia.event[j].id();

			double px = pythia.event[j].px();
			double py = pythia.event[j].py();
			double pz = pythia.event[j].pz();

			double m = pythia.event[j].m();
			double e = pythia.event[j].e();

			double eta = pythia.event[j].eta();
//			double phi = pythia.event[j].phi();
			double pabs = sqrt(pow(px,2) + pow(py,2) + pow(pz,2));

//			double eta = tanh(pz/pabs);
			std::cout<<id<< "\t "<< eta <<" \t  "<< pabs <<" \t  "<< m <<" \t  "<< e <<std::endl;

		}

	}

	return 0;

}

























