#include "TSpectrum.h"


float targetX[20000];
float targetY[20000];
float sampleX[20000];
float sampleY[20000];
float EX[20000];
float EY[20000];
float EnergyE[20000];


TGraphErrors *Rh,*Cu;
TCanvas *MyC;
TGraph *error,*errorE;



int Loader(char *fileName,float gamma,float time_step,bool newCanvas=1){
gamma=gamma/2.;
printf("%g %g\n",gamma,time_step);
gStyle->SetOptStat(0000);

ifstream f;
f.open(fileName);
if (!f.is_open()) {
	printf("FILE DOESN'T EXISTS %s!!!\nEXITING\n",fileName);
    	return(404);
	}

int lineN = 0;
char line_c_str[200];
int headerLine=2;

float ADC_Q=0;//zebrany ladunek
while (!f.eof()){


    
    stringstream line;//Musi byc tu zeby byl realokowany zawsze
    f.getline(line_c_str, 200);
    line.str(string(line_c_str));

    if(headerLine>0){headerLine--;printf("%d\n",lineN);continue;}
    //truncate white characters
    line >> ws;
    //check if comment
    if ( line.peek() == '#' ) {
    	cout<<line_c_str<<endl;
      continue;
    }
    line >> targetX[lineN];
    line >> targetY[lineN];

//    targetX[lineN]*=1e9;
//    float gamma=0;//0.1/2.;
    float t=float(lineN)*time_step;
    float a=0.0;
    float omega=sqrt(1-gamma*gamma);
    sampleX[lineN]=exp(-gamma*t)*cos(omega*t-a)-gamma*exp(-gamma*t)*sin(omega*t-a);
    sampleY[lineN]=exp(-gamma*t)*sin(omega*t-a);
    EX[lineN]=t;
    EY[lineN]=pow(targetX[lineN]-sampleX[lineN],2)+pow(targetY[lineN]-sampleY[lineN],2);
    EnergyE[lineN]=pow(targetX[lineN],2)-pow(sampleX[lineN],2)+pow(targetY[lineN],2)-pow(sampleY[lineN],2);
    lineN++;

}
lineN--;
f.close();
//.printf("\nADC_Q: %g\n",ADC_Q);
if(newCanvas){
	MyC=new TCanvas("c1","Histogramy",1596,813);
	//MyC->Divide(2,1);
	//MyC->cd(1);
}
printf("%d\n",lineN);
Cu=new TGraphErrors(lineN,targetX,targetY);
Cu->SetTitle(";Pozycja;Szybkosc");
Cu->GetXaxis()->SetLabelSize(0.05);
Cu->GetXaxis()->SetTitleSize(0.05);
Cu->GetYaxis()->SetLabelSize(0.05);
Cu->GetYaxis()->SetTitleSize(0.05);
Cu->SetLineWidth(2);

Rh=new TGraphErrors(lineN,sampleX,sampleY);
//Cu->GetXaxis()->SetRangeUser(-40e-9,0);
error=new TGraph(lineN,EX,EY);
errorE=new TGraph(lineN,EX,EnergyE);
Cu->SetMarkerStyle(10);
Cu->SetMarkerSize(0.5);
//if(newCanvas)
Cu->Draw("ALP");
//else Cu->Draw("Lsame");
Rh->SetMarkerStyle(15);
Rh->SetMarkerSize(0.5);
Rh->SetLineColor(kRed);
Rh->Draw("LPsame");
return 0;
}
