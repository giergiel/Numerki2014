#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdio>
//#include <thread>
#include <chrono>
#include <ctime>
#include <csignal>

#ifdef LINUX
#include <sys/types.h>
#include <unistd.h>
#endif

using namespace std;

template <class T> void swap ( T& a, T& b ){  T c(a); a=b; b=c;}


const int wymiar=1;//D*2
const int dlugosc[wymiar]={2};
const int wielkosc_tablicy=2;

#define prec double
//TUTAJ TYPDEF? albo przerobić świat na templaty


prec sasiedzi(prec* tablica,int pozycja,int *odleglosc=0)
{//Zwraca wskaznik do odpowiedniego elementu tabeli


//Przetlumacz pozycje na XYZ...
int vec[wymiar];

if(odleglosc==0)return tablica[pozycja];

div_t dv{};dv.quot=pozycja;
for(int i=0;i<wymiar;i++){
	dv=div(dv.quot,dlugosc[i]);
	vec[i]=dv.rem;
	//vec[i]=pozycja%dlugosc[i];//!!Zmienic na dzielenie z reszta w jednym wywolaniu
	//pozycja/=dlugosc[i];
}



//Oddalona pozycja torusowe warunki brzegowe
#ifdef debug
for(int i=0;i<wymiar;i++){
    if(abs(odleglosc[i]>=dlugosc[i])){
	printf("Ponad jeden obiegu torusa odleglosci miedzy punktami w wymiarze %d dla pozycji %d",i,pozycja);
	exit(4);
	}
    }
#endif


for(int i=0;i<wymiar;i++){
//    printf("s:    %d ",odleglosc[i]);
    vec[i]=(vec[i]+odleglosc[i]+dlugosc[i])%dlugosc[i];
    }

//Złóż int pozycja z powrotem - dla wybranego sasiada
pozycja=0;
int temp_mult=1;
for(int i=0;i<wymiar;i++){
	pozycja+=vec[i]*temp_mult;
	temp_mult*=dlugosc[i];
}
//printf("\n%d\n",pozycja);
return tablica[pozycja];
}


void output(char *filename,prec* tablica){
//FORMAT wymiar\n Dlugosci[i]\n WARTOSCI[i]\n
static FILE *out;

if(out==NULL){out = fopen(filename, "w");
    if(out==NULL){
    	fprintf(stderr, "CRITICAL: error opeining file\n");
	exit(2);}
    fprintf(out, "%d\n",wymiar);
    for(int i=0;i<wymiar;i++)fprintf(out, "%d ",dlugosc[i]);
	fprintf(out, "\n");
    }
for(int i=0;i<wielkosc_tablicy;i++){
//	printf("%d\n",i);fflush(stdout);
//	int przesun[wymiar]={0,-1};
	fprintf(out, "%.5f ",sasiedzi(tablica,i) );//!! UWAGA PAMIETAC %lf %f lub analogiczne binarne w outpucie
	}
fprintf(out, "\n");
//fclose(out);
}

prec G_gamma;
void Funkcja(prec *dane,prec *tempF){
prec gamma=G_gamma;
tempF[0]=-gamma*dane[0]-dane[1];
tempF[1]=+dane[0];
}


enum EMethod{ kRK4=1,kIMPEUL=2,kSEMIEUL};


void iteracja(prec *dane,prec *temp,prec *tempF,prec *tempYK,prec deltaT,const int method){
for(int i=0;i<wielkosc_tablicy;i++)temp[i]=dane[i];

if(method==kRK4){
//K1
Funkcja(dane,tempF);
for(int i=0;i<wielkosc_tablicy;i++){
	dane[i]+=deltaT/6.*tempF[i];
	tempYK[i]=temp[i]+deltaT/2.*tempF[i];
}
//k2
Funkcja(tempYK,tempF);
for(int i=0;i<wielkosc_tablicy;i++){
	dane[i]+=deltaT/3.*tempF[i];
	tempYK[i]=temp[i]+deltaT/2.*tempF[i];
}
//k3
Funkcja(tempYK,tempF);
for(int i=0;i<wielkosc_tablicy;i++){
	dane[i]+=deltaT/3.*tempF[i];
	tempYK[i]=temp[i]+deltaT*tempF[i];
}
//k4
Funkcja(tempYK,tempF);
for(int i=0;i<wielkosc_tablicy;i++){
	dane[i]+=deltaT/6.*tempF[i];
}
}
if(method==kIMPEUL){
for(int iter=0;iter<40000;iter++){
	Funkcja(dane,tempF);
	for(int i=0;i<wielkosc_tablicy;i++){
//		if( (iter==0) || (iter==999) )printf("%d: %g %g %g\n",iter,dane[i],temp[i]+deltaT*tempF[i],dane[i]-temp[i]-deltaT*tempF[i]);
		dane[i]=temp[i]+deltaT*tempF[i];
	}
	if( (  fabs( dane[0]-temp[0]-deltaT*tempF[0]) <fabs(temp[0])*1e-14)&&  (  fabs( dane[1]-temp[1]-deltaT*tempF[1] ) <fabs(temp[1])*1e-14))  break;
}
	if( (  fabs( dane[0]-temp[0]-deltaT*tempF[0]) >fabs(temp[0])*1e-14)||  (  fabs( dane[1]-temp[1]-deltaT*tempF[1] ) >fabs(temp[1])*1e-14))  printf("dupa - Fixed Point didnt Converge\n");

}

if(method==kSEMIEUL){
Funkcja(dane,tempF);

dane[0]=temp[0]+deltaT*tempF[0];

for(int iter=0;iter<40000;iter++){
        Funkcja(dane,tempF);
		dane[1]=temp[1]+deltaT*tempF[1];
	if(  fabs( dane[1]-temp[1]-deltaT*tempF[1] ) <fabs(temp[1])*1e-14)  break;
}
	if( fabs( dane[1]-temp[1]-deltaT*tempF[1] ) >fabs(temp[1])*1e-14)  printf("dupa - Fixed Point didnt Converge\n");


}
/*
auto secant=[tempF,dane,temp](auto point, auto direction, auto dt){

auto Val = [tempF,temp,dt](auto p){
Funkcja(p,tempF);
prec ret=0;
	for(int i=0;i<wielkosc_tablicy;i++){
		ret+=pow(temp[i]+dt*tempF[i]-p[i],2);
	}
return ret;
};
prec a[wielkosc_tablicy],b[wielkosc_tablicy],c[wielkosc_tablicy];
for(int i=0;i<wielkosc_tablicy;i++){b[i]=point[i];c[i]=b[i];c[i]+=direction[i];
        if(( fabs(b[i]-c[i])<1e-14 ) && (direction[i])  ){printf("odbyt:%g %g %g %g\n",b[i],c[i],direction[i],dt);}

}
auto scale=Val(b)/( Val(b) - Val(c));
printf("%g %g %g\n",scale,b[1],c[1]);

bool breaker=0;
for(int j=0;j<1000;j++){
    for(int i=0;i<wielkosc_tablicy;i++){
        if(( fabs(b[i]-c[i])<b[i]*1e-14 ) && (direction[i])  ){breaker=1;printf("%d\n",j);break;}
    }
    if(breaker)break;
    for(int i=0;i<wielkosc_tablicy;i++){
        a[i]=b[i]-(b[i]-c[i])*scale;
        c[i]=b[i];
        b[i]=a[i];
        }

}

return point;

};
prec dir[wielkosc_tablicy]={deltaT,deltaT};
prec *aa=secant(dane,dir,deltaT);
printf("%g %g\n",aa[0],aa[1]);
*/
/*
if(method==kIMPEUL2){
auto minimalize= [](auto point,auto direction){
return ret_point;
}

for(int iter=0;iter<40000;iter++){
 	for(int i=0;i<wielkosc_tablicy;i++)tempYK[i]=dane[i];
	Funkcja(dane,tempF);
	for(int i=0;i<wielkosc_tablicy;i++){
        minimalizacja po kierunkach pi
		dane[i]=temp[i]+deltaT*tempF[i];
	}
	wyznacz kierunek pN
	for(int i=0;i<wielkosc_tablicy;i++){
		tempYK[i]=temp[i]+deltaT*tempF[i];
	}
	MINIMALIZACJA KIERUNKU PN


	if( (  fabs( dane[0]-temp[0]-deltaT*tempF[0]) <1e-8)&&  (  fabs( dane[1]-temp[1]-deltaT*tempF[1] ) <1e-8))  break;
}
	if( (  fabs( dane[0]-temp[0]-deltaT*tempF[0]) >1e-8)||  (  fabs( dane[1]-temp[1]-deltaT*tempF[1] ) >1e-8))  printf("dupa\n");

}
*/

}

void LOL_SEGV(int signalnum){
if (signalnum == SIGSEGV) {
cerr << "LOL SEGFAULT"<<
#ifdef LINUX
" on process "<<getpid()<<
#endif
endl;
}
else {
cerr << "Unexpected signal " << signalnum << " received by LOL_SEGV\n";
}

//Here you can generate core dump

signal(signalnum, SIG_DFL);
raise(signalnum);

//You shouldn't reach here if system killed your program after running default handler of SIGSEGV
_Exit(666);
}


int main(int argc, char **argv)
{

if(argc!=5){
	printf("Usage: %s [output] [gamma-tlumienie] [delta T] [steps]\n", argv[0]);
	return -1;
	}

G_gamma=atof(argv[2]);
prec deltaT=atof(argv[3]);
int steps=atoi(argv[4]);
#ifdef LINUX
printf("PID: %d\tOutput: %s\tGamma: %g\tDelta t: %g\tSteps: %d\n",getpid(),argv[1],G_gamma,deltaT,steps);
#else
printf("Output: %s\tGamma: %g\tDelta t: %g\tSteps: %d\n",argv[1],G_gamma,deltaT,steps);
#endif
signal(SIGSEGV,LOL_SEGV);
chrono::time_point<chrono::system_clock> time_start, time_end;
time_start = chrono::system_clock::now();




prec *dane1,*dane2,*dane3,*dane4;
dane1=new prec[wielkosc_tablicy];
dane2=new prec[wielkosc_tablicy];
dane3=new prec[wielkosc_tablicy];
dane4=new prec[wielkosc_tablicy];

for(int i=0;i<wielkosc_tablicy;i++){
    dane1[i]= (i==0?1:0);
    dane2[wielkosc_tablicy-1-i]= (i==0?1:0) ;
    //printf("%f\n",sasiedzi(dane1,i));
}


while(steps--){
	if( (steps%1000)==0)
	output(argv[1],dane1);
	iteracja(dane1,dane2,dane3,dane4,deltaT,kSEMIEUL);
}
time_end = chrono::system_clock::now();
chrono::duration<double> elapsed_seconds = time_end-time_start;
time_t end_time = chrono::system_clock::to_time_t(time_end);
cout << "finished computation at "<<ctime(&end_time)<< "elapsed time: "<<elapsed_seconds.count()<< "s\n";

return 0;


}
