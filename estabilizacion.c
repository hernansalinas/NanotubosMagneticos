/*
==================================================
  Nanotubos Magnéticos-Zigzag, 
==================================================

Este programa realiza la simulación de nanotubos magnético con bordes
zigzag empleando el metodo Monte Carlo con algoritmo de metropolis y
un modelo de Heisenberg tridimensional, el hamiltoniano viene dado por
H= Edip+Eint. Se pueden realizar dos protocolos de enfriamento, rampa
lineal y "anealing". El programa permite realizar variaciones para
muestrear los espines en cono que disminuye su ángulo sólido conforme
avanza la rampa de temperatura, o el muestreo usual en la esfera de
Bloch.

---------------------------
  Instrucciones de uso.
------------------------------

++ Definir las constantes de uso del programa en los #define que se
   encuentras mas abajo.
++ Compilar el programa.      
 

Para compilar el programa use:
  $ make estabilziacion.out
Para ejecutar el programa use:
  $ GSL_RNG_SEED:numero ./estabilizacion.out radio altura altura Fmax Fmin deltaF J D >out.log
  


  Autores: 
  Hernán Salinas J.
  e-mail:hds087@gmail.com
  Oscar Iglesias.
  e-mail:oscar@ffn.ub.es
  Johans Restrepo
  jrestre@gmail.com
  Fecha Revisión :Julio 21 de 2017
===================================================
*/

//=================================================
//LIBRERIAS
//=================================================
#include<stdio.h>
#include <unistd.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h>
//==================================================
//==================================================
//CONDICIONES INICIALES
//==================================================

//Temperatura
#define Tf    50.0      
#define Ti    0.2      
#define dT    0.2      

#define alpha1  0.97   //Rampa logaritmica de temperatura
#define kf      200    //Value for anealing alpha**k     

#define pmMax 15000        //Monte Carlo Step
#define estab 13000       //stability for mean
#define kb 8.617343E-2 //Boltzman constant meV/K 

//#define atpercapa 11    //Number of atoms per layer
//#define hini 10        //H initial
//#define hfin 30        //H final
#define dh   1         //Variation de dh
//--------------------------------------------------
//Parameter of the system.
//--------------------------------------------------
#define a  0.5  //Red parameter
//#define J2  10.0 //Exchange interaction.
//#define D2  0.15  //Dipolar interaction     	    

//--------------------------------------------------
//Value for model-Heisinberg or Ising
//--------------------------------------------------
#define ps 1.0 //p=1.0 spin-up,p=0.0 spin-random, p=[0,1] spin-pseudorandom
//ps=2.0 horizontal
#define HI 1.0 //HI define el modelo, HI=1.0 Heisenberg
               //HI define el modelo, HI=0.0 Ising
//--------------------------------------------------
#define contador 0
#define dim 3
//=====================================

//=====================================
#define dircampo 1 //0 Direccion x, 1 Direccion z del Campo
#define stepH 0.1  //Paso del campo.
#define kc    20   //No se usa.


//=====================================


//Routines 
//=====================================
//#include <rutinas/geometric.h>  Tube with combination n,n
#include <rutinas/geometria.h>  
#include <rutinas/allocate.h>  
#include <rutinas/Wij.h> //Geometric Matrix for dipolar energy calculation.
#include <rutinas/mag_update_mag.h>//Magnetization update
#include <rutinas/spin_3d.h> //Spin assignment
#include <rutinas/update_campo.h>//Magnetic Field update
#include <rutinas/campo.h>     //Magnetic field calculation
#include <rutinas/energia.h>   //Energy calculation
#include <rutinas/vecinos.h>   //neighbors calculation
#include <rutinas/rotacional.h> //Rotational calculation

//=====================================



double desv_sus(double A2,double Apro,double kbt, double atomos);

int main (int argc, char *argv[])
{
  char filename0[500];
  char filename1[500];


  FILE *f2,*fs1;
  int cycle;
  cycle=contador;
  //========================================
  int lim,lim1,i,pm,f;//l;//j
  int ii;
  int count1;
  int k;

  int alt,layer,alayer;

  int numberA; //Atomos por capa 
  int Natomos; //Numero de atomos total
  int nnrows; 
  int **veci,*Nv;

  double **Rpos,*r2;
  double **Wxx,**Wyy,**Wzz,**Wxy,**Wyz,**Wxz;

  double **S;
  double sax, say, saz;
  double S0x,S0y,S0z;
  double dSx,dSy,dSz;
  double dEdip,dEint,dET;
  double Edip,Exch,ET,ET2;
  double e1,e22,Epro,Ecum22,Ecum,calor,calor2;//Epro22;

  double **hint,**hdip;  
  double *M;
  double Minty,Mintz,Mintx,Mitotal,mul,Miy,Mix,Miz;
  double m1z,m1y,m1x,m1total;
  double Magnetizacionz,susceptibilidadz;
  double Magnetizaciony,susceptibilidady;
  double Magnetizacionx,susceptibilidadx;
  double Magnetizaciontotal,Mtotal,susceptibilidadtotal;
 
  double rot1,rot2,rot3;
  double rot11,rot22,rot33;
  double rotax,parx;  
  double rotay,pary;  
  double rotaz,parz;  
  double rhox,rhoy,rhoz;
  double Rotacionaltotal,rotacionalfinal;
  double rotitotal,rottotal,Rtotal;

  double *theta;// psi[dimension],
  double q,atomos;
  double rmax,raprox; 
  double tt,pp,t,p;  //Variables para el radio promedio
  double newvalue,newvalue1;
  double rr,W,kbt;
  int Nk;


  double dEcampo,H, Hz,Hx;
  //========================================
  const gsl_rng_type * T;
  gsl_rng * r;	         
  gsl_rng_env_setup();
  T = gsl_rng_default;     //Semilla por default
  r = gsl_rng_alloc (T);   //Para variar la semilla GSL_RNG_SEED 

  long int seed;
  seed=time(NULL)*getpid();
  gsl_rng_set(r,seed);
  //========================================
  int atpercapa,hini,hfin;
  double J,D;
  double Fmax,Fmin,deltaF; //Campo magnético

  
  if(argc<6 || argc>9){
    fprintf(stdout,"Error\n");
    fprintf(stdout,"Execute ./estabilizacion.out radio h_ini hfin Fmax Fmin DF J D\n");
    exit(1);
  }
  
  atpercapa=atoi(argv[1]);
  hini=atoi(argv[2]);
  hfin=atoi(argv[3]);

  Fmax=atof(argv[4]);
  Fmin=atof(argv[5]);
  deltaF=atof(argv[6]);
  J=atof(argv[7]);
  D=atof(argv[8]);

  numberA=atpercapa;  
  for (alt=hini;alt<=hfin;alt=alt+dh){//layer number variation
    
    fprintf(stdout,"Layer: %d\n",alt);
    nnrows=alt*atpercapa;//System dimension, Total atoms

    //Initialization in zeros, all vector
    Rpos=calloc_multD(nnrows,dim); 
    S   =calloc_multD(nnrows,dim);  
    veci=calloc_multD_int(nnrows,4);
    Nv  =(int *)calloc(nnrows,sizeof(int));
    theta  =(double *)calloc(nnrows,sizeof(double));
    r2  =(double *)calloc(nnrows,sizeof(double));
    Wxx = calloc_multD(nnrows,nnrows);
    Wyy = calloc_multD(nnrows,nnrows);
    Wzz = calloc_multD(nnrows,nnrows);
    Wxy = calloc_multD(nnrows,nnrows);
    Wyz = calloc_multD(nnrows,nnrows);
    Wxz = calloc_multD(nnrows,nnrows);
    hint= calloc_multD(nnrows,3);
    hdip= calloc_multD(nnrows,3);
    M   =(double *)calloc(nnrows,sizeof(double));
    
    geometria(Rpos,&alayer,&layer,&numberA,alt,alt);//Tube Construction 
    Natomos=layer*alayer;  //take the same value the nrows

    if(Natomos!=nnrows){
      fprintf(stdout,"Error: Number of atoms!!!!");
      exit(1);
    }
   
    Asignacionspin(r,S,nnrows);
    // Asignacionspin3d(r,S,Natomos,psi,theta);
    // asignacion_vortice(r,S,atpercapa,alt,theta11);    

    vecinos(Rpos,&Natomos,veci,Nv);
    Wij(Wxx,Wyy,Wzz,Wxy,Wxz,Wyz,Rpos,r2,&lim,&lim1,&Natomos);//Matrices de distancias para el calculo de la energía dipolar
    hexch(S,hint,veci,&Natomos,Nv); //Campo intercambio    
    hdipolar(S,Wxx,Wyy,Wzz,Wxy,Wxz,Wyz,hdip,&Natomos);//Campo Dipolar    

    Energia(S,hdip,hint,&Edip,&Exch,&ET,&Natomos,J,D);    
    Magnetizacion(S,M,&Natomos);

    ET2=ET;
    newvalue=0.0;
    rmax=2.0;

    
    sprintf(filename0,"Propiedades/prop_%.3d_%.3d.dat",atpercapa,alt);
    f2=fopen(filename0,"w");


    kbt=0.113514;
    
    Nk=fabs((Fmax+Fmin))/deltaF;
    
    H=Fmax+deltaF;
    Nk=300;
    for(k=0;k<=Nk;k++){

      if(H>2.0)
	deltaF=0.5;
      if(H>-3.0 && H<=-1.0)
	deltaF=0.5;//0.05

      if(H>-3.0 && H<=-7.0)
	deltaF=0.1;//0.01
      
      if(H<-7.0)
	deltaF=0.1;//0.2

      H=H-deltaF;

      
      
      if(dircampo==1){Hz=H;Hx=0.0;}
      else{Hx=H;Hz=0.0;}
      
      
    //    for (k=1;k<=kf;k++){ //simulated annealing
    //  kbt=Tf*pow(alpha1,k); //alpha=0.97
      // for (kbt=Tf;kbt>=Ti;kbt=kbt-dT){// Linear loop  de Temperatura
      fprintf(stdout,"H=%lf  iter=%d-%d...alt=%d-%d\n",H,k,Nk,alt,hfin);
      
      Minty=0.0;
      Mintz=0.0;
      Mintx=0.0;
      m1z=0.0;
      m1y=0.0;
      m1x=0.0;
      Mtotal=0.0;
      //    my=0.0;
      //    mz=0.0;
      //    mx=0.0;
      m1total=0.0;
      Ecum=0.0;
      Ecum22=0.0;
      e1=0.0;
      e22=0.0;
      pp=0.0;
      tt=0.0;
      rot1=0.0;
      rot2=0.0;
      rot3=0.0;
      rot11=0.0;
      rot22=0.0;
      rot33=0.0;
      Rtotal=0.0;
      rottotal=0.0;
      newvalue=0.0;
      newvalue1=0.0;
      count1=0.0;
      f=0.0;
      Rtotal=0.0;
      rottotal=0.0;

      for(pm=1;pm<=pmMax;pm++){// loop  de pasos de Monte Carlo
	//fprintf(stdout,"Temperatura %lf ...step...%d Atomos %d\n",kbt,pm,Natomos);
	p=0.0;
	t=0.0;
	//      at=0;
	//      j=1;
	newvalue=0.0;
	newvalue1=0.0;
	for(ii=0;ii<Natomos;ii++){ 
	  i=ii; 
	  S0x=S[i][0];
	  S0y=S[i][1];
	  S0z=S[i][2];

	  //	  spin_3d(xr,&sax,&say,&saz);//Espín Supuesto	
	  suma_vec(r,S0x,S0y,S0z,&sax,&say,&saz,&theta[i],rmax,&raprox);
	  //	fprintf(stdout,"Rmax %lf %lf\n",rmax,raprox);

	  dSx=sax-S0x;//Cambio de espin
	  dSy=say-S0y;
	  dSz=saz-S0z;
	  dEint=-2*J*(dSx*hint[i][0]+dSy*hint[i][1]+dSz*hint[i][2]);
	  dEdip=-2*D*(dSx*hdip[i][0]+dSy*hdip[i][1]+dSz*hdip[i][2]);
	  dEcampo=-(dSx*Hx+dSz*Hz);
	  dET=dEint+dEdip+dEcampo;

	  if(dET<=0.0) {   
	    	    if(newvalue<theta[i]){
	      newvalue=theta[i];
	     }	//	fprintf(stdout,"raprox=%lf\n",raprox);
	  
	    	    if(newvalue1<=raprox){
	     newvalue1=raprox;
	     }

	    count1=count1+1;//Contador sobre los cambios de espin.
	    S[i][0]=sax;  
	    S[i][1]=say;  
	    S[i][2]=saz;  
	    t=t+theta[i];
	    p=p+raprox;
	    Edip=Edip+dEdip;
	    Exch=Exch+dEint;
	    ET2=ET2+dET;
	    updatehdip(hdip,Wxx,Wyy,Wzz,Wxy,Wxz,Wyz,dSx,dSy,dSz,i,&Natomos);
	    updatehexch(hint,veci,dSx,dSy,dSz,i,Nv);
	    updateMag(dSx,dSy,dSz,M);
	    Energia(S,hdip,hint,&Edip,&Exch,&ET,&Natomos,J,D);
	  }
	  else{
	    rr=gsl_rng_uniform(r); 
	    W=exp(-dET*1.0/kbt); //Probabilidad de Boltzman 
	    if(rr<=W){ 
	    	      if(newvalue<theta[i]){
	    		newvalue=theta[i];
	    }	//	fprintf(stdout,"raprox=%lf\n",raprox);
	    
	      if(newvalue1<=raprox){
	     	newvalue1=raprox;
	      }
	      count1=count1+1;
	      S[i][0]=sax;  
	      S[i][1]=say;  
	      S[i][2]=saz;  
	      t=t+theta[i];
	      p=p+raprox;
	      Edip=Edip+dEdip;
	      Exch=Exch+dEint;
	      ET2=ET2+dET;
	      updatehdip(hdip,Wxx,Wyy,Wzz,Wxy,Wxz,Wyz,dSx,dSy,dSz,i,&Natomos);
	      updatehexch(hint,veci,dSx,dSy,dSz,i,Nv);
	      updateMag(dSx,dSy,dSz,M); 
	      Energia(S,hdip,hint,&Edip,&Exch,&ET,&Natomos,J,D);
	    }//Final if
	  }//Final else
	  t=theta[i];
	  //	  p=psi[i];	  
	}//Final for, recorre cada espin ii;
	
	if (pm>estab){    
	  Mix=M[0];
	  Miy=M[1];
	  Miz=M[2];
	  Mitotal=sqrt(Mix*Mix+Miy*Miy+Miz*Miz); //Magnetización absoluta

	  Mintz=Miz+Mintz;    
	  m1z=Miz*Miz+m1z;

	  Minty=Miy+Minty;
	  m1y=Miy*Miy+m1y;

	  Mintx=Mix+Mintx;
	  m1x=Mix*Mix+m1x;

	  Mtotal=Mtotal+Mitotal;
	  m1total=Mitotal*Mitotal+m1total;
	
	  Ecum=ET2+Ecum;
	  Ecum22=ET2+Ecum22;

	  e1=ET*ET+e1;
	  e22=ET2*ET2+e22;
	  tt=newvalue+tt;
	  pp=newvalue1+pp;
	  //pp=raprox+pp;
	  rot(Rpos,S,&Natomos,&rhox,&rhoy,&rhoz,alt,atpercapa);
	  rotitotal=sqrt(rhox*rhox+rhoy*rhoy+rhoz*rhoz);
	  rottotal=rotitotal+rottotal;
	  Rtotal=rotitotal*rotitotal+Rtotal;
	  
	  rot1=rhox+rot1;
	  rot2=rhoy+rot2;
	  rot3=rhoz+rot3;
	  rot11=rhox*rhox+rot11;
	  rot22=rhoy*rhoy+rot22;
	  rot33=rhoz*rhoz+rot33;
	}//Final if estabilizacion
      }//FINAL PASO MONTE CARLO.


      if(cycle%2==0.0){
	
	sprintf(filename1,"Conf_espin/conf%.3d_%.3d_%.3d.dat",atpercapa,alt,cycle);
	fs1=fopen(filename1,"w");
	//  l=0.0;
	for(f=0;f<Natomos;f++){
	  fprintf(fs1,"%lf %lf %lf %lf %lf %lf\n ",
		  Rpos[f][0],Rpos[f][1],Rpos[f][2],
		  S[f][0],S[f][1],S[f][2]);
	}
	fclose(fs1); //
      }
      cycle=cycle+1;      
      
      q=pmMax-estab; //Valores para promediar,
      atomos=Natomos;  //Asignacion que se puede omitir
      mul=q*atomos; //Denominador para calculo de promedios
      
      //      fprintf(stdout,"pm=%d acep=%lf angul=%lf radio=%lf %d\n",(pm-1),count1*1.0/((pm-1)*Natomos)*100,tt/(pm-estab-1),pp/(pm-estab-1),count1/(pm-1));
      
      rmax=pp/(pm-estab-1);
      fprintf(stdout,"rmax= %lf\n",rmax);
      //==================================================
      //SUCEPTIBILIDAD MAGNÉTICA.
      //==================================================    
      rotax=rot1/mul;
      parx=desv_sus(rot11/q,rot1/q,kbt,atomos);

      rotay=rot2/mul;
      pary=desv_sus(rot22/q,rot2/q,kbt,atomos);
      
      rotaz=rot3/mul;
      parz=desv_sus(rot33/q,rot3/q,kbt,atomos);
      
      Rotacionaltotal=rottotal/mul;
      rotacionalfinal=desv_sus(Rtotal/q,rottotal/q,kbt,atomos);
            
      Magnetizacionz=Mintz/mul; //Magnetizaci\'on por sitio.
      susceptibilidadz=desv_sus(m1z/q,Mintz/q,kbt,atomos);
      
      Magnetizaciony=Minty/mul; //Magnetizaci\'on por sitio.
      susceptibilidady=desv_sus(m1y/q,Minty/q,kbt,atomos);

      Magnetizacionx=Mintx/mul; //Magnetizaci\'on por sitio.
      susceptibilidadx=desv_sus(m1x/q,Mintx/q,kbt,atomos);
      
      Magnetizaciontotal=Mtotal/mul; //Magnetizaci\'on por sitio.
      susceptibilidadtotal=desv_sus(m1total/q,Mtotal/q,kbt,atomos);

      //==================================================
      //Calor Especifico
      //==================================================
      Epro=Ecum/q;
      calor=kb*desv_sus(e1/q,Ecum/q,kbt*kbt,atomos);
      //Epro22=Ecum22/q;
      calor2=kb*desv_sus(e22/q,Ecum22/q,kbt*kbt,atomos);
      //==================================================      

      
      fprintf(f2,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",H ,(D/J),Magnetizaciontotal,Epro, Magnetizacionx,Magnetizaciony,Magnetizacionz, susceptibilidadx,susceptibilidady,susceptibilidadz,susceptibilidadtotal,calor,calor2,parx,pary,parz,rotax,rotay,rotaz,rotacionalfinal,Rotacionaltotal,Edip,Exch,rmax);     
    }//Final for campo
    fclose(f2); //

    
          
    free(Rpos);
    free(S);
    free(veci);
    free(Nv);
    free(theta);
    free(r2);
    free(Wxx);
    free(Wyy);
    free(Wzz);
    free(Wxy);
    free(Wyz);
    free(Wxz);
    free(hint);
    free(hdip);
    free(M);
  }//Final for de altura

  gsl_rng_free(r);
  return 0;
}//Final main


double desv_sus(double A2,double Apro,double kbt, double atomos){
  double A_2,prop;
  A_2=Apro*Apro;//<A>2
  prop=(A2-A_2)/(kbt*atomos);
  return prop;
} 
