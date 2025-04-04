const double pi = 3.14159265358979;
const double s_per_y = 3600.0*24.0*365.25;
const double tiny = 1.0e-16;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>


typedef struct {
  double x,y;
  double bed,topice,ice,oldice,newice,iniice,snow,topsnow,phi,psi,isostasy;
  double sedi,dHw,dHs,topsedi,Hgw,hw,dhwdt,hw_old,Pw,srate,rfail,Tair,Ts,Tb,Ta,precip;
  double bsedi,ssedi,msedi,sN10,CNprod,atten,Ka,Kq;
  double meltrate,Mb,Ms,accrate,smelt,mrate,Kgw,tprop,sfac,lee,dhw,te_new,afac,angular;
  double periglacial_erate,periglacial_erosion,fluvial_erate,fluvial_erosion;
  double hillslope_erate,hillslope_erosion,landslide_erosion,weathering,weathering_rate;
  double quarrying_rate,quarrying,abrasion_rate,abrasion,cbs;
  double dhdx,dbdx,dhdy,dbdy,dtdx,dtdy,alpha2,curv,bcurv,bslope,tslope,hslope,slidingslope;
  double exx,eyy,ezz,exy,ezz_s,hcond,hconv,kappa;
  double te2,te2i,td2,te2r,sxx,syy,sxy,szz,sxz,syz;
  double lb,pb,tn,te,te_s,ts,tbx,tby,tbz;
  double vb,vx_b,vy_b,vx_d,vy_d,sliding,deformation,vbres,vsb;
  double dsxxdx,dsyydx,dsxydx,dsxxdy,dsyydy,dsxydy;
  double S,dSdt,mw,hs,Ac,dAcdt,SLf; 
  int hydro,include,avasite;
  int np,npia; //Maxime
  double cx[4],cy[4],kc[7];
  double water,bwater,bQw,Qw,Qw_old,Qstack,qw,Wc,tau,Qt,streamload; //Maxime - bwater, bQw
  int fixflag,ri,rj,di,dj,nparent,streamorder,cnumber,inslide,margin;
  int reciever,nne,ne_i[8],ne_j[8],pa_i[8],pa_j[8];
  double wsurf,lakewater,ndist[8],rdist;
  double Vs[20];
  double Vbsedi;
  double lastslide,theta,bedStore; //Maxime
  double coolfac,bwatersurf,subglacial_erosion,water2surf; //Maxime
  double sdeformation,topicedum; //Maxime
} celltype;

typedef struct {
  double dbdx,dhdx,dtdx,dhdy,alpha2,alpha1,dalpha2dx,dalpha1dx;
  double dH,dHw,dHs,kappa,qwx,dhw,psi;
  double vx,vrx,vx_s,vx_d,vx_b,vz_b,ice,vresx,vx_di,vx_bi; 
  double dsxxdx,dsyydx,dsxydx;
  double cx[4],cy[4],kc[7],wx[10];
} hptype;

typedef struct {
  double dbdy,dhdy,dtdy,dhdx,alpha2,alpha1,dalpha2dy,dalpha1dy;
  double dH,dHw,dHs,kappa,qwy,dhw,psi;
  double vy,vry,vy_s,vy_d,vy_b,vz_b,ice,vresy,vy_di,vy_bi; 
  double dsxxdy,dsyydy,dsxydy;
  double cx[4],cy[4],kc[7],wy[10];
} vptype;

typedef struct {
  double exy,sxy;
  double Qw;
} cornertype;

typedef struct {
  int active,type; /*Maxime*/
  double x,y,z,bf,bx,by,vx,vy,vz;
  double vx_d,vx_b,vy_d,vy_b,velocity; /*Maxime*/
  double birthday,age,dl,afac;
  double sedi,N10,erate;
} particletype;


typedef struct {
  int nx,ny,nc;
  int nFiletime;
  int nmoni,gmode,slidingmode,hydromode;
  int dofluvial,doice,dosliding,doperiglacial,doglacialhydrology,doglacialsedi,doparticles,dodeposit,dodebrisablation;
  int doglacialerosion,doavalance,dohillslope,dohillslopeerosion,dosubglacialdeformation,dosediment,coldbased;
  int dolandslide,doweathering,doisostasy,docelldata,ncelldata,celldata_i[100],celldata_j[100];
  int douniformerosion, dopulse, Flaw, steadybed, gelaw, nbCores,nci; //Maxime
  long stressitt,stressitt_count,veloitt,veloitt_count,hydroitt,nice;
  double dx,dy,hmin,ct,L,H,gravity;
  double maxtime,maxdt,maxdt_ice;
  double meanice,meanele,maxspeed,maxb,maxs,meanice_b;
  double mean_periglacial_erate,mean_dHs_periglacial,mean_dHs_hillslope,mean_dHs_fluvial;
  double mean_quarrying_rate,mean_abrasion_rate,mean_fluvial_erate,mean_hillslope_erate,mean_landslide_erate;
  double meanbed,minbed,maxbed,mean_isostasy,meansedi,meanpw,meante,meanhw,meandhwdt,meanwater,meansliding,mean_weatheringrate;
  double isostasy_rate; //MAxime
  double **filetime;
} meshtype;

typedef struct {
  double gamma,gamma0,Cs,latentheat,ki,rho,cp,ksg;
  double ifac,sfac,vbfac,Ka,ap,Kq,minefac,sedifac;
  double C,L0,minbeta;
  int maxitt_v,maxitt_s;
  double maxdeformation, maxsliding; //Maxime
} iproptype;

typedef struct {
  int nTemp,nMrate_h,nMrate_T,mtype;
  double lrate,avaslope,avacurv,maxacc,dhice,T0;
  double qb,Ldebris;
  double scc; //Maxime
  double **Temp,**Mrate_h,**Mrate_T;
} mproptype;

typedef struct {
  double a2w,tscale,S0,A0,kh,h0,ds,ss0;
  double kmin,minqw,dtw,Ls,Lc,alpha;
  double dss; //MAxime
} hwproptype;

typedef struct {
  int Nc;
  double Ks,sc,Ke,Kw,Ls,gamma;
  double maxsedi;
  double kt; //Maxime
} hproptype;

typedef struct {
  int ncall;
  double Te,rho_r,rho_s,rho_a;
} flexproptype;

typedef struct {
  int npmax;
  int maxp,maxpm,minpm;
  double minice,minsedi,maxsedi,minsedim,L0,vub; /*Maxime*/
} parproptype;

typedef struct {
  double pr,rho_s,Dg,kw,tau_c,Kt,Ke,m;
  int t_sp; //Maxime
} fproptype;


typedef struct {
  int nHs,nT0;
  double rho_b,rho_s,Ke,Kt;
  double maxsedi,maxice,minslope,minsedi;
  double *Hsv,*T0v;
  double **Ci,**Tr;
} pgproptype;

typedef struct {
  int ii,jj;
  double value;
} sortarraytype;

#include "periglacial.c"
#include "hillslope.c"
#include "hydrology.c"
#include "inout.c"
#include "glacial.c"
#include "fluvial.c"
#include "flexure.c"
#include "thermal.c"
#include "particles.c"

void get_gradients(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh)
{

  int i,j;
  double meanbed = 0.0;
  double meansedi = 0.0;
  double lmaxbed = -10e3;
  double lminbed = 10e3;
  double maxbed = -10e3;
  double minbed = 10e3;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;

  double dx = (*mesh).dx;
  double dy = (*mesh).dy;

  double maxb = (*mesh).maxb;
  double maxs = (*mesh).maxs;


#pragma omp parallel shared(cells,vp,hp,meanbed,minbed,maxbed) private(i,j) firstprivate(nx,ny,dx,dy,maxs,maxb,lminbed,lmaxbed)
  {

    /*update elevations*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].topsedi = cells[i][j].bed + cells[i][j].sedi;
	cells[i][j].topice = cells[i][j].topsedi + cells[i][j].ice;
      }
    }



    /*compute gradients in x dir*/
#pragma omp for schedule(static) nowait
    for (i=0;i<ny;i++) {
      for (j=1;j<nx;j++) {
	hp[i][j].dbdx = (cells[i][j].bed-cells[i][j-1].bed)/dx;
	//hp[i][j].dbdx = (cells[i][j].topsedi-cells[i][j-1].topsedi)/dx; 
	hp[i][j].dhdx = (cells[i][j].topice-cells[i][j-1].topice)/dx;
	hp[i][j].dtdx = (cells[i][j].topsedi-cells[i][j-1].topsedi)/dx;
	hp[i][j].ice = 0.5*(cells[i][j].ice+cells[i][j-1].ice);
	if (hp[i][j].dbdx < -maxb) hp[i][j].dbdx = -maxb;
	if (hp[i][j].dbdx > maxb) hp[i][j].dbdx = maxb;
	if (hp[i][j].dhdx < -maxs) hp[i][j].dhdx = -maxs;
	if (hp[i][j].dhdx > maxs) hp[i][j].dhdx = maxs;
      }
      hp[i][0].dhdx = hp[i][1].dhdx;
      hp[i][0].dbdx = hp[i][1].dbdx;
      hp[i][0].dtdx = hp[i][1].dtdx;
      
      hp[i][nx].dhdx = hp[i][nx-1].dhdx;
      hp[i][nx].dbdx = hp[i][nx-1].dbdx;
      hp[i][nx].dtdx = hp[i][nx-1].dtdx;

    }

    /*compute gradients in y dir*/
#pragma omp for schedule(static)
    for (i=1;i<ny;i++) {
      for (j=0;j<nx;j++) {
	vp[i][j].dbdy = (cells[i][j].bed-cells[i-1][j].bed)/dy;
	//vp[i][j].dbdy = (cells[i][j].topsedi-cells[i-1][j].topsedi)/dy;
	vp[i][j].dhdy = (cells[i][j].topice-cells[i-1][j].topice)/dy;
	vp[i][j].dtdy = (cells[i][j].topsedi-cells[i-1][j].topsedi)/dy;
	vp[i][j].ice = 0.5*(cells[i][j].ice+cells[i-1][j].ice);
	if (vp[i][j].dbdy < -maxb) vp[i][j].dbdy = -maxb;
	if (vp[i][j].dbdy > maxb) vp[i][j].dbdy = maxb;
	if (vp[i][j].dhdy < -maxs) vp[i][j].dhdy = -maxs;
	if (vp[i][j].dhdy > maxs) vp[i][j].dhdy = maxs;
      }
    }
    for (j=0;j<nx;j++) {
      vp[0][j].dhdy = vp[1][j].dhdy;
      vp[0][j].dbdy = vp[1][j].dbdy;
      vp[0][j].dtdy = vp[1][j].dtdy;
      vp[ny][j].dhdy = vp[ny-1][j].dhdy;
      vp[ny][j].dbdy = vp[ny-1][j].dbdy;
      vp[ny][j].dtdy = vp[ny-1][j].dtdy;
      }


#pragma omp for schedule(static) nowait
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].dhdx = 0.5*(hp[i][j].dhdx+hp[i][j+1].dhdx);
	cells[i][j].dbdx = 0.5*(hp[i][j].dbdx+hp[i][j+1].dbdx);
	cells[i][j].dtdx = 0.5*(hp[i][j].dtdx+hp[i][j+1].dtdx);
      }
    }
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].dhdy = 0.5*(vp[i][j].dhdy+vp[i+1][j].dhdy);
	cells[i][j].dbdy = 0.5*(vp[i][j].dbdy+vp[i+1][j].dbdy);
	cells[i][j].dtdy = 0.5*(vp[i][j].dtdy+vp[i+1][j].dtdy);
	cells[i][j].bslope = sqrt(pow(cells[i][j].dbdx,2.0)+pow(cells[i][j].dbdy,2.0));
	cells[i][j].hslope = sqrt(pow(cells[i][j].dhdx,2.0)+pow(cells[i][j].dhdy,2.0));
	cells[i][j].tslope = sqrt(pow(cells[i][j].dtdx,2.0)+pow(cells[i][j].dtdy,2.0));
	cells[i][j].alpha2 = pow(cells[i][j].dhdx,2.0)+pow(cells[i][j].dhdy,2.0);
	if (cells[i][j].sliding > 1.0e-6) 
	  cells[i][j].slidingslope = -(cells[i][j].dtdx*cells[i][j].vx_b+cells[i][j].dtdy*cells[i][j].vy_b)/cells[i][j].sliding; 
	else cells[i][j].slidingslope = 0.0;
      }
    }

#pragma omp for schedule(static)
    for (i=1;i<ny-1;i++) {
      for (j=1;j<nx-1;j++) {
	cells[i][j].curv = (cells[i][j+1].topice-2.0*cells[i][j].topice+cells[i][j-1].topice)/(dx*dx)+(cells[i+1][j].topice-2.0*cells[i][j].topice+cells[i-1][j].topice)/(dy*dy);
	cells[i][j].bcurv = (cells[i][j+1].bed-2.0*cells[i][j].bed+cells[i][j-1].bed)/(dx*dx)+(cells[i+1][j].bed-2.0*cells[i][j].bed+cells[i-1][j].bed)/(dy*dy);
      }
    }


#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=1;j<nx;j++) {
	hp[i][j].dhdy = .5*(cells[i][j-1].dhdy+cells[i][j].dhdy);
	hp[i][j].alpha2 = pow(hp[i][j].dhdx,2.0)+pow(hp[i][j].dhdy,2.0); 
	hp[i][j].dalpha2dx = (cells[i][j].alpha2-cells[i][j-1].alpha2)/dx;
	hp[i][j].alpha1 =  (cells[i][j].curv+cells[i][j-1].curv)/2.0;
	hp[i][j].dalpha1dx = (cells[i][j].curv-cells[i][j-1].curv)/dx;
      }
    }

#pragma omp for schedule(static)
    for (i=1;i<ny;i++) {
      for (j=0;j<nx;j++) {
	vp[i][j].dhdx = .5*(cells[i-1][j].dhdx+cells[i][j].dhdx);
	vp[i][j].alpha2 = pow(vp[i][j].dhdx,2.0)+pow(vp[i][j].dhdy,2.0);
	vp[i][j].dalpha2dy = (cells[i][j].alpha2-cells[i-1][j].alpha2)/dy;
	vp[i][j].alpha1 = (cells[i][j].curv+cells[i-1][j].curv)/2.0;
	vp[i][j].dalpha2dy = (cells[i][j].curv-cells[i-1][j].curv)/dy;
      }
    }


  /*compute mean, min, and max metrics*/
#pragma omp for schedule(static) reduction(+:meanbed) reduction(+:meansedi)
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      meanbed += cells[i][j].bed;
      meansedi += cells[i][j].sedi;
      if (cells[i][j].bed < lminbed) lminbed = cells[i][j].bed;
      if (cells[i][j].bed > lmaxbed) lmaxbed = cells[i][j].bed; 
    }
  }
#pragma omp critical 
      { 
	if (lmaxbed > maxbed) maxbed = lmaxbed; 
	if (lminbed < minbed) minbed = lminbed; 
      }

  }/*pragma*/

  meanbed /= (double)((*mesh).nc);
  meansedi /= (double)((*mesh).nc);
  (*mesh).meanbed = meanbed;
  (*mesh).minbed = minbed;
  (*mesh).maxbed = maxbed;
  (*mesh).meansedi = meansedi;

}/*get_gradients*/




void initialize(celltype **cells,hptype **hp,vptype **vp,cornertype **cp,meshtype *mesh,hwproptype hwprop, char *path,char file[200],FILE *fp)
{

  int i,j,k;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;


  printf("nx = %d, ny = %d\n",nx,ny);

  /*read bed topography etc*/
  sprintf(file,"%s/input/meshdata.input",path); /*filename*/
  if ((fp = fopen(file,"rb")) == NULL) {
    printf("Cannot open file : %s/input/meshdata.input\n",path);
    exit(1);
  }
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].bed,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].ice,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].sedi,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].precip,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].include,sizeof(int),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].mrate,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].srate,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].fixflag,sizeof(int),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].phi,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].CNprod,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].atten,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].Ka,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].Kq,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].isostasy,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].abrasion,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].fluvial_erosion,sizeof(double),1,fp);
  for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].hillslope_erosion,sizeof(double),1,fp);
	
  fclose(fp);
  
  /*read englacial sediment concentration*/
  if ((*mesh).doglacialsedi > 0) {

    sprintf(file,"%s/input/Vsdata.input",path); /*filename*/
    if ((fp = fopen(file,"rb")) == NULL) {
      printf("Cannot open file : %s/input/Vsdata.input\n",path);
      exit(1);
    }
    for (k=0;k<20;k++) {
      for (j=0;j<nx;j++) for (i=0;i<ny;i++) fread(&cells[i][j].Vs[k],sizeof(double),1,fp);
    }/*k*/
  }/*if*/
    

  /*read cellnumbers for celldata*/
  if ((*mesh).docelldata > 0) {
    sprintf(file,"%s/input/cellnumbers.input",path); /*filename*/
    if ((fp = fopen(file,"rb")) == NULL) {
      printf("Cannot open file : %s/input/meshdata.input\n",path);
      exit(1);
    }
    fread(&(*mesh).ncelldata,sizeof(int),1,fp);
    for (i=0;i<(*mesh).ncelldata;i++) {
      fread(&(*mesh).celldata_i[i],sizeof(int),1,fp);
      fread(&(*mesh).celldata_j[i],sizeof(int),1,fp);
    }
  }

#pragma omp parallel shared(cells,vp,hp,cp,mesh) private(i,j) firstprivate(nx,ny,dx,dy)
  {

    /*initialize geometry and stress*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].te2 = 1.0;
	cells[i][j].te2r = 0.0;
	cells[i][j].td2 = 0.0;
	cells[i][j].exx = 0.0;
	cells[i][j].eyy = 0.0;
	cells[i][j].exy = 0.0;
	cells[i][j].ezz_s = 0.0;
	cells[i][j].sxx = 0.0;
	cells[i][j].syy = 0.0;
	cells[i][j].sxy = 0.0;
	cells[i][j].szz = 0.0;
	cells[i][j].sxz = 0.0;
	cells[i][j].syz = 0.0;
	cells[i][j].tn = 0.0;
	cells[i][j].pb = 0.0;
	cells[i][j].te = 0.0*cells[i][j].ice;
	cells[i][j].tbx = 0.0;
	cells[i][j].tby = 0.0;
	cells[i][j].tbz = 0.0;
	cells[i][j].ts = 0.0;
	cells[i][j].vb = 0.0;
	cells[i][j].dsxxdx = 0.0;
	cells[i][j].dsyydx = 0.0;
	cells[i][j].dsxydx = 0.0;
	cells[i][j].dsxxdy = 0.0;
	cells[i][j].dsyydy = 0.0;
	cells[i][j].dsxydy = 0.0;
	cells[i][j].cx[0] = 0.0;
	cells[i][j].cx[1] = 0.0;
	cells[i][j].cy[0] = 0.0;
	cells[i][j].cy[1] = 0.0;
	cells[i][j].kc[0] = 0.0;
	cells[i][j].kc[1] = 0.0;
	cells[i][j].kc[2] = 0.0;
	cells[i][j].x = (double)j*dx + 0.5*dx;
	cells[i][j].y = (double)i*dy + 0.5*dy;
	cells[i][j].ssedi = 0.0;
	cells[i][j].bsedi = 0.0;
	cells[i][j].msedi = 0.0;
	cells[i][j].np = 0;
	cells[i][j].npia = 0; //Maxime
	cells[i][j].topsedi = cells[i][j].bed + cells[i][j].sedi;
	cells[i][j].topice = cells[i][j].topsedi + cells[i][j].ice;
	cells[i][j].oldice = cells[i][j].ice;
	cells[i][j].newice = cells[i][j].ice;
	cells[i][j].snow = 0.0;
	cells[i][j].periglacial_erosion = 0.0;
	cells[i][j].quarrying = 0.0;
	//cells[i][j].abrasion = 0.0;
	//cells[i][j].fluvial_erosion = 0.0;
 	cells[i][j].weathering = 0.0;
	cells[i][j].weathering_rate = 0.0;
	cells[i][j].vx_d = 0.0;
	cells[i][j].vy_d = 0.0;
	cells[i][j].vx_b = 0.0;
	cells[i][j].vy_b = 0.0;
	cells[i][j].deformation = 0.0;
	cells[i][j].sliding = 0.0;
	cells[i][j].vbres = 0.0;
	cells[i][j].Hgw = cells[i][j].bed;
	cells[i][j].hw = 0.0;
	cells[i][j].Pw = 0.0;
	cells[i][j].S = 1.0e-3;
	cells[i][j].SLf = 0.0;
	cells[i][j].dSdt = 0.0;
	cells[i][j].dbdx = 0.0;
	cells[i][j].dbdy = 0.0;
	cells[i][j].dhdx = 0.0;
	cells[i][j].dhdy = 0.0;
	cells[i][j].bslope = 0.0;
	cells[i][j].hslope = 0.0;
	cells[i][j].tslope = 0.0;
	cells[i][j].cnumber = -1;
	// cells[i][j].isostasy = 0.0;
	cells[i][j].curv = 0.0;
	cells[i][j].Ts = 0.0;
	cells[i][j].Kgw = hwprop.kmin;
	cells[i][j].Mb = 0.0;
	cells[i][j].Ms = 0.0;
	cells[i][j].sfac = 0.0;
	cells[i][j].margin = 0;
	cells[i][j].alpha2 = 0.0;
	cells[i][j].Qstack = 0.0;
	cells[i][j].afac = 0.0;
	cells[i][j].cbs = 0.0;
	cells[i][j].bwater = 0; //Maxime
	cells[i][j].Vbsedi = 0.0; /*Maxime*/
	cells[i][j].iniice = cells[i][j].ice;
	cells[i][j].lastslide = 0.0; //Maxime
	cells[i][j].theta = 1.0; //Maxime
	cells[i][j].bedStore = cells[i][j].bed; //Maxime
	cells[i][j].coolfac = 0.0; //Maxime
	cells[i][j].bwatersurf = 0.0; //Maxime
	cells[i][j].subglacial_erosion = 0.0; //Maxime
	cells[i][j].water2surf = 0.0; //Maxime
	cells[i][j].sdeformation = 0.0; //Maxime
	cells[i][j].topicedum = 0.0; //Maxime

	for (k=0;k<20;k++) cells[i][j].Vs[k] = 0.0;
      }
    }

    /*initialize*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx+1;j++) {
	hp[i][j].vx = 0.0;
	hp[i][j].vx_s = 0.0;
	hp[i][j].vx_d = 0.0;
	hp[i][j].vx_b = 0.0;
	hp[i][j].vz_b = 0.0;
	hp[i][j].dhdx = 0.0;
	hp[i][j].dtdx = 0.0;
	hp[i][j].dbdx = 0.0;
	hp[i][j].dhdy = 0.0;
	hp[i][j].alpha1 = 0.0;
	hp[i][j].alpha2 = 0.0;
	hp[i][j].dalpha1dx = 0.0;
	hp[i][j].dalpha2dx = 0.0;
	hp[i][j].dsxxdx = 0.0;
	hp[i][j].dsyydx = 0.0;
	hp[i][j].dsxydx = 0.0;
	hp[i][j].cx[0] = 0.0;
	hp[i][j].cx[1] = 0.0;
	hp[i][j].cy[0] = 0.0;
	hp[i][j].cy[1] = 0.0;
	hp[i][j].kc[0] = 0.0;
	hp[i][j].kc[1] = 0.0;
	hp[i][j].kc[2] = 0.0;
	hp[i][j].wx[0] = 0.0;
	hp[i][j].wx[1] = 0.0;
	hp[i][j].wx[2] = 0.0;
	hp[i][j].wx[3] = 0.0;
	hp[i][j].dH = 0.0;
	hp[i][j].dHw = 0.0;
	hp[i][j].dHs = 0.0;
	hp[i][j].kappa = 0.0;
	hp[i][j].vrx = 0.0;
	hp[i][j].vresx = 0.0;
	hp[i][j].ice = 0.0;
      }
    }	
#pragma omp for schedule(static)
    for (i=0;i<ny+1;i++) {
      for (j=0;j<nx;j++) {
	vp[i][j].vy = 0.0;
	vp[i][j].vy_s = 0.0;
	vp[i][j].vy_d = 0.0;
	vp[i][j].vy_b = 0.0;
	vp[i][j].vz_b = 0.0;
	vp[i][j].dhdy = 0.0;
	vp[i][j].dtdy = 0.0;
	vp[i][j].dbdy = 0.0;
	vp[i][j].dhdx = 0.0;
	vp[i][j].alpha1 = 0.0;
	vp[i][j].alpha2 = 0.0;
	vp[i][j].dalpha1dy = 0.0;
	vp[i][j].dalpha2dy = 0.0;
	vp[i][j].dsxxdy = 0.0;
	vp[i][j].dsyydy = 0.0;
	vp[i][j].dsxydy = 0.0;
	vp[i][j].cx[0] = 0.0;
	vp[i][j].cx[1] = 0.0;
	vp[i][j].cy[0] = 0.0;
	vp[i][j].cy[1] = 0.0;
	vp[i][j].kc[0] = 0.0;
	vp[i][j].kc[1] = 0.0;
	vp[i][j].kc[2] = 0.0;
	vp[i][j].wy[0] = 0.0;
	vp[i][j].wy[1] = 0.0;
	vp[i][j].wy[2] = 0.0;
	vp[i][j].wy[3] = 0.0;
	vp[i][j].dH = 0.0;
	vp[i][j].dHw = 0.0;
	vp[i][j].dHs = 0.0;
	vp[i][j].kappa = 0.0;
	vp[i][j].alpha2 = 0.0;
	vp[i][j].vry = 0.0;
	vp[i][j].vresy = 0.0;
	vp[i][j].ice = 0.0;
      }
    }	
    

#pragma omp for schedule(static)
    for (i=0;i<ny+1;i++) {
      for (j=0;j<nx+1;j++) {
	cp[i][j].exy = 0.0;
	cp[i][j].sxy = 0.0;
	cp[i][j].Qw = 0.0;
      }
    }
  }/*pragma*/

  (*mesh).meanice_b = 0.0;
  (*mesh).stressitt_count = 1;
  (*mesh).veloitt_count = 1;
  (*mesh).stressitt = 0;
  (*mesh).veloitt = 0;
  (*mesh).mean_isostasy = 0.0;
  (*mesh).hydroitt = 0;
  (*mesh).isostasy_rate = 0.0; //Maxime
  (*mesh).nci = 0; //Maxime


}/*initialize*/


void get_sealevel_temperature(mproptype *mprop,double time)
{

  /*interpolate sealevel temperature*/
  int i;
  double fac;

  i = 0;
  while ((time > (*mprop).Temp[0][i])&&(i < (*mprop).nTemp-1)) i += 1;
  if (i == 0)
    (*mprop).T0 = (*mprop).Temp[1][0];
  else if (i == (*mprop).nTemp)
    (*mprop).T0 = (*mprop).Temp[1][i-1];
  else {
    fac = (time - (*mprop).Temp[0][i-1])/((*mprop).Temp[0][i]-(*mprop).Temp[0][i-1]);
    (*mprop).T0 = (1.0-fac)*(*mprop).Temp[1][i-1] + fac*(*mprop).Temp[1][i];
  }

}

void get_surface_temperature(meshtype mesh,mproptype mprop,celltype **cells)
{

  int i,j;
  int nx = mesh.nx;
  int ny = mesh.ny;
  double h;

  double lrate = mprop.lrate;
  double T0 = mprop.T0;
  double dhice = mprop.dhice;

  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {

      h = cells[i][j].topsedi + dhice*cells[i][j].ice;
      cells[i][j].Tair = T0 - lrate*h;
      if ((cells[i][j].ice > 10.0)&&(cells[i][j].Tair > 0.0)) cells[i][j].Ts = 0.0;
      else cells[i][j].Ts = cells[i][j].Tair;

    }
  }

}

void get_rain_fraction(celltype **cells,meshtype mesh,mproptype mprop)
{
	
  /*Find the fraction of the precipitation for each cell - Maxime*/
  int i,j,k;
  double fac,theta,T;

  int nx = mesh.nx;
  int ny = mesh.ny;
  
  k = 0;

  /*loop cells*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {

      T = cells[i][j].Tair;
	  
      if (T > mprop.Mrate_T[0][k]) {
	while ((T > mprop.Mrate_T[0][k])&&(k < mprop.nMrate_T-1)) k += 1;
      }
      else if (k > 0) {
	while ((T < mprop.Mrate_T[0][k-1])&&(k > 0)) k -= 1;
      }
	
      if ((k == 0)||(k == mprop.nMrate_T-1)) {
	theta = mprop.Mrate_T[3][k];
      }
      else {
	fac = (T - mprop.Mrate_T[0][k-1])/(mprop.Mrate_T[0][k]-mprop.Mrate_T[0][k-1]);
	theta = (1.0-fac)*mprop.Mrate_T[3][k-1] + fac*mprop.Mrate_T[3][k];
      }

      if (cells[i][j].include > 0) {
	  cells[i][j].theta = theta; 
      }
      else {
	cells[i][j].theta = 0.0;
      }

    }
  }
}


int main(int argc,char** argv) {

  char file[200];
  int i,j,nx,ny,nn,ndtw;
  int nnx[100],nny[100];
  double dx,dy,dt,dtw,maxHgw,maxspeed,tslope,diff,H0,Ha,Hq,Hmin,amp;
  FILE *fg;
  FILE *fp = NULL;

  long ncascade,ncatchment,npmax,npa,npia,np;
  long **catchment,**cascade;
  long **slist;
  double meanpw,wfac;
  double **Whelper,*W0helper; 

  long *indx,*pactive,*pinactive;
  double *Tz;
  double **Tcoef;
  double **al; 

  celltype **cells;
  hptype **hp;
  vptype **vp;
  cornertype **cp;
  particletype *pp;
  sortarraytype *sortarray;
  meshtype mesh;
  iproptype iprop;
  mproptype mprop;
  hwproptype hwprop;
  hproptype hprop;
  fproptype fprop;
  pgproptype pgprop; 
  flexproptype flexprop; 
  parproptype parprop;

  double time = 0.0;
  double ftime = 0.0;
  double dt_isos = 0.0; //maxime
  long step = 0;
  long fnr = 0;
  int Nt = 10; 
  int nbCores;
  double Hi = 0.0;
  double dt1 = 0.0; // Glacial
  double dt2 = 0.0; // Fluvial
  int nci = 0;

  /*Read path to working directory*/
  char *path = (char*) malloc(500*sizeof(char));
  path = argv[1];

  printf("******** This is spm-3.4.7  **********\n");

  printf("reading input\n");

  /*read input*/
  printf("  mesh\n");
  mesh = read_mesh(path,file,fp);
  
    /*set number of threads*/
  nbCores = mesh.nbCores;
  printf("Number of threads = %d\n",nbCores);
  omp_set_num_threads(nbCores);

  /*read iprop*/
  printf("  iprop\n");
  iprop = read_iprop(path,file,fp);

  /*read mprop*/
  printf("  mprop\n");
  mprop = read_mprop(path,file,fp);
 
  /*read hwprop*/
  printf("  hwprop\n");
  hwprop = read_hwprop(path,file,fp);

  /*read hprop*/
  printf("  hprop\n");
  hprop = read_hprop(path,file,fp);

  /*read fprop*/
  printf("  fprop\n");
  fprop = read_fprop(path,file,fp);

  /*read flexprop*/
  printf("  flexprop\n");
  flexprop = read_flexprop(path,file,fp);
 
  /*read parprop*/
  printf("  parprop\n");
  parprop = read_parprop(path,file,fp);

  /*read pgprop*/
  printf("  pgprop\n");
  pgprop = read_pgprop(path,file,fp);

  /*compute and extract*/
  mesh.dx = mesh.L/((double)mesh.nx);
  mesh.dy = mesh.H/((double)mesh.ny);
  nx = mesh.nx;
  ny = mesh.ny; 
  nn = nx*ny;
  dx = mesh.dx;
  dy = mesh.dy;
  dt = mesh.maxdt;

  /*Max number of particles*/
  npmax = parprop.npmax;
  np = -1; /*initialize*/
  npa = -1;
  npia = -1;

  
  /*allocate memory*/
  cells = (celltype **) malloc(ny*sizeof(celltype*)); for (i=0;i<ny;i++) cells[i] = (celltype*) malloc(nx*sizeof(celltype));
  hp = (hptype **) malloc(ny*sizeof(hptype*)); for (i=0;i<ny;i++) hp[i] = (hptype*) malloc((nx+1)*sizeof(hptype));
  vp = (vptype **) malloc((ny+1)*sizeof(vptype*)); for (i=0;i<ny+1;i++) vp[i] = (vptype*) malloc(nx*sizeof(vptype));
  cp = (cornertype **) malloc((ny+1)*sizeof(cornertype*)); for (i=0;i<ny+1;i++) cp[i] = (cornertype*) malloc((nx+1)*sizeof(cornertype));
  cascade = (long **) malloc(ny*nx*sizeof(long*)); for (i=0;i<ny*nx;i++) cascade[i] = (long*) malloc(2*sizeof(long));
  catchment = (long **) malloc(ny*nx*sizeof(long*)); for (i=0;i<ny*nx;i++) catchment[i] = (long*) malloc(2*sizeof(long));
  slist = (long **) malloc(ny*nx*sizeof(long*)); for (i=0;i<ny*nx;i++) slist[i] = (long*) malloc(2*sizeof(long));
  Whelper = (double **) malloc((ny+1)*sizeof(double*)); for (i=0;i<ny+1;i++) Whelper[i] = (double*) malloc((nx+1)*sizeof(double));
  W0helper = (double *) malloc(nx*sizeof(double));
  sortarray = (sortarraytype *) malloc(nn*sizeof(sortarraytype));
  pp = (particletype *) malloc(npmax*sizeof(particletype));
  pactive = (long *) malloc(npmax*sizeof(long));
  pinactive = (long *) malloc(npmax*sizeof(long));


  indx = (long *) malloc(Nt*sizeof(long));
  Tz = (double *) malloc(Nt*sizeof(double));
  Tcoef = (double **) malloc(Nt*sizeof(double*)); for (i=0;i<Nt;i++) Tcoef[i] = (double*) malloc(5*sizeof(double));
  al = (double **) malloc(Nt*sizeof(double*)); for (i=0;i<Nt;i++) al[i] = (double*) malloc(2*sizeof(double));

 
  printf("Initializing...");
  initialize(cells,hp,vp,cp,&mesh,hwprop,path,file,fp); 
  printf("done\n");

  /*make neighbor bookkeeping - used e.g. by avalanching*/
  initiate_drainage_info(cells,mesh);

  printf("Writing output #%d...",0);
  write_output(cells,hp,vp,mesh,path,file,fp,0);
  printf("done\n");

  /*report to status*/
  sprintf(file,"%s/status.dat",path); /*filename*/
  if ((fp = fopen(file,"w")) == NULL) {
    printf("Cannot open file : status.dat.\n");
  }
  else {
    fprintf(fp,"%2.4e %2.2f %ld\n",time,100.0*time/mesh.maxtime,fnr);
    fclose(fp);
  }

  
  /*create output file for time series*/
  sprintf(file,"%s/output/tseries.dat",path);
  if ((fp = fopen(file,"wt")) == NULL) {
    printf("could not open file for writing\n");
    exit(1);
  }
  fclose(fp);
  
  /* count cells in active area */
  for (i=0;i<ny;i++) {
	    for (j=0;j<nx;j++) {
  	if (cells[i][j].include > 0) {
		nci = nci + 1;
	}
		}
  }
  mesh.nci = nci;
  printf("threshold_sp = %d\n",fprop.t_sp);
  printf("nci = %d\n",mesh.nci);
  // printf("output times = %f\n",(mesh).filetime);

  printf("Enter time loop\n");
  while (time <= mesh.maxtime)
    {
	  
      /*climate*/
      get_sealevel_temperature(&mprop,time); 
	
      /*ice surface and bed temperatures*/
      get_surface_temperature(mesh,mprop,cells);
	
      /*compute gradients*/
      get_gradients(cells,hp,vp,&mesh);
	  
	  /*Get fraction of precipitation - Maxime*/
	  get_rain_fraction(cells,mesh,mprop); //if (mesh.doice > 0)
		  
	  /*hillslope sediment production*/
      if (mesh.doweathering > 0)
	weathering(cells,&mesh,hprop,dt);

      /*hillslope sediment production*/
      if (mesh.dohillslopeerosion > 0) {
	 hillslope_production(cells,hp,vp,&mesh,hprop,dt,time); /*Maxime 'time'*/
	 stableslope(cells,hp,vp,&mesh,hprop);
	  }

      /*hillslope sediment transport*/
      if ((mesh.dohillslope > 0)&&(mesh.dosediment > 0))
	hillslope_transport(cells,hp,vp,&mesh,hprop,dt);
	
	/*bedrock landslides*/
      if (mesh.dolandslide > 0) 
	landslide(cells,&mesh,hprop,slist,dt,time,path); //Maxime *time

      /*periglacial processes*/
      if (mesh.doperiglacial > 0)
	periglacial(cells,hp,vp,&mesh,mprop,pgprop,dt);
     
      /*fluvial transport and erosion*/
      if (mesh.dofluvial > 0) {

	/*update drainage info*/
	if ((step%10 == 0)||(step == 0))  {
	  get_dnetwork(cells,&mesh);
	  get_cascade(cells,mesh,cascade,catchment,&ncascade,&ncatchment); 
	  get_catchment(ncatchment,catchment,cells);
	  pass_water(cells,&mesh,fprop,ncascade,cascade,dt);
	}
	
	/*perform fluvial transport and erosion - Maxime*/
	if (mesh.Flaw == 1){
	fluvial_transport3(cells,&mesh,fprop,ncascade,cascade,dt);
	}
	else if (mesh.Flaw == 2) {
		fluvial_transport2(cells,&mesh,fprop,ncascade,cascade,dt);
	}
	else if (mesh.Flaw == 3) {
		fluvial_transport(cells,&mesh,fprop,ncascade,cascade,dt);
	}
	else if (mesh.Flaw == 4) {
		fluvial_transport4(cells,&mesh,fprop,ncascade,cascade,dt);
	}
	  }
      /*glacial component*/
      if (mesh.doice > 0) {

	/*update margin index*/
	get_margin(cells,&mesh);
	
	/*Compute mass balance*/
	accumulation_and_melt(cells,mesh,mprop);
	if (mesh.doavalance > 0) avalance(cells,mesh,mprop,sortarray,dt); //MAxime dt

	mass_balance(cells,&mesh,iprop,mprop,hwprop,time,dt);
	
	/*update ice temperatures*/
	if ((step%10 == 0)||(step == 0)) get_thermal(cells,mesh,iprop,mprop,Nt,indx,Tz,Tcoef,al);
	
	/*glacial hydrology*/
	if (mesh.doglacialhydrology > 0) {

	  /*Compute mean annual discharge of water*/
	  get_discharge(cells,hp,vp,&mesh,dt);
	
	  /*compute effective pressure contact ratio etc*/
	  glacial_hydrology(cells,hp,vp,&mesh,iprop,hwprop,dt);

	}
	else { /*no hydrology*/
	  meanpw = 0.0;
	  for (i=0;i<ny;i++) {
	    for (j=0;j<nx;j++) {
	      cells[i][j].Pw = 0.8*cells[i][j].ice;
	      cells[i][j].te = 0.2*cells[i][j].ice;
	      meanpw += cells[i][j].Pw;
	    }
	  }
	  mesh.meanpw = meanpw/((double)mesh.nc);

	}

	/*Compute ice velocities*/
	/*Compute velocities*/
	 if (mesh.gmode > 0) isosia(cells,hp,vp,cp,&mesh,iprop);
	 else sia(cells,hp,vp,cp,&mesh,iprop);

	/*Compute time step length*/
	dt = mesh.ct*mesh.hmin/(mesh.maxspeed+1.0); /*advection*/
	if (dt > mesh.maxdt) dt = mesh.maxdt; /*diffusion*/

	/*Update ice thicknesses*/      
      	/*get_change_MUSCL(cells,hp,vp,&mesh,dt);*/
      	get_change(cells,hp,vp,&mesh,dt);
	
	/*perform glacial erosion*/
	if (mesh.doglacialerosion > 0) glacial_erosion(cells,hp,vp,iprop,&mesh,dt,time); /*Maxime 'time'*/

    /*transport debris in ice*/
	if (mesh.doglacialsedi > 0) glacial_sediment_transport_old(cells,hp,vp,mesh,iprop,dt);

      }
	  if (mesh.douniformerosion > 0) {
			if (mesh.dopulse >0){
				for (i=0;i<ny;i++) {/*Maxime*/
					for (j=0;j<nx;j++) {
						if (cells[i][j].include >0) {
						if (cells[i][j].ice > 10) {
							cells[i][j].quarrying_rate = 1e-4;
							cells[i][j].quarrying += 1e-4;
						}
						else {
							cells[i][j].hillslope_erosion += 1e-4;//ero;
							cells[i][j].hillslope_erate = 1e-4;//ero/dt;
							
						}
							cells[i][j].sedi += 1e-4;
							
						}
					}
				}
			}
	  	}

      /*update lagrangian sediment particles*/
      if (mesh.doparticles > 0)
	update_particles(cells,hp,vp,mesh,parprop,hprop,pp,npmax,&np,&npa,&npia,pactive,pinactive,time,dt);

      /*flexural isostasy*/
      if ((mesh.doisostasy > 0)&&(step%flexprop.ncall == 0))
		  flexural_isostasy(cells,&mesh,iprop,flexprop,Whelper,dt_isos); //Maxime ftime
		  dt_isos = 0.0;
	

	/*boundary conditions - ad hoc*/      
      for (i=0;i<ny;i++) {
	//cells[i][0].bed = 0.99*cells[i][1].bed;
	cells[i][0].ice = 0.0;//cells[i][1].ice;
	cells[i][0].deformation = 0.0;
	cells[i][0].sliding = 0.0;
	// cells[i][0].bed = cells[i][1].bed;
	// cells[i][nx-1].bed = cells[i][nx-2].bed;
	//cells[i][nx-1].bed = 0.99*cells[i][nx-2].bed;
	cells[i][nx-1].ice = 0.0;//cells[i][nx-2].ice;
	cells[i][nx-1].deformation = 0.0;
	cells[i][nx-1].sliding = 0.0;
	if (cells[i][0].ice < 0.0) cells[i][0].ice = 0.0;
	if (cells[i][nx-1].ice < 0.0) cells[i][nx-1].ice = 0.0;
	  }
	for (j=1;j<nx-1;j++) {
	    //cells[0][j].bed = 0.99*cells[1][j].bed;
		cells[0][j].ice = 0.0;//cells[1][j].ice;
		cells[0][j].deformation = 0.0;
		cells[0][j].sliding = 0.0;
		// cells[0][j].bed = cells[1][j].bed;
		// cells[ny-1][j].bed = cells[ny-2][j].bed;
		//cells[ny-1][j].bed = 0.99*cells[ny-2][j].bed;
		cells[ny-1][j].ice = 0.0;//cells[ny-2][j].ice;
		cells[ny-1][j].deformation = 0.0;
		cells[ny-1][j].sliding = 0.0;
		if (cells[0][j].ice < 0.0) cells[0][j].ice = 0.0;
		if (cells[ny-1][j].ice < 0.0) cells[ny-1][j].ice = 0.0;
	  }
	  /*if ((cells[i][j].margin == 1)&&(fabs(cells[i][j].x-2870.0)<1000.0)&&(fabs(cells[i][j].y-3100.0)<1000.0)) cells[i][j].Vs[0] = 0.1;*/
	  
	  /*Uplift */
	  for (i=0;i<ny;i++) {
	    for (j=0;j<nx;j++) {
			if ((cells[i][j].fixflag <= 0) && (cells[i][j].include > 0)) cells[i][j].bed += cells[i][j].srate*dt;
			else cells[i][j].bed = cells[i][j].bedStore;
		}
	  }

	// /* handling of sediments - ad hoc */
	// for (i=0;i<ny;i++) {
	    // for (j=0;j<nx;j++) {
			// if (cells[i][j].sedi > 50.0) cells[i][j].sedi = 50.0;
		// }
	  // }
	  
	/*boundary conditions - track outflux ad hoc*/ 
	  double sed_right = 0.0; //Maxime
	  double sed_top = 0.0; //Maxime
	  double sed_left = 0.0; //Maxime
	  double sed_bottom = 0.0; //Maxime
      for (i=0;i<ny;i++) {
		sed_left += cells[i][0].sedi * (dx*dy); 
		sed_left += cells[i][0].streamload; 
		cells[i][0].sedi = 0.0; 
		sed_right += cells[i][nx-1].sedi * (dx*dy);
		sed_right += cells[i][nx-1].streamload;
		cells[i][nx-1].sedi = 0.0; 
      }
	  for (j=0;j<nx;j++) {
		sed_bottom += cells[0][j].sedi * (dx*dy);
		sed_bottom += cells[0][j].streamload;
		cells[0][j].sedi = 0.0;
		sed_top += cells[ny-1][j].sedi * (dx*dy);	
		sed_top += cells[ny-1][j].streamload;		
		cells[ny-1][j].sedi = 0.0; 
      }
	  double sed_tot = 0.0; //MAxime
	  /* Remove sediments in piedmont area */
	  for (i=0;i<ny;i++) {
			for (j=0;j<nx;j++) {
				if (cells[i][j].include <= 0) {
					sed_tot += (cells[i][j].sedi + cells[i][j].streamload);
					cells[i][j].sedi = 0.0;
					cells[i][j].streamload = 0.0;
				}
			}
	  }
	/*Check for no erosion flag*/
	if (mesh.steadybed > 0) {
		for (i=0;i<ny;i++) {
			for (j=0;j<nx;j++) {
				cells[i][j].bed = cells[i][j].bedStore;
			}
		}
	}

      step += 1;
      time += dt;
      ftime += dt;
	  dt_isos += dt; //Maxime
      
      /*output data file*/
      if (ftime >= mesh.filetime[0][fnr]) {
      /*		if (step == 1) {*/

	fnr += 1;
	fflush(stdout);
	printf("  Writing output #%ld, time = %lf kyr\n",fnr,mesh.filetime[0][fnr-1]/1e3);
	write_output(cells,hp,vp,mesh,path,file,fp,fnr);
	
	if (mesh.doparticles > 0) write_particles(pp,npa,pactive,path,file,fp,fnr);
		  
	/*printf("done\n");*/
	/*ftime = 0.0;

	/*report to status*/
	sprintf(file,"%s/status.dat",path);
	if ((fp = fopen(file,"w")) == NULL) {
	  printf("Cannot open file : status.dat.\n");
	}
	else {
	  fprintf(fp,"%2.4e %2.2f %ld\n",time,100.0*time/mesh.maxtime,fnr);
	  fclose(fp);
	}
     
      }

      if (step%mesh.nmoni == 0) {
	sprintf(file,"%s/output/tseries.dat",path);
	if ((fg = fopen(file,"at")) == NULL) {
	  printf("could not open file for writing\n");
	  exit(1);
	}
	fprintf(fg,"%4.4e ",time);
	fprintf(fg,"%4.4e ",dt);
	fprintf(fg,"%4.4e ",mesh.meanice);
	fprintf(fg,"%ld ",mesh.veloitt/mesh.veloitt_count);
	fprintf(fg,"%4.4e ",mesh.maxspeed);
	fprintf(fg,"%4.4e ",mesh.mean_periglacial_erate);
	fprintf(fg,"%4.4e ",mesh.mean_dHs_periglacial);
	fprintf(fg,"%4.4e ",mesh.mean_dHs_hillslope);
	fprintf(fg,"%4.4e ",mesh.mean_quarrying_rate);
	fprintf(fg,"%4.4e ",mprop.T0);
	fprintf(fg,"%4.4e ",dtw); 
	fprintf(fg,"%4.4e ",mesh.meanice_b);
	fprintf(fg,"%ld ",mesh.stressitt/mesh.stressitt_count);
	fprintf(fg,"%4.4e ",mesh.mean_dHs_fluvial);
	fprintf(fg,"%ld ",mesh.nice);
	fprintf(fg,"%4.4e ",mesh.meanbed);
	fprintf(fg,"%4.4e ",mesh.minbed);
	fprintf(fg,"%4.4e ",mesh.maxbed);
	fprintf(fg,"%4.4e ",mesh.mean_hillslope_erate);
	fprintf(fg,"%4.4e ",mesh.mean_landslide_erate);
	fprintf(fg,"%4.4e ",mesh.mean_abrasion_rate);
	fprintf(fg,"%4.4e ",mesh.mean_isostasy);
	fprintf(fg,"%4.4e ",mesh.meansedi);
	fprintf(fg,"%ld ",mesh.hydroitt);
	fprintf(fg,"%4.4e ",mesh.meante);
	fprintf(fg,"%4.4e ",mesh.meanhw);
	fprintf(fg,"%4.4e ",mesh.meandhwdt);
	fprintf(fg,"%4.4e ",mesh.mean_weatheringrate);
	fprintf(fg,"%ld ",np);
	fprintf(fg,"%ld ",npa);
	fprintf(fg,"%ld ",npia);
	fprintf(fg,"%4.4e ",H0);
	fprintf(fg,"%4.4e ",sed_tot);
	fprintf(fg,"%4.4e ",Hi);
	fprintf(fg,"%4.4e ",mesh.isostasy_rate);
	fprintf(fg,"%4.4e ",mesh.mean_fluvial_erate);
	fprintf(fg,"%4.4e ",sed_left);
	fprintf(fg,"%4.4e ",sed_right);
	fprintf(fg,"%4.4e ",sed_bottom);
	fprintf(fg,"%4.4e ",sed_top);
	
	fprintf(fg,"\n");
	fclose(fg);
	 if (mesh.veloitt/mesh.veloitt_count > 300) {
		write_output(cells,hp,vp,mesh,path,file,fp,fnr+500);
		//exit(1);
	}
	mesh.veloitt = 0;
	mesh.veloitt_count = 1;
	mesh.stressitt = 0;
	mesh.stressitt_count = 1;	
   

	/*output celldata*/
	if (mesh.docelldata > 0) {
	  sprintf(file,"%s/output/celldata.dat",path);
	  if ((fg = fopen(file,"at")) == NULL) {
	    printf("could not open file for writing\n");
	    exit(1);
	  }
	  fprintf(fg,"%4.4e ",time);
	  for (i=0;i<mesh.ncelldata;i++) {
	    fprintf(fg,"%4.4e ",cells[mesh.celldata_i[i]][mesh.celldata_j[i]].bed);
	    fprintf(fg,"%4.4e ",cells[mesh.celldata_i[i]][mesh.celldata_j[i]].sedi);
	    fprintf(fg,"%4.4e ",cells[mesh.celldata_i[i]][mesh.celldata_j[i]].Tair);
	    fprintf(fg,"%4.4e ",cells[mesh.celldata_i[i]][mesh.celldata_j[i]].periglacial_erate);
	    fprintf(fg,"%4.4e ",cells[mesh.celldata_i[i]][mesh.celldata_j[i]].periglacial_erosion);
	  }
	  fprintf(fg,"\n");
	  fclose(fg);
	}
	
      }
      

    }/*time*/

  for (i=0;i<ny;i++) free(cells[i]); free(cells);
  for (i=0;i<ny;i++) free(hp[i]); free(hp);
  for (i=0;i<ny+1;i++) free(vp[i]); free(vp);
  for (i=0;i<ny+1;i++) free(cp[i]); free(cp);
  for (i=0;i<ny*nx;i++) free(cascade[i]); free(cascade);
  for (i=0;i<ny*nx;i++) free(catchment[i]); free(catchment);
  for (i=0;i<ny*nx;i++) free(slist[i]); free(slist);
  for (i=0;i<ny+1;i++) free(Whelper[i]); free(Whelper);

  exit(1);


}

