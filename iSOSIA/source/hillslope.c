

void stableslope(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,hproptype hprop)
{

  int i,j;

  double sc = 0.7;

  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int dosediment = (*mesh).dosediment; //Maxime
  double dhx = dx*sc;
  double dhy = dy*sc;
  double dh;

  for (i=0;i<ny;i++) {
    for (j=1;j<nx;j++) {
      if ((cells[i][j].bed-cells[i][j-1].bed) > dhx) {
	dh = cells[i][j].bed - cells[i][j-1].bed - dhx;
	cells[i][j].bed -= dh;
	cells[i][j].landslide_erosion += dh;
	//if (dosediment > 1) cells[i][j-1].sedi += dh;
      }
      else if ((cells[i][j-1].bed-cells[i][j].bed) > dhx) {
	dh = cells[i][j-1].bed - cells[i][j].bed - dhx;
	cells[i][j-1].bed -= dh;
	cells[i][j-1].landslide_erosion += dh;
	//if (dosediment > 1) cells[i][j].sedi += dh;
      }
    }
  }


  for (i=1;i<ny;i++) {
    for (j=0;j<nx;j++) {
      if ((cells[i][j].bed-cells[i-1][j].bed) > dhy) {
	dh = cells[i][j].bed - cells[i-1][j].bed - dhy;
	cells[i][j].bed -= dh;
	cells[i][j].landslide_erosion += dh;
	//if (dosediment > 1) cells[i-1][j].sedi += dh;
      }
      else if ((cells[i-1][j].bed-cells[i][j].bed) > dhy) {
	dh = cells[i-1][j].bed - cells[i][j].bed - dhy;
	cells[i-1][j].bed -= dh;
	cells[i-1][j].landslide_erosion += dh;
	//if (dosediment > 1) cells[i][j].sedi += dh;
      }
    }
  }
}

  
void hillslope_production(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,hproptype hprop,double dt,double time)
{
  /*See Roering et al. (1999)*/
  int i,j;
  double sdiff,fac,ero,maxero;
  double meanerate = 0.0;

  double sc;
  double Ke = hprop.Ke;

  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int dosediment = (*mesh).dosediment;
  double minice = 10.0;
  double nci = (*mesh).nci;

#pragma omp parallel shared(cells,hp,vp,meanerate) private(i,j,sdiff,fac,maxero,ero) firstprivate(nx,ny,dx,dy,dt,dosediment,minice,sc,Ke)
  {

    /***************** Erosion (critical slope) *****************/
      /*loop h-points for bedrock erosion*/
#pragma omp for schedule(static)
      for (i=0;i<ny;i++) {
	for (j=1;j<nx;j++) { 
	  if (cells[i][j-1].include+cells[i][j].include > 1.0) {

	    sc = 0.5*(cells[i][j-1].phi+cells[i][j].phi);
		
	    fac = 1.0 - pow(hp[i][j].dbdx/sc,2.0); if (fac < 1.0e-4) fac = 1.0e-4;
	    sdiff = Ke/fac;
	    ero = sdiff*hp[i][j].dbdx*dt/dx*sqrt(1.0+pow(cells[i][j].bslope,2.0));
#pragma omp critical
	    {
	      if (ero < 0.0) {
		maxero = cells[i][j-1].bed-cells[i][j].bed-sc*dx;
		if (maxero < 0.0) maxero = 0.0;
		if (ero < -maxero) ero = -maxero;
		if (cells[i][j].ice > 10) ero = 0.0;
		cells[i][j-1].bed += ero;
		cells[i][j-1].hillslope_erosion -= ero;
		cells[i][j-1].hillslope_erate = -ero/dt;
		if (dosediment > 0) cells[i][j-1].sedi -= ero;
	      }
	      else {
		maxero = cells[i][j].bed-cells[i][j-1].bed-sc*dx;
		if (maxero < 0.0) maxero = 0.0;
		if (ero > maxero) ero = maxero;
		if (cells[i][j-1].ice > 10) ero = 0.0;
		cells[i][j].bed -= ero;
		cells[i][j].hillslope_erosion += ero;
		cells[i][j].hillslope_erate = ero/dt;
		if (dosediment > 0) cells[i][j].sedi += ero;
	      }
	    } 
	  }
	}
      }

    /*loop v-points*/
#pragma omp for schedule(static)
      for (j=0;j<nx;j++) {
	for (i=1;i<ny;i++) {
	  if (cells[i-1][j].include+cells[i][j].include > 1.0) {

	    sc = 0.5*(cells[i-1][j].phi+cells[i][j].phi);

	    fac = 1.0 - pow(vp[i][j].dbdy/sc,2.0); if (fac < 1.0e-4) fac = 1.0e-4;
	    sdiff = Ke/fac;
	    ero = sdiff*vp[i][j].dbdy*dt/dy*sqrt(1.0+pow(cells[i][j].bslope,2.0));
#pragma omp critical
	    {
	      if (ero < 0.0) {
		maxero = cells[i-1][j].bed-cells[i][j].bed-sc*dy;
		if (maxero < 0.0) maxero = 0.0;
		if (ero < -maxero) ero = -maxero;
		if (cells[i][j].ice > 10) ero = 0.0;
		cells[i-1][j].bed += ero;
		cells[i-1][j].hillslope_erosion -= ero;
		cells[i-1][j].hillslope_erate = -ero/dt;
		if (dosediment > 0) cells[i-1][j].sedi -= ero;
	      }
	      else {
		maxero = cells[i][j].bed-cells[i-1][j].bed-sc*dy;
		if (maxero < 0.0) maxero = 0.0;
		if (ero > maxero) ero = maxero;
		if (cells[i-1][j].ice > 10) ero = 0.0;
		cells[i][j].bed -= ero; 
	    cells[i][j].hillslope_erosion += ero;
		cells[i][j].hillslope_erate = ero/dt;
		if (dosediment > 0) cells[i][j].sedi += ero;
	      }
	    } 
	  }
	}
      }
    

  /*compute mean erosion rate and transport change*/
#pragma omp for schedule(static) reduction(+:meanerate)
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
		if (cells[i][j].include > 0.0) {
			meanerate += cells[i][j].hillslope_erate;
		}
    }
  }
#pragma omp single
  {
    meanerate /= (double)(nci);
    (*mesh).mean_hillslope_erate = meanerate;

  }

  }

}

void hillslope_transport(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,hproptype hprop,double dt)
{

  int i,j;
  double sdiff,dHs,mean_dHs=0.0,fac;

  double Ks = hprop.Ks;
  double sc = hprop.sc;

  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double minice = 10.0;


  /*********************** sediment transport *********************/

/*"Private" ensures that each thread as its own copy of the variable
  "Shared" means that variables are visible and accessible by all threads simultaneously
  */
#pragma omp parallel shared(cells,hp,vp,mean_dHs) private(i,j,sdiff,dHs,fac) firstprivate(nx,ny,dx,dy,dt,Ks,sc,minice)
  {


    /*initialize*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].dHs = 0.0;
      }
    }
    
    /*loop h-points for horizontal transport*/
#pragma omp for schedule(static) //Each thread is assigned a chunk of iterations in fixed fashion
      for (i=0;i<ny;i++) {
	for (j=1;j<nx;j++) {
		if (cells[i][j-1].include+cells[i][j].include > 1.0) {
			if (cells[i][j].ice < 20.0){ //Maxime
	  fac = 1.0 - pow(hp[i][j].dtdx/sc,2.0); if (fac < 0.001) fac = 0.001;
	  sdiff = Ks/fac;
	  dHs = -sdiff*hp[i][j].dtdx*dt/dx;
	  if (dHs < 0.0) {
		if (cells[i][j].ice > 20.0) dHs = 0.0; //Maxime - No erosion if ice
	    else if (cells[i][j].sedi <= 0.0) dHs = 0.0;
	    else if (dHs < -cells[i][j].sedi) dHs = -cells[i][j].sedi;
	  }
	  else {
		if (cells[i][j-1].ice > 20.0) dHs = 0.0; //Maxime - No erosion if ice
	    else if (cells[i][j-1].sedi <= 0.0) dHs = 0.0;
	    else if (dHs > cells[i][j-1].sedi) dHs = cells[i][j-1].sedi;
	  }
#pragma omp critical //Ensure 1 thread enters the block at a time
	  {
	    cells[i][j].dHs += dHs;
	    cells[i][j-1].dHs -= dHs;
		// if (cells[i][j].ice > 20.0) cells[i][j].Vs[0] += dHs; //Deposit on ice - Maxime
		// else cells[i][j].sedi += dHs; //Maxime
		// if (cells[i][j-1].ice > 20.0) cells[i][j].Vs[0] += dHs; //Deposit on ice - Maxime
		// else cells[i][j-1].sedi -= dHs;
	   }
			}
		}
	}
      }

      /*loop v-points*/
#pragma omp for schedule(static)
      for (j=0;j<nx;j++) {
	for (i=1;i<ny;i++) {
		if (cells[i-1][j].include+cells[i][j].include > 1.0) {
		if (cells[i][j].ice < 20.0){ //Maxime
	  fac = 1.0 - pow(vp[i][j].dtdy/sc,2.0); if (fac < 0.001) fac = 0.001;
	  sdiff = Ks/fac;
	  dHs = -sdiff*vp[i][j].dtdy*dt/dy;
	  if (dHs < 0.0) {
		if (cells[i][j].ice > 20.0) dHs = 0.0; //Maxime - No erosion if ice
	    else if (cells[i][j].sedi <= 0.0) dHs = 0.0; 
	    else if (dHs < -cells[i][j].sedi) dHs = -cells[i][j].sedi;
	  }
	  else {
		if (cells[i-1][j].ice > 20.0) dHs = 0.0; //Maxime - No erosion if ice
	    else if (cells[i-1][j].sedi <= 0.0) dHs = 0.0;
	    else if (dHs > cells[i-1][j].sedi) dHs = cells[i-1][j].sedi;
	  }
#pragma omp critical
	  {
	    cells[i][j].dHs += dHs;
	    cells[i-1][j].dHs -= dHs;
		// if (cells[i][j].ice > 20.0) cells[i][j].Vs[0] += dHs; //Deposit on ice - Maxime
		// else cells[i][j].sedi += dHs; //Maxime
		// if (cells[i-1][j].ice > 20.0) cells[i][j].Vs[0] += dHs; //Deposit on ice - Maxime
		// else cells[i-1][j].sedi -= dHs;
	  }
		}
		}
	}
      }
      
    
#pragma omp for schedule(static)
      for (i=0;i<ny;i++) {
	for (j=0;j<nx;j++) {
		if (cells[i][j].ice > 20){
			if (cells[i][j].dHs < 0.0) cells[i][j].dHs = 0.0;
			/* ad hoc - remove sediments coming from hillslope if overwhelm uppermost layer of ice*/
			if (cells[i][j].Vs[0] + cells[i][j].dHs > 0.95*(cells[i][j].ice/20)){
				cells[i][j].Vs[0] += (0.95*(cells[i][j].ice/20) - cells[i][j].Vs[0]) ;//printf("dhs = %f, Vs = %f, ice = %f",cells[i][j].dHs,cells[i][j].Vs[0],cells[i][j].ice);
			}
			cells[i][j].Vs[0] += cells[i][j].dHs; //Sediment fall onto ice surface*
			
		}
		else {
			if (cells[i][j].dHs < -cells[i][j].sedi) cells[i][j].dHs = -cells[i][j].sedi;
			cells[i][j].sedi += cells[i][j].dHs;
		}
	}
      }
    
  /*compute mean erosion rate and transport change*/
#pragma omp for schedule(static) reduction(+:mean_dHs) 
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      mean_dHs += fabs(cells[i][j].dHs);
    }
  }
#pragma omp single
  {
    mean_dHs /= (double)((*mesh).nc)*dt*2.0;
    (*mesh).mean_dHs_hillslope = mean_dHs;
  }


  }

}


void nslide(int ti,int tj,celltype **cells,long *nsc,long **slist,double x0,double y0,double h0,double *maxH,double *maxbeta)
{

  long i;
  int ni,nj;
  double hp,dist,Hc,beta; 
  double xc = cells[ti][tj].x;
  double yc = cells[ti][tj].y;
  double hc = cells[ti][tj].topsedi;

  /*printf("ti = %d, tj = %d, nne = %d\n",ti,tj,cells[ti][tj].nne);*/

  /*loop neighbours*/
  for (i=0;i<cells[ti][tj].nne;i++) {

    /*identify neighbour*/
    ni = cells[ti][tj].ne_i[i];
    nj = cells[ti][tj].ne_j[i];

    /*if not allready in slide*/
    if (cells[ni][nj].inslide == 0) {
    
      /*distance to neigbour*/
      dist = sqrt(pow(cells[ni][nj].x-xc,2.0)+pow(cells[ni][nj].y-yc,2.0));

      /*Critical hillslope height*/
      Hc = hc + cells[ni][nj].phi*dist; //Maxime
      /*Hc = hc + (double)tan(phi0)*dist;*/

      /*Height of cell*/
      hp = cells[ni][nj].topsedi;

      /*if potentially unstable*/
      if (hp > Hc) {

	/*registre cell*/
	slist[(*nsc)][0] = ni;
	slist[(*nsc)][1] = nj;
	cells[ni][nj].inslide = 1;
	cells[ni][nj].rfail = hp-Hc;

	/*add to number*/
	(*nsc) += 1;

	/*increase hillslope height*/
	if ((*maxH) < hp) (*maxH) = hp;

	/*distance to toe*/
	dist = sqrt(pow(cells[ni][nj].x-x0,2.0)+pow(cells[ni][nj].y-y0,2.0));

	/*increase hillslope angle*/
	beta = (double)atan((hp-h0)/dist);
	if ((*maxbeta) < beta) (*maxbeta) = beta;

	/*recursive call*/
	nslide(ni,nj,cells,nsc,slist,x0,y0,h0,maxH,maxbeta);
	
      }/*if*/

      /*printf("\n");*/

    }/*if*/

  }/*i*/

}




void landslide(celltype **cells,meshtype *mesh,hproptype hprop,long **slist,double dt,double time,char *path)
{

  /* See Egholm et al. (2013), Densmore (1998), Schmidt and Montgomery (1995) */
  char file[200];
  long i,j,k;
  int ti,tj,ni,nj;
  long nsc,nf;
  double HH,Hc,dH,dHmax,maxbeta,theta,phi,Pf,lfrac,ran,mdist;
  double maxH,x0,y0,h0;
  double Atot,Vtot,jgift,newsedi,Vmean;
  double ero;

  double gamma = hprop.gamma;
  int Nc = hprop.Nc;
  double kt = hprop.kt;
  FILE *fg;

  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double area = dx*dy;

  int dosediment = (*mesh).dosediment;

  /*initialize*/
  for (i=0;i<ny;i++) 
    for (j=0;j<nx;j++) {
      cells[i][j].inslide = 0;
      cells[i][j].rfail = 0.0;
    }
  
  /*mean landslide volume*/
  Vmean = 0.0;
  newsedi = 0.0;
  
  /*First check all cells to handle over-steepened slope > 45°*/
  for (i=0;i<ny;i++){ 
    for (j=0;j<nx;j++) {
		if (cells[i][j].bslope > 0.8) {
			nsc = 0;
			x0 = cells[i][j].x;
			y0 = cells[i][j].y;
			h0 = cells[i][j].topsedi;
			maxH = h0;
			maxbeta = 0.0;
			phi = (double)atan(cells[i][j].phi); //rad - Maxime

			/*recursive call*/
			nslide(i,j,cells,&nsc,slist,x0,y0,h0,&maxH,&maxbeta);

			/*if potential for slide*/
			if (nsc > 0) {
	
			/*Hillslope height*/
			HH = maxH - h0;

			/*Critical height*/
			Hc = gamma*sin(maxbeta)*cos(phi)/(1.0-cos(maxbeta-phi));

			/*Propability*/
			/*Pf = HH/Hc;*/
			Pf = 1;
			cells[i][j].lastslide = (double)time;
	
			/*pick random number between 0 and 1*/
			ran = (double)rand()/(double)RAND_MAX;

			/*if failure*/
			if (ran < Pf) {

				Atot = 0.0;
				Vtot = 0.0;

				/*loop unstable cells and collect material*/
				for (k=0;k<nsc;k++) {
	    
					/*identify cell*/
					ni = slist[k][0];
					nj = slist[k][1];

					/*Excess height*/
					dH = cells[ni][nj].rfail;

					/*erosion*/
					if (dH > 0.1) {
	      
					if (dH > cells[ni][nj].sedi) {
						newsedi = dH - cells[ni][nj].sedi;
		
						cells[ni][nj].bed -= newsedi; 
						cells[ni][nj].landslide_erosion += newsedi;
		
						if (dosediment > 0) {
							cells[ni][nj].sedi += newsedi;
						}
					}

					Atot += area;
					Vtot += area*dH;

					}/*if erosion*/

				}/*k*/
			}
			}
		}
	}
  }
			

  /*loop some random cells*/
  for (i=1;i<=Nc;i++) {

    /*random double between 0 and 1*/
    ran = (double)rand()/(double)RAND_MAX;

    /*random integer between 1 and ny*/
    ti = (long)floor(ran*(double)ny);

    /*random double between 0 and 1*/
    ran = (double)rand()/(double)RAND_MAX;

    /*random integer between 1 and nx*/
    tj = (long)floor(ran*(double)nx);

    /*if not allready in slide*/
	/*No landslide under the ice - Maxime */
    if ((ti >=0)&&(ti<ny)&&(tj>=0)&&(tj<nx)&&(cells[ti][tj].inslide == 0)&&cells[ti][tj].ice<10) {

      nsc = 0;
      x0 = cells[ti][tj].x;
      y0 = cells[ti][tj].y;
      h0 = cells[ti][tj].topsedi;
      maxH = h0;
      maxbeta = 0.0;
      phi = (double)atan(cells[ti][tj].phi); //rad - Maxime

      /*recursive call*/
      nslide(ti,tj,cells,&nsc,slist,x0,y0,h0,&maxH,&maxbeta);

      /*if potential for slide*/
      if (nsc > 0) {
	
	/*Hillslope height*/
	HH = maxH - h0;

	/*Critical height*/
	Hc = gamma*sin(maxbeta)*cos(phi)/(1.0-cos(maxbeta-phi));

	/*Propability*/
	/*Pf = HH/Hc;*/
	Pf = HH/Hc + kt*(time - cells[ti][tj].lastslide);
	cells[ti][tj].lastslide = (double)time;
	
	/*pick random number between 0 and 1*/
	ran = (double)rand()/(double)RAND_MAX;

	/*if failure*/
	if (ran < Pf) {

	  Atot = 0.0;
	  Vtot = 0.0;

	  /*loop unstable cells and collect material*/
	  for (j=0;j<nsc;j++) {
	    
	    /*identify cell*/
	    ni = slist[j][0];
	    nj = slist[j][1];

	    /*Excess height*/
	    dH = cells[ni][nj].rfail;

	    /*erosion*/
	    if (dH > 0.1) {
	      
	      if (dH > cells[ni][nj].sedi) {
		newsedi = dH - cells[ni][nj].sedi;
		
		cells[ni][nj].bed -= newsedi; 
		cells[ni][nj].landslide_erosion += newsedi;
		
		if (dosediment > 0) {
		  cells[ni][nj].sedi += newsedi;
		  /*report slide to file*/
	  /*if (cells[ni][nj].sedi <0) {
	  sprintf(file,"%s/output/sedi.dat",path);
	  if ((fg = fopen(file,"at")) == NULL) {
	      printf("yes");
	    }
		else{
	  fprintf(fg, "%2.3f  %2.3f",cells[ni][nj].sedi, newsedi);
	  fclose(fg);
		}*/
		  }
		  }

	      Atot += area;
	      Vtot += area*dH;

	    }/*if erosion*/

	  }/*j*/

	  /*
	  }
	  else {
	    fprintf(fg,"%2.5e %2.5e %2.5e %ld %ld\n",time,Atot,Vtot,tp,nsc);
	    fclose(fg);
	    }
	  */
	  /*report big landslide*/
	  /*
	  if (Vtot > 5.0e8) {

	    sprintf(file,"%s/output/bigslides.dat",path);
	    if ((fg = fopen(file,"at")) == NULL) {
	      nowrite = 1;
	    }
	    else {
	      fprintf(fg,"%2.5e %2.5e %2.5e %ld ",time,Atot,Vtot,nsc);
	      for (j=1;j<=nsc;j++) fprintf(fg,"%d ",slist[j]);
	      fprintf(fg,"\n");
	      fclose(fg);
	    }

	  }*/

	  Vmean += Vtot;


      }/*if failure*/

      }/*if nsc > 0*/

    }/*if not in slide*/

  }/*i*/

  (*mesh).mean_landslide_erate = Vmean/((*mesh).L*(*mesh).H*dt);

}/*landslide*/


void weathering(celltype **cells,meshtype *mesh,hproptype hprop,double dt)
{

  int i,j;

  double Kw = hprop.Kw;
  double Ls = hprop.Ls;
  double ero,meanwrate = 0.0;
  double minice = 10.0;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int dosediment = (*mesh).dosediment;


#pragma omp parallel shared(cells,meanwrate) private(i,j,ero) firstprivate(nx,ny,dt,Kw,Ls,minice,dosediment)
  {

#pragma omp for schedule(static) reduction(+:meanwrate)
      for (i=0;i<ny;i++) {
	for (j=0;j<nx;j++) {
	  if ((cells[i][j].ice < minice)&&(cells[i][j].include > 0)) {
	    ero = Kw*exp(-cells[i][j].sedi/Ls)*dt*sqrt(1.0+pow(cells[i][j].bslope,2.0));
	    cells[i][j].bed -= ero;
	    cells[i][j].weathering += ero;
	    cells[i][j].weathering_rate = ero/dt;
	    if (dosediment > 0) cells[i][j].sedi += ero;
	    meanwrate += ero/dt;
	  }
	}
      }
 



#pragma omp single
  {
    meanwrate /= (double)((*mesh).nc);
    (*mesh).mean_weatheringrate = meanwrate;

  }

  } 

}


