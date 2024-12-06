
double max(double val1,double val2)
{

  double maxval = val1;
  if (val2 > maxval) maxval = val2;

  return(maxval);

}


int compare_values(const void *elem1, const void *elem2)
{

  sortarraytype *i1,*i2;
  i1 = (sortarraytype*)elem1;
  i2 = (sortarraytype*)elem2;
  if (i1->value < i2->value) return 1; 
  else if (i1->value == i2->value) return 0;
  else return -1;

}


void get_discharge(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,double dt)
{
  /*Compute water discharge after massbalance has provided the melt water*/
  int i,j; 
  double vout;

  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int nn = nx*ny;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  double spy = 365.25*24.0*3600.0;
  double minice = 10.0;
  double Tfac = 0.1;  

  /*loop cells in decending order of elevation*/ 
  for (i=0;i<ny;i++) { 
    for (j=0;j<nx;j++) {
 
      if (cells[i][j].ice > minice) {
      
      vout = 0.0; 
      cells[i][j].Qstack += cells[i][j].bwater*dt;
       
      if (j > 0) if ((hp[i][j].vx < 0.0)&&(cells[i][j-1].ice > minice)) vout -= hp[i][j].vx;
      if (j < nx-1) if ((hp[i][j+1].vx > 0.0)&&(cells[i][j+1].ice > minice)) vout += hp[i][j+1].vx;
      if (i > 0) if ((vp[i][j].vy < 0.0)&&(cells[i-1][j].ice > minice)) vout -= vp[i][j].vy;
      if (i < ny-1) if ((vp[i+1][j].vy > 0.0)&&(cells[i+1][j].ice > minice)) vout += vp[i+1][j].vy;
     
      if (vout > 0.0) {
	if (j > 0) if ((hp[i][j].vx < 0.0)&&(cells[i][j-1].ice > minice)) cells[i][j-1].Qstack += cells[i][j].Qstack*(-hp[i][j].vx)/vout;
	if (j < nx-1) if ((hp[i][j+1].vx > 0.0)&&(cells[i][j+1].ice > minice)) cells[i][j+1].Qstack += cells[i][j].Qstack*(hp[i][j+1].vx)/vout;
	if (i > 0) if ((vp[i][j].vy < 0.0)&&(cells[i-1][j].ice > minice)) cells[i-1][j].Qstack += cells[i][j].Qstack*(-vp[i][j].vy)/vout;
	if (i < ny-1) if ((vp[i+1][j].vy > 0.0)&&(cells[i+1][j].ice > minice)) cells[i+1][j].Qstack += cells[i][j].Qstack*(vp[i+1][j].vy)/vout;
      } 
      else {
    
	vout = 0.0;
	if (j > 0) if (cells[i][j-1].topice < cells[i][j].topice) vout += cells[i][j].topice - cells[i][j-1].topice;
	if (j < nx-1) if (cells[i][j+1].topice < cells[i][j].topice) vout += cells[i][j].topice - cells[i][j+1].topice;
	if (i > 0) if (cells[i-1][j].topice < cells[i][j].topice) vout += cells[i][j].topice - cells[i-1][j].topice;
	if (i < ny-1) if (cells[i+1][j].topice < cells[i][j].topice) vout += cells[i][j].topice - cells[i+1][j].topice;
      
	if (vout > 0.0) {
	  if (j > 0) if (cells[i][j-1].topice < cells[i][j].topice) cells[i][j-1].Qstack += cells[i][j].Qstack*(cells[i][j].topice-cells[i][j-1].topice)/vout;
	  if (j < nx-1) if (cells[i][j+1].topice < cells[i][j].topice) cells[i][j+1].Qstack += cells[i][j].Qstack*(cells[i][j].topice-cells[i][j+1].topice)/vout;
	  if (i > 0) if (cells[i-1][j].topice < cells[i][j].topice) cells[i-1][j].Qstack += cells[i][j].Qstack*(cells[i][j].topice-cells[i-1][j].topice)/vout;
	  if (i < ny-1) if (cells[i+1][j].topice < cells[i][j].topice) cells[i+1][j].Qstack += cells[i][j].Qstack*(cells[i][j].topice-cells[i+1][j].topice)/vout;
	}	 										       
	/*else cells[i][j].Qw = -1000.0;*/

      }
 
      cells[i][j].bQw = dx*dy/(spy*dt)*cells[i][j].Qstack;  /*m3 pr. sec - basal discharge - Maxime*/ 
	  cells[i][j].bwatersurf = cells[i][j].bed + cells[i][j].sedi + cells[i][j].Qstack; //Maxime
	  cells[i][j].bwater = cells[i][j].Qstack/dt; //Maxime
      cells[i][j].Qstack = 0.0;
      /*cells[i][j].Qw = 0.1*(cells[i][j].Qw - cells[i][j].Qw_old) + cells[i][j].Qw_old;*/

    } 
      else {
	cells[i][j].bQw = 0.0;
	cells[i][j].Qstack = 0.0;
	cells[i][j].water2surf += cells[i][j].bwater; //Maxime
	cells[i][j].bwater = 0.0; //Maxime
	
      }

    }


  }/*end calculation of discharge*/ 

}


void glacial_hydrology(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,iproptype iprop,hwproptype hwprop,double dt)
{

  int i,j,k;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int nn = nx*ny;
  int nice;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;
  double dl = sqrt(dx*dy);
  double area = (*mesh).dx*(*mesh).dy; //Maxime
  double spy = 365.25*24.0*3600.0;
  double g = (*mesh).gravity;
  double rho = iprop.rho;
  double Lh = iprop.latentheat;
  double S0 = hwprop.S0;
  double A0 = hwprop.A0;
  double kh = hwprop.kh;
  double ds = hwprop.ds; 
  double ss0 = hwprop.ss0; 
  double h0 = hwprop.h0;
  double rho_w = 1000.0;
  double kmin = hwprop.kmin; 
  double Ls = hwprop.Ls; /*cavity spacing*/ 
  double Lc = hwprop.Lc; /*channel spacing*/
  double alpha = hwprop.alpha;
  double minqw = hwprop.minqw;
  double B = 73.3e6;
  double hw,S,sslope,hs;
  double meante = 0.0;
  double meanhw = 0.0;
  double meandhwdt = 0.0;
  double meanwater = 0.0;
  double meansliding = 0.0;
  double temp,trans,dhw;
  double te_s = 0.0;
  double te_c = 0.0;
  double te_ss = 0.0;
  double vout;
  double minice = 10.0;
  int dosediment = (*mesh).dosediment; //Maxime
  double sc = 1.5; // supercooling slope - Maxime
  double sflux = 0.0; //Maxime
  double dss = hwprop.dss; //Subglacial carrying capacity - Maxime
  double fp = 1.0; //porosity of englacial+subglacial - Maxime
  double sgift = 0.0; //Maxime
  double maxgift = 0.0; //Maxime
  long tn,ni,nj,ri,rj; //Maxime
  double slope,mslope,mdist; //Maxime
  double islope = 0.0; //Maxime
  double bslope = 0.0; //Maxime
  
#pragma omp parallel shared(cells,hp,vp) private(i,j,hw,S,sslope,hs,temp,trans,dhw,te_s,te_c,te_ss) firstprivate(nx,ny,dx,dy,dl,rho,g,rho_w,h0,kh,ds,ss0,B,alpha,kmin,spy,Ls,Lc,Lh,dt,S0,A0,minqw)
  {


    /*compute steady state effective pressure in cavities and channels*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) { 
	cells[i][j].qw = cells[i][j].bQw/dl; /*m2 pr. sec - computed from meltwater in get_discharge*/
	if (cells[i][j].qw < minqw) cells[i][j].qw = minqw;
	cells[i][j].psi = rho*g*(cells[i][j].hslope+1.0e-3); /*water head gradient - kg m-2 s-2 - follow ice surface gradient*/
	cells[i][j].mw = cells[i][j].qw*cells[i][j].psi/(rho*Lh); /*rate of melt water by friction - *m s-1*/
 
	if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) { 
 
	  /*cavity length and pressure*/
	  sslope = cells[i][j].slidingslope; 
	  if (sslope <= ss0) hs = kh*exp(sslope/ds)+h0; //bed step height
	  else hs = Ls*sslope;
	  cells[i][j].hs = hs;
  
	  /*cavity size and steady state effective pressure*/
	  temp = cells[i][j].S;
	  cells[i][j].S = pow(Ls*cells[i][j].qw/(kmin*sqrt(cells[i][j].psi)),0.8)/(alpha*hs);
	  if (cells[i][j].S < S0) cells[i][j].S = S0;
	  if (cells[i][j].S > Ls) cells[i][j].S = Ls-Ls*0.05; //Maxime - limit the size of the cavities
	  cells[i][j].SLf = cells[i][j].S/Ls;
	  temp = (cells[i][j].sliding*hs)/spy;
	  te_s = 1.0/(rho*g)*B*pow(16.0/(2.0*pi)*temp/pow(cells[i][j].S,2.0),1.0/3.0); 

	  /*channel cross section and steady-state effective pressure*/ 
	  temp = cells[i][j].Ac;
	  cells[i][j].Ac = pow(Lc*cells[i][j].qw/(kmin*sqrt(cells[i][j].psi)),0.8);
	  cells[i][j].dAcdt = (cells[i][j].Ac - temp)/dt;
	  temp =  Lc*cells[i][j].mw;
	  te_c = 1.0/(rho*g)*3.0*B*pow(temp/(2.0*(cells[i][j].Ac+A0)),1.0/3.0); 
 
	  /*if cavity drained*/	  
	  if (te_s > te_c) {
	    cells[i][j].te = te_s;
	    cells[i][j].hydro = -1;
	  }
	  else {
	    cells[i][j].te = te_c;
	    cells[i][j].hydro = 1;
	  }

	  /*limit effective pressure*/
	  if (cells[i][j].te > cells[i][j].tn) cells[i][j].te = cells[i][j].tn;
	  if (cells[i][j].te < 0.01*cells[i][j].tn) cells[i][j].te = 0.01*cells[i][j].tn;
 
	  /*hydrological head*/
      cells[i][j].Hgw = cells[i][j].bed + cells[i][j].tn - cells[i][j].te;

	  /*water pressure*/
	  cells[i][j].Pw = cells[i][j].tn - cells[i][j].te;
	  
	}
	else { 
	  cells[i][j].Hgw = cells[i][j].bed;
	  cells[i][j].te = 0.01*cells[i][j].tn;
	  cells[i][j].S = 0.0;
	  cells[i][j].SLf = 0.0;
	  cells[i][j].hydro = 0;
	}
      }
    } 

  }/*pragma*/

  /*compute mean effective pressure*/
  nice = 1;
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
	  /*sediment transport*/
	  if ((dosediment > 0)&&(cells[i][j].sedi > 0.0)) {	
		/* find steepest neighbor */
		tn = 0; ni = cells[i][j].ne_i[0]; nj = cells[i][j].ne_j[0];
		ri = ni; rj = nj;
		mslope = (cells[i][j].bwatersurf-cells[ni][nj].bwatersurf)/cells[i][j].ndist[0];
		for (k=1;k<cells[i][j].nne;k++) {
		  ni = cells[i][j].ne_i[k]; nj = cells[i][j].ne_j[k];
		  slope = (cells[i][j].bwatersurf-cells[ni][nj].bwatersurf)/cells[i][j].ndist[k];
		  if (slope > mslope) {
			tn = k;
			ri = ni;
			rj = nj;
			mdist = cells[i][j].ndist[k];
			mslope = slope;
		  }/*if*/
		}/*k*/
		/*bed and ice surface slopes in the direction of water flow*/
		bslope = ((cells[i][j].bed + cells[i][j].sedi) - (cells[ri][rj].bed + cells[ri][rj].sedi)) / mdist;
		islope = ((cells[i][j].bed + cells[i][j].sedi + cells[i][j].ice) - (cells[ri][rj].bed + cells[ri][rj].sedi + cells[ri][rj].ice)) / mdist;
		
		/*if meltwater supercooling - do not transport sediments*/
		if ((islope > 0.0)&&(sc >= 0.0)&&(bslope < -sc*islope)) {
			cells[i][j].coolfac = 1.0;
		}
		else { /* transport */
			
			/*coolfac[i] += wflux*bprop[1]/nneighbours[i];*/
			
			/*sediment flux per second*/
			sflux = dss*pow(cells[i][j].bQw*fp, 3.0);
			
			/*sediment package per year*/
			sgift = dt*s_per_y*sflux;
			
			/*penalties*/
			maxgift = (cells[i][j].sedi - cells[ri][rj].sedi)*area;
			
			if (maxgift < 0.0) maxgift = 0.0;
			if (sgift > maxgift) sgift = maxgift;
			
			/*pass sediment*/
			cells[i][j].sedi -= sgift/area;
			cells[ri][rj].sedi += sgift/area;
			cells[i][j].subglacial_erosion = sgift/area;
			
		}/*else*/
		
	 }/*if sedi*/

    /*if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) {*/
	meante += cells[i][j].te;
	meanwater += cells[i][j].bwater*dt;
	meansliding += cells[i][j].sliding;
	nice += 1;
	/*    }*/
    }
  }
  meante /= (double)(nice);
  meanwater /= (double)(nice);
  meansliding /= (double)(nice);
  (*mesh).meante = meante;
  (*mesh).meanwater = meanwater;
  (*mesh).meansliding = meansliding;

}/*void*/



void glacial_hydrology_Iverson_transient(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,iproptype iprop,hwproptype hwprop,double dt,double time)
{

  int i,j;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int nice;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;

  double spy = 365.25*24.0*3600.0;
  double g = (*mesh).gravity;
  double rho = iprop.rho;
  double Lh = iprop.latentheat;
  double S0 = hwprop.S0;
  double kh = hwprop.kh;
  double ds = hwprop.ds; 
  double ss0 = hwprop.ss0; 
  double h0 = hwprop.h0;
  double rho_w = 1000.0;
  double kmin = hwprop.kmin; 
  double L = hwprop.Ls;
  double alpha = hwprop.alpha;
  double maxhw = 1.0;
  double minhw = 0.1;
  double tscale = hwprop.tscale;
  double B = 73.3e6;
  double Rbed,hw,S,sslope,hs;
  double meante = 0.0;
  double meanhw = 0.0;
  double meandhwdt = 0.0;
  double meanwater = 0.0;
  double meansliding = 0.0;
  double temp,dhw,Kgw,dSdt,qwpsi;
  double te_ss = 0.0;


#pragma omp parallel shared(cells,hp,vp) private(i,j,Rbed,hw,S,sslope,hs,temp,dhw,te_ss,Kgw,qwpsi,dSdt) firstprivate(nx,ny,dx,dy,rho,g,rho_w,h0,kh,ds,ss0,minhw,maxhw,tscale,B,alpha,kmin,spy,L,Lh,dt,S0,time)
  {


    /*compute water flux qwx m2/y at time t*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=1;j<nx;j++) {
	if (cells[i][j-1].Hgw > cells[i][j].Hgw) {
	  hw = cells[i][j-1].hw; 
	}
	else {
	  hw = cells[i][j].hw;  
	}
	if (hw > 1.0) hw = 1.0; if (hw < 0.0) hw = 0.0;
	hp[i][j].psi = rho_w*g*(cells[i][j].Hgw-cells[i][j-1].Hgw)/dx;
	dhw = -dt*spy*tscale*kmin*pow(L*hw,1.2)*hp[i][j].psi/(sqrt(fabs(hp[i][j].psi)+10.0)*dx*L);
	if (dhw < 0.0) {
	  if (cells[i][j].hw <= 0.0) dhw = 0.0;
	  else if (dhw < -cells[i][j].hw) dhw = -cells[i][j].dhw;
	} 
	else {
	  if (cells[i][j-1].hw <= 0.0) dhw = 0.0;
	  else if (dhw > cells[i][j-1].hw) dhw = cells[i][j-1].hw;
	}
	hp[i][j].dhw = dhw;
      }
    }
    

    /*compute water flux qw m2/yr at time t*/
#pragma omp for schedule(static)
    for (i=1;i<ny;i++) {
      for (j=0;j<nx;j++) {
	if (cells[i-1][j].Hgw > cells[i][j].Hgw) { 
	  hw = cells[i-1][j].hw; 
	}
	else {
	  hw = cells[i][j].hw; 
	}
	if (hw < 0.0) hw = 0.0; if (hw > 1.0) hw = 1.0;
	vp[i][j].psi = rho_w*g*(cells[i][j].Hgw-cells[i-1][j].Hgw)/dy;
	dhw = -dt*spy*tscale*kmin*pow(L*hw,1.2)*vp[i][j].psi/(sqrt(fabs(vp[i][j].psi)+10.0)*dy*L);
	if (dhw < 0.0) {
	  if (cells[i][j].hw <= 0.0) dhw = 0.0; 
	  else if (dhw < -cells[i][j].hw) dhw = -cells[i][j].hw;
	}
	else {
	  if (cells[i-1][j].hw <= 0.0) dhw = 0.0;
	  else if (dhw > cells[i-1][j].hw) dhw = cells[i-1][j].hw;
	}
	vp[i][j].dhw = dhw;
      }
    }

    /*compute change in water sheet thickness and efffective pressure*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) { 
	  cells[i][j].qw = 0.5*(fabs(hp[i][j].dhw)*dx + fabs(hp[i][j+1].dhw)*dx + fabs(vp[i][j].dhw)*dy + fabs(vp[i+1][j].dhw)*dy)/dt;
	  qwpsi = 0.5*dx*(fabs(hp[i][j].dhw*hp[i][j].psi) + fabs(hp[i][j+1].dhw*hp[i][j+1].psi)) + 0.5*dy*(fabs(vp[i][j].dhw*vp[i][j].psi) + fabs(vp[i+1][j].dhw*vp[i+1][j].psi))/dt;
	  cells[i][j].mw = qwpsi/(tscale*rho*Lh);
	  cells[i][j].dhw = hp[i][j].dhw - hp[i][j+1].dhw + vp[i][j].dhw - vp[i+1][j].dhw + tscale*cells[i][j].water*dt; 
	  cells[i][j].dhwdt = cells[i][j].dhw/dt; 
	  /*if (cells[i][j].dhwdt < -1.0) cells[i][j].dhwdt = -1.0;
	    if (cells[i][j].dhwdt > 1.0) cells[i][j].dhwdt = 1.0;*/
	  cells[i][j].hw += dt*cells[i][j].dhwdt; /*at time t + dt*/
	  if (cells[i][j].hw < minhw) cells[i][j].hw = minhw;
	  if (cells[i][j].hw > maxhw) cells[i][j].hw = maxhw;
	  sslope = cells[i][j].slidingslope; 
	  if (sslope <= ss0) hs = kh*exp(sslope/ds)+h0;
	  else hs = L*sslope;
	  cells[i][j].hs = hs;
	  cells[i][j].S = L*cells[i][j].hw/(alpha*hs); 
	  cells[i][j].dSdt = L*cells[i][j].dhwdt/(alpha*hs); 
	  if (cells[i][j].S < 0.0) cells[i][j].S = 0.0;
	  /*if (cells[i][j].S > 2.0) cells[i][j].S = 2.0;*/
	  temp = (cells[i][j].sliding*hs)/spy;
	  te_ss = 1.0/(rho*g)*B*pow(16.0/(2.0*pi)*temp/pow(cells[i][j].S+S0,2.0),1.0/3.0); 
	  if (te_ss > cells[i][j].tn) te_ss = cells[i][j].tn;
	  temp += L*(cells[i][j].mw-cells[i][j].Mb)/spy;
	  if (time > 0.1) temp -= 0.0*alpha*hs*cells[i][j].dSdt/spy;
	  if (temp > 0.0) {
	    cells[i][j].te_new = 1.0/(rho*g)*B*pow(16.0/(2.0*pi)*temp/pow(cells[i][j].S+S0,2.0),1.0/3.0); 
	  }
	  else cells[i][j].te_new = 0.01*cells[i][j].tn;
	  if (cells[i][j].te_new < 0.01*cells[i][j].tn) cells[i][j].te_new = 0.01*cells[i][j].tn;
	  if (cells[i][j].te_new > cells[i][j].tn) cells[i][j].te_new = cells[i][j].tn;

	  cells[i][j].te = cells[i][j].te_new;
	  cells[i][j].Hgw = cells[i][j].bed + cells[i][j].tn - te_ss;/*cells[i][j].te;*/
	  cells[i][j].Pw = cells[i][j].tn - cells[i][j].te;
	}
	else { 
	  cells[i][j].hw = minhw;
	  cells[i][j].dhwdt = 0.0;
	  cells[i][j].Hgw = cells[i][j].bed;
	  cells[i][j].te = 0.01*cells[i][j].tn;
	  cells[i][j].S = 0.0;
	}
      }
    } 

  }/*pragma*/

  /*compute mean effective pressure*/
  nice = 1;
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      /*if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) {*/
	meante += cells[i][j].te;
	meanhw += cells[i][j].hw;
	meandhwdt += cells[i][j].dhwdt;
	meanwater += cells[i][j].water;
	meansliding += cells[i][j].sliding;
	nice += 1;
	/*    }*/
    }
  }
  meante /= (double)(nice);
  meanhw /= (double)(nice);
  meandhwdt /= (double)(nice);
  meanwater /= (double)(nice);
  meansliding /= (double)(nice);
  (*mesh).meante = meante;
  (*mesh).meanhw = meanhw;
  (*mesh).meandhwdt = meandhwdt;
  (*mesh).meanwater = meanwater;
  (*mesh).meansliding = meansliding;

}/*void*/


void glacial_hydrology_Iverson_steady(celltype **cells,meshtype *mesh,iproptype iprop,hwproptype hwprop)
{

  int i,j;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;

  double spy = 365.25*24.0*3600.0;
  double g = (*mesh).gravity;
  double rho = iprop.rho;
  double S0 = hwprop.S0;
  double kh = hwprop.kh;
  double h0 = hwprop.h0;
  double L = hwprop.Ls;
  double alpha = hwprop.alpha;
  double rho_w = 1000.0;
  double k0 = hwprop.kmin; 
  double maxhw = 1.0;
  double minhw = 0.1;
  double tscale = hwprop.tscale;
  double B = 73.3e6;
  double tb0 = 0.001;
  double Rbed,meante,hw,S,sslope,hs,psi;

#pragma omp parallel shared(cells) private(i,j,Rbed,hw,S,sslope,hs,psi) firstprivate(nx,ny,dx,dy,k0,rho,g,rho_w,h0,kh,minhw,maxhw,tscale,B,spy,alpha,L,tb0)
  {


    /*compute change in water sheet thickness at time t and efffective pressure at time t+dt/2*/
#pragma omp for schedule(static)
    for (i=1;i<(ny-1);i++) {
      for (j=1;j<(nx-1);j++) {
	if (cells[i][j].ice > 10.0) {
	  /*cells[i][j].qw = qfac*cells[i][j].ice*(cells[i][j].deformation+cells[i][j].sliding);*/
	  /*psi = rho*g*(cells[i][j].tslope + tb0); 
	    cells[i][j].hw = pow(cells[i][j].qw*vis_w/(k0*psi),1.0/3.0);*/ 
	  cells[i][j].hw = 0.1*cells[i][j].water;
	  if (cells[i][j].hw < minhw) cells[i][j].hw = minhw;
	  if (cells[i][j].hw > maxhw) cells[i][j].hw = maxhw;
	  sslope = cells[i][j].slidingslope; if (sslope < 0.0) sslope = 0.0; 
	  hs = kh*sslope + h0; 
	  Rbed = hs/L;
	  S = cells[i][j].hw/(alpha*Rbed); if (S < 1.0e-3) S = 1.0e-3;
	  cells[i][j].te = 1.0/(rho*g)*B*pow((16.0*cells[i][j].sliding/spy*hs/(2.0*pi))/(S*S),1.0/3.0); /*at time t + dt*/ 
	  if (cells[i][j].te < 0.1*cells[i][j].tn) cells[i][j].te = 0.1*cells[i][j].tn;
	  if (cells[i][j].te > cells[i][j].tn) cells[i][j].te = cells[i][j].tn;
	  cells[i][j].Hgw = cells[i][j].bed + cells[i][j].tn - cells[i][j].te;
	}
	else { 
	  cells[i][j].hw = minhw;
	  cells[i][j].dhwdt = 0.0;
	  cells[i][j].Hgw = cells[i][j].bed;
	  cells[i][j].te = 0.5*cells[i][j].tn;
	}
      }
    } 

#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      cells[i][0].te = cells[i][1].te;
      cells[i][nx-1].te = cells[i][nx-2].te;
    }

#pragma omp for schedule(static)
    for (j=0;j<nx;j++) {
      cells[0][j].te = cells[1][j].te;
      cells[ny-1][j].te = cells[ny-2][j].te;
    }

  }/*pragma*/

  /*compute mean effective pressure*/
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      meante += cells[i][j].te;
    }
  }
  meante /= (double)((*mesh).nc);
  (*mesh).meante = meante;

}/*void*/


void glacial_hydrology_transient(celltype **cells,hptype **hp,vptype **vp,meshtype *mesh,iproptype iprop,hwproptype hwprop,double dt,double time)
{

  int i,j;
  int nx = (*mesh).nx;
  int ny = (*mesh).ny;
  int nice;
  double dx = (*mesh).dx;
  double dy = (*mesh).dy;

  double spy = 365.25*24.0*3600.0;
  double g = (*mesh).gravity;
  double rho = iprop.rho;
  double Lh = iprop.latentheat;
  double S0 = hwprop.S0;
  double A0 = hwprop.A0;
  double kh = hwprop.kh;
  double ds = hwprop.ds; 
  double ss0 = hwprop.ss0; 
  double h0 = hwprop.h0;
  double rho_w = 1000.0;
  double kmin = hwprop.kmin; 
  double Ls = hwprop.Ls; /*cavity spacing*/ 
  double Lc = hwprop.Lc; /*channel spacing*/
  double alpha = hwprop.alpha;
  double maxhw = 1.0;
  double minhw = 0.1;
  double tscale = hwprop.tscale;
  double B = 73.3e6;
  double Rbed,hw,S,sslope,hs;
  double meante = 0.0;
  double meanhw = 0.0;
  double meandhwdt = 0.0;
  double meanwater = 0.0;
  double meansliding = 0.0;
  double temp,dhw,Kgw,dSdt,qwpsi;
  double te_s = 0.0;
  double te_c = 0.0;
  double te_ss = 0.0;

#pragma omp parallel shared(cells,hp,vp) private(i,j,Rbed,hw,S,sslope,hs,temp,dhw,te_s,te_c,te_ss,Kgw,qwpsi,dSdt) firstprivate(nx,ny,dx,dy,rho,g,rho_w,h0,kh,ds,ss0,minhw,maxhw,tscale,B,alpha,kmin,spy,Ls,Lc,Lh,dt,S0,A0,time)
  {


    /*compute water flux qwx m2/y at time t*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=1;j<nx;j++) {
	if (cells[i][j-1].Hgw > cells[i][j].Hgw) {
	  hw = cells[i][j-1].hw; 
	}
	else {
	  hw = cells[i][j].hw;  
	}
	if (hw > 1.0) hw = 1.0; if (hw < 0.0) hw = 0.0;
	hp[i][j].psi = rho_w*g*(cells[i][j].Hgw-cells[i][j-1].Hgw)/dx;
	dhw = -dt*spy*tscale*kmin*pow(Ls*hw,1.2)*hp[i][j].psi/(sqrt(fabs(hp[i][j].psi)+10.0)*dx*Ls);
	if (dhw < 0.0) {
	  if (cells[i][j].hw <= 0.0) dhw = 0.0;
	  else if (dhw < -cells[i][j].hw) dhw = -cells[i][j].dhw;
	} 
	else {
	  if (cells[i][j-1].hw <= 0.0) dhw = 0.0;
	  else if (dhw > cells[i][j-1].hw) dhw = cells[i][j-1].hw;
	}
	hp[i][j].dhw = dhw;
      }
    }
    

    /*compute water flux qw m2/yr at time t*/
#pragma omp for schedule(static)
    for (i=1;i<ny;i++) {
      for (j=0;j<nx;j++) {
	if (cells[i-1][j].Hgw > cells[i][j].Hgw) { 
	  hw = cells[i-1][j].hw; 
	}
	else {
	  hw = cells[i][j].hw; 
	}
	if (hw < 0.0) hw = 0.0; if (hw > 1.0) hw = 1.0;
	vp[i][j].psi = rho_w*g*(cells[i][j].Hgw-cells[i-1][j].Hgw)/dy;
	dhw = -dt*spy*tscale*kmin*pow(Ls*hw,1.2)*vp[i][j].psi/(sqrt(fabs(vp[i][j].psi)+10.0)*dy*Ls);
	if (dhw < 0.0) {
	  if (cells[i][j].hw <= 0.0) dhw = 0.0; 
	  else if (dhw < -cells[i][j].hw) dhw = -cells[i][j].hw;
	}
	else {
	  if (cells[i-1][j].hw <= 0.0) dhw = 0.0;
	  else if (dhw > cells[i-1][j].hw) dhw = cells[i-1][j].hw;
	}
	vp[i][j].dhw = dhw;
      }
    }

    /*compute change in water sheet thickness and efffective pressure*/
#pragma omp for schedule(static)
    for (i=0;i<ny;i++) {
      for (j=0;j<nx;j++) {
	cells[i][j].qw = 0.5*(fabs(hp[i][j].dhw + hp[i][j+1].dhw)*dx + fabs(vp[i][j].dhw + vp[i+1][j].dhw)*dy)/dt;
	qwpsi = 0.5*dx*(fabs(hp[i][j].dhw*hp[i][j].psi) + fabs(hp[i][j+1].dhw*hp[i][j+1].psi)) + 0.5*dy*(fabs(vp[i][j].dhw*vp[i][j].psi) + fabs(vp[i+1][j].dhw*vp[i+1][j].psi))/dt;
	cells[i][j].mw = qwpsi/(tscale*rho*Lh);
	cells[i][j].dhw = hp[i][j].dhw - hp[i][j+1].dhw + vp[i][j].dhw - vp[i+1][j].dhw + tscale*cells[i][j].water*dt; 
	cells[i][j].dhwdt = cells[i][j].dhw/dt; 
	/*if (cells[i][j].dhwdt < -1.0) cells[i][j].dhwdt = -1.0;
	  if (cells[i][j].dhwdt > 1.0) cells[i][j].dhwdt = 1.0;*/
	cells[i][j].hw += dt*cells[i][j].dhwdt; /*at time t + dt*/
	/*if (cells[i][j].hw < minhw) cells[i][j].hw = minhw;*/
	/*if (cells[i][j].hw > maxhw) cells[i][j].hw = maxhw;*/
 
	if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) { 
 
	  /*cavity length and pressure*/
	  sslope = cells[i][j].slidingslope; 
	  if (sslope <= ss0) hs = kh*exp(sslope/ds)+h0;
	  else hs = Ls*sslope;
	  cells[i][j].hs = hs;
	  cells[i][j].S = Ls*cells[i][j].hw/(alpha*hs); 
	  cells[i][j].dSdt = Ls*cells[i][j].dhwdt/(alpha*hs); 
	  if (cells[i][j].S < 0.0) cells[i][j].S = 0.0;
	  /*if (cells[i][j].S > 2.0) cells[i][j].S = 2.0;*/
	  temp = (cells[i][j].sliding*hs)/spy;
	  te_ss = 1.0/(rho*g)*B*pow(16.0/(2.0*pi)*temp/pow(cells[i][j].S+S0,2.0),1.0/3.0); 
	  temp -= alpha*hs*cells[i][j].dSdt/spy;
	  te_s = 1.0/(rho*g)*B*pow(16.0/(2.0*pi)*temp/pow(cells[i][j].S+S0,2.0),1.0/3.0); 

	  /*channel cross section and pressure*/
	  cells[i][j].Ac = Lc*cells[i][j].hw;
	  cells[i][j].dAcdt = Lc*cells[i][j].dhwdt;
	  temp =  Lc*(cells[i][j].mw-cells[i][j].Mb)/spy;
	  /*temp -= cells[i][j].dAcdt/spy;*/
	  te_c = 1.0/(rho*g)*3.0*B*pow(temp/(2.0*(cells[i][j].Ac+A0)),1.0/3.0); 
 
	  /*if cavity drained*/	  
	  if (te_s > te_c) {
	    cells[i][j].te = te_s;
	    cells[i][j].hydro = -1;
	  }
	  else {
	    cells[i][j].te = te_c;
	    cells[i][j].hydro = 1;
	  }

	  if (te_ss > cells[i][j].tn) te_ss = cells[i][j].tn;
	  if (te_ss < 0.01*cells[i][j].tn) te_ss = 0.01*cells[i][j].tn;


	  if (cells[i][j].te > cells[i][j].tn) cells[i][j].te = cells[i][j].tn;
	  if (cells[i][j].te < 0.01*cells[i][j].tn) cells[i][j].te = 0.01*cells[i][j].tn;
 
	  cells[i][j].te_s = te_ss;
	  cells[i][j].Hgw = cells[i][j].bed + cells[i][j].tn - te_ss;
          /*cells[i][j].Hgw = cells[i][j].bed + cells[i][j].tn - cells[i][j].te;*/
	  cells[i][j].Pw = cells[i][j].tn - cells[i][j].te;
	}
	else { 
	  /*cells[i][j].hw = minhw;*/
	  /*cells[i][j].dhwdt = 0.0;*/
	  cells[i][j].Hgw = cells[i][j].bed + 1.0e5*cells[i][j].hw;
	  cells[i][j].te = 0.01*cells[i][j].tn;
	  cells[i][j].S = 0.0; 
	  cells[i][j].hydro = 0;
	}
      }
    } 

  }/*pragma*/

  /*compute mean effective pressure*/
  nice = 1;
  for (i=0;i<ny;i++) {
    for (j=0;j<nx;j++) {
      /*if ((cells[i][j].ice > 10.0)&&(cells[i][j].margin < 0)) {*/
	meante += cells[i][j].te;
	meanhw += cells[i][j].hw;
	meandhwdt += cells[i][j].dhwdt;
	meanwater += cells[i][j].water;
	meansliding += cells[i][j].sliding;
	nice += 1;
	/*    }*/
    }
  }
  meante /= (double)(nice);
  meanhw /= (double)(nice);
  meandhwdt /= (double)(nice);
  meanwater /= (double)(nice);
  meansliding /= (double)(nice);
  (*mesh).meante = meante;
  (*mesh).meanhw = meanhw;
  (*mesh).meandhwdt = meandhwdt;
  (*mesh).meanwater = meanwater;
  (*mesh).meansliding = meansliding;

}/*void*/
